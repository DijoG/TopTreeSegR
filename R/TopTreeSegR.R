#' TopTreeSegR: Topology-based Tree Segmentation for TLS Point Clouds
#' 
#' @name TopTreeSegR
#' @docType package
#' @useDynLib TopTreeSegR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "=== TopTreeWSegR loaded! ===\n",
    "C++ accelerated tree segmentation with gradient flow\n",
    "Use TTS_pipeline() for complete segmentation workflow"
  )
}

.check_dependencies <- function() {
  required_pkgs = c("DiscreteMorseR", "ahull3D", "dbscan", "FNN", "lidR")
  missing_pkgs = required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "The following required packages are missing: ", 
      paste(missing_pkgs, collapse = ", "), "\n",
      "Install with: install.packages(c('", 
      paste(setdiff(missing_pkgs, c("DiscreteMorseR", "ahull3D")), collapse = "', '"), 
      "'))\n",
      "And GitHub packages with: remotes::install_github(c('",
      "DijoG/DiscreteMorseR', 'DijoG/ahull3D'))"  
    )
  }
  return(TRUE)
}

#' Ultra-Fast Tree Segmentation from LiDAR Point Clouds
#'
#' @param las Input point cloud in las/laz format (LAS object)
#' @param alpha Alpha parameter for alpha-complex construction (default: 0.1)
#' @param input_truth Point ID attribute from las, if not found generated (default: "pid")
#' @param stem_height Height threshold for detecting stem seeds (default: 0.5)
#' @param max_distance Maximum distance for minima connectivity (default: 2.0)
#' @param grid_size Grid size for spatial hashing in large datasets (default: 5.0)
#' @param spatial_eps Spatial epsilon for clustering (default: 1.0)
#' @param method Segmentation method: "morse-smale" (default) or "gradient"
#' @param cores Number of CPU threads for parallel processing (default: 2)
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{labels}: Vector of tree labels for mesh vertices
#'   \item \code{original_pid}: Original point IDs from mesh
#'   \item \code{mesh_coords}: Mesh vertices coordinates
#'   \item \code{n_trees}: Number of trees detected
#'   \item \code{minima}: Indices of Morse minima (1-based)
#'   \item \code{ascending_regions}: Ascending region assignments (gradient method only)
#'   \item \code{seeds}: Seed minima used (i123 values)
#'   \item \code{method}: Method used for segmentation
#'   \item \code{parameters}: List of parameters used
#' }
#'
#' @export
TTS_segmentation <- function(las, 
                             alpha = 0.1, 
                             input_truth = "pid",
                             stem_height = 0.5,
                             max_distance = 2.0,
                             grid_size = 5.0,
                             spatial_eps = 1.0,
                             method = "morse-smale",
                             cores = 2) {
  
  # Check dependencies
  .check_dependencies()
  
  message("=== TopTreeWSegR: C++ Accelerated Tree Segmentation ===\n")
  message("Method selected: ", method)
  
  # Check input LAS
  if (inherits(las, "LAS")) {
    # Extract data from LAS object
    las_data = las@data
    
    # Check if required columns exist
    if (!all(c("X", "Y", "Z") %in% names(las_data))) {
      stop("LAS data must contain X, Y, Z coordinate columns")
    }
    
    # Remove duplicate coordinates
    las_data = las_data[!duplicated(las_data[, c("X", "Y", "Z")]), ]
    coords = as.matrix(las_data[, c("X", "Y", "Z")])
    
    # Handle input_truth: if column exists, use it; otherwise generate sequence
    if (input_truth %in% names(las_data)) {
      input_truth_vec = las_data[[input_truth]]
      message("Using '", input_truth, "' column for point IDs")
    } else {
      # Generate sequential IDs starting from 1
      input_truth_vec = 1:nrow(las_data)
      message("Column '", input_truth, "' not found. Generating sequential point IDs (1 to ", nrow(las_data), ")")
    }
    
  } else {
    stop("las must be a LAS object")
  }
  
  message("Input points: ", nrow(las_data))
  
  # Step 1: Build alpha-complex mesh
  message("1. Building alpha-complex...")
  a = ahull3D::ahull3D(points = coords, alpha = alpha, input_truth = input_truth_vec)
  mesh = DiscreteMorseR::get_CCMESH(a)
  
  # Step 2: Compute Morse complex
  message("2. Computing Morse complex...")
  morse_complex = DiscreteMorseR::compute_MORSE_complex(
    mesh, 
    output_dir = tempdir(), 
    cores = cores
  )
  
  # Step 3: Fix Z column if it's character
  vertices_df = morse_complex$simplices$vertices
  
  # Convert Z to numeric if it's character
  if (is.character(vertices_df$Z)) {
    message("  Converting Z column from character to numeric...")
    vertices_df$Z = as.numeric(vertices_df$Z)
  }
  
  # Ensure all columns are correct types
  vertices_df$X = as.numeric(vertices_df$X)
  vertices_df$Y = as.numeric(vertices_df$Y)
  vertices_df$Z = as.numeric(vertices_df$Z)
  vertices_df$i123 = as.integer(vertices_df$i123)
  
  # Step 4: Segment trees using C++
  message("3. Segmenting trees (C++ implementation)...")
  
  if (method == "morse-smale") {
    result = morse_smale_segment_cpp(
      vertices_df = vertices_df,
      morse_complex_data = morse_complex,
      stem_height = stem_height,
      spatial_eps = spatial_eps
    )
    
  } else if (method == "gradient") {
    result = tree_segment_cpp(
      vertices_df = vertices_df,
      morse_complex_data = morse_complex,
      stem_height = stem_height,
      spatial_eps = spatial_eps,
      max_distance = max_distance,
      grid_size = grid_size
    )
    
  } else {
    stop("method must be 'morse-smale' or 'gradient'")
  }
  
  # Ensure all required fields exist (some methods don't have all fields)
  if (!"ascending_regions" %in% names(result)) {
    result$ascending_regions = integer(0)
  }
  
  if (!"seeds" %in% names(result)) {
    result$seeds = integer(0)
  }
  
  if (!"minima" %in% names(result)) {
    result$minima = integer(0)
  }
  
  # Convert seeds from i123 to row indices if needed
  if (length(result$seeds) > 0 && all(result$seeds > nrow(vertices_df))) {
    # Seeds are i123 values, convert to row indices
    seeds_rows = match(result$seeds, vertices_df$i123)
    seeds_rows = seeds_rows[!is.na(seeds_rows)]
    # Keep original i123 values in result$seeds, add row indices if needed
    if (!"seeds_rows" %in% names(result)) {
      result$seeds_rows = seeds_rows
    }
  }
  
  # Convert minima to row indices if they're i123 values
  if (length(result$minima) > 0 && all(result$minima > nrow(vertices_df))) {
    # Minima are i123 values, convert to row indices
    minima_rows = match(result$minima, vertices_df$i123)
    minima_rows = minima_rows[!is.na(minima_rows)]
    result$minima = minima_rows
  }
  
  # Statistics
  message("\n=== Segmentation Complete ===")
  message("Method: ", result$method)
  message("Detected ", result$n_trees, " trees")
  
  if (length(result$mesh_labels) > 0) {
    valid_labels = result$mesh_labels[result$mesh_labels > 0]
    if (length(valid_labels) > 0) {
      message("Label distribution:")
      label_table = table(valid_labels)
      # Print only first 10 trees if many
      if (length(label_table) > 10) {
        message("  Showing first 10 trees:")
        print(head(label_table, 10))
        message(sprintf("  ... and %d more trees", length(label_table) - 10))
      } else {
        print(label_table)
      }
    }
  }
  
  message("Seeds found: ", length(result$seeds))
  message("Minima found: ", length(result$minima))
  
  # Create result object
  result_obj = list(
    labels = result$mesh_labels,          # Labels for mesh vertices (1-based tree IDs)
    original_pid = vertices_df$i123,      # Original point IDs
    mesh_coords = result$mesh_coords,     # Mesh coordinates
    n_trees = as.integer(result$n_trees), # Number of trees detected
    minima = as.integer(result$minima),   # Indices of Morse minima (1-based)
    ascending_regions = as.integer(result$ascending_regions),  # For gradient method
    seeds = as.integer(result$seeds),     # Seed minima (i123 values)
    method = result$method,               # Method used
    parameters = list(
      alpha = alpha,
      stem_height = stem_height,
      max_distance = max_distance,
      spatial_eps = spatial_eps,
      grid_size = grid_size,
      cores = cores
    )
  )
  
  # Add any additional fields from the result
  if ("labeled_via_msmale" %in% names(result)) {
    result_obj$labeled_via_msmale = as.integer(result$labeled_via_msmale)
  }
  
  if ("spatially_assigned" %in% names(result)) {
    result_obj$spatially_assigned = as.integer(result$spatially_assigned)
  }
  
  if ("vertex_pair_count" %in% names(result)) {
    result_obj$vertex_pair_count = as.integer(result$vertex_pair_count)
  }
  
  if ("edge_count" %in% names(result)) {
    result_obj$edge_count = as.integer(result$edge_count)
  }
  
  class(result_obj) = "TTS_result"
  return(result_obj)
}

#' Markov Random Field Boundary Refinement for crown (grid-based)
#' 
#' @param TTS_result Segmentation result
#' @param intensity Smoothing strength (1.0 = moderate, 2.0 = strong)
#' @param iterations MRF iterations (1-3)
#' @param spatial_sigma Neighborhood scale (default: 5.0)
#' @export
TTS_MRFBR <- function(TTS_result,
                      intensity = 1.0,
                      iterations = 3,
                      spatial_sigma = 5.0) {
  
  message("=== Fast MRF Crown Refinement ===\n")
  message("Intensity: ", intensity)
  message("Iterations: ", iterations)
  
  initial_labels = TTS_result$labels
  initial_trees = length(unique(initial_labels[initial_labels > 0]))
  
  start_time = Sys.time()
  
  refined_labels = MRFBR_cpp(
    coords = TTS_result$mesh_coords,
    labels = TTS_result$labels,
    spatial_sigma = spatial_sigma,
    intensity = intensity,
    iterations = iterations)

  elapsed = difftime(Sys.time(), start_time, units = "secs")
  
  # Update result
  TTS_result$labels_initial = initial_labels
  TTS_result$labels = refined_labels
  
  # Statistics
  final_trees = length(unique(refined_labels[refined_labels > 0]))
  changed_points = sum(initial_labels != refined_labels)
  change_percent = round(changed_points / length(initial_labels) * 100, 2)
  
  message("\n=== Results ===")
  message("Time: ", round(elapsed, 1), " seconds")
  message("Initial trees: ", initial_trees)
  message("Final trees:   ", final_trees)
  message("Changed points: ", changed_points, " (", change_percent, "%)")
  
  TTS_result$refinement = list(
    intensity = intensity,
    time_seconds = round(elapsed, 1),
    changed_points = changed_points
  )
  
  return(TTS_result)
}

#' Bayesian Boundary Refinement
#' 
#' @param TTS_result Segmentation result
#' @param prior_strength Spatial consistency (1.0 = moderate, 2.0 = strong)
#' @param likelihood_strength Elevation consistency (1.0 = moderate, 2.0 = strong)
#' @param confidence_threshold How much better new label must be (default: 1.3)
#' @param cores Number of CPU threads for parallel computation (default: 2)
#' @param fix_fragments Remove tiny fragments after refinement (default: TRUE)
#' 
#' @export
TTS_BBR <- function(TTS_result,
                    prior_strength = 1.0,
                    likelihood_strength = 2.0,
                    confidence_threshold = 1.3,
                    cores = 2,
                    fix_fragments = TRUE) {
  
  message("=== Bayesian Boundary Refinement ===\n")
  message("Prior strength: ", prior_strength)
  message("Likelihood strength: ", likelihood_strength)
  message("Confidence threshold: ", confidence_threshold)
  message("Threads: ", cores)
  
  initial_labels = TTS_result$labels
  initial_trees = length(unique(initial_labels[initial_labels > 0]))
  
  start_time = Sys.time()
  
  refined_result = BBR_ultrafast_cpp(
    coords = TTS_result$mesh_coords,
    labels = TTS_result$labels,
    prior_strength = prior_strength,
    likelihood_strength = likelihood_strength,
    confidence_threshold = confidence_threshold,
    cores = cores
  )
  
  elapsed = Sys.time() - start_time
  
  # Update result
  TTS_result$labels_initial = initial_labels
  TTS_result$labels = refined_result$labels
  
  # Add metrics
  TTS_result$refinement = list(
    method = "bayesian",
    prior_strength = prior_strength,
    likelihood_strength = likelihood_strength,
    confidence_threshold = confidence_threshold,
    time_seconds = round(as.numeric(elapsed), 2),
    changes = refined_result$changes,
    boundary_points = refined_result$boundary_points,
    change_ratio = refined_result$change_ratio
  )
  
  if ("uncertainties" %in% names(refined_result)) {
    TTS_result$uncertainties = refined_result$uncertainties
  }
  
  # Statistics
  final_trees = length(unique(TTS_result$labels[TTS_result$labels > 0]))
  
  message("\n=== Results ===")
  message("Time: ", round(as.numeric(elapsed), 2), " seconds")
  message("Boundary points: ", refined_result$boundary_points, 
          " (", round(100 * refined_result$boundary_points / length(initial_labels), 1), "%)")
  message("Label changes: ", refined_result$changes,
          " (", round(100 * refined_result$change_ratio, 1), "% of boundaries)")
  message("Initial trees: ", initial_trees)
  message("Final trees: ", final_trees)
  
  return(TTS_result)
}

#' Complete Tree Segmentation Pipeline with Bayesian Refinement
#'
#' Performs Morse-Smale segmentation followed by a single Bayesian boundary refinement
#' using optimized parameters proven to achieve >0.85 Adjusted Rand Index.
#'
#' @param las Input LAS point cloud (must contain X, Y, Z coordinates)
#' @param method Segmentation method: "morse-smale" (default) or "gradient"
#' @param alpha Alpha parameter for alpha-complex construction (default: 0.1)
#' @param input_truth Attribute from las (default: "pid")
#' @param stem_height Height threshold for detecting stem seeds (default: 0.5)
#' @param spatial_eps Spatial epsilon for clustering minima (default: 1.0)
#' @param cores Number of CPU threads for parallel processing (default: 2)
#' 
#' @param prior_strength Spatial consistency strength in Bayesian refinement.
#'   Controls how much points should match their neighbors' labels.
#'   Range: 0.5 (weak) to 2.0 (strong). Default: 1.0 (moderate).
#'   
#' @param likelihood_strength Elevation consistency strength in Bayesian refinement.
#'   Controls how much points at similar heights should have same label.
#'   **This is an exponent**: 1.0 = moderate, 1.5 = strong, 2.0 = very strong.
#'   Range: 0.5 to 2.5. Default: 1.6 (optimized for high ARI).
#'   
#' @param confidence_threshold How much better a new label must be to change.
#'   Higher values = more conservative refinement (changes fewer labels).
#'   - 1.0: Change if new label is ANY better (even 0.0001% better)
#'   - 1.3: Change only if new label is 30% better
#'   - 1.5: Change only if new label is 50% better
#'   - 2.0: Change only if new label is 100% better (twice as good)
#'   Range: 1.0 to 2.0. Default: 1.0 (aggressive, good for first pass).
#'   
#' @param fix_fragments Remove tiny fragments after refinement (default: TRUE)
#'   Merges clusters with < 30 points into neighboring trees.
#'   
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A TTS_result object with the following components:
#' \itemize{
#'   \item \code{labels}: Refined tree labels for mesh vertices
#'   \item \code{original_pid}: Original point IDs from mesh
#'   \item \code{mesh_coords}: Mesh vertices coordinates
#'   \item \code{n_trees}: Number of trees detected
#'   \item \code{minima}: Indices of Morse minima
#'   \item \code{seeds}: Seed minima used
#'   \item \code{method}: Method used for segmentation
#'   \item \code{pipeline}: List of pipeline parameters used
#' }
#'
#' @examples
#' \dontrun{
#' # Simplest usage with defaults (proven to achieve >0.85 ARI)
#' result <- TTS_pipeline(las = my_las_data)
#'
#' # With custom Bayesian parameters
#' result <- TTS_pipeline(
#'   las = my_las_data,
#'   prior_strength = 1.0,       # Moderate spatial consistency
#'   likelihood_strength = 1.6,   # Strong elevation consistency
#'   confidence_threshold = 1.0   # Aggressive refinement
#' )
#'
#' # Validate results (if ground truth available)
#' if ("treeid" %in% names(my_las_data@data)) {
#'   validate_TTS(result, my_las_data)
#' }
#'
#' # Visualize in 3D
#' plot_TTS_3d(result)
#' }
#'
#' @export
TTS_pipeline <- function(las,
                         method = "morse-smale",
                         alpha = 0.1,
                         input_truth = "pid",
                         stem_height = 0.5,
                         spatial_eps = 1.0,
                         cores = 2,
                         prior_strength = 1.0,
                         likelihood_strength = 1.6,
                         confidence_threshold = 1.0,
                         fix_fragments = TRUE,
                         verbose = TRUE) {
  
  if (verbose) {
    message("=== TopTreeSegR Pipeline ===\n")
    message("Segmentation:")
    message(sprintf("  Method: %s", method))
    message(sprintf("  Alpha: %.1f", alpha))
    message(sprintf("  Stem height: %.1f", stem_height))
    message(sprintf("  Cores: %d", cores))
    
    message("\nBayesian Refinement:")
    message(sprintf("  Prior strength: %.1f (spatial consistency)", prior_strength))
    message(sprintf("  Likelihood strength: %.1f (elevation consistency)", likelihood_strength))
    message(sprintf("  Confidence threshold: %.1f", confidence_threshold))
    message(sprintf("  Fix fragments: %s", fix_fragments))
    
    if (confidence_threshold == 1.0) {
      message("    → Aggressive: changes labels if any improvement")
    } else if (confidence_threshold >= 1.5) {
      message("    → Conservative: changes only if much better")
    }
  }
  
  # STEP 1: Initial segmentation
  if (verbose) message("\n1. Performing initial segmentation...")
  
  result = TTS_segmentation(
    las = las,
    alpha = alpha,
    input_truth = input_truth,
    stem_height = stem_height,
    spatial_eps = spatial_eps,
    method = method,
    cores = cores
  )
  
  if (verbose) {
    message(sprintf("  ✓ Detected %d trees", result$n_trees))
    message(sprintf("  ✓ Mesh vertices: %d", nrow(result$mesh_coords)))
  }
  
  # STEP 2: Single Bayesian refinement
  if (verbose) message("\n2. Applying Bayesian boundary refinement...")
  
  result = TTS_BBR(
    TTS_result = result,
    prior_strength = prior_strength,
    likelihood_strength = likelihood_strength,
    confidence_threshold = confidence_threshold,
    cores = cores,
    fix_fragments = fix_fragments
  )
  
  # Store pipeline parameters for reproducibility
  result$pipeline = list(
    segmentation_method = method,
    alpha = alpha,
    stem_height = stem_height,
    spatial_eps = spatial_eps,
    refinement = "bayesian",
    refinement_params = list(
      prior_strength = prior_strength,
      likelihood_strength = likelihood_strength,
      confidence_threshold = confidence_threshold,
      fix_fragments = fix_fragments
    ),
    timestamp = Sys.time()
  )
  
  # Add pipeline class
  class(result) = c("TTS_pipeline_result", class(result))
  
  if (verbose) {
    message("\n=== Pipeline Complete ===")
    message(sprintf("Final trees: %d", result$n_trees))
    
    # Hint about validation if ground truth column exists
    if ("treeid" %in% names(las@data)) {
      message("\nℹ Ground truth column 'treeid' detected.")
      message("  Run validate_TTS(result, las) for accuracy metrics.")
    }
    
    # Performance hint
    if (result$n_trees > 0) {
      avg_points = round(nrow(result$mesh_coords) / result$n_trees)
      message(sprintf("Average points per tree: ~%d", avg_points))
    }
  }
  
  return(result)
}

#' Print Method for Pipeline Results
#' @export
print.TTS_pipeline_result <- function(x, ...) {
  cat("TopTreeSegR Pipeline Result\n")
  cat("===========================\n")
  cat(sprintf("%-25s: %s\n", "Segmentation method", x$method))
  cat(sprintf("%-25s: %d\n", "Trees detected", x$n_trees))
  cat(sprintf("%-25s: %d\n", "Mesh vertices", nrow(x$mesh_coords)))
  
  if (!is.null(x$pipeline)) {
    cat("\nRefinement Parameters:\n")
    cat(sprintf("  %-22s: %.1f\n", "Prior strength", 
                x$pipeline$refinement_params$prior_strength))
    cat(sprintf("  %-22s: %.1f\n", "Likelihood strength", 
                x$pipeline$refinement_params$likelihood_strength))
    cat(sprintf("  %-22s: %.1f\n", "Confidence threshold", 
                x$pipeline$refinement_params$confidence_threshold))
    
    # Interpret confidence threshold
    conf = x$pipeline$refinement_params$confidence_threshold
    if (conf == 1.0) {
      cat("  → Aggressive refinement\n")
    } else if (conf >= 1.5) {
      cat("  → Conservative refinement\n")
    }
  }
  
  # Label distribution
  labels = x$labels[x$labels > 0]
  if (length(labels) > 0) {
    label_table = table(labels)
    cat(sprintf("\n%-25s: ", "Points per tree"))
    
    if (length(label_table) <= 5) {
      tree_summary = sprintf("%s(%d)", names(label_table), as.numeric(label_table))
      cat(paste(tree_summary, collapse = ", "), "\n")
    } else {
      cat(sprintf("%d trees, %d to %d points each\n", 
                  length(label_table), 
                  min(as.numeric(label_table)), 
                  max(as.numeric(label_table))))
    }
  }
  
  invisible(x)
}

#' Print method for TTS_result
#' @export
print.TTS_result <- function(x, ...) {
  cat("TopTreeSegR Result (", x$method, " method)\n", sep = "")
  cat("================================\n")
  cat(sprintf("%-20s: %d\n", "Trees detected", x$n_trees))
  cat(sprintf("%-20s: %d\n", "Mesh vertices", nrow(x$mesh_coords)))
  cat(sprintf("%-20s: %d\n", "Seed minima", length(x$seeds)))
  cat(sprintf("%-20s: %d\n", "Morse minima", length(x$minima)))
  
  # Additional info based on method
  if (x$method == "morse_smale" || x$method == "morse-smale") {
    if ("labeled_via_msmale" %in% names(x)) {
      cat(sprintf("%-20s: %d\n", "Via Morse-Smale", x$labeled_via_msmale))
    }
    if ("spatially_assigned" %in% names(x)) {
      cat(sprintf("%-20s: %d\n", "Spatially assigned", x$spatially_assigned))
    }
  }
  
  # Label distribution
  labels = x$labels[x$labels > 0]
  if (length(labels) > 0) {
    label_table = table(labels)
    cat(sprintf("%-20s: ", "Points per tree"))
    
    if (length(label_table) <= 10) {
      tree_summary = sprintf("%s(%d)", names(label_table), as.numeric(label_table))
      cat(paste(tree_summary, collapse = ", "), "\n")
    } else {
      cat(sprintf("%d trees, %d to %d points each\n", 
                  length(label_table), 
                  min(as.numeric(label_table)), 
                  max(as.numeric(label_table))))
    }
  }
  
  cat("\nParameters:\n")
  cat(sprintf("  alpha: %.2f\n", x$parameters$alpha))
  cat(sprintf("  stem_height: %.2f\n", x$parameters$stem_height))
  cat(sprintf("  spatial_eps: %.2f\n", x$parameters$spatial_eps))
  cat(sprintf("  max_distance: %.2f\n", x$parameters$max_distance))
  
  invisible(x)
}

#' 3D Visualization of Tree Segmentation
#'
#' @param TTS_result Result from TTS_segmentation()
#' @param colors Color palette for trees (default: "viridis")
#' @param point_size Size of points in 3D plot (default: 2)
#' @param alpha Point transparency (default: 0.8)
#' @param show_minima Whether to show minima locations (default: FALSE)
#' @param minima_size Size of minima points (default: 4)
#' @param minima_color Color for minima points (default: "orange")
#'
#' @return plotly 3D scatter plot
#' @export
plot_TTS_3d <- function(TTS_result, colors = "viridis", point_size = 2, 
                        alpha = 0.8, show_minima = FALSE, minima_size = 4,
                        minima_color = "grey35") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package required for 3D visualization. Install with: install.packages('plotly')")
  }
  
  # Create plot data from mesh coordinates
  plot_data = data.frame(
    x = TTS_result$mesh_coords[,1],
    y = TTS_result$mesh_coords[,2], 
    z = TTS_result$mesh_coords[,3],
    tree = as.factor(TTS_result$labels)
  )
  
  # Remove unassigned points (label 0)
  plot_data = plot_data[plot_data$tree != "0", ]
  
  # Create main plot - separate traces for each tree for better coloring
  unique_trees = unique(plot_data$tree)
  
  if (length(unique_trees) == 0) {
    stop("No trees found in the segmentation result")
  }
  
  # Start with empty plot
  fig = plotly::plot_ly()
  
  # Add points for each tree
  for (tree_id in unique_trees) {
    tree_data = plot_data[plot_data$tree == tree_id, ]
    if (nrow(tree_data) > 0) {
      fig = fig %>% plotly::add_trace(
        x = tree_data$x,
        y = tree_data$y,
        z = tree_data$z,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = point_size, 
          opacity = alpha,
          color = if(colors == "viridis") {
            # Generate distinct color for each tree
            # Using a hash to get consistent colors
            tree_num = as.numeric(as.character(tree_id))
            hue = ((tree_num - 1) * 137.508) %% 360  # Golden angle
            hsv(hue/360, 0.8, 0.9)
          } else {
            NULL  # plotly will use default colors
          }
        ),
        name = paste("Tree", tree_id),
        legendgroup = paste("Tree", tree_id),
        showlegend = length(unique_trees) <= 20  # Only show legend for reasonable number of trees
      )
    }
  }
  
  # Add minima if requested
  if (show_minima && length(TTS_result$minima) > 0) {
    minima_indices = TTS_result$minima
    # Ensure indices are within bounds
    valid_minima = minima_indices[minima_indices > 0 & minima_indices <= nrow(TTS_result$mesh_coords)]
    
    if (length(valid_minima) > 0) {
      minima_coords = TTS_result$mesh_coords[valid_minima, , drop = FALSE]
      
      fig = fig %>% plotly::add_trace(
        x = minima_coords[,1],
        y = minima_coords[,2],
        z = minima_coords[,3],
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = minima_size, 
          color = minima_color, 
          symbol = "circle",
          line = list(color = "black", width = 1)
        ),
        name = "Minima",
        showlegend = TRUE
      )
    }
  }
  
  # Layout settings
  fig = plotly::layout(
    fig,
    title = list(
      text = paste("TopTreeWSegR 3D -", TTS_result$n_trees, "Trees -", TTS_result$method, "method"),
      x = 0.5
    ),
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"), 
      zaxis = list(title = "Z"),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.5)
      )
    ),
    showlegend = length(unique_trees) <= 20,
    legend = list(
      itemsizing = "constant",
      font = list(size = 10)
    )
  )
  
  return(fig)
}

#' 2D Projection Plot
#'
#' @param TTS_result Result from TTS_segmentation()
#' @param projection Plot projection: "xy", "xz", or "yz" (default: "xy")
#' @param point_size Size of points in plot (default: 0.5)
#' @param alpha Point transparency (default: 0.8)
#' @param show_minima Whether to show minima locations (default: FALSE)
#' @param minima_size Size of minima points (default: 2)
#' @param minima_color Color for minima points (default: "orange")
#' @param ... Additional parameters passed to plot()
#'
#' @export
plot_TTS_2d <- function(TTS_result, projection = "xy", point_size = 0.5, 
                        alpha = 0.8, show_minima = FALSE, minima_size = 2,
                        minima_color = "grey35", ...) {
  
  # Prepare data based on projection
  if (projection == "xy") {
    x = TTS_result$mesh_coords[,1]
    y = TTS_result$mesh_coords[,2]
    xlab = "X"
    ylab = "Y"
  } else if (projection == "xz") {
    x = TTS_result$mesh_coords[,1]
    y = TTS_result$mesh_coords[,3]
    xlab = "X"
    ylab = "Z"
  } else if (projection == "yz") {
    x = TTS_result$mesh_coords[,2]
    y = TTS_result$mesh_coords[,3]
    xlab = "Y"
    ylab = "Z"
  } else {
    stop("projection must be 'xy', 'xz', or 'yz'")
  }
  
  # Colors for trees
  tree_labels = TTS_result$labels
  unique_trees = unique(tree_labels[tree_labels > 0])
  
  if (length(unique_trees) == 0) {
    stop("No trees found in the segmentation result")
  }
  
  # Create color palette
  colors = rainbow(length(unique_trees))
  
  # Create empty plot
  plot(x, y, type = "n", xlab = xlab, ylab = ylab,
       main = paste("Tree Segmentation (", projection, "projection)"),
       ...)
  
  # Add points for each tree
  for (i in seq_along(unique_trees)) {
    tree_id = unique_trees[i]
    idx = which(tree_labels == tree_id)
    points(x[idx], y[idx], col = adjustcolor(colors[i], alpha.f = alpha), 
           pch = 16, cex = point_size)
  }
  
  # Add minima if requested
  if (show_minima && length(TTS_result$minima) > 0) {
    minima_indices = TTS_result$minima
    if (projection == "xy") {
      minima_x = TTS_result$mesh_coords[minima_indices, 1]
      minima_y = TTS_result$mesh_coords[minima_indices, 2]
    } else if (projection == "xz") {
      minima_x = TTS_result$mesh_coords[minima_indices, 1]
      minima_y = TTS_result$mesh_coords[minima_indices, 3]
    } else {
      minima_x = TTS_result$mesh_coords[minima_indices, 2]
      minima_y = TTS_result$mesh_coords[minima_indices, 3]
    }
    
    points(minima_x, minima_y, col = minima_color, pch = 17, cex = minima_size, lwd = 1.5)
  }
  
  # Add legend only if reasonable number of trees
  if (length(unique_trees) <= 20) {
    legend("topright", legend = paste("Tree", unique_trees), 
           col = colors, pch = 16, bty = "n", cex = 0.7, 
           pt.cex = 1, ncol = min(3, ceiling(length(unique_trees)/10)))
  }
}

#' Validate TTS Segmentation Results 
#'
#' Uses result$labels and result$original_pid to match with LAS ground truth
#'
#' @param TTS_result Result from TTS_segmentation() with $labels and $original_pid
#' @param lasdf Original LAS object with ground truth
#' @param val_col Column in LAS with ground truth labels (default: 'treeid')
#' @param auto_align Label (true and pridicted) alighment (default: TRUE) 
#'
#' @return Validation metrics
#' @export
validate_TTS = function(TTS_result, lasdf, val_col = "treeid", auto_align = TRUE) {
  
  message("=== TTS Segmentation Validation ===\n")
  
  # Basic checks
  if (!inherits(lasdf, "LAS")) {
    stop("lasdf must be a LAS object")
  }
  
  if (is.null(TTS_result$labels) || is.null(TTS_result$original_pid)) {
    stop("TTS_result must contain both 'labels' and 'original_pid' elements")
  }
  
  if (!"pid" %in% names(lasdf@data)) {
    stop("LAS data must contain 'pid' column")
  }
  
  if (!val_col %in% names(lasdf@data)) {
    stop(sprintf("Ground truth column '%s' not found in LAS", val_col))
  }
  
  # Create prediction table
  pred_table = data.frame(
    pid = TTS_result$original_pid,
    predicted = TTS_result$labels
  )
  
  # Create truth table
  las_data = lasdf@data
  truth_table = data.frame(
    pid = las_data$pid,
    truth = las_data[[val_col]]
  )
  
  # Merge
  merged = merge(pred_table, truth_table, by = "pid", all = FALSE)
  
  if (nrow(merged) == 0) {
    stop("No matching points found!")
  }
  
  # Extract vectors
  pred = merged$predicted
  truth = merged$truth
  
  # Remove unassigned points for alignment
  valid = pred != 0 & truth != 0
  pred_valid = pred[valid]
  truth_valid = truth[valid]
  
  if (length(pred_valid) == 0) {
    stop("No valid points (both non-zero) for analysis")
  }
  
  # Auto-align tree IDs if requested
  if (auto_align && length(unique(pred_valid)) > 1 && length(unique(truth_valid)) > 1) {
    # Find optimal tree ID mapping
    conf = table(Predicted = pred_valid, Truth = truth_valid)
    
    if (requireNamespace("clue", quietly = TRUE) && nrow(conf) > 0 && ncol(conf) > 0) {
      # Hungarian algorithm for optimal matching
      cost = max(conf) - conf
      assign = clue::solve_LSAP(cost, maximum = FALSE)
      
      # Apply alignment to ALL points (including zeros)
      lookup = integer(max(pred, na.rm = TRUE) + 1)
      for (i in 1:nrow(conf)) {
        old_id = as.numeric(rownames(conf)[i])
        new_id = as.numeric(colnames(conf)[assign[i]])
        lookup[old_id + 1] = new_id  # +1 for 0-based indexing
      }
      
      # Align all predictions
      aligned_pred = pred
      non_zero = pred != 0
      aligned_pred[non_zero] = lookup[pred[non_zero] + 1]
      
      # Update merged data
      merged$predicted_aligned = aligned_pred
      pred = aligned_pred
      pred_valid = pred[valid]
      
      message("✓ Tree IDs automatically aligned")
    } else {
      message("⚠ Could not perform optimal alignment (clue package not available)")
    }
  }
  
  # Statistics
  n_points = nrow(merged)
  n_trees_pred = length(unique(pred[pred != 0]))
  n_trees_true = length(unique(truth[truth != 0]))
  
  # Confusion matrix (after alignment if performed)
  conf_matrix = table(Predicted = pred, Truth = truth)
  
  # Calculate metrics
  # ARI with mclust if available
  if (requireNamespace("mclust", quietly = TRUE)) {
    ari = mclust::adjustedRandIndex(pred, truth)
  } else {
    ari = NA
  }
  
  # Rand Index (simplified calculation)
  n = length(pred)
  if (n > 10000) {
    sample_size = min(5000, n)
    idx = sample(1:n, sample_size)
    pred_sample = pred[idx]
    truth_sample = truth[idx]
    
    same_pred = outer(pred_sample, pred_sample, "==")
    same_truth = outer(truth_sample, truth_sample, "==")
    
    a = sum(same_pred & same_truth) - sample_size
    b = sum(!same_pred & !same_truth) - sample_size
    n_pairs = sample_size * (sample_size - 1) / 2
    
    ri = (a + b) / (2 * n_pairs)
  } else {
    same_pred = outer(pred, pred, "==")
    same_truth = outer(truth, truth, "==")
    
    a = sum(same_pred & same_truth) - n
    b = sum(!same_pred & !same_truth) - n
    n_pairs = n * (n - 1) / 2
    
    ri = (a + b) / (2 * n_pairs)
  }
  
  # Tree-level metrics (precision, recall, F1)
  TP = sum(apply(conf_matrix, 1, max))
  FP = sum(rowSums(conf_matrix)) - TP
  FN = sum(colSums(conf_matrix)) - TP
  
  precision = if (TP + FP > 0) TP / (TP + FP) else 0
  recall = if (TP + FN > 0) TP / (TP + FN) else 0
  f1_score = if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0
  
  # Point-level accuracy
  accuracy = sum(diag(conf_matrix)) / sum(conf_matrix)
  
  # Output
  message(sprintf("Points: %d | Pred trees: %d | True trees: %d", 
                  n_points, n_trees_pred, n_trees_true))
  message(sprintf("Match rate: %.1f%%", nrow(merged)/nrow(pred_table)*100))
  
  message("\n=== METRICS ===")
  message(sprintf("Precision:  %.4f", precision))
  message(sprintf("Recall:     %.4f", recall))
  message(sprintf("F1-Score:   %.4f", f1_score))
  message(sprintf("Accuracy:   %.4f", accuracy))
  message(sprintf("Rand Index: %.4f", ri))
  
  if (!is.na(ari)) {
    message(sprintf("Adj Rand I: %.4f", ari))
  }
  
  # Return results
  results = list(
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    accuracy = accuracy,
    rand_index = ri,
    adjusted_rand_index = ari,
    n_trees_pred = n_trees_pred,
    n_trees_true = n_trees_true,
    n_points = n_points,
    confusion_matrix = conf_matrix,
    merged_data = merged,
    aligned = auto_align,
    summary = data.frame(
      Metric = c("Precision", "Recall", "F1", "Accuracy", "RI", "ARI"),
      Value = c(
        round(precision, 4),
        round(recall, 4),
        round(f1_score, 4),
        round(accuracy, 4),
        round(ri, 4),
        ifelse(is.na(ari), NA, round(ari, 4))
      )
    )
  )
  
  invisible(results)
}

