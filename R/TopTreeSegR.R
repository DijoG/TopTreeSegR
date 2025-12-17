#' Topological Tree Segmentation using Discrete Morse Theory
#'
#' @name TopTreeSegR
#' @description Ultra-fast topology-based tree segmentation from LiDAR point clouds 
#' using Discrete Morse Theory and RcppArmadillo for maximum performance.
#' @docType package
#' @useDynLib TopTreeSegR, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
NULL

.onAttach = function(libname, pkgname) {
  packageStartupMessage(
    "------ TopTreeSegR loaded! Ultra-fast tree segmentation with Discrete Morse Theory -----\n",
    "------ GitHub dependencies: DiscreteMorseR & AlphaHull3D -----\n",
    "------ Use TTS_segmentation() for blazing-fast segmentation! ------"
  )
}

.check_dependencies = function() {
  required_pkgs = c("DiscreteMorseR", "AlphaHull3D", "dbscan", "ggplot2", "FNN")
  missing_pkgs = required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "The following required packages are missing: ", 
      paste(missing_pkgs, collapse = ", "), "\n",
      "Install with: install.packages(c('", 
      paste(setdiff(missing_pkgs, c("DiscreteMorseR", "AlphaHull3D")), collapse = "', '"), 
      "'))\n",
      "And GitHub packages with: remotes::install_github(c('",
      "DijoG/DiscreteMorseR', 'stla/AlphaHull3D'))"  
    )
  }
  return(TRUE)
}

#' Ultra-Fast Tree Segmentation from LiDAR Point Clouds
#'
#' @param lasdf Input point cloud in las/laz or matrix/data.frame format
#' @param alpha Alpha parameter for alpha-complex construction (default: 0.1)
#' @param input_truth Attribute from lasd (default: 'pid')
#' @param clip_height Height threshold for detecting stem seeds (default: 0.5)
#' @param max_distance Maximum distance for minima connectivity (default: 2.0)
#' @param grid_size Grid size for spatial hashing in large datasets (default: 5.0)
#' @param cores Number of cores to use for Morse complex computation (default: 12)
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{segmentation}: Vector of tree labels for each point
#'   \item \code{points}: Original point coordinates
#'   \item \code{n_trees}: Number of trees detected
#'   \item \code{minima}: Indices of Morse minima
#'   \item \code{seeds}: Indices of seed minima used for tree initialization
#'   \item \code{gradient_flow}: Complete gradient flow information
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Load your point cloud data
#' # lasdf = read.csv("point_cloud.csv")
#' 
#' # Run ultra-fast segmentation
#' result = TTS_segmentation(lasdf)
#' 
#' # Plot results
#' plot(result)
#' }
TTS_segmentation = function(lasdf, alpha = 0.1, 
                            input_truth = "pid", 
                            clip_height = 0.5, 
                            max_distance = 2.0, 
                            grid_size = 5.0, 
                            cores = 12) {
  
  # Check dependencies
  .check_dependencies()
  
  message("=== TopTreeSegR: FAST TTS with DiscreteMorseR and RcppArmadillo ===\n")
  
  # Handle different input types (LAS object, data frame, or matrix)
  if (inherits(lasdf, "LAS")) {
    points = as.matrix(unique(as.data.frame(lasdf@data[,1:3])))
  } else if (is.data.frame(lasdf) || is.matrix(lasdf)) {
    points = as.matrix(unique(as.data.frame(lasdf[,1:3])))
  } else {
    stop("lasdf must be a LAS object, data frame, or matrix")
  }
  
  points = as.matrix(unique(as.data.frame(lasdf@data[,1:3])))
  
  message("Input points: ", nrow(lasdf), "\n")
  
  # Step 1: Build alpha-complex mesh
  message("1. Building alpha-complex...\n")
  a = ahull3D::ahull3D(points, alpha = alpha, input_truth = input_truth)
  mesh = DiscreteMorseR::get_CCMESH(a)
  
  # Step 2: Compute Morse complex
  message("2. Computing Morse complex...\n")
  morse_complex = DiscreteMorseR::compute_MORSE_complex(
    mesh, 
    output_dir = tempdir(), 
    cores = cores
  )
  
  # Step 3: ULTRA-FAST gradient flow with Armadillo
  message("3. Building ULTRA-FAST gradient flow...\n")
  gradient_flow = build_gradient_flow(morse_complex, mesh$vertices, 
                                      max_distance, grid_size)
  
  # Step 4: Detect tree seeds
  message("4. Detecting tree seeds...\n")
  seeds = detect_tree_seeds_proper(gradient_flow$minima, mesh$vertices, clip_height)
  
  # Step 5: Propagate labels
  message("5. Propagating labels...\n")
  labeled_data = propagate_labels_simple_fixed(gradient_flow, seeds, mesh$vertices)
  
  # Step 6: Create segmentation
  message("6. Creating final segmentation...\n")
  segmentation = create_segmentation_proper(mesh$vertices, gradient_flow, labeled_data)
  
  message("ULTRA-FAST implementation complete. Found ", length(unique(segmentation)), "trees\n")
  
  structure(list(
    segmentation = segmentation,
    input_truth = mesh$vertices$i123,
    points = mesh$vertices,
    n_trees = length(unique(segmentation)),
    minima = gradient_flow$minima,
    seeds = seeds$seed_minima,
    gradient_flow = gradient_flow
  ), class = "TTS_result")
}

#' Build Ultra-Fast Gradient Flow
#' 
#' @keywords internal
build_gradient_flow = function(morse_complex, vertices, max_distance = 2.0, grid_size = 5.0) {
  message("  Building ULTRA-FAST gradient flow with Armadillo...\n")
  
  # Extract minima using C++ implementation
  minima = extract_minima_corrected_cpp(morse_complex, vertices, tempdir())
  message("  Found ", length(minima), "minima\n")
  
  # Convert to Armadillo types for maximum speed
  minima_arma = as.numeric(minima) - 1  # Convert to 0-based for C++
  vertices_arma = as.matrix(vertices)
  
  # ULTRA-FAST: Parse gradient pairs
  message("  [Armadillo] Parsing gradient pairs...\n")
  parse_time = system.time({
    gradient_network = parse_gradient_network_fast_cpp(morse_complex$vector_field, 
                                                       vertices[,3],     # Add elevations
                                                       nrow(vertices))   # Add n_vertices
  })
  
  message(sprintf("  [Armadillo] Found %d vertex->edge pairs, %d edge->face pairs (%.2fs)\n",
                  gradient_network$vertex_pair_count, gradient_network$edge_count, parse_time[3]))
  
  # ULTRA-FAST: Compute ascending regions
  message("  [Armadillo] Computing ascending regions...\n")
  ascend_time = system.time({
    ascending_regions = compute_ascending_regions_fast_cpp(gradient_network, minima_arma, nrow(vertices))
  })
  
  # ULTRA-FAST: Build minima connectivity with spatial hashing
  message("  [Armadillo] Building minima connectivity...\n")
  connect_time = system.time({
    if (length(minima) > 1000) {
      # Use spatial hashing for large datasets
      minima_connectivity = build_minima_connectivity_spatial_cpp(minima_arma, vertices_arma, 
                                                              max_distance, grid_size)
    } else {
      # Use direct method for smaller datasets
      minima_connectivity = build_minima_connectivity_fast_cpp(minima_arma, vertices_arma, max_distance)
    }
  })
  
  avg_degree = mean(sapply(minima_connectivity, length))
  message(sprintf("  [Armadillo] Connectivity: %d minima, avg degree: %.2f (%.2fs)\n", 
                  length(minima), avg_degree, connect_time[3]))
  
  list(
    minima = minima,
    ascending_regions = as.integer(ascending_regions),
    gradient_network = gradient_network,
    minima_connectivity = minima_connectivity
  )
}

#' Detect Tree Seeds from Morse Minima
#' 
#' @keywords internal
detect_tree_seeds_proper = function(minima, vertices, clip_height) {
  message("  Detecting tree seeds from stem minima...\n")
  
  # Find minima in stem region (low elevation)
  low_minima = minima[vertices[minima, 3] < clip_height]
  
  if (length(low_minima) == 0) {
    message("  No low minima found, using all minima\n")
    low_minima = minima
  }
  
  # Cluster stem minima spatially
  stem_coords = vertices[low_minima, 1:2, drop = FALSE]
  dbscan_result = dbscan::dbscan(stem_coords, eps = 1, minPts = 1)
  
  clusters = unique(dbscan_result$cluster)
  clusters = clusters[clusters > 0]
  
  seed_minima = integer(0)
  seed_labels = integer(0)
  
  for (cluster_id in clusters) {
    cluster_mask = dbscan_result$cluster == cluster_id
    cluster_minima = low_minima[cluster_mask]
    
    if (length(cluster_minima) > 0) {
      # Take the lowest minimum in each cluster
      cluster_elevations = vertices[cluster_minima, 3]
      lowest_min_idx = which.min(cluster_elevations)
      seed_minima = c(seed_minima, cluster_minima[lowest_min_idx])
      seed_labels = c(seed_labels, cluster_id)
    }
  }
  
  # Ensure unique labels
  seed_labels = as.numeric(as.factor(seed_labels))
  
  message("  Detected ", length(seed_minima), "seed minima\n")
  list(seed_minima = seed_minima, seed_labels = seed_labels)
}

#' Propagate Labels via Gradient Flow
#' 
#' @keywords internal
propagate_labels_simple_fixed = function(gradient_flow, seeds, vertices) {
  minima = gradient_flow$minima
  n_minima = length(minima)
  labels = integer(n_minima)
  
  # Mark seed minima
  seed_indices = match(seeds$seed_minima, minima)
  valid_seeds = which(!is.na(seed_indices))
  
  if (length(valid_seeds) == 0) {
    message("  No seed matches, using simple spatial assignment\n")
    # Simple spatial assignment
    seed_coords = vertices[seeds$seed_minima, 1:2, drop = FALSE]
    minima_coords = vertices[minima, 1:2, drop = FALSE]
    
    # Assign each minimum to nearest seed
    nn_idx = FNN::get.knnx(seed_coords, minima_coords, k = 1)$nn.index
    labels = seeds$seed_labels[nn_idx]
  } else {
    message("  Found ", length(valid_seeds), "seed matches\n")
    seed_indices = seed_indices[valid_seeds]
    labels[seed_indices] = seeds$seed_labels[valid_seeds]
    
    # Simple spatial propagation for remaining
    unlabeled = which(labels == 0)
    if (length(unlabeled) > 0) {
      labeled_coords = vertices[minima[labels > 0], 1:2, drop = FALSE]
      unlabeled_coords = vertices[minima[unlabeled], 1:2, drop = FALSE]
      nn_idx = FNN::get.knnx(labeled_coords, unlabeled_coords, k = 1)$nn.index
      labels[unlabeled] = labels[labels > 0][nn_idx]
    }
  }
  
  message("  Simple label distribution:\n")
  print(table(labels))
  list(minima = minima, labels = labels)
}

#' Create Final Segmentation
#' 
#' @keywords internal
create_segmentation_proper = function(vertices, gradient_flow, labeled_data) {
  n_vertices = nrow(vertices)
  segmentation = integer(n_vertices)
  
  # Assign based on ascending regions
  for (i in 1:n_vertices) {
    region_id = gradient_flow$ascending_regions[i]
    if (region_id > 0 && region_id <= length(labeled_data$minima)) {
      segmentation[i] = labeled_data$labels[region_id]
    }
  }
  
  # Handle any unassigned vertices
  unassigned = which(segmentation == 0)
  if (length(unassigned) > 0) {
    assigned = which(segmentation > 0)
    if (length(assigned) > 0) {
      nn_idx = FNN::get.knnx(
        vertices[assigned, , drop = FALSE],
        vertices[unassigned, , drop = FALSE],
        k = 1
      )$nn.index
      segmentation[unassigned] = segmentation[assigned[nn_idx]]
    }
  }
  
  segmentation
}

#' 2D Visualization of Tree Segmentation 
#'
#' @param x TTS_result object from TTS_segmentation()
#' @param point_size Size of points in plot (default: 0.5)
#' @param alpha Point transparency (default: 0.8)
#' @param show_seeds Whether to show seed locations (default: TRUE)
#' @param projection Coordinate projection to plot: "XY", "XZ", or "YZ" (default: "XY")
#' @param ... Additional plotting parameters
#'
#' @return ggplot object
#' @export
#' @method plot TTS_result
plot_TTS_2d = function(x, point_size = 0.5, alpha = 0.8, 
                        show_seeds = TRUE, projection = "XY", ...) {
  
  # Validate projection argument
  if (!projection %in% c("XY", "XZ", "YZ")) {
    stop('projection must be one of: "XY", "XZ", or "YZ"')
  }
  
  # Create plot data based on projection
  if (projection == "XY") {
    plot_data = data.frame(
      x = x$points[,1],
      y = x$points[,2],
      tree = as.factor(x$segmentation)
    )
    x_lab = "X"
    y_lab = "Y"
  } else if (projection == "XZ") {
    plot_data = data.frame(
      x = x$points[,1],
      y = x$points[,3],
      tree = as.factor(x$segmentation)
    )
    x_lab = "X"
    y_lab = "Z"
  } else if (projection == "YZ") {
    plot_data = data.frame(
      x = x$points[,2],
      y = x$points[,3],
      tree = as.factor(x$segmentation)
    )
    x_lab = "Y"
    y_lab = "Z"
  }
  
  p = ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = tree)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(
      title = paste("TopTreeSegR -", x$n_trees, "Trees (", projection, "projection)"),
      x = x_lab, 
      y = y_lab
    ) +
    ggplot2::theme_minimal() +
    ggplot2::coord_equal() +
    ggplot2::theme(legend.position = "none")
  
  # Add seed points if available
  if (show_seeds && !is.null(x$seeds)) {
    if (projection == "XY") {
      seed_data = data.frame(
        x = x$points[x$seeds, 1],
        y = x$points[x$seeds, 2]
      )
    } else if (projection == "XZ") {
      seed_data = data.frame(
        x = x$points[x$seeds, 1],
        y = x$points[x$seeds, 3]
      )
    } else if (projection == "YZ") {
      seed_data = data.frame(
        x = x$points[x$seeds, 2],
        y = x$points[x$seeds, 3]
      )
    }
    
    p = p + ggplot2::geom_point(
      data = seed_data, 
      ggplot2::aes(x = x, y = y), 
      color = "black", size = 4, shape = "x"
    )
  }
  
  return(p)
}

#' 3D Visualization of Tree Segmentation
#'
#' @param TTS_result Result from TTS_segmentation()
#' @param colors Color of labeled points (default: "viridis")
#' @param point_size Size of points in 3D plot (default: 2)
#' @param alpha Point transparency (default: 0.8)
#'
#' @return plotly 3D scatter plot
#' @export
plot_TTS_3d = function(TTS_result, colors = "viridis", point_size = 2, alpha = 0.8) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package required for 3D visualization. Install with: install.packages('plotly')")
  }
  
  plot_data = data.frame(
    x = TTS_result$points[,1],
    y = TTS_result$points[,2], 
    z = TTS_result$points[,3],
    tree = as.factor(TTS_result$segmentation)
  )
  
  fig = plotly::plot_ly(plot_data, 
                        x = ~x, y = ~y, z = ~z, 
                        color = ~tree,
                        colors = colors,
                        type = "scatter3d", 
                        mode = "markers",
                        marker = list(size = point_size, opacity = alpha))
  
  fig = plotly::layout(
    fig,
    title = paste("TopTreeSegR 3D -", TTS_result$n_trees, "Trees"),
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"), 
      zaxis = list(title = "Z")
    )
  )
  
  return(fig)
}

#' Benchmark Segmentation Performance
#'
#' @param lasdf Point cloud data
#' @param ... Parameters passed to TTS_segmentation()
#'
#' @return Benchmark results with timing information
#' @export
benchmark_TTS = function(lasdf, ...) {
  message("=== TopTreeSegR BENCHMARK ===\n")
  
  total_time = system.time({
    result = TTS_segmentation(lasdf, ...)
  })
  
  message("\n=== BENCHMARK RESULTS ===\n")
  message(sprintf("Total time: %.2f seconds\n", total_time[3]))
  message(sprintf("Points processed: %d\n", nrow(lasdf)))
  message(sprintf("Trees detected: %d\n", result$n_trees))
  message(sprintf("Minima found: %d\n", length(result$minima)))
  message(sprintf("Seeds detected: %d\n", length(result$seeds)))
  message(sprintf("Processing rate: %.2f points/second\n", nrow(lasdf)/total_time[3]))
  
  invisible(result)
}

#' ' Validate Segmentation Quality
#'
#' @param TTS_result Result from TTS_segmentation()
#' @param lasdf Input point cloud las/laz format
#' @param input_truth Attribute from lasdf (default: 'pid')
#' @param val_col Column/field in the LAS data storing true labels
#'
#' @return Validation metrics including precision, recall, F1, accuracy, RI, ARI
#' @export
validate_TTS = function(TTS_result, lasdf, input_truth = "pid", val_col = "treeid") {
  message("=== TopTreeSegR VALIDATION ===")
  
  # 1. Create prediction table
  pred_table = data.frame(
    mergeon = TTS_result$input_truth,  # This should contain pids or another identifier
    predicted = TTS_result$segmentation
  )
  
  # 2. Check if lasdf is a LAS object
  if (!inherits(lasdf, "LAS")) {
    stop("Argument 'lasdf' must be a LAS object.")
  }
  
  las_data = lasdf@data
  
  # 3. Check if required columns exist
  if (!input_truth %in% names(las_data)) {
    stop(sprintf("Column '%s' (input_truth) not found in LAS data.", input_truth))
  }
  
  if (!val_col %in% names(las_data)) {
    stop(sprintf("Column '%s' (val_col) not found in LAS data.", val_col))
  }
  
  # 4. Create truth table
  truth_table = data.frame(
    mergeon = las_data[[input_truth]],
    truth = las_data[[val_col]]
  )
  
  # 5. Merge prediction and truth tables
  # Note: Use the actual column name "mergeon" in the merge
  merged = merge(pred_table, truth_table, by = "mergeon", all = FALSE)
  
  if (nrow(merged) == 0) {
    stop("No matching points found. Check that values in TTS_result$input_truth match those in the '", 
         input_truth, "' column of lasdf.")
  }
  
  # 6. Extract aligned vectors
  segmentation = merged$predicted
  truth = merged$truth
  
  # Basic statistics
  n_trees_pred = length(unique(segmentation[segmentation != 0]))
  n_trees_true = length(unique(truth[truth != 0]))
  n_matched = nrow(merged)
  
  message(sprintf("Matched points: %d", n_matched))
  message(sprintf("Predicted trees: %d | True trees: %d", n_trees_pred, n_trees_true))
  
  # Check for unassigned points
  unassigned_pred = sum(segmentation == 0)
  unassigned_true = sum(truth == 0)
  
  if (unassigned_pred > 0) {
    message(sprintf("Unassigned (pred): %d (%.1f%%)", 
                    unassigned_pred, unassigned_pred/n_matched*100))
  }
  
  # Get unique cluster IDs
  pred_clusters = unique(segmentation)
  true_clusters = unique(truth)
  
  pred_clusters_nozero = pred_clusters[pred_clusters != 0]
  true_clusters_nozero = true_clusters[true_clusters != 0]
  
  # Initialize confusion matrix
  confusion_matrix = matrix(0, nrow = length(pred_clusters_nozero), 
                            ncol = length(true_clusters_nozero))
  rownames(confusion_matrix) = pred_clusters_nozero
  colnames(confusion_matrix) = true_clusters_nozero
  
  # Fill confusion matrix
  for (i in seq_along(pred_clusters_nozero)) {
    for (j in seq_along(true_clusters_nozero)) {
      pred_mask = segmentation == pred_clusters_nozero[i]
      true_mask = truth == true_clusters_nozero[j]
      confusion_matrix[i, j] = sum(pred_mask & true_mask)
    }
  }
  
  # Calculate metrics from confusion matrix
  TP = sum(apply(confusion_matrix, 1, max))
  FP = sum(rowSums(confusion_matrix)) - TP
  FN = sum(colSums(confusion_matrix)) - TP
  
  # Precision, Recall, F1
  precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall = ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1_score = ifelse((precision + recall) > 0, 
                    2 * precision * recall / (precision + recall), 0)
  
  # Overall accuracy
  correct_assignments = sum(segmentation == truth)
  accuracy = correct_assignments / n_matched
  
  # Rand Index (RI)
  n = n_matched
  tp_pairs = 0
  tn_pairs = 0
  fp_pairs = 0
  fn_pairs = 0
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      same_pred = segmentation[i] == segmentation[j]
      same_true = truth[i] == truth[j]
      
      if (same_pred && same_true) {
        tp_pairs = tp_pairs + 1
      } else if (!same_pred && !same_true) {
        tn_pairs = tn_pairs + 1
      } else if (same_pred && !same_true) {
        fp_pairs = fp_pairs + 1
      } else {
        fn_pairs = fn_pairs + 1
      }
    }
  }
  
  total_pairs = n * (n-1) / 2
  ri = (tp_pairs + tn_pairs) / total_pairs
  
  # Adjusted Rand Index (ARI)
  valid_points = truth != 0 & segmentation != 0
  if (sum(valid_points) > 0) {
    pred_valid = segmentation[valid_points]
    true_valid = truth[valid_points]
    
    contingency = table(pred_valid, true_valid)
    
    sum_row = rowSums(contingency)
    sum_col = colSums(contingency)
    n_valid = sum(contingency)
    
    index = sum(choose(contingency, 2))
    expected_index = sum(choose(sum_row, 2)) * sum(choose(sum_col, 2)) / choose(n_valid, 2)
    max_index = 0.5 * (sum(choose(sum_row, 2)) + sum(choose(sum_col, 2)))
    
    if (max_index - expected_index != 0) {
      ari = (index - expected_index) / (max_index - expected_index)
    } else {
      ari = 0
    }
  } else {
    ari = NA
  }
  
  # Output metrics
  message("\n=== METRICS ===")
  message(sprintf("Precision:  %.4f", precision))
  message(sprintf("Recall:     %.4f", recall))
  message(sprintf("F1-Score:   %.4f", f1_score))
  message(sprintf("Accuracy:   %.4f", accuracy))
  message(sprintf("Rand Index: %.4f", ri))
  
  if (!is.na(ari)) {
    message(sprintf("Adj Rand I: %.4f", ari))
  } else {
    message("Adj Rand I: NA (no valid points)")
  }
  
  # Return all metrics
  invisible(list(
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    accuracy = accuracy,
    rand_index = ri,
    adjusted_rand_index = ari,
    n_trees_pred = n_trees_pred,
    n_trees_true = n_trees_true,
    n_matched_points = n_matched,
    unassigned_pred = unassigned_pred,
    unassigned_true = unassigned_true,
    confusion_matrix = confusion_matrix,
    metrics_summary = data.frame(
      Metric = c("Trees Predicted", "Trees True", "Precision", "Recall", 
                 "F1-Score", "Accuracy", "RI", "ARI"),
      Value = c(n_trees_pred, n_trees_true, precision, recall, 
                f1_score, accuracy, ri, ifelse(is.na(ari), NA, ari))
    )
  ))
}

#' Print Method for TTS_result
#'
#' @param x TTS_result object
#' @param ... Additional parameters
#'
#' @export
#' @method print TTS_result
print_TTS_result <- function(x, ...) {
  cat("TopTreeSegR Segmentation Result\n")
  cat("================================\n")
  cat(sprintf("%-20s: %d\n", "Trees detected", length(unique(x$segmentation[x$segmentation != 0]))))
  
  # Add the true tree count line
  if (!is.null(x$input_truth)) {
    cat(sprintf("%-20s: %d\n", "Trees true", length(unique(x$input_truth[x$input_truth != 0]))))
  } else {
    cat(sprintf("%-20s: %s\n", "Trees true", "No ground truth provided"))
  }
  
  cat(sprintf("%-20s: %d\n", "Total points", nrow(x$points)))
  cat(sprintf("%-20s: %d\n", "Morse minima", length(x$minima)))
  cat(sprintf("%-20s: %d\n", "Seed minima", length(x$seeds)))
  
  # Show a sample of segmentation IDs
  seg_ids = unique(x$segmentation)
  if (length(seg_ids) <= 10) {
    cat(sprintf("%-20s: %s\n", "Segmentation IDs", 
                paste(seg_ids, collapse = ", ")))
  } else {
    cat(sprintf("%-20s: %s ... (%d total)\n", "Segmentation IDs", 
                paste(seg_ids[1:10], collapse = ", "), length(seg_ids)))
  }
  invisible(x)
}

