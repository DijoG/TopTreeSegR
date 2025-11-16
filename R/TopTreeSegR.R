#' Discrete Morse Theory Analysis
  #' 
  #' @description Compute Morse vector field and critical simplices from 3D meshes.
  #' @docType package
  #' @useDynLib DiscreteMorseR, .registration = TRUE
  #' @importFrom Rcpp evalCpp
  #' @import dplyr purrr stringr data.table gtools readr future furrr
  NULL

#' Add decimal formatting
#' 
#' @param x Numeric vector
#' @param k Number of decimal places
#' @return Formatted character vector
#' @keywords internal
add_DECIMAL = function(x, k) format(round(x, k), nsmall = k)

#' Compute lexicographically sorted simplex labels
#'
#' @param df Data frame with label and idlabel columns
#' @return Data table with lexicographically sorted labels
#' @keywords internal
get_lexIDLAB = function(df) {
  # Trim labels
  LA = stringr::str_trim(df$label, "both")
  IDLA = stringr::str_trim(df$idlabel, "both")
  
  # Split and get sorting indices
  LA_split = stringr::str_split(LA, " ")
  m1 = purrr::map(LA_split, ~get_MIXEDSORT_cpp(add_DECIMAL(as.numeric(.), k = 9)))
  
  # Apply sorting to both label types
  mlab = purrr::map2(LA_split, m1, ~paste(.x[.y], collapse = " "))
  mid = purrr::map2(stringr::str_split(IDLA, " "), m1, ~paste(.x[.y], collapse = " "))
  
  # Create result data.table
  result = data.table::data.table(
    lexi_label = unlist(mlab),
    lexi_id = unlist(mid)
  )
  
  return(result)
}

#' Ultra-fast mesh preparation with connected components
#' @param vertices Vertex data from alphahull 
#' @param faces Face data from alphahull
#' @param select_largest If TRUE, returns only largest connected component (default: TRUE).
#'                      If FALSE, returns list of all connected components.
#' @return Single mesh (if select_largest=TRUE) or list of meshes (if select_largest=FALSE)
#' @export
get_CCMESH = function(vertices, faces, select_largest = TRUE) {
  
  vertices = as.matrix(vertices)
  faces = as.matrix(faces)
  
  # Handle alpha hull format
  if(nrow(vertices) == 4) vertices = t(vertices)
  if(nrow(faces) == 3) faces = t(faces)
  if(ncol(vertices) == 4) vertices = vertices[, 1:3]
  if(min(faces) == 0) faces = faces + 1
  
  # Input validation
  if (nrow(vertices) == 0) stop("Vertex matrix is empty")
  if (nrow(faces) == 0) stop("Face matrix is empty")
  
  message("Input: ", nrow(vertices), " vertices, ", nrow(faces), " faces")
  
  # Call C++ function
  result = get_CCMESH_cpp(vertices, faces, select_largest)
  
  # Handle different return types
  if ("components" %in% names(result)) {
    # Multiple components returned (select_largest = FALSE)
    components = result$components
    n_components = result$n_components
    
    # Convert each component to proper mesh format
    meshes = lapply(components, function(comp) {
      mesh = list(
        vertices = as.matrix(comp$vertices),
        faces = as.matrix(comp$faces),
        edges = as.data.frame(comp$edges)
      )
      
      # Set column names
      colnames(mesh$vertices) = c("X", "Y", "Z")
      colnames(mesh$faces) = c("i1", "i2", "i3")
      colnames(mesh$edges) = c("i1", "i2")
      
      return(mesh)
    })
    
    # Add component info as attributes
    attr(meshes, "n_components") = n_components
    attr(meshes, "component_sizes") = sapply(meshes, function(m) nrow(m$faces))
    
    message("Returning ", n_components, " connected components")
    
    return(meshes)
    
  } else {
    # Single mesh returned (select_largest = TRUE)
    mesh = list(
      vertices = as.matrix(result$vertices),
      faces = as.matrix(result$faces),
      edges = as.data.frame(result$edges)
    )
    
    # Set column names
    colnames(mesh$vertices) = c("X", "Y", "Z")
    colnames(mesh$faces) = c("i1", "i2", "i3")
    colnames(mesh$edges) = c("i1", "i2")
    
    message("Largest mesh: ", nrow(mesh$vertices), " vertices, ",
            nrow(mesh$faces), " faces, ", nrow(mesh$edges), " edges")
    
    return(mesh)
  }
}

#' Extract and process simplices from mesh
#'
#' @param mesh Mesh object from get_CCMESH()
#' @param txt_dirout Directory for output files (optional)
#' @return List of processed simplices
#' @keywords internal
get_SIMPLICES = function(mesh, txt_dirout = "") {
  
  # Vertices (0-simplex) - Safe for parallel
  mesh_ver = as.data.frame(mesh$vertices)
  colnames(mesh_ver) = c("X", "Y", "Z")
  mesh_ver$i123 = 1:nrow(mesh_ver)
  mesh_ver = mesh_ver[!duplicated(mesh_ver), ]
  mesh_ver$Z = as.character(mesh_ver$Z)
  
  if (txt_dirout != "") {
    readr::write_tsv(mesh_ver, stringr::str_c(txt_dirout, "/vertices.txt"))
  }
  
  # Edges (1-simplex) - Safe for parallel
  edges_df = as.data.frame(mesh$edges)
  colnames(edges_df) = c("i1", "i2")
  
  # First join - i1 to vertices
  mesh_edge = dplyr::left_join(edges_df, mesh_ver, by = c("i1" = "i123")) %>%
    dplyr::rename(ii1X = X, ii1Y = Y, ii1Z = Z)
  
  # Second join - i2 to vertices  
  mesh_edge = dplyr::left_join(mesh_edge, mesh_ver, by = c("i2" = "i123")) %>%
    dplyr::rename(ii2X = X, ii2Y = Y, ii2Z = Z) %>%
    dplyr::mutate(
      label = stringr::str_c(as.character(ii1Z), as.character(ii2Z), sep = " "),
      idlabel = stringr::str_c(as.character(i1), as.character(i2), sep = " ")
    )
  
  out = get_lexIDLAB(mesh_edge)
  
  # Use cbind! (Restored from original)
  mesh_edge_final = cbind(mesh_edge, out)
  
  if (txt_dirout != "") {
    readr::write_tsv(mesh_edge_final, stringr::str_c(txt_dirout, "/edges.txt"))
  }
  
  # Faces (2-simplex) - Already has correct column names!
  faces_df = as.data.frame(mesh$faces)  # Keeps i1, i2, i3 from matrix
  
  # Three explicit joins for faces
  mesh_f = dplyr::left_join(faces_df, mesh_ver, by = c("i1" = "i123")) %>%
    dplyr::rename(ii1X = X, ii1Y = Y, ii1Z = Z)
  
  mesh_f = dplyr::left_join(mesh_f, mesh_ver, by = c("i2" = "i123")) %>%
    dplyr::rename(ii2X = X, ii2Y = Y, ii2Z = Z)
  
  mesh_face = dplyr::left_join(mesh_f, mesh_ver, by = c("i3" = "i123")) %>%
    dplyr::rename(ii3X = X, ii3Y = Y, ii3Z = Z) %>%
    dplyr::mutate(
      label = stringr::str_c(as.character(ii1Z), as.character(ii2Z), as.character(ii3Z), sep = " "),
      idlabel = stringr::str_c(as.character(i1), as.character(i2), as.character(i3), sep = " ")
    )
  
  out = get_lexIDLAB(mesh_face)
  
  # Use cbind! (Restored from original)
  mesh_face_final = cbind(mesh_face, out)
  
  if (txt_dirout != "") {
    readr::write_tsv(mesh_face_final, stringr::str_c(txt_dirout, "/faces.txt"))
  }
  
  return(list(
    vertices = mesh_ver,
    edges = mesh_edge_final,
    faces = mesh_face_final
  ))
}

#' Compute lower star filtration for vertices
#'
#' @param vertex Vertex data
#' @param edge Edge data  
#' @param face Face data
#' @param dirout Output directory (optional)
#' @param cores Number of cores (for consistency)
#' @return List of lower star sets
#' @keywords internal
get_lowerSTAR = function(vertex, edge, face, dirout = NULL, cores = 1) {
  
  # Pre-compute connections
  all_connections = get_vertTO_cpp(vertex, edge, face)
  lexi_labels = all_connections$lexi_label
  lexi_ids = all_connections$lexi_id
  
  # Pre-compute everything upfront
  first_verts = sapply(strsplit(lexi_labels, " "), function(x) as.numeric(x[1]))
  simplex_vertices = strsplit(lexi_ids, " ")
  
  # Build reverse index: vertex_id -> simplex_indices
  vertex_to_simplices = list()
  for (i in seq_along(simplex_vertices)) {
    for (v_id in simplex_vertices[[i]]) {
      if (is.null(vertex_to_simplices[[v_id]])) {
        vertex_to_simplices[[v_id]] = integer(0)
      }
      vertex_to_simplices[[v_id]] = c(vertex_to_simplices[[v_id]], i)
    }
  }
  
  lowerSTAR = vector("list", nrow(vertex))
  
  for (i in 1:nrow(vertex)) {
    v_id = as.character(vertex$i123[i])
    v_z = as.numeric(vertex$Z[i])
    
    simplex_indices = vertex_to_simplices[[v_id]]
    
    if (!is.null(simplex_indices) && length(simplex_indices) > 0) {
      valid_indices = simplex_indices[first_verts[simplex_indices] <= v_z]
      
      if (length(valid_indices) > 0) {
        valid_simplices = all_connections[valid_indices, ]
        
        order_idx = gtools::mixedorder(valid_simplices$lexi_label, decreasing = FALSE)
        sorted_simplices = valid_simplices[order_idx, ]
        
        # Create data.frame WITHOUT factors
        result = data.frame(
          id = rep(as.integer(v_id), nrow(sorted_simplices)),
          lexi_label = sorted_simplices$lexi_label,
          lexi_id = sorted_simplices$lexi_id,
          stringsAsFactors = FALSE
        )
        
        lowerSTAR[[i]] = result
        
        if (!is.null(dirout)) {
          output_lines = apply(result, 1, function(row) {
            paste(row, collapse = "\t")
          })
          
          # Append to file with cat() for exact format matching
          cat(output_lines, file = file.path(dirout, "lowerSTAR.txt"), 
              sep = "\n", append = TRUE)
        }
      }
    }
  }
  
  lowerSTAR = lowerSTAR[!sapply(lowerSTAR, is.null)]
  return(lowerSTAR)
}

#' Process lower star to get Morse pairings
#'
#' @param list_lowerSTAR Lower star data
#' @param vertex Vertex data
#' @return List with vector field and critical simplices
#' @keywords internal
proc_lowerSTAR <- function(list_lowerSTAR, vertex) {
  result <- proc_lowerSTAR_cpp(list_lowerSTAR, vertex)
  return(list(
    VE_ = result$VE_,
    CR_ = gtools::mixedsort(result$CR_)
  ))
}

#' Compute lower star in parallel 
#'
#' @param vertex Vertex data
#' @param edge Edge data
#' @param face Face data
#' @param output_dir Output directory
#' @param cores Number of cores (default: available cores-1)
#' @param batch_size Number of vertices per batch
#' @return Combined lower star results
#' @export
compute_lowerSTAR_parallel = function(vertex, edge, face, output_dir = NULL, 
                                      cores = NULL, batch_size = NULL) {
  
  if (!requireNamespace("clustermq", quietly = TRUE)) {
    stop("Package 'clustermq' required. Install with: install.packages('clustermq')")
  }
  
  # Validate inputs
  if (is.null(vertex) || nrow(vertex) == 0) {
    stop("Vertex data is empty or NULL")
  }
  if (is.null(edge) || nrow(edge) == 0) {
    stop("Edge data is empty or NULL")
  }
  if (is.null(face) || nrow(face) == 0) {
    stop("Face data is empty or NULL")
  }
  
  # Save original options to restore later
  original_scheduler = getOption("clustermq.scheduler")
  original_timeout = getOption("clustermq.worker.timeout")
  
  # Set appropriate scheduler with cleanup
  on.exit({
    options(clustermq.scheduler = original_scheduler)
    options(clustermq.worker.timeout = original_timeout)
  })
  
  if (.Platform$OS.type == "windows") {
    options(clustermq.scheduler = "multiprocess")
    message("🔧 Windows: Using 'multiprocess' scheduler")
  } else {
    options(clustermq.scheduler = "multicore")
    if (Sys.info()["sysname"] == "Darwin") {
      message("🔧 macOS: Using 'multicore' scheduler")
    } else {
      message("🔧 Linux/Unix: Using 'multicore' scheduler")
    }
  }
  
  # Set reasonable timeout (10 mins)
  options(clustermq.worker.timeout = 600)  
  
  # Cap number of cores to 32
  if (is.null(cores)) {
    cores = min(parallel::detectCores() - 1, 32)  
  } else {
    cores = min(cores, parallel::detectCores(), 32) 
  }
  
  n_vertex = nrow(vertex)
  
  # Optimal batch sizing with validation
  if (is.null(batch_size)) {
    batch_size = optimal_BATCH_size(n_vertex, cores)
  }
  batch_size = max(1, min(batch_size, n_vertex)) # Ensure valid batch size
  
  batches = split(1:n_vertex, ceiling(seq_along(1:n_vertex) / batch_size))
  total_batches = length(batches)
  
  message("🚀 ROCKET-PARALLEL clustermq: ", n_vertex, " vertices, ", 
          total_batches, " batches, ", cores, " cores")
  
  # Pre-compute connections - let C++ errors propagate naturally
  message("1 > Pre-computing vertex connections...")
  all_connections = get_vertTO_cpp(vertex, edge, face)
  
  message("2 > Building optimized vertex-simplex index...")
  precomputed_data = get_PRECOMPUTEDvert_cpp(
    all_connections$lexi_id, 
    all_connections$lexi_label
  )
  
  vertex_to_simplices = precomputed_data$vertex_index
  first_verts_z = precomputed_data$first_verts_z
  
  # CRITICAL FIX: Always returns PROPER results (for no data ~ 0-simplex)
  worker_function = function(batch_indices, vertex_i123, vertex_Z, 
                             connections, vertex_index, first_z, output_path) {
    
    batch_results = list()
    
    for (i in batch_indices) {
      v_id = vertex_i123[i]
      v_z = as.numeric(vertex_Z[i])
      
      # Ultra-fast lookup using pre-built index
      simplex_indices = vertex_index[[as.character(v_id)]]
      
      result = NULL
      
      if (!is.null(simplex_indices) && length(simplex_indices) > 0) {
        # Get valid simplices using pre-computed Z values
        valid_indices = simplex_indices[first_z[simplex_indices] <= v_z]
        
        if (length(valid_indices) > 0) {
          # Extract valid simplices
          valid_simplices = connections[valid_indices, ]
          
          # Fast sorting (only if needed for this vertex)
          if (nrow(valid_simplices) > 1) {
            order_idx = order(gtools::mixedorder(valid_simplices$lexi_label, decreasing = FALSE))
            valid_simplices = valid_simplices[order_idx, ]
          }
          
          # Create results efficiently
          result = data.frame(
            id = rep(v_id, nrow(valid_simplices)),
            lexi_label = valid_simplices$lexi_label,
            lexi_id = valid_simplices$lexi_id,
            stringsAsFactors = FALSE
          )
          
          # Write to file ONLY for non-empty results
          if (!is.null(output_path)) {
            write.table(result, file = output_path,
                        sep = "\t", append = TRUE, 
                        col.names = FALSE, row.names = FALSE,
                        quote = FALSE)
          }
        }
      }
      
      # CRITICAL FIX: Create PROPER empty data frame
      if (is.null(result)) {
        # Create properly structured empty data frame
        result = data.frame(
          id = integer(),
          lexi_label = character(),
          lexi_id = character(),
          stringsAsFactors = FALSE
        )
      }
      
      batch_results[[length(batch_results) + 1]] = result
    }
    
    return(batch_results)
  } 
  
  # Prepare data for workers
  vertex_i123 = vertex$i123
  vertex_Z = vertex$Z
  output_path = if (!is.null(output_dir)) file.path(output_dir, "lowerSTAR.txt")
  
  # Initialize output file
  if (!is.null(output_path)) {
    if (file.exists(output_path)) {
      file.remove(output_path)
    }
    # Create directory if it doesn't exist
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  }
  
  # Diagnostic message
  message("3 > Starting parallel workers with scheduler: '", 
          getOption("clustermq.scheduler"), "'")
  
  # Run with clustermq - let errors propagate naturally
  results = clustermq::Q(
    fun = worker_function,
    batch_indices = batches,
    const = list(
      vertex_i123 = vertex_i123,
      vertex_Z = vertex_Z,
      connections = all_connections,
      vertex_index = vertex_to_simplices,
      first_z = first_verts_z,
      output_path = output_path
    ),
    n_jobs = min(cores, total_batches), # Don't create more jobs than batches
    template = list(),
    export = list(gtools = "gtools"),
    chunk_size = 1
  )
  
  # Combine results
  final_results = unlist(results, recursive = FALSE)
  
  # Calculate success rate (should now be 100%!)
  success_rate = round(length(final_results) / n_vertex * 100, 1)
  message("✅ PARALLEL complete: ", length(final_results), 
          " lower star sets (", success_rate, "%)")
  
  if (success_rate < 100) {
    warning("Some vertices may not have been processed correctly. Success rate: ", success_rate, "%")
  }
  
  return(final_results)
}

#' @keywords internal
optimal_BATCH_size = function(n_vertex, cores) {
  if (n_vertex > 200000) return(10000)
  if (n_vertex > 100000) return(5000)
  return(2000)
}

#' Compute Morse complex from mesh
#'
#' @param mesh Mesh object from MeshesOperations
#' @param output_dir Optional directory to save results. If provided, writes:
#'   - vertices.txt, edges.txt, faces.txt: Mesh simplices
#'   - vector_field.txt: Gradient vector field pairs  
#'   - critical_simplices.txt: Critical simplices
#'   - lowerSTAR.txt: Lower star filtration data
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param cores Number of cores for parallel processing (default: 4)
#' @param batch_size Number of vertices per batch in parallel processing
#' @return List with Morse vector field and critical simplices
#' @export
compute_MORSE_complex = function(mesh, output_dir = NULL, parallel = TRUE, 
                                 cores = 4, batch_size = NULL) {
  
  message("___ Step 1: Computing simplices")
  simplices = get_SIMPLICES(mesh, output_dir)
  
  message("___ Step 2: Computing lower star filtration")
  if (parallel) {
    lower_star = compute_lowerSTAR_parallel(
      simplices$vertices, simplices$edges, simplices$faces, 
      output_dir, cores = cores, batch_size = batch_size
    )
  } else {
    lower_star = get_lowerSTAR(simplices$vertices, simplices$edges, simplices$faces, output_dir)
  }
  
  message("___ Step 3: Processing Morse pairings")
  morse_complex = proc_lowerSTAR(lower_star, simplices$vertices)
  
  # Write final results to files if output_dir provided
  if (!is.null(output_dir)) {
    message("___ Step 4: Writing final results to files")
    
    # Write vector field (one pair per line)
    writeLines(morse_complex$VE_, file.path(output_dir, "vector_field.txt"))
    
    # Write critical simplices (one simplex per line)
    writeLines(morse_complex$CR_, file.path(output_dir, "critical_simplices.txt"))
    
    # Write simplices as tab-separated text files (no row limits)
    write.table(simplices$vertices, file.path(output_dir, "vertices.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(simplices$edges, file.path(output_dir, "edges.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(simplices$faces, file.path(output_dir, "faces.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    message("___ Final results written to:")
    message("  - ", file.path(output_dir, "vector_field.txt"))
    message("  - ", file.path(output_dir, "critical_simplices.txt"))
    message("  - ", file.path(output_dir, "vertices.txt"))
    message("  - ", file.path(output_dir, "edges.txt")) 
    message("  - ", file.path(output_dir, "faces.txt"))
    message("  - ", file.path(output_dir, "lowerSTAR.txt"))
  }
  
  message("___ Discrete Morse Theory analysis complete!")
  return(list(
    vector_field = morse_complex$VE_,    # Gradient vector field
    critical = morse_complex$CR_,        # Critical simplices
    simplices = simplices
  ))
}

#' Fast simplex center calculation 
#' 
#' @param simplex Simplex identifier string
#' @param vertices_matrix Vertex coordinates as numeric matrix (columns: X, Y, Z)
#' @return Numeric vector of coordinates
#' @keywords internal
get_simplexCENTER <- function(simplex, vertices_matrix) {
  get_simplexCENTER_cpp(simplex, vertices_matrix)
}

#' Fast Morse complex visualization using get_simplexCENTER()
#'
#' @param morse_complex Output from compute_MORSE_complex()
#' @param projection Projection plane: "XY", "XZ", or "YZ" (default: "XZ")
#' @param point_alpha Point transparency (default: 0.6)
#' @param point_size Point size (default: 1)
#' @param max_points Maximum points to plot per category (default: 30000)
#' @param plot_gradient Whether to plot gradient arrows (default: TRUE)
#' @param plot_critical Whether to plot critical points (default: TRUE)
#' @return ggplot2 object
#' @export
visualize_MORSE_2d <- function(morse_complex, 
                               projection = "XZ",
                               point_alpha = 0.6, 
                               point_size = 1, 
                               max_points = 30000,
                               plot_gradient = TRUE,
                               plot_critical = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }
  
  vertices = morse_complex$simplices$vertices
  vector_field = morse_complex$vector_field
  critical = morse_complex$critical
  
  # Convert vertices to matrix for C++ (FASTER)
  vertices_matrix = as.matrix(data.frame(
    X = as.numeric(vertices$X),
    Y = as.numeric(vertices$Y),
    Z = as.numeric(vertices$Z)  
  ))
  
  # Separate critical points by type for proper layering
  crit_types <- sapply(strsplit(critical, " "), length)
  crit_vertices = critical[crit_types == 1]
  crit_edges = critical[crit_types == 2] 
  crit_faces = critical[crit_types == 3]
  
  # Apply sampling per type to maintain proportions
  if (plot_critical) {
    if (length(crit_edges) > max_points) {
      crit_edges = sample(crit_edges, max_points)
      message("Sampled ", max_points, " critical edges")
    }
    if (length(crit_faces) > max_points) {
      crit_faces = sample(crit_faces, max_points)
      message("Sampled ", max_points, " critical faces")
    }
    # Always show all vertices (they're the most important)
    if (length(crit_vertices) > max_points) {
      crit_vertices = sample(crit_vertices, max_points)
      message("Sampled ", max_points, " critical vertices")
    }
  }
  
  # Pre-sample gradient arrows if needed
  if (plot_gradient && length(vector_field) > max_points) {
    vector_field = sample(vector_field, max_points)
    message("Sampled ", max_points, " gradient arrows")
  }
  
  plot_data = list()
  
  # Fast gradient processing with C++
  if (plot_gradient && length(vector_field) > 0) {
    gradient_list = vector("list", length(vector_field))
    valid_count = 0
    
    for (i in seq_along(vector_field)) {
      pair = vector_field[i]
      parts = strsplit(pair, ":", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        # Use C++ version
        from_coords = get_simplexCENTER(parts[1], vertices_matrix)
        to_coords = get_simplexCENTER(parts[2], vertices_matrix)
        
        if (!any(is.na(from_coords)) && !any(is.na(to_coords))) {
          if (projection == "XY") {
            x_from = from_coords[1]; y_from = from_coords[2]
            x_to = to_coords[1]; y_to = to_coords[2]
          } else if (projection == "XZ") {
            x_from = from_coords[1]; y_from = from_coords[3]
            x_to = to_coords[1]; y_to = to_coords[3]
          } else {
            x_from = from_coords[2]; y_from = from_coords[3]
            x_to = to_coords[2]; y_to = to_coords[3]
          }
          
          if (!any(is.na(c(x_from, y_from, x_to, y_to)))) {
            gradient_list[[valid_count + 1]] = data.frame(
              x_from = x_from, y_from = y_from,
              x_to = x_to, y_to = y_to,
              type = "gradient"
            )
            valid_count = valid_count + 1
          }
        }
      }
    }
    
    if (valid_count > 0) {
      plot_data$gradient = do.call(rbind, gradient_list[1:valid_count])
    }
    message("Valid gradient arrows: ", valid_count, "/", length(vector_field))
  }
  
  # Process critical points by type for proper layering
  if (plot_critical) {
    # Process edges first (bottom layer)
    if (length(crit_edges) > 0) {
      edge_list = vector("list", length(crit_edges))
      valid_count = 0
      
      for (i in seq_along(crit_edges)) {
        crit = crit_edges[i]
        coords = get_simplexCENTER(crit, vertices_matrix)
        if (!any(is.na(coords))) {
          if (projection == "XY") {
            x_coord = coords[1]; y_coord = coords[2]
          } else if (projection == "XZ") {
            x_coord = coords[1]; y_coord = coords[3]
          } else {
            x_coord = coords[2]; y_coord = coords[3]
          }
          
          if (!any(is.na(c(x_coord, y_coord)))) {
            edge_list[[valid_count + 1]] = data.frame(
              x = x_coord, y = y_coord, 
              type = "edge",
              shape = "edge"  # NEW: shape identifier
            )
            valid_count = valid_count + 1
          }
        }
      }
      
      if (valid_count > 0) {
        plot_data$edges = do.call(rbind, edge_list[1:valid_count])
      }
      message("Valid critical edges: ", valid_count, "/", length(crit_edges))
    }
    
    # Process faces second (middle layer)
    if (length(crit_faces) > 0) {
      face_list = vector("list", length(crit_faces))
      valid_count = 0
      
      for (i in seq_along(crit_faces)) {
        crit = crit_faces[i]
        coords = get_simplexCENTER(crit, vertices_matrix)
        if (!any(is.na(coords))) {
          if (projection == "XY") {
            x_coord = coords[1]; y_coord = coords[2]
          } else if (projection == "XZ") {
            x_coord = coords[1]; y_coord = coords[3]
          } else {
            x_coord = coords[2]; y_coord = coords[3]
          }
          
          if (!any(is.na(c(x_coord, y_coord)))) {
            face_list[[valid_count + 1]] = data.frame(
              x = x_coord, y = y_coord, 
              type = "face",
              shape = "face"  # NEW: shape identifier
            )
            valid_count = valid_count + 1
          }
        }
      }
      
      if (valid_count > 0) {
        plot_data$faces = do.call(rbind, face_list[1:valid_count])
      }
      message("Valid critical faces: ", valid_count, "/", length(crit_faces))
    }
    
    # Process vertices last (top layer)
    if (length(crit_vertices) > 0) {
      vertex_list = vector("list", length(crit_vertices))
      valid_count = 0
      
      for (i in seq_along(crit_vertices)) {
        crit = crit_vertices[i]
        coords = get_simplexCENTER(crit, vertices_matrix)
        if (!any(is.na(coords))) {
          if (projection == "XY") {
            x_coord = coords[1]; y_coord = coords[2]
          } else if (projection == "XZ") {
            x_coord = coords[1]; y_coord = coords[3]
          } else {
            x_coord = coords[2]; y_coord = coords[3]
          }
          
          if (!any(is.na(c(x_coord, y_coord)))) {
            vertex_list[[valid_count + 1]] = data.frame(
              x = x_coord, y = y_coord, 
              type = "vertex",
              shape = "vertex"  # NEW: shape identifier
            )
            valid_count = valid_count + 1
          }
        }
      }
      
      if (valid_count > 0) {
        plot_data$vertices = do.call(rbind, vertex_list[1:valid_count])
      }
      message("Valid critical vertices: ", valid_count, "/", length(crit_vertices))
    }
  }
  
  # Create the plot with proper layering
  p = ggplot2::ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Morse Complex -", projection, "Projection"),
      x = ifelse(projection == "XY", "X", ifelse(projection == "XZ", "X", "Y")),
      y = ifelse(projection == "XY", "Y", ifelse(projection == "XZ", "Z", "Z"))
    )
  
  # Add gradient arrows first (bottom layer)
  if (!is.null(plot_data$gradient) && nrow(plot_data$gradient) > 0) {
    p = p + ggplot2::geom_segment(
      data = plot_data$gradient,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      color = "grey60", alpha = point_alpha * 0.7, 
      arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "cm")),
      linewidth = 0.2
    )
  }
  
  # Add critical points in reverse order (vertices on top) with different shapes
  if (!is.null(plot_data$edges) && nrow(plot_data$edges) > 0) {
    p = p + ggplot2::geom_point(
      data = plot_data$edges,
      ggplot2::aes(x = x, y = y, color = type),
      alpha = point_alpha, 
      size = point_size,
      shape = 3  # '+' shape for edges
    )
  }
  
  if (!is.null(plot_data$faces) && nrow(plot_data$faces) > 0) {
    p = p + ggplot2::geom_point(
      data = plot_data$faces,
      ggplot2::aes(x = x, y = y, color = type),
      alpha = point_alpha, 
      size = point_size,
      shape = 18  # Diamond shape for faces
    )
  }
  
  if (!is.null(plot_data$vertices) && nrow(plot_data$vertices) > 0) {
    p = p + ggplot2::geom_point(
      data = plot_data$vertices,
      ggplot2::aes(x = x, y = y, color = type),
      alpha = point_alpha, 
      size = point_size * 1.5,  # Larger for emphasis
      shape = 16  # Solid circle for vertices
    )
  }
  
  # Color scale
  p = p + ggplot2::scale_color_manual(
    values = c(
      vertex = "firebrick2", 
      edge = "steelblue2", 
      face = "forestgreen"
    ),
    name = "Critical Simplices",
    breaks = c("vertex", "edge", "face"),
    labels = c("Vertices", "Edges", "Faces")
  ) 
  
  message("Ultra-fast 2D visualization complete!")
  return(p)
}

#' Create multiple 2D projection plots
#'
#' @param morse_complex Output from compute_MORSE_complex()
#' @param point_alpha Point transparency (default: 0.6)
#' @param point_size Point size (default: 1)
#' @param max_points Maximum points to plot per category (default: 30000)
#' @param plot_gradient Whether to plot gradient arrows (default: TRUE)
#' @param plot_critical Whether to plot critical points (default: TRUE)
#' @return List of ggplot2 objects
#' @export
visualize_MORSE_2d_panel <- function(morse_complex, 
                                     point_alpha = 0.6, 
                                     point_size = 1, 
                                     max_points = 30000,
                                     plot_gradient = TRUE,
                                     plot_critical = TRUE) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required for panel plots. Install with: install.packages('patchwork')")
  }
  
  plots = list()
  projections = c("XY", "XZ", "YZ")
  
  for (proj in projections) {
    plots[[proj]] = visualize_MORSE_2d(
      morse_complex, 
      projection = proj,
      point_alpha = point_alpha, 
      point_size = point_size, 
      max_points = max_points,
      plot_gradient = plot_gradient,
      plot_critical = plot_critical
    )
  }
  
  combined_plot = plots$XY + plots$XZ + plots$YZ + 
    patchwork::plot_layout(ncol = 2, guides = "collect") & 
    ggplot2::theme(legend.position = "bottom")
  
  return(combined_plot)
}

#' Save 2D visualization to file
#'
#' @param morse_complex Output from compute_MORSE_complex() 
#' @param filename Output file name
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
#' @param panel_2d If TRUE, saves multi-panel plot. If FALSE, saves single projection plot.
#' @param ... Additional arguments to visualize_MORSE_2d() or visualize_MORSE_2d_panel()
#' @export
save_MORSE_2d <- function(morse_complex, filename, 
                          width = 10, height = 8, dpi = 300, 
                          panel_2d = FALSE, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required.")
  }
  
  if (panel_2d) {
    # Save multi-panel plot
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' required for panel plots.")
    }
    
    plot = visualize_MORSE_2d_panel(morse_complex, ...)
    
    # Adjust default size for panel plots (wider)
    if (missing(width)) width <- 12
    if (missing(height)) height <- 10
    
  } else {
    # Save single projection plot
    plot = visualize_MORSE_2d(morse_complex, ...)
  }
  
  # FIX: Explicitly set background and grid colors for saving
  plot = plot + 
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_line(color = "grey90")
    )
  
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi,
                  bg = "white")  # Critical: set background color
  
  message("Plot saved to: ", filename)
}   