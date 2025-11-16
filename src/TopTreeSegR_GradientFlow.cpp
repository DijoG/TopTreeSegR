// TopTreeSegR_GradientFlow.cpp
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// C++ implementation of parse_simplex_vertices
// [[Rcpp::export]]
std::vector<arma::uword> parse_simplex_vertices_cpp(const std::string& simplex_str) {
  std::vector<arma::uword> vertices;
  
  // Remove dimension info and parentheses (like the R version)
  std::string clean_str = simplex_str;
  
  // Remove "dim=[0-9]+"
  size_t dim_pos = clean_str.find("dim=");
  while (dim_pos != std::string::npos) {
    size_t end_pos = clean_str.find(" ", dim_pos);
    if (end_pos == std::string::npos) {
      clean_str.erase(dim_pos);
      break;
    }
    clean_str.erase(dim_pos, end_pos - dim_pos + 1);
    dim_pos = clean_str.find("dim=");
  }
  
  // Remove parentheses
  clean_str.erase(std::remove(clean_str.begin(), clean_str.end(), '('), clean_str.end());
  clean_str.erase(std::remove(clean_str.begin(), clean_str.end(), ')'), clean_str.end());
  
  // Trim whitespace
  clean_str.erase(0, clean_str.find_first_not_of(" \t\n\r"));
  clean_str.erase(clean_str.find_last_not_of(" \t\n\r") + 1);
  
  if (clean_str.empty()) {
    return vertices;
  }
  
  // Split by commas and spaces
  const char* cstr = clean_str.c_str();
  size_t len = clean_str.length();
  size_t i = 0;
  
  while (i < len) {
    // Skip non-digits
    while (i < len && !std::isdigit(cstr[i])) i++;
    
    if (i >= len) break;
    
    // Parse number
    arma::uword num = 0;
    while (i < len && std::isdigit(cstr[i])) {
      num = num * 10 + (cstr[i] - '0');
      i++;
    }
    
    if (num > 0) {
      vertices.push_back(num);
    }
  }
  
  return vertices;
}

// C++ implementation of extract_minima_corrected
// [[Rcpp::export]]
arma::uvec extract_minima_corrected_cpp(const List& morse_complex, 
                                        const arma::mat& vertices,
                                        const std::string& temp_dir = "") {
  std::vector<arma::uword> minima_vec;
  
  // First try: Read from critical_simplices.txt file
  std::string critical_file;
  if (temp_dir.empty()) {
    critical_file = "critical_simplices.txt";
  } else {
    critical_file = temp_dir + "/critical_simplices.txt";
  }
  
  std::ifstream file(critical_file);
  if (file.is_open()) {
    Rcpp::Rcout << "  [C++] Reading minima from critical simplices file..." << std::endl;
    std::string line;
    arma::uword file_minima_count = 0;
    
    while (std::getline(file, line)) {
      // Trim the line
      line.erase(0, line.find_first_not_of(" \t\n\r"));
      line.erase(line.find_last_not_of(" \t\n\r") + 1);
      
      if (line.empty()) continue;
      
      // Parse vertices from the line
      std::vector<arma::uword> verts = parse_simplex_vertices_cpp(line);
      
      // Single vertex = minimum (0-simplex)
      if (verts.size() == 1) {
        minima_vec.push_back(verts[0]);
        file_minima_count++;
      }
    }
    file.close();
    
    Rcpp::Rcout << "  [C++] Extracted " << file_minima_count << " minima from critical simplices file" << std::endl;
    
    if (file_minima_count > 0) {
      return arma::uvec(minima_vec);
    }
  }
  
  // Fallback: Use morse_complex object (like R version)
  Rcpp::Rcout << "  [C++] Using fallback: extracting minima from morse_complex object" << std::endl;
  
  CharacterVector critical = morse_complex["critical"];
  arma::uword morse_minima_count = 0;
  
  for (int i = 0; i < critical.size(); i++) {
    std::string crit_str = Rcpp::as<std::string>(critical[i]);
    std::vector<arma::uword> verts = parse_simplex_vertices_cpp(crit_str);
    
    // Single vertex = minimum
    if (verts.size() == 1) {
      minima_vec.push_back(verts[0]);
      morse_minima_count++;
    }
  }
  
  Rcpp::Rcout << "  [C++] Found " << morse_minima_count << " minima from morse_complex object" << std::endl;
  
  return arma::uvec(minima_vec);
}

// Ultra-fast parsing using manual parsing
std::vector<arma::uword> parse_simplex_vertices_fast(const std::string& str) {
  std::vector<arma::uword> vertices;
  const char* cstr = str.c_str();
  size_t len = str.length();
  size_t i = 0;
  
  while (i < len) {
    // Skip non-digits
    while (i < len && !std::isdigit(cstr[i])) i++;
    
    if (i >= len) break;
    
    // Parse number
    arma::uword num = 0;
    while (i < len && std::isdigit(cstr[i])) {
      num = num * 10 + (cstr[i] - '0');
      i++;
    }
    
    if (num > 0) {
      vertices.push_back(num);
    }
  }
  
  return vertices;
}

// Optimized gradient network parsing WITH ELEVATION SUPPORT
// [[Rcpp::export]]
List parse_gradient_network_fast(const std::vector<std::string>& vector_field, 
                                 const arma::vec& elevations,
                                 arma::uword n_vertices) {
  // Pre-allocate for performance
  std::vector<std::vector<arma::uvec>> vertex_to_edge(n_vertices + 1);
  arma::uvec vertex_flow = zeros<arma::uvec>(n_vertices + 1);
  
  arma::uword vertex_pair_count = 0;
  arma::uword edge_count = 0;
  
  for (const auto& line : vector_field) {
    size_t colon_pos = line.find(':');
    if (colon_pos == std::string::npos) continue;
    
    std::string left = line.substr(0, colon_pos);
    std::string right = line.substr(colon_pos + 1);
    
    // Fast parsing
    std::vector<arma::uword> left_verts = parse_simplex_vertices_fast(left);
    std::vector<arma::uword> right_verts = parse_simplex_vertices_fast(right);
    
    // Vertex -> Edge pairs
    if (left_verts.size() == 1 && right_verts.size() == 2) {
      arma::uword vertex_id = left_verts[0];
      if (vertex_id > 0 && vertex_id <= n_vertices) {
        arma::uvec edge_verts = { right_verts[0], right_verts[1] };
        vertex_to_edge[vertex_id].push_back(edge_verts);
        
        // FIX: Use elevation to determine correct flow direction
        if (vertex_flow(vertex_id) == 0) {
          // Get elevations for both potential targets
          double elev1 = elevations(right_verts[0] - 1);  // 0-based indexing
          double elev2 = elevations(right_verts[1] - 1);
          
          // Flow to the vertex with LOWER elevation (downhill)
          if (elev1 < elev2) {
            vertex_flow(vertex_id) = right_verts[0];
          } else {
            vertex_flow(vertex_id) = right_verts[1];
          }
        }
        
        vertex_pair_count++;
      }
    }
    // Edge -> Face pairs
    else if (left_verts.size() == 2 && right_verts.size() == 3) {
      edge_count++;
    }
  }
  
  // Build reverse flow efficiently in single pass
  std::vector<std::vector<arma::uword>> reverse_flow(n_vertices + 1);
  for (arma::uword i = 1; i <= n_vertices; i++) {
    arma::uword target = vertex_flow(i);
    if (target > 0 && target <= n_vertices) {
      reverse_flow[target].push_back(i);
    }
  }
  
  // Convert to R efficiently
  List vertex_to_edge_r(n_vertices + 1);
  for (arma::uword i = 1; i <= n_vertices; i++) {
    if (!vertex_to_edge[i].empty()) {
      arma::umat edges(vertex_to_edge[i].size(), 2);
      for (size_t j = 0; j < vertex_to_edge[i].size(); j++) {
        edges(j, 0) = vertex_to_edge[i][j](0);
        edges(j, 1) = vertex_to_edge[i][j](1);
      }
      vertex_to_edge_r[i] = edges;
    } else {
      vertex_to_edge_r[i] = arma::umat(0, 2);
    }
  }
  
  List reverse_flow_r(n_vertices + 1);
  for (arma::uword i = 1; i <= n_vertices; i++) {
    reverse_flow_r[i] = arma::uvec(reverse_flow[i]);
  }
  
  return List::create(
    _["vertex_to_edge"] = vertex_to_edge_r,
    _["vertex_flow"] = vertex_flow,
    _["reverse_flow"] = reverse_flow_r,
    _["vertex_pair_count"] = vertex_pair_count,
    _["edge_count"] = edge_count
  );
}

// PROPER ascending regions with NO BASIN MERGING
// [[Rcpp::export]]
arma::uvec compute_ascending_regions_fast(const List& gradient_network, 
                                          const arma::uvec& minima, 
                                          arma::uword n_vertices) {
  arma::uvec ascending_regions = zeros<arma::uvec>(n_vertices);
  arma::uvec vertex_flow = gradient_network["vertex_flow"];
  
  // Extract reverse_flow list once
  List reverse_flow_list = gradient_network["reverse_flow"];
  
  // Convert to 1-based and assign minima in one pass
  for (arma::uword i = 0; i < minima.n_elem; i++) {
    ascending_regions(minima(i)) = i + 1; // Already 0-based
  }
  
  // Use queue-based BFS for better cache performance
  std::queue<arma::uword> to_process;
  for (arma::uword i = 0; i < minima.n_elem; i++) {
    to_process.push(minima(i));
  }
  
  while (!to_process.empty()) {
    arma::uword current = to_process.front();
    to_process.pop();
    arma::uword current_region = ascending_regions(current);
    
    // Process all vertices that flow to current
    // Get the reverse flow for current vertex (1-based indexing in R list)
    arma::uvec sources = reverse_flow_list[current + 1];
    for (arma::uword j = 0; j < sources.n_elem; j++) {
      arma::uword source = sources(j) - 1; // Convert to 0-based
      if (ascending_regions(source) == 0) {
        ascending_regions(source) = current_region;
        to_process.push(source);
      }
    }
  }
  
  return ascending_regions;
}

// Optimized minima connectivity
// [[Rcpp::export]]
List build_minima_connectivity_fast(const arma::uvec& minima, 
                                    const arma::mat& vertices, 
                                    double max_distance = 2.0) {
  arma::uword n_minima = minima.n_elem;
  
  // Extract coordinates for minima
  arma::mat minima_coords(n_minima, 2);
  for (arma::uword i = 0; i < n_minima; i++) {
    arma::uword idx = minima(i);
    minima_coords(i, 0) = vertices(idx, 0);  // x
    minima_coords(i, 1) = vertices(idx, 1);  // y
  }
  
  // Use Armadillo's efficient distance calculations
  List adjacency(n_minima);
  double max_dist_sq = max_distance * max_distance;
  
  for (arma::uword i = 0; i < n_minima; i++) {
    std::vector<arma::uword> neighbors_vec;
    
    for (arma::uword j = 0; j < n_minima; j++) {
      if (i == j) continue;
      
      double dx = minima_coords(i, 0) - minima_coords(j, 0);
      double dy = minima_coords(i, 1) - minima_coords(j, 1);
      double dist_sq = dx * dx + dy * dy;
      
      if (dist_sq <= max_dist_sq) {
        neighbors_vec.push_back(j);
      }
    }
    
    adjacency[i] = arma::uvec(neighbors_vec);
  }
  
  return adjacency;
}

// Spatial hashing for large datasets
// [[Rcpp::export]]
List build_minima_connectivity_spatial(const arma::uvec& minima, 
                                       const arma::mat& vertices, 
                                       double max_distance = 2.0,
                                       double grid_size = 5.0) {
  arma::uword n_minima = minima.n_elem;
  arma::mat minima_coords(n_minima, 2);
  
  // Extract coordinates
  for (arma::uword i = 0; i < n_minima; i++) {
    arma::uword idx = minima(i);
    minima_coords(i, 0) = vertices(idx, 0);
    minima_coords(i, 1) = vertices(idx, 1);
  }
  
  // Spatial hashing for O(n) performance
  std::unordered_map<arma::uword, std::vector<arma::uword>> grid;
  double max_dist_sq = max_distance * max_distance;
  
  // Assign points to grid cells
  for (arma::uword i = 0; i < n_minima; i++) {
    arma::uword grid_x = floor(minima_coords(i, 0) / grid_size);
    arma::uword grid_y = floor(minima_coords(i, 1) / grid_size);
    arma::uword cell_id = grid_x * 10000 + grid_y;  // Simple hash
    
    grid[cell_id].push_back(i);
  }
  
  List adjacency(n_minima);
  
  // Check neighbors in adjacent cells only
  for (arma::uword i = 0; i < n_minima; i++) {
    std::vector<arma::uword> neighbors_vec;
    arma::uword grid_x = floor(minima_coords(i, 0) / grid_size);
    arma::uword grid_y = floor(minima_coords(i, 1) / grid_size);
    
    // Check 3x3 neighborhood
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        arma::uword neighbor_cell = (grid_x + dx) * 10000 + (grid_y + dy);
        auto it = grid.find(neighbor_cell);
        if (it != grid.end()) {
          for (arma::uword j : it->second) {
            if (i == j) continue;
            
            double dx_val = minima_coords(i, 0) - minima_coords(j, 0);
            double dy_val = minima_coords(i, 1) - minima_coords(j, 1);
            double dist_sq = dx_val * dx_val + dy_val * dy_val;
            
            if (dist_sq <= max_dist_sq) {
              neighbors_vec.push_back(j);
            }
          }
        }
      }
    }
    
    adjacency[i] = arma::uvec(neighbors_vec);
  }
  
  return adjacency;
}

// Spatial-constrained region to tree assignment 
// [[Rcpp::export]]
arma::uvec assign_regions_to_trees(const arma::uvec& ascending_regions,
                                   const arma::uvec& seed_minima,
                                   const arma::uvec& seed_labels, 
                                   const arma::mat& points,
                                   double spatial_threshold = 3.0) {
  
  arma::uword n_points = points.n_rows;
  arma::uvec segmentation = zeros<arma::uvec>(n_points);
  
  if (seed_minima.n_elem == 0) return segmentation;
  
  // Precompute seed coordinates once
  arma::mat seed_coords = points.submat(seed_minima, arma::uvec{0, 1});
  double threshold_sq = spatial_threshold * spatial_threshold; // Use squared distance
  
  // Step 1: Process regions efficiently
  std::vector<arma::uword> region_ids;
  std::vector<arma::rowvec> region_centers;
  
  // Reserve memory to avoid reallocations
  region_ids.reserve(10000);
  region_centers.reserve(10000);
  
  for (arma::uword reg_id = 1; reg_id <= arma::max(ascending_regions); reg_id++) {
    arma::uvec region_points = find(ascending_regions == reg_id);
    
    if (region_points.n_elem > 0) {
      region_ids.push_back(reg_id);
      // Fast center calculation using column-wise mean
      arma::mat region_xy = points.submat(region_points, arma::uvec{0, 1});
      region_centers.push_back(mean(region_xy, 0));
    }
  }
  
  // Step 2: Efficient distance calculation
  arma::uvec region_to_tree = zeros<arma::uvec>(arma::max(ascending_regions) + 1);
  
  for (size_t i = 0; i < region_centers.size(); i++) {
    double min_dist_sq = std::numeric_limits<double>::max();
    arma::uword best_seed = 0;
    
    // Fast distance computation with early termination
    for (arma::uword j = 0; j < seed_coords.n_rows; j++) {
      double dx = region_centers[i](0) - seed_coords(j, 0);
      double dy = region_centers[i](1) - seed_coords(j, 1);
      double dist_sq = dx * dx + dy * dy;
      
      if (dist_sq < min_dist_sq) {
        min_dist_sq = dist_sq;
        best_seed = j;
        // Early termination if we find a very close seed
        if (min_dist_sq < 0.1) break;
      }
    }
    
    if (min_dist_sq <= threshold_sq) {
      region_to_tree(region_ids[i]) = seed_labels(best_seed);
    }
  }
  
  // Step 3: Fast segmentation assignment
  for (arma::uword i = 0; i < n_points; i++) {
    arma::uword region_id = ascending_regions(i);
    if (region_id > 0) {
      segmentation(i) = region_to_tree(region_id);
    }
  }
  
  // Step 4: Optimized unassigned point handling
  arma::uvec unassigned = find(segmentation == 0);
  if (unassigned.n_elem > 0) {
    arma::uvec assigned = find(segmentation > 0);
    arma::mat unassigned_xy = points.submat(unassigned, arma::uvec{0, 1});
    arma::mat assigned_xy = points.submat(assigned, arma::uvec{0, 1});
    
    arma::uvec nn_indices(unassigned.n_elem);
    
    // Efficient brute force with squared distances
    for (arma::uword i = 0; i < unassigned.n_elem; i++) {
      double min_dist_sq = std::numeric_limits<double>::max();
      arma::uword best_idx = 0;
      
      for (arma::uword j = 0; j < assigned.n_elem; j++) {
        double dx = unassigned_xy(i, 0) - assigned_xy(j, 0);
        double dy = unassigned_xy(i, 1) - assigned_xy(j, 1);
        double dist_sq = dx * dx + dy * dy;
        
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best_idx = j;
          // Early termination for very close points
          if (min_dist_sq < 0.01) break;
        }
      }
      nn_indices(i) = best_idx;
    }
    
    segmentation(unassigned) = segmentation(assigned(nn_indices));
  }
  
  return segmentation;
}