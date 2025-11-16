// TopTreeSegR_GradientFlow.cpp - CORRECTED VERSION with consistent _cpp suffix
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

// Optimized gradient network parsing WITH ELEVATION SUPPORT
// [[Rcpp::export]]
List parse_gradient_network_fast_cpp(const std::vector<std::string>& vector_field, 
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
    std::vector<arma::uword> left_verts = parse_simplex_vertices_cpp(left);
    std::vector<arma::uword> right_verts = parse_simplex_vertices_cpp(right);
    
    // Vertex -> Edge pairs
    if (left_verts.size() == 1 && right_verts.size() == 2) {
      arma::uword vertex_id = left_verts[0];
      
      // CRITICAL FIX: Add bounds checking
      if (vertex_id > 0 && vertex_id <= n_vertices && 
          right_verts[0] > 0 && right_verts[0] <= n_vertices &&
          right_verts[1] > 0 && right_verts[1] <= n_vertices) {
        
        arma::uvec edge_verts = { right_verts[0], right_verts[1] };
        vertex_to_edge[vertex_id].push_back(edge_verts);
        
        // FIXED: Use elevation to determine correct flow direction with bounds checking
        if (vertex_flow(vertex_id) == 0) {
          // Get elevations for both potential targets (convert to 0-based)
          double elev1 = elevations(right_verts[0] - 1);
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
arma::uvec compute_ascending_regions_fast_cpp(const List& gradient_network, 
                                              const arma::uvec& minima, 
                                              arma::uword n_vertices) {
  arma::uvec ascending_regions = zeros<arma::uvec>(n_vertices);
  arma::uvec vertex_flow = gradient_network["vertex_flow"];
  
  // Extract reverse_flow list once
  List reverse_flow_list = gradient_network["reverse_flow"];
  
  // Convert to 1-based and assign minima in one pass
  for (arma::uword i = 0; i < minima.n_elem; i++) {
    // CRITICAL FIX: Check bounds
    if (minima(i) < n_vertices) {
      ascending_regions(minima(i)) = i + 1;
    }
  }
  
  // Use queue-based BFS for better cache performance
  std::queue<arma::uword> to_process;
  for (arma::uword i = 0; i < minima.n_elem; i++) {
    if (minima(i) < n_vertices) {
      to_process.push(minima(i));
    }
  }
  
  while (!to_process.empty()) {
    arma::uword current = to_process.front();
    to_process.pop();
    arma::uword current_region = ascending_regions(current);
    
    // Process all vertices that flow to current
    // Get the reverse flow for current vertex (1-based indexing in R list)
    if ((current + 1) < reverse_flow_list.size()) {
      arma::uvec sources = reverse_flow_list[current + 1];
      for (arma::uword j = 0; j < sources.n_elem; j++) {
        arma::uword source = sources(j) - 1; // Convert to 0-based
        if (source < n_vertices && ascending_regions(source) == 0) {
          ascending_regions(source) = current_region;
          to_process.push(source);
        }
      }
    }
  }
  
  return ascending_regions;
}

// Optimized minima connectivity
// [[Rcpp::export]]
List build_minima_connectivity_fast_cpp(const arma::uvec& minima, 
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
List build_minima_connectivity_spatial_cpp(const arma::uvec& minima, 
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

// NEW 
// Find previous vertex in gradient flow (helper function)
arma::uword find_previous_vertex_cpp(const std::vector<std::string>& vector_field, 
                                     arma::uword current_vertex) {
  for (const auto& line : vector_field) {
    size_t colon_pos = line.find(':');
    if (colon_pos == std::string::npos) continue;
    
    std::string left = line.substr(0, colon_pos);
    std::string right = line.substr(colon_pos + 1);
    
    std::vector<arma::uword> left_verts = parse_simplex_vertices_cpp(left);
    std::vector<arma::uword> right_verts = parse_simplex_vertices_cpp(right);
    
    // Look for vertex->edge pairs where edge contains current_vertex
    if (left_verts.size() == 1 && right_verts.size() == 2) {
      // Check if current_vertex is in the edge
      if (right_verts[0] == current_vertex || right_verts[1] == current_vertex) {
        return left_verts[0];  // Return the source vertex
      }
    }
  }
  return 0;  // Not found
}

// Build minima connectivity from 1-saddles with progress reporting
// [[Rcpp::export]]
List build_minima_connectivity_from_saddles_cpp(const std::vector<std::string>& critical_simplices,
                                                const std::vector<std::string>& vector_field,
                                                const arma::uvec& minima) {
  arma::uword n_minima = minima.n_elem;
  std::vector<std::vector<arma::uword>> adjacency(n_minima);
  
  // Step 1: Extract 1-saddles (critical edges)
  std::vector<std::vector<arma::uword>> saddles;
  
  Rcpp::Rcout << "  [C++] Processing " << critical_simplices.size() << " critical simplices..." << std::endl;
  
  for (size_t i = 0; i < critical_simplices.size(); i++) {
    std::vector<arma::uword> vertices = parse_simplex_vertices_cpp(critical_simplices[i]);
    
    // 1-saddles are edges (2 vertices)
    if (vertices.size() == 2) {
      saddles.push_back(vertices);
    }
    
    // Progress reporting every 10,000 simplices
    if (i % 10000 == 0 && i > 0) {
      Rcpp::Rcout << "  [C++] Processed " << i << "/" << critical_simplices.size() 
                  << " critical simplices, found " << saddles.size() << " 1-saddles" << std::endl;
    }
  }
  
  Rcpp::Rcout << "  [C++] Found " << saddles.size() << " 1-saddles total" << std::endl;
  
  // Step 2: For each saddle, find connected minima via V-paths
  Rcpp::Rcout << "  [C++] Building minima connectivity graph..." << std::endl;
  
  for (size_t s = 0; s < saddles.size(); s++) {
    const auto& saddle = saddles[s];
    std::vector<arma::uword> connected_minima_indices;
    
    // Progress reporting every 1000 saddles
    if (s % 1000 == 0 && s > 0) {
      Rcpp::Rcout << "  [C++] Processed " << s << "/" << saddles.size() << " 1-saddles" << std::endl;
    }
    
    // For each vertex in the saddle, follow V-path to find minima
    for (arma::uword saddle_vertex : saddle) {
      arma::uword current_vertex = saddle_vertex;
      std::unordered_set<arma::uword> visited;
      
      while (true) {
        // Check if current vertex is a minimum
        arma::uvec min_match = find(minima == current_vertex);
        if (min_match.n_elem > 0) {
          connected_minima_indices.push_back(min_match(0));
          break;
        }
        
        // Avoid infinite loops
        if (visited.find(current_vertex) != visited.end()) break;
        visited.insert(current_vertex);
        
        // Find the vertex that flows TO this one (reverse gradient)
        arma::uword next_vertex = find_previous_vertex_cpp(vector_field, current_vertex);
        if (next_vertex == 0) break;  // Not found
        
        current_vertex = next_vertex;
      }
    }
    
    // Remove duplicates and connect minima found through this saddle
    std::sort(connected_minima_indices.begin(), connected_minima_indices.end());
    connected_minima_indices.erase(
      std::unique(connected_minima_indices.begin(), connected_minima_indices.end()),
      connected_minima_indices.end()
    );
    
    // Connect all minima found through this saddle (if at least 2)
    if (connected_minima_indices.size() >= 2) {
      for (size_t i = 0; i < connected_minima_indices.size(); i++) {
        for (size_t j = i + 1; j < connected_minima_indices.size(); j++) {
          arma::uword idx1 = connected_minima_indices[i];
          arma::uword idx2 = connected_minima_indices[j];
          
          // Add bidirectional connection
          if (std::find(adjacency[idx1].begin(), adjacency[idx1].end(), idx2) == adjacency[idx1].end()) {
            adjacency[idx1].push_back(idx2);
          }
          if (std::find(adjacency[idx2].begin(), adjacency[idx2].end(), idx1) == adjacency[idx2].end()) {
            adjacency[idx2].push_back(idx1);
          }
        }
      }
    }
  }
  
  Rcpp::Rcout << "  [C++] Minima connectivity graph complete" << std::endl;
  
  // Convert to R list
  List adjacency_r(n_minima);
  for (arma::uword i = 0; i < n_minima; i++) {
    adjacency_r[i] = arma::uvec(adjacency[i]);
  }
  
  return adjacency_r;
}

// Optimized label propagation with priority queue 
// [[Rcpp::export]]
arma::uvec propagate_labels_minima_graph_cpp(const List& minima_graph,
                                             const arma::uvec& seed_indices,
                                             const arma::uvec& seed_labels,
                                             const arma::vec& elevations,
                                             arma::uword n_minima) {
  arma::uvec labels = zeros<arma::uvec>(n_minima);
  
  // Initialize seed labels
  for (arma::uword i = 0; i < seed_indices.n_elem; i++) {
    if (seed_indices(i) < n_minima) {
      labels(seed_indices(i)) = seed_labels(i);
    }
  }
  
  // Priority queue: (elevation, minima_index) - lowest elevation first
  std::priority_queue<std::pair<double, arma::uword>, 
                      std::vector<std::pair<double, arma::uword>>,
                      std::greater<std::pair<double, arma::uword>>> pq;
  
  // Initialize queue with seeds
  for (arma::uword i = 0; i < seed_indices.n_elem; i++) {
    if (seed_indices(i) < n_minima) {
      double elev = elevations(seed_indices(i));
      pq.push(std::make_pair(elev, seed_indices(i)));
    }
  }
  
  // Process queue
  while (!pq.empty()) {
    auto current = pq.top();
    pq.pop();
    
    // Remove the unused variable - just use the index
    arma::uword current_idx = current.second;
    arma::uword current_label = labels(current_idx);
    
    // Get neighbors from minima graph
    arma::uvec neighbors = minima_graph[current_idx];
    
    for (arma::uword j = 0; j < neighbors.n_elem; j++) {
      arma::uword neighbor_idx = neighbors(j);
      
      if (labels(neighbor_idx) == 0) {  // Unlabeled
        labels(neighbor_idx) = current_label;
        double neighbor_elev = elevations(neighbor_idx);
        pq.push(std::make_pair(neighbor_elev, neighbor_idx));
      }
    }
  }
  
  return labels;
}

// [[Rcpp::export]]
List build_minima_connectivity_optimized_vpath_cpp(const std::vector<std::string>& critical_simplices,
                                                   const std::vector<std::string>& vector_field,
                                                   const arma::uvec& minima) {
  arma::uword n_minima = minima.n_elem;
  std::vector<std::vector<arma::uword>> adjacency(n_minima);
  
  Rcpp::Rcout << "  [C++] Building minima connectivity (BIDIRECTIONAL METHOD)..." << std::endl;
  
  // Step 1: Extract 1-saddles (critical edges with 2 vertices)
  std::vector<std::vector<arma::uword>> saddle_edges;
  for (const auto& simplex_str : critical_simplices) {
    std::vector<arma::uword> vertices = parse_simplex_vertices_cpp(simplex_str);
    if (vertices.size() == 2) {
      saddle_edges.push_back(vertices);
    }
  }
  
  Rcpp::Rcout << "  [C++] Found " << saddle_edges.size() << " 1-saddles" << std::endl;
  
  // Step 2: Build fast lookup maps
  // Map: vertex → edges it flows to (forward gradient)
  std::unordered_map<arma::uword, std::vector<arma::uword>> forward_flow;
  // Map: edge vertex → source vertices that flow to it (reverse gradient)  
  std::unordered_map<arma::uword, std::vector<arma::uword>> reverse_flow;
  
  for (const auto& line : vector_field) {
    size_t colon_pos = line.find(':');
    if (colon_pos == std::string::npos) continue;
    
    std::string left = line.substr(0, colon_pos);
    std::string right = line.substr(colon_pos + 1);
    
    std::vector<arma::uword> left_verts = parse_simplex_vertices_cpp(left);
    std::vector<arma::uword> right_verts = parse_simplex_vertices_cpp(right);
    
    if (left_verts.size() == 1 && right_verts.size() == 2) {
      arma::uword source = left_verts[0];
      arma::uword edge_v1 = right_verts[0];
      arma::uword edge_v2 = right_verts[1];
      
      forward_flow[source].push_back(edge_v1);
      forward_flow[source].push_back(edge_v2);
      
      reverse_flow[edge_v1].push_back(source);
      reverse_flow[edge_v2].push_back(source);
    }
  }
  
  Rcpp::Rcout << "  [C++] Built flow maps" << std::endl;
  
  // Step 3: For each 1-saddle, find connected minima via BOTH directions
  arma::uword connections_found = 0;
  
  for (const auto& saddle_edge : saddle_edges) {
    std::unordered_set<arma::uword> connected_minima;
    
    // For each vertex in the saddle edge
    for (arma::uword saddle_vertex : saddle_edge) {
      // Direction 1: From saddle vertex FOLLOW FORWARD flow to find minima
      std::queue<arma::uword> forward_queue;
      std::unordered_set<arma::uword> forward_visited;
      
      forward_queue.push(saddle_vertex);
      forward_visited.insert(saddle_vertex);
      
      while (!forward_queue.empty()) {
        arma::uword current = forward_queue.front();
        forward_queue.pop();
        
        // Check if current is a minimum
        arma::uvec min_match = find(minima == current);
        if (min_match.n_elem > 0) {
          connected_minima.insert(min_match(0));
        }
        
        // Follow forward flow
        if (forward_flow.find(current) != forward_flow.end()) {
          for (arma::uword next_vertex : forward_flow[current]) {
            if (forward_visited.find(next_vertex) == forward_visited.end()) {
              forward_visited.insert(next_vertex);
              forward_queue.push(next_vertex);
            }
          }
        }
      }
      
      // Direction 2: From saddle vertex FOLLOW REVERSE flow to find minima  
      std::queue<arma::uword> reverse_queue;
      std::unordered_set<arma::uword> reverse_visited;
      
      reverse_queue.push(saddle_vertex);
      reverse_visited.insert(saddle_vertex);
      
      while (!reverse_queue.empty()) {
        arma::uword current = reverse_queue.front();
        reverse_queue.pop();
        
        // Check if current is a minimum
        arma::uvec min_match = find(minima == current);
        if (min_match.n_elem > 0) {
          connected_minima.insert(min_match(0));
        }
        
        // Follow reverse flow
        if (reverse_flow.find(current) != reverse_flow.end()) {
          for (arma::uword prev_vertex : reverse_flow[current]) {
            if (reverse_visited.find(prev_vertex) == reverse_visited.end()) {
              reverse_visited.insert(prev_vertex);
              reverse_queue.push(prev_vertex);
            }
          }
        }
      }
    }
    
    // Connect all minima found through this saddle
    std::vector<arma::uword> minima_list(connected_minima.begin(), connected_minima.end());
    
    if (minima_list.size() >= 2) {
      connections_found++;
      for (size_t i = 0; i < minima_list.size(); i++) {
        for (size_t j = i + 1; j < minima_list.size(); j++) {
          arma::uword idx1 = minima_list[i];
          arma::uword idx2 = minima_list[j];
          
          if (std::find(adjacency[idx1].begin(), adjacency[idx1].end(), idx2) == adjacency[idx1].end()) {
            adjacency[idx1].push_back(idx2);
          }
          if (std::find(adjacency[idx2].begin(), adjacency[idx2].end(), idx1) == adjacency[idx2].end()) {
            adjacency[idx2].push_back(idx1);
          }
        }
      }
    }
  }
  
  Rcpp::Rcout << "  [C++] Bidirectional method complete - " << connections_found 
              << "/" << saddle_edges.size() << " saddles connected minima" << std::endl;
  
  // Convert to R list
  List adjacency_r(n_minima);
  for (arma::uword i = 0; i < n_minima; i++) {
    adjacency_r[i] = arma::uvec(adjacency[i]);
  }
  
  return adjacency_r;
}