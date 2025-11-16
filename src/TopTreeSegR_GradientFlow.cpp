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

// Optimized gradient network parsing
// [[Rcpp::export]]
List parse_gradient_network_fast(const std::vector<std::string>& vector_field, arma::uword n_vertices) {
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
        
        if (vertex_flow(vertex_id) == 0) {
          vertex_flow(vertex_id) = right_verts[0];
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

// Ultra-fast ascending regions - FIXED VERSION
// [[Rcpp::export]]
arma::uvec compute_ascending_regions_fast(const List& gradient_network, 
                                          const arma::uvec& minima, 
                                          arma::uword n_vertices) {
  arma::uvec ascending_regions = zeros<arma::uvec>(n_vertices);
  arma::uvec vertex_flow = gradient_network["vertex_flow"];
  List reverse_flow_r = gradient_network["reverse_flow"];
  
  // Convert minima to 1-based indexing for C++
  arma::uvec minima_1based = minima + 1;
  
  // Use unique region IDs starting from 1
  arma::uword region_counter = 1;
  
  for (arma::uword min_idx = 0; min_idx < minima_1based.n_elem; min_idx++) {
    arma::uword min_vertex = minima_1based(min_idx);
    
    // Only process if this minimum has not been assigned yet
    if (ascending_regions(min_vertex - 1) == 0) {
      std::queue<arma::uword> queue;
      queue.push(min_vertex);
      ascending_regions(min_vertex - 1) = region_counter;
      
      while (!queue.empty()) {
        arma::uword current = queue.front();
        queue.pop();
        
        arma::uvec sources = reverse_flow_r[current];
        for (arma::uword j = 0; j < sources.n_elem; j++) {
          arma::uword source = sources(j);
          if (ascending_regions(source - 1) == 0) {
            ascending_regions(source - 1) = region_counter;
            queue.push(source);
          }
        }
      }
      
      region_counter++;  // Increment for next unique region
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