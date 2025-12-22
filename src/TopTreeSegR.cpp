// TopTreeSegR.cpp - C++ implementation with gradient flow and Morse-Smale segmentation
#include <RcppArmadillo.h>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Data structures
struct Point3D {
  double x, y, z;
  int id;
  int original_id;
  int label;
};

struct Edge {
  uword v1, v2;
};

struct GradientNetwork {
  std::vector<std::vector<Edge>> vertex_to_edge;
  uvec vertex_flow;
  std::vector<std::vector<uword>> reverse_flow;
  uword vertex_pair_count;
  uword edge_count;
};

// Parsing functions
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

// Function to map mesh labels to original points
// [[Rcpp::export]]
IntegerVector map_mesh_to_points_cpp(const arma::mat& mesh_coords,
                                     const IntegerVector& mesh_labels,
                                     const arma::mat& point_coords,
                                     double max_distance = 2.0) {
  
  uword n_points = point_coords.n_rows;
  uword n_mesh = mesh_coords.n_rows;
  IntegerVector point_labels(n_points, 0);
  
  for (uword i = 0; i < n_points; i++) {
    double min_dist = std::numeric_limits<double>::max();
    int best_label = 0;
    
    for (uword j = 0; j < n_mesh; j++) {
      double dx = point_coords(i, 0) - mesh_coords(j, 0);
      double dy = point_coords(i, 1) - mesh_coords(j, 1);
      double dz = point_coords(i, 2) - mesh_coords(j, 2);
      
      double dist = dx * dx + dy * dy + dz * dz;
      
      if (dist < min_dist) {
        min_dist = dist;
        best_label = mesh_labels[j];
      }
    }
    
    if (std::sqrt(min_dist) <= max_distance) {
      point_labels[i] = best_label;
    }
  }
  
  return point_labels;
}

// Main segmentation class for gradient flow method
class TreeSegmenter {
private:
  // Data
  std::vector<Point3D> vertices;
  std::vector<uword> minima;
  GradientNetwork gradient_net;
  uvec ascending_regions;
  
  // Parameters
  double stem_threshold;
  double spatial_eps;
  double max_distance;
  double grid_size;
  
public:
  TreeSegmenter(double stem_ht = 0.5, double eps = 1.0, 
                double max_dist = 2.0, double grid_sz = 5.0) :
  stem_threshold(stem_ht), spatial_eps(eps), 
  max_distance(max_dist), grid_size(grid_sz) {}
  
  // Main segmentation function for gradient method
  List segmentTrees(DataFrame vertices_df, List morse_complex_data) {
    Rcout << "\n=== Tree Segmentation with Gradient Flow ===\n";
    
    // 1. Load vertices
    loadVertices(vertices_df);
    
    // 2. Extract Morse complex data
    extractMorseData(morse_complex_data);
    
    // 3. Build gradient network
    buildGradientNetwork(morse_complex_data);
    
    // 4. Find seed minima
    std::vector<uword> seeds = findSeedMinima();
    
    // 5. Propagate labels
    std::vector<int> labels = propagateLabels(seeds);
    
    // 6. Create output
    return createOutput(labels, seeds);
  }
  
private:
  void loadVertices(DataFrame vertices_df) {
    NumericVector x = vertices_df["X"];
    NumericVector y = vertices_df["Y"];
    NumericVector z = vertices_df["Z"];
    IntegerVector original_ids = vertices_df["i123"];
    
    vertices.resize(x.size());
    for (int i = 0; i < x.size(); i++) {
      vertices[i] = {x[i], y[i], z[i], i, original_ids[i], 0};
    }
    
    Rcout << "Loaded " << vertices.size() << " mesh vertices\n";
  }
  
  void extractMorseData(List morse_complex_data) {
    minima.clear();
    
    if (morse_complex_data.containsElementNamed("critical")) {
      CharacterVector crit = morse_complex_data["critical"];
      
      for (int i = 0; i < crit.size(); i++) {
        std::vector<uword> verts = parse_simplex_vertices_cpp(as<std::string>(crit[i]));
        
        // Single vertex = minimum (0-simplex)
        if (verts.size() == 1) {
          uword vertex_id = verts[0];
          if (vertex_id > 0 && vertex_id <= vertices.size()) {
            minima.push_back(vertex_id);
          }
        }
      }
    }
    
    Rcout << "Found " << minima.size() << " minima\n";
  }
  
  void buildGradientNetwork(List morse_complex_data) {
    Rcout << "Building gradient network...\n";
    
    CharacterVector vector_field = morse_complex_data["vector_field"];
    uword n_vertices = vertices.size();
    
    // Initialize gradient network
    gradient_net.vertex_to_edge.resize(n_vertices + 1);
    gradient_net.vertex_flow = zeros<uvec>(n_vertices + 1);
    gradient_net.reverse_flow.resize(n_vertices + 1);
    gradient_net.vertex_pair_count = 0;
    gradient_net.edge_count = 0;
    
    // Extract elevations
    vec elevations(n_vertices);
    for (uword i = 0; i < n_vertices; i++) {
      elevations(i) = vertices[i].z;
    }
    
    // Parse vector field
    for (int i = 0; i < vector_field.size(); i++) {
      std::string line = as<std::string>(vector_field[i]);
      size_t colon_pos = line.find(':');
      if (colon_pos == std::string::npos) continue;
      
      std::string left = line.substr(0, colon_pos);
      std::string right = line.substr(colon_pos + 1);
      
      std::vector<uword> left_verts = parse_simplex_vertices_cpp(left);
      std::vector<uword> right_verts = parse_simplex_vertices_cpp(right);
      
      // Vertex -> Edge pairs
      if (left_verts.size() == 1 && right_verts.size() == 2) {
        uword vertex_id = left_verts[0];
        
        if (vertex_id > 0 && vertex_id <= n_vertices && 
            right_verts[0] > 0 && right_verts[0] <= n_vertices &&
            right_verts[1] > 0 && right_verts[1] <= n_vertices) {
          
          // Store edge as struct
          Edge edge = {right_verts[0], right_verts[1]};
          gradient_net.vertex_to_edge[vertex_id].push_back(edge);
          
          // Determine flow direction (downhill)
          if (gradient_net.vertex_flow(vertex_id) == 0) {
            double elev1 = elevations(right_verts[0] - 1);
            double elev2 = elevations(right_verts[1] - 1);
            
            if (elev1 < elev2) {
              gradient_net.vertex_flow(vertex_id) = right_verts[0];
            } else {
              gradient_net.vertex_flow(vertex_id) = right_verts[1];
            }
          }
          
          gradient_net.vertex_pair_count++;
        }
      }
      // Edge -> Face pairs
      else if (left_verts.size() == 2 && right_verts.size() == 3) {
        gradient_net.edge_count++;
      }
    }
    
    // Build reverse flow (for ascending manifolds)
    for (uword i = 1; i <= n_vertices; i++) {
      uword target = gradient_net.vertex_flow(i);
      if (target > 0 && target <= n_vertices) {
        gradient_net.reverse_flow[target].push_back(i);
      }
    }
    
    Rcout << "Gradient network: " << gradient_net.vertex_pair_count 
          << " vertex-edge pairs, " << gradient_net.edge_count << " edge-face pairs\n";
  }
  
  std::vector<uword> findSeedMinima() {
    std::vector<uword> seeds;
    
    // Find minima in stem region
    std::vector<std::pair<double, uword>> low_minima;
    for (uword min_id : minima) {
      if (vertices[min_id - 1].z < stem_threshold) {
        low_minima.push_back({vertices[min_id - 1].z, min_id});
      }
    }
    
    if (low_minima.empty()) {
      Rcout << "No low minima found, using all minima\n";
      for (uword min_id : minima) {
        low_minima.push_back({vertices[min_id - 1].z, min_id});
      }
    }
    
    // Sort by elevation
    std::sort(low_minima.begin(), low_minima.end());
    
    // Spatial clustering for seeds
    std::unordered_map<uword, int> cluster_map;
    int cluster_id = 1;
    
    for (size_t i = 0; i < low_minima.size(); i++) {
      uword v1 = low_minima[i].second;
      const Point3D& p1 = vertices[v1 - 1];
      
      if (cluster_map.find(v1) == cluster_map.end()) {
        cluster_map[v1] = cluster_id;
        
        // Find nearby minima
        for (size_t j = i + 1; j < low_minima.size(); j++) {
          uword v2 = low_minima[j].second;
          
          if (cluster_map.find(v2) == cluster_map.end()) {
            const Point3D& p2 = vertices[v2 - 1];
            
            double dist = std::sqrt(std::pow(p1.x - p2.x, 2) + 
                                    std::pow(p1.y - p2.y, 2));
            
            if (dist < spatial_eps) {
              cluster_map[v2] = cluster_id;
            }
          }
        }
        
        // Select lowest minimum in cluster as seed
        uword lowest_in_cluster = v1;
        double lowest_elevation = p1.z;
        
        for (const auto& entry : cluster_map) {
          if (entry.second == cluster_id) {
            if (vertices[entry.first - 1].z < lowest_elevation) {
              lowest_in_cluster = entry.first;
              lowest_elevation = vertices[entry.first - 1].z;
            }
          }
        }
        
        seeds.push_back(lowest_in_cluster);
        cluster_id++;
      }
    }
    
    Rcout << "Selected " << seeds.size() << " seed minima\n";
    return seeds;
  }
  
  std::vector<int> propagateLabels(const std::vector<uword>& seeds) {
    std::vector<int> labels(vertices.size(), 0);
    
    // First, assign ascending regions
    computeAscendingRegions();
    
    // Create mapping from minima to tree IDs
    std::unordered_map<uword, int> minima_to_tree;
    for (size_t i = 0; i < seeds.size(); i++) {
      minima_to_tree[seeds[i]] = i + 1;
    }
    
    // Assign labels to vertices based on ascending regions
    for (size_t i = 0; i < vertices.size(); i++) {
      uword region_id = ascending_regions(i);
      if (region_id > 0 && region_id <= minima.size()) {
        uword min_id = minima[region_id - 1];
        auto it = minima_to_tree.find(min_id);
        if (it != minima_to_tree.end()) {
          labels[i] = it->second;
        }
      }
    }
    
    // Handle unassigned vertices with spatial propagation
    handleUnassignedVertices(labels, seeds);
    
    return labels;
  }
  
  void computeAscendingRegions() {
    uword n_vertices = vertices.size();
    ascending_regions = zeros<uvec>(n_vertices);
    
    // Assign minima to ascending regions
    for (uword i = 0; i < minima.size(); i++) {
      uword min_id = minima[i];
      if (min_id > 0 && min_id <= n_vertices) {
        ascending_regions(min_id - 1) = i + 1;
      }
    }
    
    // BFS propagation through reverse flow (ascending manifolds)
    std::queue<uword> to_process;
    for (uword min_id : minima) {
      if (min_id > 0 && min_id <= n_vertices) {
        to_process.push(min_id - 1);
      }
    }
    
    while (!to_process.empty()) {
      uword current = to_process.front();
      to_process.pop();
      uword current_region = ascending_regions(current);
      
      // Process vertices that flow to current (uphill)
      uword current_1based = current + 1;
      for (uword source : gradient_net.reverse_flow[current_1based]) {
        uword source_idx = source - 1;
        if (source_idx < n_vertices && ascending_regions(source_idx) == 0) {
          ascending_regions(source_idx) = current_region;
          to_process.push(source_idx);
        }
      }
    }
    
    Rcout << "Ascending regions computed: " 
          << sum(ascending_regions > 0) << " vertices assigned\n";
  }
  
  void handleUnassignedVertices(std::vector<int>& labels, 
                                const std::vector<uword>& seeds) {
    // Find unassigned vertices
    std::vector<size_t> unassigned;
    for (size_t i = 0; i < labels.size(); i++) {
      if (labels[i] == 0) {
        unassigned.push_back(i);
      }
    }
    
    if (unassigned.empty()) {
      Rcout << "All vertices assigned via gradient flow\n";
      return;
    }
    
    Rcout << "Assigning " << unassigned.size() << " unassigned vertices...\n";
    
    // Assign to nearest seed spatially
    for (size_t idx : unassigned) {
      const Point3D& p = vertices[idx];
      double min_dist = std::numeric_limits<double>::max();
      int best_label = 0;
      
      for (size_t j = 0; j < seeds.size(); j++) {
        const Point3D& seed_p = vertices[seeds[j] - 1];
        double dist = std::sqrt(std::pow(p.x - seed_p.x, 2) + 
                                std::pow(p.y - seed_p.y, 2) + 
                                std::pow(p.z - seed_p.z, 2));
        
        if (dist < min_dist) {
          min_dist = dist;
          best_label = j + 1;
        }
      }
      
      if (best_label > 0) {
        labels[idx] = best_label;
      }
    }
    
    Rcout << "Completed spatial assignment\n";
  }
  
  List createOutput(const std::vector<int>& labels, 
                    const std::vector<uword>& seeds) {
    // Create mesh coordinates matrix
    arma::mat mesh_coords(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); i++) {
      mesh_coords(i, 0) = vertices[i].x;
      mesh_coords(i, 1) = vertices[i].y;
      mesh_coords(i, 2) = vertices[i].z;
    }
    
    // Create mesh labels vector
    IntegerVector mesh_labels(labels.size());
    for (size_t i = 0; i < labels.size(); i++) {
      mesh_labels[i] = labels[i];
    }
    
    // Create seeds vector
    IntegerVector seeds_vec(seeds.size());
    for (size_t i = 0; i < seeds.size(); i++) {
      seeds_vec[i] = seeds[i];
    }
    
    // Count trees
    std::unordered_set<int> unique_labels;
    for (int label : labels) {
      if (label > 0) unique_labels.insert(label);
    }
    
    Rcout << "\n=== Gradient Flow Segmentation Complete ===\n";
    Rcout << "Detected " << unique_labels.size() << " trees in mesh\n";
    Rcout << "Seeds: " << seeds.size() << "\n";
    Rcout << "Minima: " << minima.size() << "\n";
    
    return List::create(
      _["mesh_labels"] = mesh_labels,
      _["mesh_coords"] = mesh_coords,
      _["minima"] = minima,
      _["ascending_regions"] = ascending_regions,
      _["seeds"] = seeds_vec,
      _["n_trees"] = unique_labels.size(),
      _["method"] = "gradient"
    );
  }
};

// Rcpp interface for gradient flow segmentation
// [[Rcpp::export]]
List tree_segment_cpp(DataFrame vertices_df, 
                      List morse_complex_data,
                      double stem_height = 0.5,
                      double spatial_eps = 1.0,
                      double max_distance = 2.0,
                      double grid_size = 5.0) {
  
  try {
    TreeSegmenter segmenter(stem_height, spatial_eps, max_distance, grid_size);
    return segmenter.segmentTrees(vertices_df, morse_complex_data);
  } catch (const std::exception& e) {
    Rcpp::stop("Segmentation error: %s", e.what());
  } catch (...) {
    Rcpp::stop("Unknown segmentation error");
  }
}

// Mathematically correct Morse-Smale segmentation using ascending manifolds
// [[Rcpp::export]]
List morse_smale_segment_cpp(DataFrame vertices_df,
                             List morse_complex_data,
                             double stem_height = 0.5,
                             double spatial_eps = 1.0) {
  
  try {
    Rcout << "\n=== Morse-Smale Complex Segmentation (Ascending Manifolds) ===\n";
    
    // 1. Load vertices
    NumericVector x = vertices_df["X"];
    NumericVector y = vertices_df["Y"];
    NumericVector z = vertices_df["Z"];
    IntegerVector i123 = vertices_df["i123"];
    
    int n_vertices = x.size();
    
    // 2. Build i123 to index mapping
    std::unordered_map<int, int> i123_to_idx;
    for (int i = 0; i < n_vertices; i++) {
      i123_to_idx[i123[i]] = i;
    }
    
    Rcout << "Loaded " << n_vertices << " vertices\n";
    
    // 3. Find minima (tree roots) below stem_height
    std::vector<int> potential_seeds;
    
    if (morse_complex_data.containsElementNamed("critical")) {
      CharacterVector crit = morse_complex_data["critical"];
      
      for (int i = 0; i < crit.size(); i++) {
        std::vector<arma::uword> verts = parse_simplex_vertices_cpp(as<std::string>(crit[i]));
        
        // Single vertex = minimum (0-simplex)
        if (verts.size() == 1) {
          int vertex_id = verts[0];
          auto it = i123_to_idx.find(vertex_id);
          if (it != i123_to_idx.end()) {
            int idx = it->second;
            if (z[idx] < stem_height) {
              potential_seeds.push_back(vertex_id);
            }
          }
        }
      }
    }
    
    if (potential_seeds.empty()) {
      Rcout << "No minima below stem_height, using all low points\n";
      for (int i = 0; i < n_vertices; i++) {
        if (z[i] < stem_height) {
          potential_seeds.push_back(i123[i]);
        }
      }
    }
    
    Rcout << "Found " << potential_seeds.size() << " potential seed minima\n";
    
    // 4. Spatial clustering for seeds
    std::unordered_map<int, int> seed_to_cluster;
    std::vector<int> seeds;
    int cluster_id = 1;
    
    for (size_t i = 0; i < potential_seeds.size(); i++) {
      int seed_i123 = potential_seeds[i];
      auto it = i123_to_idx.find(seed_i123);
      if (it == i123_to_idx.end()) continue;
      
      if (seed_to_cluster.find(seed_i123) == seed_to_cluster.end()) {
        seed_to_cluster[seed_i123] = cluster_id;
        int seed_idx = it->second;
        
        // Find nearby minima
        for (size_t j = i + 1; j < potential_seeds.size(); j++) {
          int other_i123 = potential_seeds[j];
          auto it_other = i123_to_idx.find(other_i123);
          if (it_other == i123_to_idx.end()) continue;
          
          if (seed_to_cluster.find(other_i123) == seed_to_cluster.end()) {
            int other_idx = it_other->second;
            double dx = x[seed_idx] - x[other_idx];
            double dy = y[seed_idx] - y[other_idx];
            double dist = std::sqrt(dx*dx + dy*dy);
            
            if (dist < spatial_eps) {
              seed_to_cluster[other_i123] = cluster_id;
            }
          }
        }
        
        // Select lowest minimum in cluster as seed
        int lowest_i123 = seed_i123;
        double lowest_z = z[seed_idx];
        
        for (const auto& entry : seed_to_cluster) {
          if (entry.second == cluster_id) {
            int idx = i123_to_idx[entry.first];
            if (z[idx] < lowest_z) {
              lowest_i123 = entry.first;
              lowest_z = z[idx];
            }
          }
        }
        
        seeds.push_back(lowest_i123);
        cluster_id++;
      }
    }
    
    Rcout << "Selected " << seeds.size() << " seed minima after spatial clustering\n";
    
    // 5. Parse gradient field to build flow maps
    CharacterVector vector_field = morse_complex_data["vector_field"];
    
    // Primary flow: vertex -> lower vertex (downhill)
    std::unordered_map<int, int> flow_downhill;
    // Reverse flow: lower vertex -> vertices that flow to it (for ascending manifolds)
    std::unordered_map<int, std::vector<int>> flows_to_this_vertex;
    
    int gradient_pairs = 0;
    
    for (int i = 0; i < vector_field.size(); i++) {
      std::string line = as<std::string>(vector_field[i]);
      size_t colon_pos = line.find(':');
      if (colon_pos == std::string::npos) continue;
      
      std::string left = line.substr(0, colon_pos);
      std::string right = line.substr(colon_pos + 1);
      
      std::vector<arma::uword> left_verts = parse_simplex_vertices_cpp(left);
      std::vector<arma::uword> right_verts = parse_simplex_vertices_cpp(right);
      
      // Vertex -> Edge pair
      if (left_verts.size() == 1 && right_verts.size() == 2) {
        int vertex_i123 = left_verts[0];
        int edge_v1 = right_verts[0];
        int edge_v2 = right_verts[1];
        
        if (i123_to_idx.find(vertex_i123) != i123_to_idx.end() &&
            i123_to_idx.find(edge_v1) != i123_to_idx.end() &&
            i123_to_idx.find(edge_v2) != i123_to_idx.end()) {
          
          int v1_idx = i123_to_idx[edge_v1];
          int v2_idx = i123_to_idx[edge_v2];
          
          double z_v1 = z[v1_idx];
          double z_v2 = z[v2_idx];
          
          // Gradient flows FROM vertex TO lower vertex
          int lower_vertex;
          if (z_v1 < z_v2) {
            lower_vertex = edge_v1;
          } else {
            lower_vertex = edge_v2;
          }
          
          // Store downhill flow
          flow_downhill[vertex_i123] = lower_vertex;
          
          // Store reverse flow (for ascending manifolds)
          flows_to_this_vertex[lower_vertex].push_back(vertex_i123);
          
          gradient_pairs++;
        }
      }
    }
    
    Rcout << "Built gradient flow graph: " << gradient_pairs << " gradient pairs\n";
    
    // 6. Map seeds to tree IDs
    std::unordered_map<int, int> seed_to_tree;
    std::queue<int> bfs_queue;
    
    for (size_t i = 0; i < seeds.size(); i++) {
      int seed_i123 = seeds[i];
      if (i123_to_idx.find(seed_i123) != i123_to_idx.end()) {
        seed_to_tree[seed_i123] = i + 1;  // Tree IDs start at 1
        bfs_queue.push(seed_i123);
        int seed_idx = i123_to_idx[seed_i123];
        Rcout << "Seed " << seed_i123 << " -> Tree " << (i + 1) 
              << " (z=" << z[seed_idx] << ")\n";
      }
    }
    
    if (seed_to_tree.empty()) {
      Rcpp::stop("No valid seeds found!");
    }
    
    // 7. Compute ASCENDING MANIFOLDS via BFS following REVERSE flow
    // This correctly computes all vertices whose downhill flow reaches each minimum
    IntegerVector labels(n_vertices, 0);
    
    // Assign seeds first
    for (const auto& pair : seed_to_tree) {
      int seed_i123 = pair.first;
      int tree_id = pair.second;
      int seed_idx = i123_to_idx[seed_i123];
      labels[seed_idx] = tree_id;
    }
    
    int labeled_count = 0;
    
    while (!bfs_queue.empty()) {
      int current_i123 = bfs_queue.front();
      bfs_queue.pop();
      
      int current_idx = i123_to_idx[current_i123];
      int current_tree = labels[current_idx];
      
      // Find vertices that flow TO current (uphill direction)
      // These are in current's ascending manifold
      auto reverse_flow_it = flows_to_this_vertex.find(current_i123);
      if (reverse_flow_it != flows_to_this_vertex.end()) {
        for (int uphill_i123 : reverse_flow_it->second) {
          auto uphill_it = i123_to_idx.find(uphill_i123);
          if (uphill_it != i123_to_idx.end()) {
            int uphill_idx = uphill_it->second;
            if (labels[uphill_idx] == 0) {
              labels[uphill_idx] = current_tree;
              bfs_queue.push(uphill_i123);
              labeled_count++;
            }
          }
        }
      }
    }
    
    Rcout << "Computed ascending manifolds: " << labeled_count 
          << " vertices assigned via Morse-Smale flow\n";
    
    // 8. Handle any unassigned vertices (spatial fallback)
    int unassigned = 0;
    for (int i = 0; i < n_vertices; i++) {
      if (labels[i] == 0) unassigned++;
    }
    
    if (unassigned > 0) {
      Rcout << "Spatially assigning " << unassigned << " unassigned vertices\n";
      
      std::vector<int> seed_idxs;
      std::vector<int> tree_ids;
      for (const auto& pair : seed_to_tree) {
        seed_idxs.push_back(i123_to_idx[pair.first]);
        tree_ids.push_back(pair.second);
      }
      
      for (int i = 0; i < n_vertices; i++) {
        if (labels[i] == 0) {
          double min_dist = 1e100;
          int best_tree = 0;
          
          for (size_t j = 0; j < seed_idxs.size(); j++) {
            double dx = x[i] - x[seed_idxs[j]];
            double dy = y[i] - y[seed_idxs[j]];
            double dz = z[i] - z[seed_idxs[j]];
            double dist = dx*dx + dy*dy + dz*dz;
            
            if (dist < min_dist) {
              min_dist = dist;
              best_tree = tree_ids[j];
            }
          }
          
          labels[i] = best_tree;
        }
      }
    }
    
    // 9. Prepare output
    arma::mat mesh_coords(n_vertices, 3);
    for (int i = 0; i < n_vertices; i++) {
      mesh_coords(i, 0) = x[i];
      mesh_coords(i, 1) = y[i];
      mesh_coords(i, 2) = z[i];
    }
    
    // Count trees
    std::unordered_set<int> unique_trees;
    for (int label : labels) {
      if (label > 0) unique_trees.insert(label);
    }
    
    // Prepare seeds vector
    IntegerVector seeds_vec(seeds.size());
    for (size_t i = 0; i < seeds.size(); i++) {
      seeds_vec[i] = seeds[i];
    }
    
    // Extract minima
    std::vector<int> minima_vec;
    if (morse_complex_data.containsElementNamed("critical")) {
      CharacterVector crit = morse_complex_data["critical"];
      
      for (int i = 0; i < crit.size(); i++) {
        std::vector<arma::uword> verts = parse_simplex_vertices_cpp(as<std::string>(crit[i]));
        
        if (verts.size() == 1) {
          int vertex_id = verts[0];
          auto it = i123_to_idx.find(vertex_id);
          if (it != i123_to_idx.end()) {
            minima_vec.push_back(it->second + 1);  // 1-based index
          }
        }
      }
    }
    
    IntegerVector minima(minima_vec.size());
    for (size_t i = 0; i < minima_vec.size(); i++) {
      minima[i] = minima_vec[i];
    }
    
    Rcout << "\n=== Morse-Smale Segmentation Complete ===\n";
    Rcout << "Detected " << unique_trees.size() << " trees\n";
    
    std::unordered_map<int, int> tree_counts;
    for (int label : labels) {
      if (label > 0) tree_counts[label]++;
    }
    
    for (const auto& pair : tree_counts) {
      Rcout << "  Tree " << pair.first << ": " << pair.second << " vertices\n";
    }
    
    return List::create(
      _["mesh_labels"] = labels,
      _["mesh_coords"] = mesh_coords,
      _["seeds"] = seeds_vec,
      _["minima"] = minima,
      _["n_trees"] = (int)unique_trees.size(),
      _["labeled_via_msmale"] = labeled_count,
      _["spatially_assigned"] = unassigned,
      _["method"] = "morse_smale",
      _["ascending_regions"] = IntegerVector(0)  // Not used in Morse-Smale
    );
    
  } catch (const std::exception& e) {
    Rcpp::stop("Morse-Smale segmentation error: %s", e.what());
  } catch (...) {
    Rcpp::stop("Unknown Morse-Smale segmentation error");
  }
}

// Helper function to extract mesh vertices from Morse complex
// [[Rcpp::export]]
arma::mat get_mesh_coords_cpp(DataFrame vertices_df) {
  NumericVector x = vertices_df["X"];
  NumericVector y = vertices_df["Y"];
  NumericVector z = vertices_df["Z"];
  
  arma::mat coords(x.size(), 3);
  for (int i = 0; i < x.size(); i++) {
    coords(i, 0) = x[i];
    coords(i, 1) = y[i];
    coords(i, 2) = z[i];
  }
  
  return coords;
}

// Markov Random Field Boundary Refinement (MRFBR), smoothing with spatial grid acceleration
// [[Rcpp::export]]
IntegerVector MRFBR_cpp(const arma::mat& coords,
                        const IntegerVector& labels,
                        double spatial_sigma = 1.0,
                        double intensity = 1.0,
                        int iterations = 2) {
  
  int n = coords.n_rows;
  IntegerVector current_labels = clone(labels);
  
  Rcout << "Fast MRF: " << n << " points\n";
  
  // 1. Build spatial grid for fast neighbor search
  double grid_size = 2.0; // Should be > max_dist
  double max_dist = 1.5;  // Neighborhood radius
  
  // Find bounds
  double min_x = arma::min(coords.col(0));
  double max_x = arma::max(coords.col(0));
  double min_y = arma::min(coords.col(1));
  double max_y = arma::max(coords.col(1));
  
  int grid_x = std::ceil((max_x - min_x) / grid_size) + 1;
  int grid_y = std::ceil((max_y - min_y) / grid_size) + 1;
  
  std::vector<std::vector<int>> grid(grid_x * grid_y);
  
  // Assign points to grid cells
  for (int i = 0; i < n; i++) {
    if (current_labels[i] == 0) continue;
    
    int gx = (int)((coords(i, 0) - min_x) / grid_size);
    int gy = (int)((coords(i, 1) - min_y) / grid_size);
    int cell_idx = gy * grid_x + gx;
    
    if (cell_idx >= 0 && cell_idx < grid_x * grid_y) {
      grid[cell_idx].push_back(i);
    }
  }
  
  // 2. Get unique labels
  std::unordered_set<int> unique_labels_set;
  for (int i = 0; i < n; i++) {
    if (current_labels[i] > 0) {
      unique_labels_set.insert(current_labels[i]);
    }
  }
  std::vector<int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
  
  // 3. MRF iterations (only on boundary points)
  for (int iter = 0; iter < iterations; iter++) {
    IntegerVector new_labels = clone(current_labels);
    int changes = 0;
    
    // Pre-compute which points are boundary points
    std::vector<bool> is_boundary(n, false);
    int boundary_count = 0;
    
#pragma omp parallel for reduction(+:boundary_count) schedule(dynamic, 100)
    for (int i = 0; i < n; i++) {
      if (current_labels[i] == 0) continue;
      
      int gx = (int)((coords(i, 0) - min_x) / grid_size);
      int gy = (int)((coords(i, 1) - min_y) / grid_size);
      int current_label = current_labels[i];
      
      // Check neighboring grid cells
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          int ngx = gx + dx;
          int ngy = gy + dy;
          
          if (ngx < 0 || ngx >= grid_x || ngy < 0 || ngy >= grid_y) continue;
          
          int cell_idx = ngy * grid_x + ngx;
          
          for (int j : grid[cell_idx]) {
            if (i == j) continue;
            
            // Quick distance check
            double dx_val = coords(i, 0) - coords(j, 0);
            double dy_val = coords(i, 1) - coords(j, 1);
            
            if (dx_val*dx_val + dy_val*dy_val < max_dist*max_dist) {
              if (current_labels[j] != current_label) {
                is_boundary[i] = true;
                boundary_count++;
                // Break out of loops
                dx = 2; dy = 2; // Exit both loops
                break;
              }
            }
          }
          if (is_boundary[i]) break;
        }
        if (is_boundary[i]) break;
      }
    }
    
    Rcout << "  Iteration " << iter + 1 << ": " << boundary_count << " boundary points\n";
    
    // Only process boundary points
#pragma omp parallel for reduction(+:changes) schedule(dynamic, 50)
    for (int i = 0; i < n; i++) {
      if (!is_boundary[i] || current_labels[i] == 0) continue;
      
      int gx = (int)((coords(i, 0) - min_x) / grid_size);
      int gy = (int)((coords(i, 1) - min_y) / grid_size);
      int current_label = current_labels[i];
      
      // Count neighbor labels
      std::unordered_map<int, double> label_weights;
      label_weights[current_label] = 1.0; // Self weight
      
      // Check neighboring grid cells
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          int ngx = gx + dx;
          int ngy = gy + dy;
          
          if (ngx < 0 || ngx >= grid_x || ngy < 0 || ngy >= grid_y) continue;
          
          int cell_idx = ngy * grid_x + ngx;
          
          for (int j : grid[cell_idx]) {
            if (i == j) continue;
            
            double dx_val = coords(i, 0) - coords(j, 0);
            double dy_val = coords(i, 1) - coords(j, 1);
            double dz_val = coords(i, 2) - coords(j, 2);
            double dist_sq = dx_val*dx_val + dy_val*dy_val;
            
            if (dist_sq < max_dist*max_dist) {
              double dist = std::sqrt(dist_sq);
              double weight = std::exp(-dist / spatial_sigma) * 
                std::exp(-dz_val*dz_val / (2.0 * 0.5 * 0.5));
              
              label_weights[current_labels[j]] += weight;
            }
          }
        }
      }
      
      // Find best label
      int best_label = current_label;
      double best_weight = 0.0;
      
      for (const auto& pair : label_weights) {
        if (pair.second > best_weight) {
          best_weight = pair.second;
          best_label = pair.first;
        }
      }
      
      if (best_label != current_label) {
        new_labels[i] = best_label;
        changes++;
      }
    }
    
    current_labels = new_labels;
    
    Rcout << "  Changed " << changes << " labels\n";
    
    if (changes < boundary_count * 0.01) { // Less than 1% of boundaries changed
      Rcout << "  Convergence reached\n";
      break;
    }
  }
  
  return current_labels;
}

// Bayesian Boundary Refinement (BBR) with smart decisions
// [[Rcpp::export]]
List BBR_ultrafast_cpp(const arma::mat& coords,
                       const IntegerVector& labels,
                       double prior_strength = 1.0,
                       double likelihood_strength = 2.0,
                       double confidence_threshold = 1.3,
                       int cores = 4) {
  
  int n = coords.n_rows;
  IntegerVector refined = clone(labels);
  
  Rcout << "ULTRA-FAST Bayesian Boundary Refinement: " << n << " points\n";
  
  // 1. SUPER FAST GRID BUILDING
  const double grid_size = 4.0;
  const double boundary_dist = 2.0;
  const double boundary_dist_sq = boundary_dist * boundary_dist;
  const double neighbor_dist = 4.0;
  const double neighbor_dist_sq = neighbor_dist * neighbor_dist;
  
  double min_x = arma::min(coords.col(0));
  double max_x = arma::max(coords.col(0));
  double min_y = arma::min(coords.col(1));
  double max_y = arma::max(coords.col(1));
  
  int grid_x = std::max(1, (int)std::ceil((max_x - min_x) / grid_size));
  int grid_y = std::max(1, (int)std::ceil((max_y - min_y) / grid_size));
  int grid_cells = grid_x * grid_y;
  
  // Use arrays for speed
  std::vector<std::vector<int>> grid(grid_cells);
  std::vector<std::vector<float>> grid_z(grid_cells);
  
  // Pre-allocate
#pragma omp parallel for num_threads(cores)
  for (int i = 0; i < grid_cells; i++) {
    grid[i].reserve(128);
    grid_z[i].reserve(128);
  }
  
  // Fill grid
  for (int i = 0; i < n; i++) {
    if (labels[i] == 0) continue;
    
    int gx = (int)((coords(i, 0) - min_x) / grid_size);
    int gy = (int)((coords(i, 1) - min_y) / grid_size);
    gx = std::max(0, std::min(gx, grid_x - 1));
    gy = std::max(0, std::min(gy, grid_y - 1));
    
    int cell_idx = gy * grid_x + gx;
    grid[cell_idx].push_back(i);
    grid_z[cell_idx].push_back(coords(i, 2));
  }
  
  // 2. IDENTIFY BOUNDARY POINTS ONLY (FAST PARALLEL)
  std::vector<char> is_boundary(n, 0); // Use char for cache efficiency
  std::vector<int> boundary_indices;
  boundary_indices.reserve(n / 2); // Reserve ~50%
  
#pragma omp parallel num_threads(cores)
{
  std::vector<int> local_boundary;
  
#pragma omp for nowait schedule(static, 256)
  for (int i = 0; i < n; i++) {
    if (labels[i] == 0) continue;
    
    int current_label = labels[i];
    double xi = coords(i, 0);
    double yi = coords(i, 1);
    
    int gx = (int)((xi - min_x) / grid_size);
    int gy = (int)((yi - min_y) / grid_size);
    
    bool boundary_found = false;
    
    // Check only immediate neighbors (4 cells) for speed
    for (int dx = -1; dx <= 1 && !boundary_found; dx++) {
      for (int dy = -1; dy <= 1 && !boundary_found; dy++) {
        int ngx = gx + dx;
        int ngy = gy + dy;
        
        if (ngx < 0 || ngx >= grid_x || ngy < 0 || ngy >= grid_y) continue;
        
        int cell_idx = ngy * grid_x + ngx;
        const auto& cell_points = grid[cell_idx];
        
        for (int j : cell_points) {
          if (i == j) continue;
          
          double dx_val = xi - coords(j, 0);
          double dy_val = yi - coords(j, 1);
          double dist_sq = dx_val*dx_val + dy_val*dy_val;
          
          if (dist_sq < boundary_dist_sq && labels[j] != current_label) {
            boundary_found = true;
            local_boundary.push_back(i);
            break;
          }
        }
      }
    }
  }
  
#pragma omp critical
{
  for (int idx : local_boundary) {
    is_boundary[idx] = 1;
    boundary_indices.push_back(idx);
  }
}
}

int boundary_count = boundary_indices.size();
Rcout << "Boundary points: " << boundary_count << " (" 
      << 100.0 * boundary_count / n << "%)\n";

// 3. PRE-COMPUTE REGION STATISTICS (Z distribution per label)
std::unordered_map<int, double> label_sum_z;
std::unordered_map<int, double> label_sum_z2;
std::unordered_map<int, int> label_count;

for (int i = 0; i < n; i++) {
  int label = labels[i];
  if (label == 0) continue;
  
  double z = coords(i, 2);
  label_sum_z[label] += z;
  label_sum_z2[label] += z * z;
  label_count[label]++;
}

// Compute means and variances
std::unordered_map<int, double> label_mean_z;
std::unordered_map<int, double> label_var_z;

for (const auto& pair : label_count) {
  int label = pair.first;
  double count = pair.second;
  double mean_z = label_sum_z[label] / count;
  double var_z = std::max(0.1, (label_sum_z2[label] / count) - mean_z * mean_z);
  
  label_mean_z[label] = mean_z;
  label_var_z[label] = var_z;
}

// 4. ULTRA-FAST BAYESIAN REFINEMENT (ONLY BOUNDARIES)
std::vector<double> uncertainties(n, 0.0);
int changes = 0;

#pragma omp parallel num_threads(cores) reduction(+:changes)
{
  // Thread-local random number generator for tie-breaking
  unsigned int seed = omp_get_thread_num() + 1;
  
#pragma omp for schedule(dynamic, 64)
  for (int idx = 0; idx < boundary_count; idx++) {
    int i = boundary_indices[idx];
    
    int current_label = labels[i];
    double current_z = coords(i, 2);
    double xi = coords(i, 0);
    double yi = coords(i, 1);
    
    int gx = (int)((xi - min_x) / grid_size);
    int gy = (int)((yi - min_y) / grid_size);
    
    // SMARTER NEIGHBOR COLLECTION: Only collect unique labels with weights
    std::unordered_map<int, double> label_weights;
    label_weights[current_label] = 1.0; // Self-prior
    double total_weight = 1.0;
    
    // Check 3x3 grid neighborhood
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        int ngx = gx + dx;
        int ngy = gy + dy;
        
        if (ngx < 0 || ngx >= grid_x || ngy < 0 || ngy >= grid_y) continue;
        
        int cell_idx = ngy * grid_x + ngx;
        const auto& cell_points = grid[cell_idx];
        const auto& cell_zs = grid_z[cell_idx];
        
        for (size_t j_idx = 0; j_idx < cell_points.size(); j_idx++) {
          int j = cell_points[j_idx];
          if (i == j) continue;
          
          double dx_val = xi - coords(j, 0);
          double dy_val = yi - coords(j, 1);
          double dist_sq = dx_val*dx_val + dy_val*dy_val;
          
          if (dist_sq < neighbor_dist_sq) {
            int neighbor_label = labels[j];
            double weight = std::exp(-std::sqrt(dist_sq) / 2.0) * prior_strength;
            
            label_weights[neighbor_label] += weight;
            total_weight += weight;
          }
        }
      }
    }
    
    // BAYESIAN INFERENCE WITH PRE-COMPUTED STATISTICS
    double max_posterior = -1.0;
    int best_label = current_label;
    
    for (const auto& pair : label_weights) {
      int label = pair.first;
      double prior = pair.second / total_weight;
      
      // Use pre-computed region statistics
      double mean_z = label_mean_z[label];
      double var_z = label_var_z[label];
      
      // Fast Gaussian likelihood
      double z_diff = current_z - mean_z;
      double z_score_sq = z_diff * z_diff / (2.0 * var_z);
      double likelihood = std::exp(-z_score_sq) / std::sqrt(var_z * 6.283185);
      
      // Posterior = Prior * Likelihood^strength
      double posterior = prior * std::pow(likelihood, likelihood_strength);
      
      if (posterior > max_posterior) {
        max_posterior = posterior;
        best_label = label;
      } else if (std::abs(posterior - max_posterior) < 1e-6) {
        // Tie-breaking: prefer current label or smaller label
        if (label == current_label) {
          best_label = current_label;
        } else if (label < best_label) {
          // For deterministic tie-breaking
          best_label = label;
        }
      }
    }
    
    // CONFIDENCE-BASED DECISION
    double current_posterior = label_weights[current_label] / total_weight *
      std::pow(std::exp(-(current_z - label_mean_z[current_label])*
      (current_z - label_mean_z[current_label]) / 
      (2.0 * label_var_z[current_label])) / 
      std::sqrt(label_var_z[current_label] * 6.283185), 
      likelihood_strength);
    
    // Only change if MUCH better (reduces noise)
    if (best_label != current_label && max_posterior > current_posterior * confidence_threshold) {
      refined[i] = best_label;
      changes++;
    }
    
    // Fast uncertainty computation
    double entropy = 0.0;
    for (const auto& pair : label_weights) {
      double p = pair.second / total_weight;
      if (p > 1e-6) {
        entropy -= p * std::log(p);
      }
    }
    uncertainties[i] = entropy;
  }
}

Rcout << "Changed " << changes << " labels (" 
      << 100.0 * changes / boundary_count << "% of boundaries)\n";

// 5. POST-PROCESSING: Remove tiny fragments (optional but fast)
if (changes > 0) {
  // Count new label sizes
  std::unordered_map<int, int> new_label_counts;
  for (int i = 0; i < n; i++) {
    if (refined[i] > 0) {
      new_label_counts[refined[i]]++;
    }
  }
  
  // Merge very small regions (< 30 points) back to original
  int fragment_fixes = 0;
#pragma omp parallel for reduction(+:fragment_fixes)
  for (int i = 0; i < n; i++) {
    if (refined[i] > 0 && new_label_counts[refined[i]] < 30) {
      refined[i] = labels[i]; // Revert to original
      fragment_fixes++;
    }
  }
  
  if (fragment_fixes > 0) {
    Rcout << "Fixed " << fragment_fixes << " tiny fragments\n";
    changes -= fragment_fixes;
  }
}

return List::create(
  _["labels"] = refined,
  _["uncertainties"] = uncertainties,
  _["changes"] = changes,
  _["boundary_points"] = boundary_count,
  _["change_ratio"] = (double)changes / boundary_count
);
}

