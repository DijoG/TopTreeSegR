# TopTreeSegR - Ultra-Fast Topological Tree Segmentation

![R](https://img.shields.io/badge/-%2764?style=for-the-badge&logo=r&logoColor=grey)
![LiDAR](https://img.shields.io/badge/LiDAR-green?style=for-the-badge)
![Topology](https://img.shields.io/badge/Topology-Discrete%20Morse%20Theory-blue?style=for-the-badge)
![C++](https://img.shields.io/badge/C++-RcppArmadillo-blue?style=for-the-badge&logo=cplusplus)
![Development](https://img.shields.io/badge/Development-brightgreen?style=for-the-badge)

**Blazing-fast individual tree segmentation from terrestrial LiDAR point clouds using Discrete Morse Theory with Bayesian refinement achieving >0.85 Adjusted Rand Index (ARI).**

`TopTreeSegR` combines **Morse-Smale complex analysis** with **Bayesian boundary optimization** to achieve state-of-the-art tree segmentation accuracy. This R package delivers a production-ready pipeline that transforms raw LiDAR points into clean, validated tree segments in under two minutes for a typical plot.

> **⚠️ Important:** TopTreeSegR segments the **alpha-complex mesh** (Delaunay triangulation), **not** the original raw point cloud. The mesh captures the essential tree structure while reducing point count by ~90% for **ultra-fast processing**.

🎯 Core Features: >0.85 ARI | Single TTS_pipeline() | Bayesian refinement | Full validation suite | Interactive 3D viz

## 🏆 Key Achievement

**Proven >0.85 Adjusted Rand Index (ARI)** with the 2-pass Bayesian Boundary Refinement pipeline on benchmark TLS datasets.

## 🚀 Quick Start

### Installation
```r
# Important: ahull3D, DiscreteMorseR, and lidR are not on CRAN
# Install them from GitHub first:

remotes::install_github("DijoG/ahull3D")
remotes::install_github("DijoG/DiscreteMorseR")
remotes::install_github("r-lidar/lidR")

# Install TopTreeSegR
remotes::install_github("DijoG/TopTreeSegR")

# Additional dependencies from CRAN
install.packages(c("dbscan", "FNN", "plotly", "tictoc"))

library(TopTreeSegR)
```

### One-Line Pipeline (Recommended)
```r
# Read LAS file with ground truth
trees <- lidR::readLAS("your_forest.las")

# Add "pid" field to LAS attribute
pid <- 1:lidR::npoints(trees)  
trees_filtered <- lidR::add_lasattribute(trees, pid, "pid", "Unique point ID")

# Complete pipeline: segmentation + Bayesian refinement
result <- TopTreeSegR::TTS_pipeline(
  las = trees
  cores = 20                  
)

# Boom! You've got individual trees 🎉
print(result)  

# Visualize the topological segmentation (3D)
TopTreeSegR::plot_TTS_3d(result)

# Visualize the topological segmentation (2D)
TopTreeSegR::plot_TTS_2d(result, projection = "XZ")

# Validate 
TopTreeSegR::validate_TTS(result, trees)
```
## ⚡ Performance

```r
# Complete pipeline benchmark
tictoc::tic()
result <- TTS_pipeline(trees, cores = 8)
tictoc::toc()
# 118.68 seconds (~2 minutes)

validate_TTS(result, trees)
```
```text
=== TTS Segmentation Validation ===

Tree IDs automatically aligned
Points: 232420 | Pred trees: 5 | True trees: 5
Match rate: 100.0%

=== METRICS ===
Precision:  0.9051
Recall:     0.9051
F1-Score:   0.9051
Accuracy:   0.9051
Rand Index: 0.9416
Adj Rand I: 0.8403
```

## 🛠️ Advanced Usage

### Complete Pipeline with 2-Pass Bayesian Boundary Refinement (BBR)

```r
# TTS_pipeline includes the first BBR pass
res <- TTS_pipeline(
  las = trees,
  method = "morse-smale",
  input_truth = "pid"              # Las attribute of point ids
  alpha = 0.1,                     # Alpha value for alpha hull 
  stem_height = 0.5,               # Find seeds below 
  prior_strength = 1.0,            # Spatial consistency
  likelihood_strength = 1.6,       # Elevation consistency (key!)
  confidence_threshold = 1.0,      # 1) Aggressive refinement (conf=1.0)
  cores = 16,
  fix_fragments = TRUE
) # ~2 minutes

validate_TTS(res, trees)  # ARI: 0.8408

# Second BBR pass
tictoc::tic()
res2 <- TTS_BBR(
  res1,
  prior_strength = 1.0,
  likelihood_strength = 1.6,
  confidence_threshold = 1.9,      # 2) Conservative cleanup, conf=1.9 ~ only very (> 0.9) probable improvements
  cores = 16) 
tictoc::toc()# ~83 seconds

validate_TTS(res2, trees) # ARI: 0.8501
```

### Results
```r
# 2D projections
plot_TTS_2d(res2, projection = "XY")  # Top-down view
plot_TTS_2d(res2, projection = "XZ")  # Side view

# 3D interactive plot
plot_TTS_3d(res2)
```
| Input Point Cloud | Segmented Trees |
|:---:|:---:|
| <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSRms_01.png" width="400"> | <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSRms_03.png" width="400"> |
| <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSRms_02.png" width="400"> | <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSRms_04.png" width="400"> |

## 🏗️ Architecture

```text
RAW POINTS → ALPHA-COMPLEX → MORSE COMPLEX → SEGMENTATION → BAYESIAN REFINEMENT → OUTPUT
    ↓              ↓              ↓               ↓                 ↓               ↓
 1.2M pts      134K pts       Critical        Descending        Boundary         Individual
                              simplices       manifolds         optimization     trees
```
## 🔬 How It Works

`TopTreeSegR` uses Discrete Morse Theory to segment trees by analyzing the topological structure of the point cloud:

```text
🔺 1. Delaunay Triangulation
    - Convert discrete points to a continuous triangular mesh
    - Preserves topological connectivity of the forest structure

🔺 2. Alpha-Complex Construction
    - Simplify triangulation using alpha-shape filtration
    - Removes noisy edges and isolates tree structures
    - Output: **134K mesh vertices** (from 1.2M input points)

🎯 3. Morse Complex Computation
    - Compute discrete gradient vector field on the mesh
    - Identify critical points (minima, maxima, saddles)
    - Tree trunks correspond to Morse minima (lowest points in each tree)

🌊 4. Morse-Smale Complex Analysis
    - Partition mesh into ascending/descending manifolds
    - Each minimum's ascending manifold = influence region of a tree
    - Points flow downhill to their corresponding trunk minimum

🔍 5. Seed Detection & Initial Segmentation
    - Cluster low minima near ground as tree trunks
    - Assign each mesh vertex to the trunk whose ascending manifold it belongs to
    - Create initial tree segments based on gradient flow structure

🧠 6. Bayesian Boundary Refinement (Key Innovation!)
    - Identify boundary points between trees
    - Compute posterior probabilities using:
        • Prior: Neighborhood label consistency
        • Likelihood: Elevation coherence (Gaussian model)
    - Perform Maximum A Posteriori (MAP) estimation
    - Fix boundary errors using elevation consistency

✅ 7. Output
    - Tree labels assigned to mesh vertices
    - Bayesian refinement improves boundary accuracy
    - Clean individual tree segments (mesh vertices labeled)
```
## Key Features

```text
⚡ Ultra-Fast: RcppArmadillo + OpenMP parallel processing
🧠 Bayesian-Optimized: Elevation-consistent boundary refinement (>0.85 ARI)
🎯 Topology-Based: Morse-Smale complex analysis for robust segmentation
🌳 Structure-Aware: Leverages tree height coherence for accurate boundaries
📊 Validation-Ready: Built-in ARI, accuracs, precision, recall, F1-score metrics
🔍 Uncertainty-Aware: Bayesian framework quantifies segmentation confidence
🚀 Production-Ready: Single TTS_pipeline() function for end-to-end workflow
📈 Benchmark-Proven: Validated on TLS datasets with ground truth
```
## Mathematical Foundation

```text
- Discrete Morse Theory: Combinatorial framework for topological analysis of discrete data
- Morse-Smale Complex: Partition of space into ascending/descending manifolds from critical points
- Alpha Complex: Topologically correct subset of Delaunay triangulation for point cloud simplification
- Forman Gradient: Discrete vector field representing flow direction on the mesh
- Bayesian Inference: Posterior = Prior × Likelihood, with MAP estimation for optimal labeling
- Gaussian Likelihood: Models elevation consistency within each tree segment
- Markovian Prior: Captures spatial coherence through neighborhood relationships
```