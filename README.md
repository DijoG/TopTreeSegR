# TopTreeSegR - Ultra-Fast Topological Tree Segmentation

![R](https://img.shields.io/badge/R-%2764?style=for-the-badge&logo=r&logoColor=grey)
![LiDAR](https://img.shields.io/badge/LiDAR-Tree%20Segmentation-green?style=for-the-badge)
![Topology](https://img.shields.io/badge/Topology-Discrete%20Morse%20Theory-blue?style=for-the-badge)
![C++](https://img.shields.io/badge/C++-RcppArmadillo-blue?style=for-the-badge&logo=cplusplus)
![Development](https://img.shields.io/badge/Development-Active-brightgreen?style=for-the-badge)

**Blazing-fast individual tree segmentation from terrestrial LiDAR point clouds using Discrete Morse Theory and RcppArmadillo.**

TopTreeSegR produces clean, structural point clouds representing tree architecture using alpha-complex simplification. This removes redundant points while perfectly preserving topological structure, enabling fast, accurate ecological analysis without foliage noise.


## 🎯 Quick Start

### Installation
```r
# Install from GitHub
remotes::install_github("DijoG/DiscreteMorseR")
remotes::install_github("stla/AlphaHull3D") 
remotes::install_github("DijoG/TopTreeSegR")

# Load the package
library(TopTreeSegR)
```

### One-Line Segmentation
```r
library(lidR)
trees <- lidR::readLAS(trees.las)

# One-line segmentation on structural points
result <- TopTreeSegR::TTS_segmentation(trees)

# # Boom! You've got individual trees 🎉
print(result)  

# Visualize the topological segmentation
TopTreeSegR::plot_TTS_3d(result)

# Access individual trees
trees_list <- split(result$points, result$segmentation)
```
## ⚡ Performance

```r
TopTreeSegR::benchmark_TTS(result)
```
```text
=== BENCHMARK RESULTS ===

Total time: 164.92 seconds

Points processed: 1281618

Trees detected: 5

Minima found: 1055

Seeds detected: 5

Processing rate: 7771.15 points/second
```

## 🛠️ Advanced Usage

### Custom Parameters

```r
# Fine-tune for your forest type
tts <- TopTreeSegR::TTS_segmentation(
  lasdf = trees,
  alpha = 0.1,          # Alpha-complex parameter
  clip_height = 1.5,    # Stem detection height
  max_distance = 2.0,   # Minima connectivity
  cores = 12            # Parallel processing
)
```

### Validation & Quality Control
```r
# Check segmentation quality
validation <- TopTreeSegR::validate_segmentation(tts)

# Extract individual trees
tree_points <- split(tts$points, tts$segmentation)
tree_1 <- tree_points[[1]]  # First tree point cloud
```

### Results
```r
# 2D projections
TopTreeSegR::plot_TTS_2d(result, projection = "XY")  # Top-down view
TopTreeSegR::plot_TTS_2d(result, projection = "XZ")  # Side view

# 3D interactive plot
TopTreeSegR::plot_TTS_3d(tts)
```
| Input Point Cloud | Segmented Trees |
|:---:|:---:|
| <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR01.png" width="400"> | <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR03.png" width="400"> |
| <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR02.png" width="400"> | <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR04.png" width="400"> |

## 🔬 How It Works

### The Magic of Discrete Morse Theory

`TopTreeSegR` uses mathematical topology to understand tree structure:

```text
1. 🔄 Alpha-Complex Construction - Convert points to topo logical structure
2. 🎯 Morse Complex Computation - Find critical points (minima, saddles)
3. 🌳 Gradient Flow Analysis - Understand how points connect
4. 🔍 Seed Detection - Automatically find tree bottoms
5. 🚀 Label Propagation - Grow trees from seeds using gradient flow
```

## 🏗️ Architecture

```text
Raw Points → Alpha-Complex → Morse Theory → Tree Segmentation
    ↓              ↓              ↓              ↓
 1.2M pts     134K pts      Critical       Clean tree
                          points only     structure
```

## Key Features

```text
⚡ Ultra-Fast: RcppArmadillo backend with parallel processing
🎯 Topology-Preserving: Mathematical foundation ensures accuracy
🌳 Structure-Aware: Focuses on trunks and branches, ignores foliage noise
🚀 Automated: Zero parameter tuning required
📊 Validation-Ready: Built-in quality assessment metrics
```