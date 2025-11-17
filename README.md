# TopTreeSegR - Ultra-Fast Topological Tree Segmentation

![R](https://img.shields.io/badge/-%2764?style=for-the-badge&logo=r&logoColor=grey)
![LiDAR](https://img.shields.io/badge/LiDAR-green?style=for-the-badge)
![Topology](https://img.shields.io/badge/Topology-Discrete%20Morse%20Theory-blue?style=for-the-badge)
![C++](https://img.shields.io/badge/C++-RcppArmadillo-blue?style=for-the-badge&logo=cplusplus)
![Development](https://img.shields.io/badge/Development-brightgreen?style=for-the-badge)

**Blazing-fast individual tree segmentation from terrestrial LiDAR point clouds using Discrete Morse Theory and Gradient Flow Watershed.**

`TopTreeSegR` produces clean, structural point clouds representing tree architecture using alpha-complex simplification and gradient flow analysis. This revolutionary approach treats tree segmentation as **reverse hydrological modeling** - where labels flow upstream from trunk minima to capture complete tree watersheds naturally.


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
libary(dplyr) # Required by DiscreteMorseR
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
TopTreeSegR::benchmark_TTS(trees)
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
  alpha = .1,           # Alpha-complex parameter
  clip_height = .5,     # Max stem detection height
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
| <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR_01.png" width="400"> | <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR_03.png" width="400"> |
| <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR_02.png" width="400"> | <img src="https://raw.githubusercontent.com/DijoG/storage/main/TTSR/TTSR_04.png" width="400"> |

## 🔬 How It Works

### The Magic of Gradient Flow Watersheds

`TopTreeSegR` uses mathematical topology to to understand tree structure through **natural drainage patterns**:

```text
🔄 Alpha-Complex Construction - Convert points to topological mesh
🎯 Morse Complex Computation - Find critical points (tree trunks = minima)
🌊 Gradient Flow Analysis - Map natural "downhill flow" on tree surfaces 
🔍 Seed Detection - Automatically find tree bottoms (drainage points)
🚀 Reverse Upstream Propagation - Start at trunks, capture ALL points that flow to them
🌳️ Complete Watershed Capture - Each tree = everything connected to its trunk
```
## 🏗️ Architecture

```text
Raw Points → Alpha-Complex → Morse Theory → Gradient Flow → Tree Watersheds
    ↓              ↓              ↓               ↓               ↓
 1.2M pts      134K pts       Critical       Natural flow    Complete tree
                               points          analysis       structures
```
## Key Features

```text
⚡ Ultra-Fast: RcppArmadillo backend with parallel processing
🌊 Watershed-Based: Revolutionary gradient flow approach
🎯 Topology-Preserving: Mathematical foundation ensures accuracy
🌳 Structure-Aware: Focuses on trunks and branches, ignores foliage noise
🚀 Automated: Zero parameter tuning required
📊 Validation-Ready: Built-in quality assessment metrics
🔄 Natural Boundaries: Trees segmented as connected watersheds, not clusters
🔗 Perfect Connectivity: Every tree guaranteed to be one connected component
```
## Mathematical Foundation

- Discrete Morse Theory: Combinatorial alternative to smooth Morse theory
- Alpha Complexes: Subcomplex of Delaunay triangulation for topological simplification
- Forman Gradient: Discrete representation of gradient flow
- Ascending Regions: Influence regions of minima via gradient paths