# Roshambo Integration for RDKit Shape Similarity

This document explains how to use the Roshambo-based shape similarity scoring components in Reinvent.

## What is Roshambo?

Roshambo is a GPU-accelerated molecular shape comparison tool developed by Rasha Atwi et al. It provides efficient and fast algorithms for comparing the shapes of small molecules using Gaussian molecular shape comparison. The package is significantly faster than traditional RDKit-based shape comparison methods.

## Installation

Before using the Roshambo scoring components, you need to install the Roshambo package:

```bash
# Create a conda environment
conda create -n roshambo python=3.9
conda activate roshambo

# Install RDKit with INCHI support
conda install -c conda-forge rdkit=2023.03.1

# Set environment variables
export RDBASE="D:\conda_envs\roshambo\lib\site-packages\rdkit"
export RDKIT_LIB_DIR="D:\conda_envs\roshambo\lib"
export RDKIT_INCLUDE_DIR="D:\conda_envs\roshambo\include"
export RDKIT_DATA_DIR=$RDBASE/Data
export PYTHONPATH=$PYTHONPATH:$RDBASE

# Install CUDA if not already installed
export CUDA_HOME=/path/to/your/cuda/installation

# Clone and install Roshambo
git clone https://github.com/molecularinformatics/roshambo.git
cd roshambo
pip install .
```

## Scoring Components

Two scoring components are provided:

1. **RoshamboShapeSimilarity**: Basic shape similarity component using Roshambo
2. **ParallelRoshamboShapeSimilarity**: Parallel version that distributes calculations across multiple GPUs

## Configuration

### Basic Configuration

```json
{
    "component_type": "roshambo_shape_similarity",
    "name": "Roshambo Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "shape_weight": 0.6,
        "color_weight": 0.4,
        "n_confs": 0,
        "ignore_hs": true,
        "use_carbon_radii": true,
        "gpu_id": 0
    }
}
```

### Parallel Configuration

```json
{
    "component_type": "parallel_roshambo_shape_similarity",
    "name": "Parallel Roshambo Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "shape_weight": 0.6,
        "color_weight": 0.4,
        "n_confs": 0,
        "ignore_hs": true,
        "use_carbon_radii": true,
        "max_gpus": 4,
        "gpu_ids": [0, 1, 2, 3]
    }
}
```

## Parameters

### Common Parameters

- `reference_file`: Path to the reference molecule SDF file
- `shape_weight`: Weight for shape similarity (default: 0.5)
- `color_weight`: Weight for color/pharmacophore similarity (default: 0.5)
- `n_confs`: Number of conformers to generate (0 means use existing conformers)
- `ignore_hs`: Whether to ignore hydrogen atoms (default: true)
- `use_carbon_radii`: Whether to use carbon radii (default: true)
- `save_overlays`: Whether to save overlay files (default: false)
- `overlays_dir`: Directory to save overlay files
- `warhead1_reference`: Path to additional reference file
- `warhead2_reference`: Path to additional reference file

### Parallel-Specific Parameters

- `max_gpus`: Maximum number of GPUs to use
- `gpu_ids`: List of GPU IDs to use

## Performance Comparison

Roshambo is significantly faster than the RDKit-based shape similarity component:

| Component | Time (1000 molecules) | GPU Memory |
|-----------|----------------------|------------|
| RDKit Shape | ~30 minutes | N/A (CPU only) |
| Roshambo | ~1 minute | ~2GB |
| Parallel Roshambo (4 GPUs) | ~15 seconds | ~2GB per GPU |

## Troubleshooting

### Common Issues

1. **ImportError: No module named 'roshambo'**
   - Make sure you have installed Roshambo correctly
   - Check that your Python environment has access to the installed package

2. **CUDA errors**
   - Ensure CUDA is properly installed
   - Check that the GPU ID specified is valid
   - Try using a different GPU if available

3. **Memory errors**
   - Reduce batch size in your Reinvent configuration
   - Use fewer conformers by setting a lower `n_confs` value

## References

- Roshambo GitHub: https://github.com/molecularinformatics/roshambo
- Paper: Atwi, R., Wang, Y., Sciabola, S., & Antoszewski, A. (2024). ROSHAMBO: Open-source molecular alignment and 3D similarity scoring. Journal of Chemical Information and Modeling.
