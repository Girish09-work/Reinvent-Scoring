# Enhanced Roshambo with Dynamic Environment Configuration

This guide covers the enhanced Roshambo shape similarity components with dynamic environment configuration, enabling flexible deployment across different user setups without hardcoded paths.

## Overview

The enhanced Roshambo components provide GPU-accelerated molecular shape similarity calculations with automatic environment discovery and configuration. Key features include:

- **🔧 Dynamic Environment Configuration**: Automatically discover and configure conda environments and RDKit installations
- **📁 Flexible Path Management**: Support for custom RDBASE, conda, and CUDA paths
- **🧪 RDKit Conformer Generation**: Use RDKit instead of relying on Roshambo's internal conformer generation
- **🐍 Smart Environment Discovery**: Auto-discover conda installations and environments
- **📊 Enhanced Parameter Management**: Dedicated parameter enum for better organization
- **🧬 Multiple Reference Files**: Support for PROTAC design with multiple warhead references
- **🐛 Improved Error Handling**: Better debugging with detailed diagnostics and suggestions
- **✅ Environment Validation**: Comprehensive validation with helpful error messages

## Components

### 1. RoshamboShapeSimilarity (Enhanced)
Basic shape similarity component with dynamic environment configuration.

### 2. ParallelRoshamboShapeSimilarity (Enhanced)
Parallel version that distributes calculations across multiple GPUs.

### 3. RoshamboDynamicEnvironmentManager
New component that handles automatic discovery and configuration of:
- Conda installations and environments
- RDKit build directories (RDBASE)
- CUDA installations
- Environment variable setup

## Dynamic Configuration Parameters

### Core Environment Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `rdbase_path` | string | **YES** | None | Path to RDKit build directory |
| `conda_env_name` | string | **YES** | None | Name of the conda environment with Roshambo |
| `conda_base_path` | string | No | auto-discover | Custom path to conda installation |
| `cuda_home_path` | string | No | "/usr/local/cuda" | Path to CUDA installation |
| `auto_setup_env` | bool | No | true | Automatically set environment variables |

### Legacy Support
| Parameter | Type | Description |
|-----------|------|-------------|
| `environment_path` | string | Legacy environment script path (still supported) |

## Configuration Examples

### 1. Minimal Required Configuration
```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Roshambo Minimal Required",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "rdbase_path": "/home/user/Desktop/roshambo/rdkit",
    "conda_env_name": "roshambo",
    "debug": true
  }
}
```

### 2. Custom Paths Configuration
```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Roshambo Custom Paths",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "rdbase_path": "/custom/path/to/rdkit",
    "conda_env_name": "my_roshambo_env",
    "conda_base_path": "/home/user/miniconda3",
    "cuda_home_path": "/usr/local/cuda-12.0",
    "auto_setup_env": true,
    "debug": true
  }
}
```

### 3. PROTAC Design with Dynamic Environment
```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "PROTAC Roshambo",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "protac_reference.sdf",
    "warhead1_reference": "warhead1.sdf",
    "warhead2_reference": "warhead2.sdf",
    "rdbase_path": "/home/user/Desktop/roshambo/rdkit",
    "conda_env_name": "roshambo",
    "auto_setup_env": true,
    "n_confs": 20,
    "save_overlays": true,
    "debug": true
  }
}
```

### 4. Multi-GPU Parallel Configuration
```json
{
  "component_type": "parallel_roshambo_shape_similarity",
  "name": "Parallel Roshambo",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "reference.sdf",
    "rdbase_path": "/opt/rdkit",
    "conda_env_name": "roshambo",
    "max_gpus": 2,
    "gpu_ids": [0, 1],
    "auto_setup_env": true,
    "debug": true
  }
}
```

## Environment Discovery Process

The dynamic environment manager follows this discovery process:

### 1. Conda Discovery
- Searches for conda in PATH
- Checks common installation locations:
  - `/opt/conda/bin/conda`
  - `/home/*/miniconda3/bin/conda`
  - `/home/*/anaconda3/bin/conda`
- Uses `conda_base_path` if specified

### 2. Environment Discovery
- Lists conda environments
- Finds environment matching `conda_env_name`
- Validates environment exists and is accessible

### 3. RDBASE Discovery
- Uses `rdbase_path` if specified and valid
- Searches common RDKit build locations:
  - `/home/*/Desktop/roshambo/rdkit`
  - `/home/*/rdkit`
  - `/opt/rdkit`
  - `/usr/local/rdkit`
- Validates required components (lib, Code, Data directories)

### 4. Environment Variable Setup
Automatically sets:
```bash
export RDBASE=<discovered_or_specified_path>
export RDKIT_LIB_DIR=$RDBASE/lib
export RDKIT_INCLUDE_DIR=$RDBASE/Code
export RDKIT_DATA_DIR=$RDBASE/Data
export PYTHONPATH=$RDBASE:$CONDA_PREFIX/lib/python*/site-packages
export LD_LIBRARY_PATH=$RDBASE/lib:$CONDA_PREFIX/lib:$CUDA_HOME/lib64
export CUDA_HOME=<cuda_home_path>
```

## Error Handling and Validation

### Validation Checks
- ✅ Conda executable exists and is functional
- ✅ Specified conda environment exists
- ✅ RDBASE path contains required RDKit components
- ✅ CUDA installation is accessible
- ✅ Python packages (RDKit, Roshambo) can be imported

### Error Messages
The system provides detailed error messages with suggestions:

```
Environment setup failed:
  - Conda environment 'roshambo' not found
  - No valid RDBASE path found

Environment Summary:
  Conda executable: /opt/conda/bin/conda
  Conda environment: Not found
  RDBASE path: Not found

Suggestions:
  1. Ensure conda environment 'roshambo' exists
  2. Specify 'rdbase_path' in configuration if RDKit is in a custom location
  3. Set 'conda_base_path' if conda is installed in a non-standard location
```

## Migration from Legacy Configuration

### Old Configuration (Hardcoded Paths)
```json
{
  "environment_path": "source /home/user/setup_env.sh &&"
}
```

### New Configuration (Dynamic)
```json
{
  "rdbase_path": "/home/user/Desktop/roshambo/rdkit",
  "conda_env_name": "roshambo",
  "auto_setup_env": true,
  "debug": true
}
```

## Benefits of Dynamic Configuration

### 1. **Flexibility**
- Works with any user's RDKit/Roshambo installation
- No hardcoded paths in configuration
- Supports different conda installations

### 2. **Reliability**
- Validates environment before execution
- Clear error messages with suggestions
- Automatic path discovery reduces configuration errors

### 3. **Maintainability**
- Easier to deploy across different systems
- Reduced configuration complexity
- Better debugging capabilities

### 4. **User Experience**
- Minimal configuration required for standard setups
- Detailed feedback during setup
- Helpful suggestions for fixing issues

## Troubleshooting

### Common Issues and Solutions

#### Issue: "Conda environment not found"
**Solution**:
```bash
conda create -n roshambo python=3.8
conda activate roshambo
pip install rdkit roshambo
```

#### Issue: "No valid RDBASE path found"
**Solution**: Specify custom path in configuration:
```json
{
  "rdbase_path": "/path/to/your/rdkit/build"
}
```

#### Issue: "RDKit import failed"
**Solution**: Check RDKit installation in conda environment:
```bash
conda activate roshambo
python -c "from rdkit import Chem; print('RDKit OK')"
```

### Debug Mode
Enable debug mode for detailed diagnostics:
```json
{
  "debug": true
}
```

This provides:
- Environment discovery details
- Validation results
- Environment variable settings
- Import test results

## Performance Considerations

- **Auto-discovery overhead**: Minimal, runs once during initialization
- **Environment validation**: Fast, uses cached results
- **Dynamic script generation**: Only when needed for subprocess execution
- **Memory usage**: Negligible additional overhead

The dynamic configuration system adds minimal overhead while providing significant flexibility and reliability improvements.
