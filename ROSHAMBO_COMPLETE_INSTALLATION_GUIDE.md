# Complete Roshambo Installation Guide for DGX

## Overview
This guide provides step-by-step instructions for installing Roshambo with RDKit on DGX systems, including fixes for common Boost library issues.

## Prerequisites
- DGX system with CUDA support
- Conda/Miniconda installed
- Git installed
- CMake installed

## Step 1: Clone Roshambo Repository
```bash
cd /home/protacinvent/Desktop
git clone https://github.com/molecularinformatics/roshambo.git
cd roshambo
```

## Step 2: Create and Configure Conda Environment
```bash
# Create conda environment with Python 3.9
conda create -n roshambo python=3.9 -y

# Activate the environment
conda activate roshambo

# Install essential packages with specific versions to avoid conflicts
conda install -c conda-forge \
    boost-cpp=1.85.0 \
    boost=1.85.0 \
    cairo \
    pandas \
    pillow \
    freetype \
    cmake \
    numpy \
    eigen \
    matplotlib \
    pytest \
    -y

# Update matplotlib
pip install -U matplotlib
```

## Step 3: Clone and Build RDKit
```bash
# Clone RDKit repository
git clone https://github.com/rdkit/rdkit.git
cd rdkit

# Create build directory
mkdir build
cd build

# Configure CMake with proper Boost settings
cmake \
    -DPy_ENABLE_SHARED=1 \
    -DRDK_INSTALL_INTREE=ON \
    -DRDK_BUILD_AVALON_SUPPORT=ON \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DRDK_INSTALL_STATIC_LIBS=OFF \
    -DRDK_BUILD_CPP_TESTS=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=ON \
    -DRDK_BUILD_YAEHMOP_SUPPORT=ON \
    -DBoost_NO_SYSTEM_PATHS=ON \
    -DRDK_BUILD_CAIRO_SUPPORT=ON \
    -DRDK_BUILD_XYZ2MOL_SUPPORT=ON \
    -DRDK_BUILD_FREESASA_SUPPORT=ON \
    -DPYTHON_NUMPY_INCLUDE_PATH="$(python -c 'import numpy; print(numpy.get_include())')" \
    -DBOOST_ROOT="$CONDA_PREFIX" \
    -DBoost_INCLUDE_DIR="$CONDA_PREFIX/include" \
    -DBoost_LIBRARY_DIRS="$CONDA_PREFIX/lib" \
    -DRDK_INSTALL_COMIC_FONTS=OFF \
    -DINCHI_URL=https://github.com/IUPAC-InChI/InChI/releases/download/v1.07.3/INCHI-1-SRC.zip \
    ..

# Build RDKit (this may take 30-60 minutes)
make -j4 install

# Go back to rdkit root directory
cd ..

# Set environment variables
export RDBASE=$(pwd)
export PYTHONPATH=$RDBASE
export LD_LIBRARY_PATH=$RDBASE/lib:$CONDA_PREFIX/lib

# Run tests to verify installation
cd build
ctest -j4 --output-on-failure

# Go back to roshambo directory
cd ../../
```

## Step 4: Configure Environment Variables
```bash
# Set RDKit environment variables
export RDKIT_LIB_DIR=$RDBASE/lib
export RDKIT_INCLUDE_DIR=$RDBASE/Code
export RDKIT_DATA_DIR=$RDBASE/Data
export PYTHONPATH=$PYTHONPATH:$RDBASE

# Set CUDA environment (adjust path if needed)
export CUDA_HOME=/usr/local/cuda
```

## Step 5: Install Roshambo
```bash
# Install Roshambo in development mode
pip install -e .

# Install additional dependencies
pip install IPython pandas jupyter
```

## Step 6: Fix Boost Library Issues (If Needed)
If you encounter `libboost_python39.so.1.85.0: cannot open shared object file` error:

```bash
# Navigate to conda lib directory
cd $CONDA_PREFIX/lib

# Check available boost libraries
ls -la libboost_python*

# Create symlink if needed (replace X.Y.Z with your actual version)
ln -sf libboost_python39.so.X.Y.Z libboost_python39.so.1.85.0

# Or reinstall boost with specific version
conda install -c conda-forge boost-cpp=1.85.0 boost=1.85.0 --force-reinstall
```

## Step 7: Verification
Test your installation:

```bash
# Test RDKit import
python -c "from rdkit import Chem; print('RDKit version:', Chem.rdBase.rdkitVersion)"

# Test Roshambo import
python -c "import roshambo; from roshambo.api import get_similarity_scores; print('Roshambo imported successfully')"
```

## Common Issues and Solutions

### Issue 1: Boost Library Version Mismatch
**Error**: `libboost_python39.so.1.85.0: cannot open shared object file`

**Solution**:
1. Check available boost libraries: `ls $CONDA_PREFIX/lib/libboost_python*`
2. Create symlink or reinstall boost with matching version
3. Use the updated `roshambo.sh` script which handles this automatically

### Issue 2: CMake Configuration Fails
**Error**: CMake cannot find Boost

**Solution**:
1. Ensure `CONDA_PREFIX` is set: `echo $CONDA_PREFIX`
2. Reinstall boost: `conda install -c conda-forge boost-cpp boost`
3. Clear CMake cache: `rm -rf build/*` and reconfigure

### Issue 3: RDKit Tests Fail
**Error**: Some RDKit tests fail during `ctest`

**Solution**:
1. This is often normal - check if core functionality works
2. Test basic RDKit import in Python
3. If critical tests fail, rebuild with different options

## Environment Setup Script
Use the provided `roshambo.sh` script for easy environment setup:

```bash
# Make executable
chmod +x roshambo.sh

# Source to set up environment
source roshambo.sh
```

## Final Verification
Run a complete test:

```bash
# Run the test script
bash run_roshambo_test.sh
```

This should complete without errors and generate similarity scores in CSV format.

## Notes
- Build time for RDKit can be 30-60 minutes depending on system
- Boost version compatibility is critical - use version 1.85.0 for best results
- The environment setup script handles most common issues automatically
- Always activate the conda environment before using Roshambo
