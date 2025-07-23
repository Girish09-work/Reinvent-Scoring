#!/bin/bash

# =============================================================================
# Roshambo Environment Setup Script for DGX
# =============================================================================
# This script sets up the environment variables for Roshambo with RDKit on DGX
# Usage: source roshambo.sh
# =============================================================================

echo "Setting up Roshambo environment for DGX..."

# =============================================================================
# Conda Environment Configuration for DGX
# =============================================================================
echo "Activating conda environment..."

# Activate the roshambo conda environment (DGX path)
conda activate /home/protacinvent/.conda/envs/roshambo

# Verify conda environment activation
if [ "$CONDA_DEFAULT_ENV" != "roshambo" ] && [ "$CONDA_DEFAULT_ENV" != "/home/protacinvent/.conda/envs/roshambo" ]; then
    echo "Warning: Conda environment may not be properly activated"
    echo "Current environment: $CONDA_DEFAULT_ENV"
fi

# =============================================================================
# RDKit Configuration for DGX
# =============================================================================
echo "Setting up RDKit environment variables..."

# Navigate to RDKit directory and set RDBASE
cd /home/protacinvent/Desktop/roshambo/rdkit
export RDBASE=$(pwd)
export RDKIT_LIB_DIR="$RDBASE/lib"
export RDKIT_INCLUDE_DIR="$RDBASE/Code"
export RDKIT_DATA_DIR="$RDBASE/Data"

# =============================================================================
# Python Path Configuration (DGX specific)
# =============================================================================
# Set PYTHONPATH to include RDKit
export PYTHONPATH="$PYTHONPATH:$RDBASE"

# =============================================================================
# Library Path Configuration for DGX (Fix Boost library issue)
# =============================================================================
# Set library paths for DGX environment with Boost fix
export LD_LIBRARY_PATH="$RDBASE/lib:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# Fix for Boost library version mismatch
echo "Checking and fixing Boost library paths..."
if [ -d "$CONDA_PREFIX/lib" ]; then
    # Find available boost_python libraries
    BOOST_PYTHON_LIBS=$(find "$CONDA_PREFIX/lib" -name "libboost_python*.so*" 2>/dev/null | head -5)
    if [ -n "$BOOST_PYTHON_LIBS" ]; then
        echo "Available Boost Python libraries:"
        echo "$BOOST_PYTHON_LIBS"

        # Create symlinks for missing boost libraries if needed
        cd "$CONDA_PREFIX/lib"

        # Check if the required boost_python39.so.1.85.0 exists
        if [ ! -f "libboost_python39.so.1.85.0" ]; then
            # Find the actual boost_python library
            ACTUAL_BOOST=$(find . -name "libboost_python*.so*" | grep -E "libboost_python[0-9]*\.so" | head -1)
            if [ -n "$ACTUAL_BOOST" ]; then
                echo "Creating symlink: libboost_python39.so.1.85.0 -> $ACTUAL_BOOST"
                ln -sf "$ACTUAL_BOOST" libboost_python39.so.1.85.0 2>/dev/null || echo "Note: Could not create symlink (may need sudo)"
            fi
        fi

        cd - > /dev/null
    else
        echo "Warning: No Boost Python libraries found in conda environment"
        echo "You may need to reinstall boost: conda install -c conda-forge boost-cpp"
    fi
fi

# =============================================================================
# CUDA Configuration for WSL
# =============================================================================
export CUDA_HOME="/usr/local/cuda"
export CUDA_ROOT="/usr/local/cuda"
export PATH="$CUDA_HOME/bin:$PATH"

# =============================================================================
# Roshambo Specific Configuration
# =============================================================================
export ROSHAMBO_ROOT="/mnt/d/roshambo"
export ROSHAMBO_DATA_DIR="$ROSHAMBO_ROOT/data"

# =============================================================================
# GPU Configuration
# =============================================================================
export CUDA_VISIBLE_DEVICES="0"  # Use GPU 0 by default
export NVIDIA_VISIBLE_DEVICES="0"

# Navigate back to the main roshambo directory
cd /mnt/d/roshambo

# =============================================================================
# Environment Verification
# =============================================================================
echo ""
echo "Environment variables set:"
echo "  CONDA_PREFIX: $CONDA_PREFIX"
echo "  CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"
echo "  RDBASE: $RDBASE"
echo "  PYTHONPATH: $PYTHONPATH"
echo "  LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "  CUDA_HOME: $CUDA_HOME"
echo "  ROSHAMBO_ROOT: $ROSHAMBO_ROOT"
echo "  Current directory: $(pwd)"

# =============================================================================
# Package Import Testing
# =============================================================================
echo ""
echo "Testing RDKit import..."
python -c "
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    print('✓ RDKit imported successfully')
    print('  RDKit version:', Chem.rdBase.rdkitVersion)
except ImportError as e:
    print('✗ RDKit import failed:', e)
    print('  Check RDBASE and PYTHONPATH settings')
"

echo ""
echo "Testing Roshambo import..."
python -c "
try:
    import roshambo
    from roshambo.api import get_similarity_scores
    print('✓ Roshambo imported successfully')
except ImportError as e:
    print('✗ Roshambo import failed:', e)
    print('  Make sure Roshambo is installed in the conda environment')
"

# =============================================================================
# CUDA Testing
# =============================================================================
echo ""
echo "Testing CUDA availability..."
python -c "
try:
    import torch
    if torch.cuda.is_available():
        print('✓ CUDA available')
        print('  CUDA version:', torch.version.cuda)
        print('  GPU count:', torch.cuda.device_count())
        for i in range(torch.cuda.device_count()):
            print(f'  GPU {i}: {torch.cuda.get_device_name(i)}')
    else:
        print('⚠ CUDA not available')
except ImportError:
    print('⚠ PyTorch not installed, cannot check CUDA')
except Exception as e:
    print('⚠ CUDA check failed:', e)
"

# =============================================================================
# Boost Library Check (WSL specific issue from memory)
# =============================================================================
echo ""
echo "Checking Boost library configuration..."
python -c "
import os
conda_prefix = os.environ.get('CONDA_PREFIX', '')
if conda_prefix:
    boost_lib = os.path.join(conda_prefix, 'lib')
    boost_include = os.path.join(conda_prefix, 'include')
    print(f'Conda Boost lib path: {boost_lib}')
    print(f'Conda Boost include path: {boost_include}')
    if os.path.exists(boost_lib):
        print('✓ Boost library directory found in conda environment')
    else:
        print('⚠ Boost library directory not found in conda environment')
else:
    print('⚠ CONDA_PREFIX not set')
"

echo ""
echo "==============================================================================="
echo "Roshambo WSL environment setup complete!"
echo "==============================================================================="
echo ""
echo "Usage instructions:"
echo "  1. To activate this environment: source roshambo.sh"
echo "  2. To run Roshambo API tests: python simple_roshambo_test.py"
echo "  3. To run full test suite: bash run_roshambo_test.sh"
echo ""
echo "If you encounter Boost version issues:"
echo "  conda install -c conda-forge boost=1.81.0"
echo ""
echo "Environment is ready for Roshambo operations!"
