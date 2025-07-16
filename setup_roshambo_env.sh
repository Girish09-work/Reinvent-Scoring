#!/bin/bash

# =============================================================================
# Roshambo Environment Setup Script
# =============================================================================
# This script sets up the environment variables for Roshambo with RDKit
# Usage: source setup_roshambo_env.sh
# =============================================================================

echo "Setting up Roshambo environment..."

# =============================================================================
# Conda Environment Configuration
# =============================================================================
export CONDA_PREFIX="/home/protacinvent/.conda/envs/roshambo"

# Activate the roshambo conda environment
echo "Activating roshambo conda environment..."

# Initialize conda properly
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source /opt/conda/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/miniconda3/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/miniconda3/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/.conda/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/.conda/etc/profile.d/conda.sh
else
    echo "Warning: Could not find conda.sh, trying direct activation"
fi

conda activate roshambo

# =============================================================================
# RDKit Configuration
# =============================================================================
export RDBASE="/home/protacinvent/Desktop/roshambo/rdkit"
export RDKIT_LIB_DIR="$RDBASE/lib"
export RDKIT_INCLUDE_DIR="$RDBASE/Code"
export RDKIT_DATA_DIR="$RDBASE/Data"

# =============================================================================
# Python Path Configuration (Based on Working Setup)
# =============================================================================
# Set PYTHONPATH exactly as in your working configuration
export PYTHONPATH="$RDBASE"
export PYTHONPATH="$PYTHONPATH:$RDBASE"

# =============================================================================
# Library Path Configuration
# =============================================================================
# Update CONDA_PREFIX to actual activated environment if needed
if [ -z "$CONDA_PREFIX" ]; then
    export CONDA_PREFIX="$(conda info --base)/envs/roshambo"
fi

export LD_LIBRARY_PATH="$RDBASE/lib:$CONDA_PREFIX/lib"

# =============================================================================
# CUDA Configuration
# =============================================================================
export CUDA_HOME="/usr/local/cuda"
export CUDA_ROOT="/usr/local/cuda"
export PATH="$CUDA_HOME/bin:$PATH"

# =============================================================================
# Roshambo Specific Configuration
# =============================================================================
export ROSHAMBO_ROOT="/home/protacinvent/Desktop/roshambo"
export ROSHAMBO_DATA_DIR="$ROSHAMBO_ROOT/data"

# =============================================================================
# GPU Configuration
# =============================================================================
export CUDA_VISIBLE_DEVICES="0"  # Use GPU 0 by default
export NVIDIA_VISIBLE_DEVICES="0"

# =============================================================================
# Verification and Testing
# =============================================================================
echo "Environment variables set:"
echo "  CONDA_PREFIX: $CONDA_PREFIX"
echo "  RDBASE: $RDBASE"
echo "  PYTHONPATH: $PYTHONPATH"
echo "  LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "  CUDA_HOME: $CUDA_HOME"
echo "  ROSHAMBO_ROOT: $ROSHAMBO_ROOT"

# Test RDKit import
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
"

# Test Roshambo import
echo ""
echo "Testing Roshambo import..."
python -c "
try:
    import roshambo
    from roshambo.api import get_similarity_scores
    print('✓ Roshambo imported successfully')
except ImportError as e:
    print('✗ Roshambo import failed:', e)
    print('  Make sure Roshambo is installed: pip install git+https://github.com/molecularinformatics/roshambo.git')
"

# Test CUDA availability
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
"

# =============================================================================
# Build and Test RDKit (if needed)
# =============================================================================
build_rdkit() {
    echo ""
    echo "Building and testing RDKit..."
    cd "$RDBASE/build" || {
        echo "Error: RDKit build directory not found at $RDBASE/build"
        return 1
    }

    echo "Running RDKit tests..."
    ctest -j4 --output-on-failure

    cd - > /dev/null
}

# Uncomment the line below if you need to build/test RDKit
# build_rdkit

echo ""
echo "Roshambo environment setup complete!"
echo ""
echo "To use this environment in your reinvent-scoring configuration:"
echo "  \"environment_path\": \"source /path/to/setup_roshambo_env.sh &&\""
echo ""
echo "Or to activate manually:"
echo "  source setup_roshambo_env.sh"
