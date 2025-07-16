#!/bin/bash

# =============================================================================
# Roshambo Environment for Reinvent-Scoring
# =============================================================================
# Minimal environment setup for use with reinvent-scoring
# This script is designed to be called from the environment_path parameter
# =============================================================================

# Activate conda environment
# Try different conda initialization paths
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source /opt/conda/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/miniconda3/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/miniconda3/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/.conda/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/.conda/etc/profile.d/conda.sh
fi

conda activate roshambo

# Set essential environment variables
export RDBASE="/home/protacinvent/Desktop/roshambo/rdkit"
export PYTHONPATH="$RDBASE:$CONDA_PREFIX/lib/python3.8/site-packages"
export LD_LIBRARY_PATH="$RDBASE/lib:$CONDA_PREFIX/lib:/usr/local/cuda/lib64:$LD_LIBRARY_PATH"
export CUDA_HOME="/usr/local/cuda"
export CUDA_VISIBLE_DEVICES="0"

# Verify critical imports (silent)
python -c "
import sys
try:
    from rdkit import Chem
    import roshambo
    # Silent success
except ImportError as e:
    print(f'Environment setup error: {e}', file=sys.stderr)
    sys.exit(1)
" 2>/dev/null || {
    echo "Error: Failed to import required packages" >&2
    exit 1
}
