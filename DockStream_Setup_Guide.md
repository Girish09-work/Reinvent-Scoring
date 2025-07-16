# DockStream Setup Guide: AutoDock Vina with Open Source Tools

## Overview

This guide provides complete instructions for setting up and running DockStream with AutoDock Vina using only open-source tools. No proprietary software is required.

## Requirements and Dependencies

### 1. Core Software Requirements

**AutoDock Vina**:
- **Version**: 1.2.3 (included in DockStream)
- **Location**: Already included in `DockStreamCommunity/AutoDock-Vina-1.2.3/`
- **Status**: ✅ Open Source

**OpenBabel**:
- **Purpose**: File format conversions (PDB ↔ PDBQT ↔ SDF)
- **Installation**: Via conda/pip
- **Status**: ✅ Open Source

**RDKit**:
- **Purpose**: Molecule handling and conformer generation
- **Installation**: Via conda/pip
- **Status**: ✅ Open Source

### 2. Python Dependencies

**Core Python Packages**:
```bash
# Essential packages
conda install -c conda-forge rdkit
conda install -c conda-forge openbabel
conda install pydantic
conda install multiprocessing-logging

# Additional dependencies
pip install typing-extensions
pip install numpy pandas
```

**Complete Requirements List**:
```
rdkit>=2022.03.1
openbabel>=3.1.1
pydantic>=1.8.0
typing-extensions>=4.0.0
numpy>=1.21.0
pandas>=1.3.0
multiprocessing-logging>=0.3.0
```

### 3. System Requirements

**Operating System**: Windows, Linux, macOS
**Python Version**: 3.7-3.10
**Memory**: Minimum 4GB RAM (8GB+ recommended)
**CPU**: Multi-core recommended for parallelization

## Installation Steps

### Step 1: Create Conda Environment

```bash
# Create new environment
conda create -n DockStream python=3.9
conda activate DockStream

# Install core dependencies
conda install -c conda-forge rdkit openbabel
pip install pydantic typing-extensions
```

### Step 2: Verify DockStream Structure

Ensure your DockStream directory has this structure:
```
DockStream-master/
├── docker.py                          # Main entry point
├── dockstream/                         # Core package
│   ├── core/
│   │   └── AutodockVina/
│   │       └── AutodockVina_docker.py
│   └── utils/
│       └── parallelization/
├── DockStreamCommunity/
│   └── AutoDock-Vina-1.2.3/          # Vina binaries
└── dockstream_config.json             # Configuration file
```

### Step 3: Fix Import Issues

The import error you encountered is due to missing `Parallelization` class. I've already fixed this by updating the `__init__.py` file.

### Step 4: Configure AutoDock Vina Binary

**Option A: Use Included Binary**
```bash
# Make the included Vina binary executable (Linux/macOS)
chmod +x DockStreamCommunity/AutoDock-Vina-1.2.3/bin/vina

# For Windows, use the .exe file directly
```

**Option B: Install Vina Separately**
```bash
# Via conda
conda install -c conda-forge autodock-vina

# Or download from official source
# https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

## Configuration

### 1. Basic Configuration File

Create or modify `dockstream_config.json`:

```json
{
    "ligand_preparation": {
        "type": "RDkit",
        "target_ph": 7.4,
        "ph_tolerance": 1.0
    },
    "docking_runs": [{
        "backend": "AutoDockVina",
        "run_id": "AutoDockVina_run",
        "input_pools": ["RDkit_pool"],
        "parameters": {
            "binary_location": "vina",
            "parallelization": {
                "number_cores": 4
            },
            "seed": 42,
            "receptor_pdbqt_path": ["path/to/your/receptor.pdbqt"],
            "number_poses": 2,
            "search_space": {
                "--center_x": 0.0,
                "--center_y": 0.0,
                "--center_z": 0.0,
                "--size_x": 20.0,
                "--size_y": 20.0,
                "--size_z": 20.0
            }
        },
        "output": {
            "poses": {
                "poses_path": "output/docked_poses.sdf",
                "overwrite": true
            },
            "scores": {
                "scores_path": "output/docking_scores.csv",
                "overwrite": true
            }
        }
    }]
}
```

### 2. Prepare Receptor File

**Convert PDB to PDBQT**:
```bash
# Using OpenBabel
obabel receptor.pdb -opdbqt -O receptor.pdbqt -xr

# Or using AutoDockTools (if available)
prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
```

### 3. Prepare Input Ligands

**SMILES File Format** (`ligands.smi`):
```
CCO ethanol
CCN ethylamine
c1ccccc1 benzene
```

**SDF File**: Use any standard SDF file with 3D coordinates

## Running DockStream

### 1. Command Line Execution

**Basic Usage**:
```bash
cd DockStream-master
python docker.py -conf dockstream_config.json
```

**With SMILES Input**:
```bash
python docker.py -conf dockstream_config.json -smiles "CCO;CCN;c1ccccc1"
```

**With Output Prefix**:
```bash
python docker.py -conf dockstream_config.json -output_prefix "my_docking_run"
```

**Print Scores**:
```bash
python docker.py -conf dockstream_config.json -print_scores
```

### 2. Python Script Usage

```python
from dockstream.core.AutodockVina.AutodockVina_docker import AutodockVina
from dockstream.core.ligand.ligand import Ligand

# Configuration
config = {
    "backend": "AutoDockVina",
    "parameters": {
        "binary_location": "vina",
        "receptor_pdbqt_path": ["receptor.pdbqt"],
        "search_space": {
            "--center_x": 0.0,
            "--center_y": 0.0,
            "--center_z": 0.0,
            "--size_x": 20.0,
            "--size_y": 20.0,
            "--size_z": 20.0
        },
        "number_poses": 2,
        "seed": 42
    }
}

# Initialize docker
docker = AutodockVina(**config)

# Add ligands
smiles_list = ["CCO", "CCN", "c1ccccc1"]
ligands = [Ligand(smi, None, None) for smi in smiles_list]
docker.add_molecules(ligands)

# Run docking
docker.dock()

# Save results
docker.write_docked_ligands("output.sdf")
```

## Troubleshooting

### 1. Import Errors

**Error**: `ImportError: cannot import name 'Parallelization'`
**Solution**: The `__init__.py` file has been updated to fix this issue.

**Error**: `ModuleNotFoundError: No module named 'pydantic'`
**Solution**: 
```bash
pip install pydantic
```

### 2. Binary Execution Errors

**Error**: `vina: command not found`
**Solutions**:
1. Use full path to vina binary in configuration
2. Add vina to system PATH
3. Install vina via conda: `conda install -c conda-forge autodock-vina`

**Error**: Permission denied
**Solution** (Linux/macOS):
```bash
chmod +x /path/to/vina/binary
```

### 3. File Format Issues

**Error**: `Receptor PDBQT file not found`
**Solution**: Ensure receptor is properly converted to PDBQT format:
```bash
obabel receptor.pdb -opdbqt -O receptor.pdbqt -xr
```

**Error**: `Invalid SMILES`
**Solution**: Validate SMILES using RDKit:
```python
from rdkit import Chem
mol = Chem.MolFromSmiles("your_smiles")
if mol is None:
    print("Invalid SMILES")
```

### 4. Performance Issues

**Slow execution**: 
- Increase `number_cores` in configuration
- Reduce `number_poses` for faster execution
- Use smaller search space

**Memory issues**:
- Reduce batch size via `max_compounds_per_subjob`
- Process ligands in smaller chunks

## Output Files

### 1. Poses Output (SDF)

Contains 3D coordinates of docked poses:
- Multiple conformations per ligand
- Binding scores in SDF properties
- Ranked by binding affinity

### 2. Scores Output (CSV)

Contains numerical scores:
```csv
ligand_id,pose_id,binding_score,rmsd
mol_001,1,-8.5,0.0
mol_001,2,-7.2,1.5
mol_002,1,-9.1,0.0
```

### 3. Log Files

Detailed execution logs including:
- Processing steps
- Error messages
- Performance metrics
- File paths used

## Performance Optimization

### 1. Parallelization Settings

```json
"parallelization": {
    "number_cores": 8,              // Use available CPU cores
    "max_compounds_per_subjob": 50  // Batch size per process
}
```

### 2. Search Space Optimization

- Use smaller search boxes for faster execution
- Center search space on known binding site
- Typical box size: 15-25 Å per dimension

### 3. Pose Generation

```json
"number_poses": 1,    // Fewer poses = faster execution
"local_only": true,   // Local optimization only
"score_only": false   // Full docking vs scoring only
```

This setup provides a complete, open-source molecular docking pipeline using AutoDock Vina without any proprietary dependencies.
