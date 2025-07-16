# Reinvent Scoring RDKit - Usage Guide

This guide explains how to install and use the `reinvent-scoring-rdkit` package, which provides RDKit-based shape similarity components for the Reinvent molecular generation platform.

## Installation

### 1. Install the package from PyPI

```bash
pip install reinvent-scoring-rdkit
```

### 2. Verify installation

```python
import reinvent_scoring
print(reinvent_scoring.__version__)
```

## Key Differences from Original Package

This package is a fork of the original `reinvent-scoring` package with these key differences:

1. **Open-Source Shape Similarity**: Uses RDKit instead of OpenEye's ROCS for 3D shape comparison
2. **No License Required**: All components work without commercial licenses
3. **Import Compatibility**: Uses the same import structure as the original package

## Using RDKit Shape Components

The package provides two main components:
- `rdkit_shape_similarity`: Basic RDKit shape similarity scoring
- `parallel_rdkit_shape_similarity`: Multi-threaded version for faster processing

### Configuration in JSON

#### Basic RDKit Shape Similarity

```json
{
    "component_type": "rdkit_shape_similarity",
    "name": "RDKit Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "/path/to/your/reference.sdf",
        "method": "usrcat",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "max_confs": 50,
        "ewindow": 10,
        "max_stereo": 0
    }
}
```

#### Parallel RDKit Shape Similarity

```json
{
    "component_type": "parallel_rdkit_shape_similarity",
    "name": "Parallel RDKit Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "/path/to/your/reference.sdf",
        "method": "o3a",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "max_confs": 50,
        "max_num_cpus": 8,
        "save_overlays": true,
        "overlays_dir": "/path/to/save/overlays",
        "overlay_prefix": "mol_"
    }
}
```

### Parameter Explanation

#### Common Parameters

- `reference_file`: (Required) Full path to the SDF file containing your reference molecule(s)
- `method`: (Required) Shape comparison method - either "usrcat" or "o3a"
- `shape_weight`: Weight given to shape similarity (default: 0.5)
- `color_weight`: Weight given to pharmacophore/color similarity (default: 0.5)
- `max_confs`: Maximum number of conformers to generate (default: 50)
- `ewindow`: Energy window for conformer selection in kcal/mol (default: 10)
- `max_stereo`: Maximum number of stereoisomers to consider (default: 0, meaning use input stereochemistry)

#### Parallel-Specific Parameters

- `max_num_cpus`: Maximum number of CPU cores to use (default: 4)
- `save_overlays`: Whether to save overlay structures (default: false)
- `overlays_dir`: Directory to save overlay structures (default: "overlays")
- `overlay_prefix`: Prefix for overlay filenames (default: "mol_")

## Complete Scoring Function Example

Here's a complete example of a scoring function that combines RDKit shape similarity with other components:

```json
{
    "name": "custom_product",
    "parallel": true,
    "parameters": [
        {
            "component_type": "parallel_rdkit_shape_similarity",
            "name": "Shape Similarity",
            "weight": 1.0,
            "specific_parameters": {
                "reference_file": "/path/to/your/reference.sdf",
                "method": "usrcat",
                "shape_weight": 0.5,
                "color_weight": 0.5,
                "max_confs": 50,
                "max_num_cpus": 8
            }
        },
        {
            "component_type": "qed_score",
            "name": "QED Score",
            "weight": 1.0,
            "specific_parameters": {
                "transformation": {
                    "transformation_type": "no_transformation"
                }
            }
        },
        {
            "component_type": "molecular_weight",
            "name": "Molecular Weight",
            "weight": 1.0,
            "specific_parameters": {
                "transformation": {
                    "transformation_type": "double_sigmoid",
                    "high": 500,
                    "low": 200,
                    "coef_div": 500,
                    "coef_si": 20,
                    "coef_se": 20
                }
            }
        }
    ]
}
```

## Using in Reinvent

To use this scoring function in Reinvent, save your configuration to a JSON file and reference it in your Reinvent configuration:

```json
{
    "version": 3,
    "run_type": "reinforcement_learning",
    "model_type": "default",
    "parameters": {
        "scoring_function": {
            "name": "custom_product",
            "parallel": true,
            "parameters": [
                {
                    "component_type": "parallel_rdkit_shape_similarity",
                    "name": "Shape Similarity",
                    "weight": 1.0,
                    "specific_parameters": {
                        "reference_file": "/path/to/your/reference.sdf",
                        "method": "usrcat",
                        "shape_weight": 0.5,
                        "color_weight": 0.5,
                        "max_confs": 50,
                        "max_num_cpus": 8
                    }
                },
                // Other components...
            ]
        },
        // Other Reinvent parameters...
    }
}
```

## Custom Configuration (Optional)

If you need to customize the package configuration:

1. Create a JSON configuration file (e.g., `config.json`):

```json
{
  "DEVELOPMENT_ENVIRONMENT": false,
  "MAIN_TEST_PATH": "tmp_test_folder",
  "COMPONENT_SPECIFIC": {
    "AZDOCK": {
      "AZDOCK_DOCKER_SCRIPT_PATH": "",
      "AZDOCK_ENV_PATH": "",
      "AZDOCK_DEBUG": false
    },
    "DOCKSTREAM": {
      "DOCKSTREAM_DOCKER_SCRIPT_PATH": "",
      "DOCKSTREAM_ENV_PATH": "",
      "DOCKSTREAM_DEBUG": false
    },
    "ICOLOS": {
      "ICOLOS_EXECUTOR_PATH": "",
      "ICOLOS_DEBUG": false
    }
  },
  "ENVIRONMENTAL_VARIABLES": {
    "PIP_URL": "",
    "PIP_KEY": "",
    "PIP_GET_RESULTS": ""
  }
}
```

2. When running your scripts, specify the config file location:

```python
import sys
sys.argv.append('--base_config')
sys.argv.append('/path/to/your/config.json')

# Then import reinvent_scoring
import reinvent_scoring
```

## Preparing Reference Files

The reference SDF file should contain 3D coordinates for your reference molecule(s). If you have a SMILES string, you can convert it to a 3D SDF file using RDKit:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Create molecule from SMILES
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
mol = Chem.MolFromSmiles(smiles)

# Add hydrogens
mol = Chem.AddHs(mol)

# Generate 3D coordinates
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# Optimize the geometry
AllChem.MMFFOptimizeMolecule(mol)

# Write to SDF file
writer = Chem.SDWriter("/path/to/your/reference.sdf")
writer.write(mol)
writer.close()
```

## Troubleshooting

### Common Issues

1. **Missing reference file**: Ensure the path to your reference SDF file is absolute and the file exists
2. **Invalid SDF format**: Verify your SDF file contains valid 3D coordinates
3. **Method selection**: "usrcat" is faster but less accurate, "o3a" is slower but more accurate

### Debugging

To debug issues with the shape similarity calculation, you can:

1. Set `save_overlays` to `true` to inspect the molecular alignments
2. Try with a simpler reference molecule first
3. Reduce `max_confs` to speed up calculations during testing

## Advanced Usage

### Custom Transformation Functions

You can apply transformations to the raw similarity scores:

```json
"specific_parameters": {
    "reference_file": "/path/to/your/reference.sdf",
    "method": "usrcat",
    "shape_weight": 0.5,
    "color_weight": 0.5,
    "transformation": {
        "transformation_type": "sigmoid",
        "high": 1.0,
        "low": 0.0,
        "k": 0.5
    }
}
```

### Combining with Docking

For more advanced workflows, you can combine RDKit shape similarity with docking:

1. First use shape similarity as a fast pre-filter
2. Then dock the most promising candidates

This can be achieved by using both components in your scoring function with appropriate weights.
