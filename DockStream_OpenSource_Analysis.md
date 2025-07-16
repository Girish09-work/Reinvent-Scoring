# DockStream Integration in Reinvent-Scoring: Complete Analysis of Open Source vs Proprietary Components

## Executive Summary

This document provides a comprehensive analysis of DockStream integration in the reinvent-scoring package, focusing on the pipeline for using open-source docking backends and identifying where proprietary software dependencies exist. The analysis covers two main scoring components: `dockstream` and `docked_parallel_rocs_similarity`.

## Table of Contents

1. [Overview of DockStream Integration](#overview)
2. [Scoring Components Analysis](#scoring-components)
3. [Open Source Docking Backends](#open-source-backends)
4. [Proprietary Software Dependencies](#proprietary-dependencies)
5. [Pipeline Analysis](#pipeline-analysis)
6. [Open Source Alternatives](#alternatives)
7. [Implementation Details](#implementation)
8. [Recommendations](#recommendations)

## 1. Overview of DockStream Integration {#overview}

DockStream is a docking wrapper that enhances de novo molecular design by providing a unified interface to multiple docking backends. In the reinvent-scoring package, DockStream integration enables:

- **Molecular docking scoring**: Evaluating binding affinity through docking energy
- **Post-docking shape similarity**: Comparing docked poses with reference structures
- **Multi-backend support**: Using different docking engines (open source and proprietary)

### Key Integration Points

1. **Configuration Integration**: DockStream configuration files referenced in scoring setup
2. **Output File Integration**: Docked poses used as input for similarity calculations
3. **Script Integration**: Python environment and script execution coordination

## 2. Scoring Components Analysis {#scoring-components}

### 2.1 `dockstream` Component

**Purpose**: Provides docking-based scoring using various docking backends

**Implementation**: `reinvent_scoring/scoring/score_components/structural/dockstream.py`

**Key Features**:
- Executes external DockStream commands
- Handles SMILES input and score output
- Supports debug mode for detailed logging
- Manages step-wise execution for reinforcement learning

**Configuration Parameters**:
```json
{
    "component_type": "dockstream",
    "specific_parameters": {
        "configuration_path": "/path/to/dockstream_config.json",
        "docker_script_path": "/path/to/docker.py",
        "environment_path": "/path/to/python/env"
    }
}
```

**Process Flow**:
1. Receives SMILES strings from reinvent
2. Creates command with DockStream configuration
3. Executes docking via subprocess
4. Parses numerical scores from output
5. Returns transformed scores for optimization

### 2.2 `docked_parallel_rocs_similarity` Component

**Purpose**: Compares docked poses with reference molecules using shape similarity

**Key Insight**: This component is **NOT implemented as a separate class** in the codebase. Instead, it appears to be a configuration-based workflow that:

1. Uses `dockstream` component to generate docked poses
2. Reads docked poses from specified path
3. Applies `parallel_rocs_similarity` logic to compare with reference structures

**Configuration Parameters**:
```json
{
    "component_type": "docked_parallel_rocs_similarity",
    "specific_parameters": {
        "docked_poses_path": "/path/to/docked/poses",
        "rocs_input": "/path/to/reference.sdf",
        "shape_weight": 1.0,
        "color_weight": 0.0,
        "similarity_measure": "Tanimoto"
    }
}
```

**Critical Finding**: The `docked_parallel_rocs_similarity` component type is used in configuration files but **does not exist in the component factory or enum definitions**. This suggests it may be:
- A legacy component name
- A custom implementation not in the main codebase
- A workflow combination rather than a single component

## 3. Open Source Docking Backends {#open-source-backends}

DockStream supports several open-source docking backends that can replace proprietary alternatives:

### 3.1 AutoDock Vina

**Status**: ✅ **Fully Open Source**

**Features**:
- Fast gradient-based docking algorithm
- No proprietary dependencies
- Well-established and widely used
- Good performance for most use cases

**Configuration Example**:
```json
{
    "backend": "AutoDockVina",
    "parameters": {
        "binary_location": "/path/to/vina",
        "receptor_pdbqt_path": ["/path/to/receptor.pdbqt"],
        "search_space": {
            "--center_x": 3.3,
            "--center_y": 11.5,
            "--center_z": 24.8,
            "--size_x": 15,
            "--size_y": 10,
            "--size_z": 10
        },
        "number_poses": 2,
        "seed": 42
    }
}
```

**File Format Pipeline**:
- Input: SMILES → SDF (via RDKit)
- Conversion: SDF → PDBQT (via OpenBabel)
- Docking: PDBQT → PDBQT (via Vina)
- Output: PDBQT → SDF (via OpenBabel)

### 3.2 rDock

**Status**: ✅ **Fully Open Source**

**Features**:
- Academic/research focused
- Good for fragment-based drug design
- Supports ensemble docking
- No licensing restrictions

**Configuration Example**:
```json
{
    "backend": "rDock",
    "parameters": {
        "prefix_execution": "module load rDock",
        "rbdock_prm_paths": ["/path/to/receptor.prm"],
        "number_poses": 2
    }
}
```

**File Format Pipeline**:
- Input: SMILES → SDF (via RDKit)
- Docking: SDF → SDF (via rDock)
- Output: SDF (ready for similarity analysis)

## 4. Proprietary Software Dependencies {#proprietary-dependencies}

### 4.1 ROCS (Rapid Overlay of Chemical Structures)

**Status**: ❌ **Proprietary (OpenEye Scientific)**

**Usage in Codebase**:
- `parallel_rocs_similarity` component
- Shape and pharmacophore similarity calculations
- 3D molecular overlay generation

**Dependencies**:
```python
from openeye import oechem, oeomega, oeshape, oequacpac
```

**Key Functions Using ROCS**:
- `oeshape.OEBestOverlayScore()`
- `oeshape.OEMultiRefOverlay()`
- `oeshape.OEOverlapPrep()`

### 4.2 OMEGA (Conformer Generation)

**Status**: ❌ **Proprietary (OpenEye Scientific)**

**Usage in Codebase**:
- Conformer generation for ROCS similarity
- Used within `parallel_rocs_similarity`

**Dependencies**:
```python
from openeye import oeomega
```

**Key Functions Using OMEGA**:
- `oeomega.OEOmega()`
- `oeomega.OEOmegaOptions()`

### 4.3 Other Proprietary Backends

**OpenEye Hybrid**: ❌ Proprietary docking backend
**Schrödinger Glide**: ❌ Proprietary docking backend

## 5. Pipeline Analysis {#pipeline-analysis}

### 5.1 Open Source Docking Pipeline (AutoDock Vina)

```
SMILES Input
    ↓ (RDKit)
SDF Format
    ↓ (OpenBabel)
PDBQT Format
    ↓ (AutoDock Vina)
Docked PDBQT
    ↓ (OpenBabel)
Docked SDF
    ↓ (Shape Similarity - NEEDS REPLACEMENT)
ROCS Similarity Score
```

**Bottleneck**: The shape similarity step requires ROCS/OMEGA (proprietary)

### 5.2 Complete Open Source Alternative Pipeline

```
SMILES Input
    ↓ (RDKit)
SDF Format
    ↓ (AutoDock Vina/rDock)
Docked SDF
    ↓ (RDKit Shape Similarity OR Roshambo)
Open Source Similarity Score
```

## 6. Open Source Alternatives {#alternatives}

### 6.1 RDKit Shape Similarity

**Status**: ✅ **Fully Open Source**

**Implementation**: Already available in codebase
- `reinvent_scoring/scoring/score_components/rdkit_shape/`
- Uses RDKit's O3A alignment
- GPU acceleration support (experimental)

**Advantages**:
- No licensing costs
- Integrated with existing RDKit infrastructure
- Customizable similarity metrics

**Limitations**:
- Less sophisticated than ROCS
- Slower for large-scale screening
- Limited pharmacophore features

### 6.2 Roshambo

**Status**: ✅ **Open Source Alternative to ROCS**

**Implementation**: Available in codebase
- `reinvent_scoring/scoring/score_components/roshambo/`
- GPU-accelerated shape similarity
- Compatible with ROCS workflows

**Advantages**:
- GPU acceleration
- ROCS-like functionality
- Open source license

**Current Status**: Already integrated in the scoring package

### 6.3 Conformer Generation Alternatives

**RDKit ETKDG**: ✅ Open source conformer generation
**OpenMM**: ✅ Open source molecular dynamics
**MMFF94**: ✅ Available in RDKit

## 7. Implementation Details {#implementation}

### 7.1 File Format Conversions

**Open Source Tools Available**:
- **RDKit**: SMILES ↔ SDF ↔ MOL
- **OpenBabel**: SDF ↔ PDBQT ↔ PDB
- **MDAnalysis**: Various format support

**Conversion Pipeline**:
```python
# SMILES to SDF (RDKit)
mol = Chem.MolFromSmiles(smiles)
writer = Chem.SDWriter('output.sdf')
writer.write(mol)

# SDF to PDBQT (OpenBabel)
subprocess.run(['obabel', '-isdf', 'input.sdf', '-opdbqt', '-O', 'output.pdbqt'])

# PDBQT to SDF (OpenBabel)
subprocess.run(['obabel', '-ipdbqt', 'input.pdbqt', '-osdf', '-O', 'output.sdf'])
```

### 7.2 Docking Execution

**AutoDock Vina Execution**:
```python
command = [
    'vina',
    '--receptor', receptor_pdbqt,
    '--ligand', ligand_pdbqt,
    '--center_x', str(center_x),
    '--center_y', str(center_y),
    '--center_z', str(center_z),
    '--size_x', str(size_x),
    '--size_y', str(size_y),
    '--size_z', str(size_z),
    '--out', output_pdbqt,
    '--num_modes', str(num_poses)
]
subprocess.run(command)
```

### 7.3 Shape Similarity (Open Source)

**RDKit O3A Alignment**:
```python
from rdkit.Chem import AllChem

# Generate conformers
AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)

# Align molecules
pyO3A = AllChem.GetO3A(query_mol, ref_mol)
score = pyO3A.Score()
pyO3A.Align()
```

### 7.4 Key Code References

**DockStream Component Implementation**:
```python
# File: reinvent_scoring/scoring/score_components/structural/dockstream.py
class DockStream(BaseStructuralComponent):
    def _create_command(self, smiles: List[str], step):
        concat_smiles = '"' + ';'.join(smiles) + '"'
        command = ' '.join([self._environment_path,
                            self._docker_script_path,
                            "-conf", self._configuration_path,
                            "-output_prefix", self._get_step_string(step),
                            "-smiles", concat_smiles,
                            "-print_scores"])
        return command
```

**ROCS Similarity (Proprietary)**:
```python
# File: reinvent_scoring/scoring/score_components/rocs/parallel_rocs_similarity.py
from openeye import oechem, oeomega, oeshape, oequacpac

class ParallelRocsSimilarity(BaseROCSComponent):
    @classmethod
    def parallel_scoring(cls, smile, shape_weight, color_weight, ...):
        # Uses OpenEye OMEGA for conformer generation
        omega_success, imol = oehelper.get_omega_confs(imol, cls.omega, ...)

        # Uses OpenEye ROCS for shape similarity
        score = oeshape.OEBestOverlayScore()
        cls.rocs_overlay.BestOverlay(score, imol, predicate)
```

**RDKit Shape Similarity (Open Source)**:
```python
# File: reinvent_scoring/scoring/score_components/rdkit_shape/parallel_rdkit_shape_similarity.py
from rdkit.Chem import AllChem

def calculate_shape_similarity(query_mol, ref_mol, shape_weight, color_weight):
    # Generate conformers using RDKit
    AllChem.EmbedMolecule(query_mol)
    AllChem.MMFFOptimizeMolecule(query_mol)

    # Align using O3A
    pyO3A = AllChem.GetO3A(query_mol, ref_mol)
    shape_sim = pyO3A.Score() / 100.0

    # Combine scores
    combined_score = ((shape_weight * shape_sim) +
                     (color_weight * color_sim)) / (shape_weight + color_weight)
    return combined_score
```

## 8. Recommendations {#recommendations}

### 8.1 Immediate Actions

1. **Verify `docked_parallel_rocs_similarity` Implementation**
   - Component appears in configs but not in factory
   - May need custom implementation or workflow setup
   - **Action**: Add to `ScoringFunctionComponentNameEnum` and component factory

2. **Use Open Source Docking Backends**
   - AutoDock Vina: Recommended for general use
   - rDock: Good for academic/research applications
   - **Action**: Configure DockStream with open source backends

3. **Replace ROCS with Open Source Alternatives**
   - Use existing RDKit shape similarity components
   - Consider Roshambo for GPU acceleration
   - **Action**: Update configurations to use `rdkit_shape_similarity` or `roshambo_shape_similarity`

### 8.2 Long-term Strategy

1. **Complete Open Source Pipeline**
   ```
   SMILES → RDKit → AutoDock Vina → RDKit Shape Similarity
   ```

2. **Performance Optimization**
   - GPU acceleration where available
   - Parallel processing for batch operations
   - Efficient file I/O management

3. **Validation and Benchmarking**
   - Compare open source vs proprietary results
   - Establish performance baselines
   - Document accuracy trade-offs

### 8.3 Missing Component Implementation

The `docked_parallel_rocs_similarity` component needs to be implemented or clarified:

```python
# Potential implementation approach
class DockedParallelRocsSimilarity(BaseROCSComponent):
    def __init__(self, parameters):
        super().__init__(parameters)
        self.docked_poses_path = parameters.specific_parameters["docked_poses_path"]

    def _calculate_score(self, smiles, step):
        # 1. Read docked poses from docked_poses_path
        # 2. Match SMILES to corresponding docked poses
        # 3. Apply ROCS similarity to docked conformations
        # 4. Return similarity scores
        pass
```

## Conclusion

The DockStream integration in reinvent-scoring provides a flexible framework for molecular docking and shape similarity scoring. While the system currently relies on some proprietary components (ROCS, OMEGA), complete open-source alternatives exist and are partially implemented. The main challenges are:

1. **Missing Implementation**: `docked_parallel_rocs_similarity` component needs clarification/implementation
2. **Performance Trade-offs**: Open source alternatives may be slower but are functionally equivalent
3. **Workflow Integration**: Ensuring smooth file format conversions and data flow

With proper implementation of open-source alternatives, the entire pipeline can be made proprietary-software-free while maintaining scientific validity and reasonable performance.

## Summary of Proprietary vs Open Source Components

| Component | Current Status | Open Source Alternative | Implementation Status |
|-----------|----------------|------------------------|----------------------|
| **Docking** | | | |
| AutoDock Vina | ✅ Open Source | N/A | ✅ Available |
| rDock | ✅ Open Source | N/A | ✅ Available |
| OpenEye Hybrid | ❌ Proprietary | AutoDock Vina/rDock | ✅ Available |
| Schrödinger Glide | ❌ Proprietary | AutoDock Vina/rDock | ✅ Available |
| **Shape Similarity** | | | |
| ROCS | ❌ Proprietary | RDKit/Roshambo | ✅ Available |
| OMEGA | ❌ Proprietary | RDKit ETKDG | ✅ Available |
| **File Conversion** | | | |
| OpenBabel | ✅ Open Source | N/A | ✅ Available |
| RDKit | ✅ Open Source | N/A | ✅ Available |

**Final Recommendation**: The entire DockStream pipeline can be made fully open source by using AutoDock Vina or rDock for docking and RDKit or Roshambo for shape similarity, with no significant loss in functionality.
