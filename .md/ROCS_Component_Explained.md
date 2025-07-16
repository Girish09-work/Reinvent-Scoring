# ROCS in Reinvent Scoring

This document provides a detailed explanation of ROCS (Rapid Overlay of Chemical Structures) in the context of Reinvent's scoring components.

## 1. What is ROCS?

ROCS is a shape-based similarity search tool developed by OpenEye Scientific that:

- Evaluates 3D shape and pharmacophore similarity between molecules
- Aligns molecules to reference structures based on shape and chemical features
- Calculates similarity scores that guide the reinforcement learning process
- Helps identify molecules with similar binding modes to known active compounds
- Enables scaffold hopping by finding molecules with similar 3D properties but different 2D structures

## 2. ROCS Components in Reinvent

Reinvent implements ROCS functionality through two main components:

### 2.1 `rocs_similarity`

**What it is:** A scoring component that evaluates the 3D shape and pharmacophore (color) similarity between generated molecules and reference structures.

**How it works:**
- Uses OpenEye's ROCS technology to compare 3D shapes
- Aligns molecules to reference structures
- Calculates a combined score based on shape and color (pharmacophore) similarity

**Key parameters:**
- `shape_weight`: Weight given to shape similarity (typically 0.5)
- `color_weight`: Weight given to pharmacophore similarity (typically 0.5)
- `sim_measure`: Similarity measure to use (options: "Tanimoto", "RefTversky", "FitTversky")
- `rocs_input`: Path to reference structure file (SDF or shape query)
- `input_type`: Type of reference file ("sdf" or "shape_query")

**Example configuration:**
```json
{
    "component_type": "rocs_similarity",
    "name": "ROCS Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "rocs_input": "path/to/reference.sdf",
        "input_type": "sdf",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "sim_measure": "tanimoto"
    }
}
```

### 2.2 `parallel_rocs_similarity`

**What it is:** A parallel version of the ROCS similarity component that distributes calculations across multiple CPU cores.

**How it works:**
- Same as `rocs_similarity` but uses Python's multiprocessing to parallelize calculations
- Significantly faster for large batches of molecules
- Supports additional features like saving overlay files and negative volume constraints

**Additional parameters:**
- `max_num_cpus`: Maximum number of CPU cores to use (default: 4)
- `save_rocs_overlays`: Whether to save overlay files (default: false)
- `rocs_overlays_dir`: Directory to save overlays
- `rocs_overlays_prefix`: Prefix for overlay filenames
- `negative_volume`: Whether to use negative volume constraints (default: false)
- `protein_neg_vol_file`: Protein file for negative volume
- `ligand_neg_vol_file`: Ligand file for negative volume

## 3. ROCS Implementation Details

### 3.1 Core ROCS Similarity Calculation

The core ROCS similarity calculation is implemented in the `_calculate_omega_score` method:

```python
def _calculate_omega_score(self, smiles, step=-1) -> np.array:
    scores = []
    predicate = getattr(oeshape, self.sim_func_name_set.predicate)()
    for smile in smiles:
        imol = oechem.OEMol()
        best_score = 0.0
        if oechem.OESmilesToMol(imol, smile):
            if self.omega(imol):
                self.prep.Prep(imol)
                score = oeshape.OEBestOverlayScore()
                self.overlay.BestOverlay(score, imol, predicate)
                best_score_shape = getattr(score, self.sim_func_name_set.shape)()
                best_score_color = getattr(score, self.sim_func_name_set.color)()
                best_score_color = correct_color_score(best_score_color)
                best_score = ((self.shape_weight * best_score_shape) + (
                        self.color_weight * best_score_color)) / (self.shape_weight + self.color_weight)
        scores.append(best_score)
    return np.array(scores)
```

### 3.2 Similarity Measures

ROCS supports different similarity measures for comparing molecular shapes:

1. **Tanimoto** (default): Measures the overlap volume divided by the total volume of both molecules
2. **RefTversky**: Biases the similarity toward the reference molecule
3. **FitTversky**: Biases the similarity toward the fit molecule

These are implemented through a collection of similarity functions:

```python
def __similarity_collection(self, sim_measure_type):
    _SIM_FUNC = namedtuple('sim_func', ['shape', 'color', 'predicate'])
    _SIM_DEF_DICT = {
        self.sim_measure_enum.TANIMOTO: _SIM_FUNC('GetTanimoto', 'GetColorTanimoto', 'OEHighestTanimotoCombo'),
        self.sim_measure_enum.REF_TVERSKY: _SIM_FUNC('GetRefTversky', 'GetRefColorTversky',
                                                     'OEHighestRefTverskyCombo'),
        self.sim_measure_enum.FIT_TVERSKY: _SIM_FUNC('GetFitTversky', 'GetFitColorTversky',
                                                     'OEHighestFitTverskyCombo'),
    }
    return _SIM_DEF_DICT.get(sim_measure_type)
```

## 4. ROCS in Protac-INVENT and Link-INVENT

In Protac-INVENT and Link-INVENT, ROCS serves critical functions for 3D molecular shape analysis:

- Evaluates 3D shape similarity between generated molecules and reference structures
- Aligns molecules to reference compounds based on shape and chemical features
- Calculates similarity scores that guide the reinforcement learning process
- Helps identify molecules with similar binding modes to known active compounds
- Enables scaffold hopping by finding molecules with similar 3D properties but different 2D structures

These capabilities are particularly important in Link-INVENT, where the 3D arrangement of the linker between two binding warheads is critical for proper protein-protein interaction.

## 5. Open-Source Alternatives to ROCS

Several open-source alternatives to OpenEye's ROCS are available through RDKit:

### 5.1 USR and USRCAT

Ultrafast Shape Recognition (USR) and its extension USRCAT (USR with Chemical Annotation Types) are efficient alternatives to ROCS for shape-based similarity calculations.

```python
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT

def calculate_shape_similarity(query_mol, ref_mol):
    # Generate USR descriptors
    query_descriptor = GetUSRCAT(query_mol)
    ref_descriptor = GetUSRCAT(ref_mol)
    # Calculate similarity (returns value between 0-1)
    similarity = GetUSRScore(query_descriptor, ref_descriptor)
    return similarity
```

### 5.2 RDKit's Shape Alignment

For more ROCS-like explicit alignment and scoring:

```python
from rdkit.Chem import AllChem

def align_molecules(query_mol, ref_mol):
    # Align query molecule to reference
    pyO3A = AllChem.GetO3A(query_mol, ref_mol)
    score = pyO3A.Score()
    pyO3A.Align()
    return score, query_mol
```

## 6. Implementing RDKit-based ROCS Alternative

To implement an RDKit-based alternative to ROCS in Reinvent, you would create a new scoring component:

```python
class RDKitShapeSimilarity(BaseROCSComponent):
    """RDKit-based shape similarity calculator to replace OpenEye's ROCS."""
    
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        
        # Extract parameters
        self.shape_weight = self.parameters.specific_parameters.get("shape_weight", 0.5)
        self.color_weight = self.parameters.specific_parameters.get("color_weight", 0.5)
        self.method = self.parameters.specific_parameters.get("method", "usrcat")
        
        # Load reference molecule
        ref_file = self.parameters.specific_parameters.get("rocs_input")
        self.reference_mol = self._load_reference(ref_file)
        
        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=self.parameters.specific_parameters.get("max_confs", 200),
            energy_window=self.parameters.specific_parameters.get("ewindow", 10),
            max_stereo=self.parameters.specific_parameters.get("max_stereo", 0)
        )
```

## 7. Conclusion

ROCS is a powerful tool for 3D shape-based similarity calculations in Reinvent's scoring system. It enables the generation of molecules with similar 3D properties to known active compounds, which is crucial for structure-based drug design. While it requires an OpenEye license, open-source alternatives based on RDKit are available and can be integrated into the Reinvent framework.
