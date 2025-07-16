# ROCS in Protac-INVENT: Functionality, Implementation, and Open-Source Alternatives

## 1. What ROCS Does in Protac-INVENT

In Protac-INVENT, ROCS (Rapid Overlay of Chemical Structures) serves several critical functions that enable the design of effective PROTACs (Proteolysis Targeting Chimeras):

### 1.1 Shape-Based Similarity Scoring

ROCS evaluates the 3D shape similarity between generated linker molecules and reference structures. This is crucial because:

- **Spatial Arrangement**: PROTACs require precise spatial positioning of the two warheads (E3 ligase ligand and target protein ligand) to facilitate protein-protein interaction
- **Linker Conformation**: The linker's 3D shape determines whether the two proteins can be brought into proximity for ubiquitination
- **Bioactive Conformation**: Ensures the generated molecules can adopt conformations similar to known effective PROTACs

### 1.2 Pharmacophore Feature Matching

Beyond pure shape, ROCS also evaluates "color" similarity, which represents pharmacophore features:

- **H-bond Donors/Acceptors**: Identifies potential hydrogen bonding interactions
- **Hydrophobic Regions**: Maps areas that contribute to hydrophobic interactions
- **Charged Groups**: Locates positively and negatively charged regions
- **Aromatic Systems**: Identifies π-stacking opportunities

### 1.3 Guiding the Generative Model

Within the reinforcement learning framework of Protac-INVENT:

- ROCS similarity scores serve as rewards that guide the generative model
- Higher scores are given to molecules with shapes similar to known effective linkers
- This steers the model toward generating linkers with appropriate 3D properties
- Enables "scaffold hopping" by finding molecules with similar 3D properties but different 2D structures

### 1.4 Conformational Analysis

ROCS (together with OMEGA) helps analyze the conformational preferences of generated linkers:

- Evaluates flexibility and rigidity of proposed linkers
- Identifies preferred conformations that might facilitate protein-protein interactions
- Helps predict whether the linker will allow the two warheads to adopt productive binding orientations

## 2. How ROCS Works in Protac-INVENT

### 2.1 Overall Workflow

The implementation of ROCS in Protac-INVENT follows these steps:

1. Generate 3D conformers for candidate molecules using OMEGA
2. Align these conformers to reference structure(s)
3. Calculate shape and pharmacophore similarity scores
4. Combine scores into a final similarity metric
5. Feed this metric into the reinforcement learning scoring function

### 2.2 Code Implementation

The core ROCS functionality is implemented in the `RocsSimilarity` and `ParallelRocsSimilarity` classes. Here's how the scoring works:

```python
def _calculate_omega_score(self, smiles, step=-1) -> np.array:
    scores = []
    predicate = getattr(oeshape, self.sim_func_name_set.predicate)()
    for smile in smiles:
        imol = oechem.OEMol()
        best_score = 0.0
        if oechem.OESmilesToMol(imol, smile):
            if self.omega(imol):  # Generate conformers with OMEGA
                self.prep.Prep(imol)  # Prepare molecule for ROCS
                score = oeshape.OEBestOverlayScore()
                self.overlay.BestOverlay(score, imol, predicate)  # Perform alignment
                
                # Get shape and color scores
                best_score_shape = getattr(score, self.sim_func_name_set.shape)()
                best_score_color = getattr(score, self.sim_func_name_set.color)()
                best_score_color = correct_color_score(best_score_color)
                
                # Calculate weighted average
                best_score = ((self.shape_weight * best_score_shape) + 
                             (self.color_weight * best_score_color)) / 
                             (self.shape_weight + self.color_weight)
        scores.append(best_score)
    return np.array(scores)
```

### 2.3 Key Parameters

The behavior of ROCS in Protac-INVENT can be customized through several parameters:

- **Reference Structure**: SDF file or shape query containing the reference linker conformation
- **Shape Weight**: Importance of pure shape similarity (typically 0.5)
- **Color Weight**: Importance of pharmacophore feature matching (typically 0.5)
- **Similarity Measure**: Method for calculating similarity (Tanimoto or Tversky)
- **Conformer Settings**: Maximum number of conformers, energy window, etc.

### 2.4 Integration with Reinforcement Learning

ROCS similarity scores are incorporated into the reinforcement learning process:

```python
# Example configuration in scoring function
{
    "component_type": "parallel_rocs_similarity",
    "name": "Linker Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "rocs_input": "path/to/reference_linker.sdf",
        "input_type": "sdf",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "sim_measure": "tanimoto"
    }
}
```

This component is then used as part of the scoring function that guides the generative model toward producing linkers with appropriate 3D properties.

## 3. Open-Source Alternatives to ROCS

Several open-source alternatives can replace ROCS functionality in Protac-INVENT:

### 3.1 USR and USRCAT (Ultrafast Shape Recognition)

USR is a method for rapid shape similarity calculations without explicit alignment:

```python
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT

def calculate_usr_similarity(query_mol, ref_mol):
    # Generate USR descriptors
    query_descriptor = GetUSRCAT(query_mol)
    ref_descriptor = GetUSRCAT(ref_mol)
    
    # Calculate similarity (returns value between 0-1)
    similarity = GetUSRScore(query_descriptor, ref_descriptor)
    return similarity
```

**Advantages:**
- Very fast (orders of magnitude faster than ROCS)
- No alignment required
- USRCAT extends USR with chemical feature information
- Fully integrated in RDKit

**Limitations:**
- Less accurate than explicit shape overlay methods
- No visualization of the alignment
- Less sophisticated pharmacophore handling

### 3.2 RDKit Shape Alignment with O3A

RDKit's Open3DAlign (O3A) provides more ROCS-like functionality:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def calculate_shape_similarity_o3a(query_mol, ref_mol):
    # Align molecules using Open3DAlign
    pyO3A = AllChem.GetO3A(query_mol, ref_mol)
    
    # Get alignment score
    score = pyO3A.Score() / 100.0  # Normalize to 0-1 range
    
    # Perform the alignment
    pyO3A.Align()
    
    return score, query_mol  # Return score and aligned molecule
```

**Advantages:**
- Explicit alignment like ROCS
- Can visualize the alignment
- Considers both shape and chemical features
- Free and open-source

**Limitations:**
- Slower than USR
- Less mature than ROCS
- Requires more parameter tuning

### 3.3 Shape-it

Shape-it is an open-source shape alignment tool from Silicos-it:

```python
# Example of calling Shape-it from Python
import subprocess

def run_shape_it(query_mol_file, ref_mol_file, output_file):
    cmd = [
        "shape-it",
        "--reference", ref_mol_file,
        "--dbase", query_mol_file,
        "--out", output_file,
        "--best", "1"
    ]
    subprocess.run(cmd, check=True)
    
    # Parse results
    with open(output_file, 'r') as f:
        lines = f.readlines()
        # Extract score from output
        score = float(lines[1].split()[1])
    
    return score
```

**Advantages:**
- Specifically designed for shape comparison
- Relatively fast
- Command-line tool that can be integrated into workflows

**Limitations:**
- Less feature-rich than ROCS
- Less integration with Python
- Less active development

### 3.4 Implementing a Complete ROCS Alternative

A complete RDKit-based alternative to ROCS for Protac-INVENT would combine conformer generation with shape comparison:

```python
class RDKitShapeSimilarity:
    """RDKit-based shape similarity calculator to replace OpenEye's ROCS."""
    
    def __init__(self, reference_file, shape_weight=0.5, feature_weight=0.5, method="usrcat"):
        self.shape_weight = shape_weight
        self.feature_weight = feature_weight
        self.method = method
        
        # Load reference molecule
        self.reference_mol = self._load_reference(reference_file)
        
        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=50,
            energy_window=10
        )
    
    def calculate_similarity(self, smiles):
        # Generate conformers for query molecule
        query_mol = self.conformer_generator.generate_conformers(smiles)
        if query_mol is None:
            return 0.0
        
        # Calculate similarity based on method
        if self.method == "usrcat":
            return self._calculate_usrcat_similarity(query_mol)
        elif self.method == "o3a":
            return self._calculate_o3a_similarity(query_mol)
        else:
            return 0.0
    
    def _calculate_usrcat_similarity(self, query_mol):
        # Implementation using USRCAT
        best_score = 0.0
        
        # Calculate for all conformer pairs
        for q_conf_id in range(query_mol.GetNumConformers()):
            # Get USRCAT descriptors
            query_descriptor = GetUSRCAT(query_mol, confId=q_conf_id)
            ref_descriptor = GetUSRCAT(self.reference_mol)
            
            # Calculate shape similarity
            shape_sim = GetUSRScore(query_descriptor, ref_descriptor)
            
            # For feature similarity, we could use feature-based USRCAT components
            feature_sim = self._calculate_feature_similarity(query_mol, self.reference_mol)
            
            # Combine scores
            combined_score = ((self.shape_weight * shape_sim) + 
                             (self.feature_weight * feature_sim)) / (self.shape_weight + self.feature_weight)
            
            if combined_score > best_score:
                best_score = combined_score
                
        return best_score
```

## 4. Conclusion

ROCS plays a crucial role in Protac-INVENT by evaluating the 3D shape and pharmacophore similarity of generated linkers to reference structures. This guides the generative model toward producing molecules with appropriate spatial properties for bringing two proteins into proximity.

While ROCS is a powerful commercial tool, several open-source alternatives are available through RDKit and other packages. These alternatives, particularly USR/USRCAT and O3A, can provide similar functionality without requiring commercial licenses, making them viable options for implementing shape-based scoring in Protac-INVENT.

By replacing ROCS with these open-source alternatives, researchers can develop and deploy Protac-INVENT models without the licensing constraints of proprietary software while still benefiting from 3D shape-based molecular design.
