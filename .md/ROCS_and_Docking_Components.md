# Understanding ROCS and Docking Components in Reinvent Scoring


This document provides a simplified explanation of the following scoring components used in Reinvent's reinforcement learning:

- `rocs_similarity`
- `parallel_rocs_similarity`
- `dockstream`
- `docked_parallel_rocs_similarity`

## 1. ROCS Similarity Components

### 1.1 `rocs_similarity`

**What it is:** A scoring component that evaluates the 3D shape and pharmacophore (color) similarity between generated molecules and reference structures.

**How it works:**
- Uses OpenEye's ROCS (Rapid Overlay of Chemical Structures) technology
- Generates 3D conformers of molecules using OpenEye's OMEGA
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

### 1.2 `parallel_rocs_similarity`

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

**Example configuration:**
```json
{
    "component_type": "parallel_rocs_similarity",
    "name": "Parallel ROCS Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "rocs_input": "path/to/reference.sdf",
        "input_type": "sdf",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "sim_measure": "tanimoto",
        "max_num_cpus": 8,
        "save_rocs_overlays": true,
        "rocs_overlays_dir": "overlays",
        "rocs_overlays_prefix": "mol_"
    }
}
```

## 2. Docking Components

### 2.1 `dockstream`

**What it is:** A scoring component that evaluates molecules based on their docking scores against a protein target.

**How it works:**
- Interfaces with the DockStream package, which provides a unified API for multiple docking engines
- Generates 3D conformers of molecules
- Docks them into protein binding sites
- Returns docking scores (lower scores typically indicate better binding)

**Key parameters:**
- `configuration_path`: Path to DockStream configuration file
- `docker_script_path`: Path to DockStream docker script
- `environment_path`: Path to Python environment for DockStream
- `debug`: Whether to run in debug mode (default: false)
- `transformation`: Score transformation parameters

**Example configuration:**
```json
{
    "component_type": "dockstream",
    "name": "DockStream",
    "weight": 1.0,
    "specific_parameters": {
        "configuration_path": "configs/dockstream_config.json",
        "docker_script_path": "path/to/docker_script.py",
        "environment_path": "path/to/environment",
        "debug": false,
        "transformation": {
            "transformation_type": "reverse_sigmoid",
            "low": -12,
            "high": -6,
            "k": 0.5
        }
    }
}
```

### 2.2 `docked_parallel_rocs_similarity`

**Note:** Based on the codebase analysis, this component does not appear to be explicitly implemented as a separate class. It likely refers to a combination of docking and ROCS similarity approaches.

**What it likely represents:**
- A workflow that first docks molecules to a protein target
- Then performs ROCS similarity calculations on the docked poses
- This would allow for evaluating both binding affinity and 3D similarity simultaneously

**How it might be used:**
- First dock molecules using `dockstream`
- Then use the docked poses as input for `parallel_rocs_similarity`
- Combine the scores to favor molecules that both dock well and match a reference shape

## 3. Usage in Reinforcement Learning

These components are typically used in the scoring function of a reinforcement learning run to guide the model toward generating molecules with desired 3D properties:

1. **Define the scoring function** in your configuration file:
   ```json
   "scoring_function": {
       "name": "custom_product",
       "parallel": false,
       "parameters": [
           {
               "component_type": "parallel_rocs_similarity",
               "name": "Shape Similarity",
               "weight": 1.0,
               "specific_parameters": {...}
           },
           {
               "component_type": "dockstream",
               "name": "Docking Score",
               "weight": 1.0,
               "specific_parameters": {...}
           }
       ]
   }
   ```

2. **Run reinforcement learning** to optimize molecules against these criteria
3. **Analyze results** to find molecules with good 3D shape matching and docking scores

## 4. Requirements

- OpenEye license for ROCS and OMEGA (commercial software)
- DockStream installation for docking components
- Appropriate configuration files for each component

## 5. Alternatives

If you don't have access to OpenEye tools, consider:
- RDKit's shape alignment and scoring for 3D similarity
  - `rdkit_shape_similarity`: A component that uses RDKit for 3D shape similarity
  - `parallel_rdkit_shape_similarity`: A parallel version for faster processing

### 5.1 `rdkit_shape_similarity`

**What it is:** An open-source alternative to ROCS similarity that uses RDKit for 3D shape comparison.

**How it works:**
- Uses RDKit's conformer generation and shape alignment capabilities
- Supports multiple methods for shape comparison (USRCAT and O3A)
- Calculates similarity scores between generated molecules and reference structures

**Key parameters:**
- `reference_file`: Path to reference structure file (SDF)
- `method`: Shape comparison method to use ("usrcat" or "o3a")
- `shape_weight`: Weight given to shape similarity (typically 0.5)
- `color_weight`: Weight given to pharmacophore similarity (typically 0.5)
- `max_confs`: Maximum number of conformers to generate (default: 50)
- `ewindow`: Energy window for conformer selection (default: 10)
- `max_stereo`: Maximum number of stereoisomers to consider (default: 0)

**Example configuration:**
```json
{
    "component_type": "rdkit_shape_similarity",
    "name": "RDKit Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "method": "usrcat",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "max_confs": 50
    }
}
```

### 5.2 `parallel_rdkit_shape_similarity`

**What it is:** A parallel version of the RDKit shape similarity component for faster processing of large datasets.

**How it works:**
- Same as `rdkit_shape_similarity` but distributes calculations across multiple CPU cores
- Can optionally save overlay structures for visualization

**Additional parameters:**
- `max_num_cpus`: Maximum number of CPU cores to use (default: 4)
- `save_overlays`: Whether to save overlay structures (default: false)
- `overlays_dir`: Directory to save overlay structures (default: "overlays")
- `overlay_prefix`: Prefix for overlay filenames (default: "mol_")

**Example configuration:**
```json
{
    "component_type": "parallel_rdkit_shape_similarity",
    "name": "Parallel RDKit Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "method": "o3a",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "max_confs": 50,
        "max_num_cpus": 8,
        "save_overlays": true,
        "overlays_dir": "my_overlays",
        "overlay_prefix": "compound_"
    }
}
```

## 6. Implementation Details

### 6.1 RDKit-based Conformer Generation

The code below shows how RDKit is used to generate conformers in the `RDKitConformerGenerator` class:

```python
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers

class RDKitConformerGenerator:
    """RDKit-based conformer generator to replace OpenEye's OMEGA."""

    def __init__(self, max_confs=200, energy_window=10, max_stereo=0):
        self.max_confs = max_confs
        self.energy_window = energy_window
        self.max_stereo = max_stereo

    def generate_conformers(self, smiles):
        """Generate conformers for a molecule from SMILES."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Handle stereochemistry if needed
        if self.max_stereo > 0:
            isomers = self._enumerate_stereoisomers(mol)
            if not isomers:
                return None

            # Generate conformers for each stereoisomer
            all_confs_mol = None
            for i, isomer in enumerate(isomers):
                confs = self._generate_confs_for_mol(isomer)
                if confs and confs.GetNumConformers() > 0:
                    if all_confs_mol is None:
                        all_confs_mol = confs
                    else:
                        # Add conformers from this isomer to the main molecule
                        for conf in confs.GetConformers():
                            all_confs_mol.AddConformer(conf, assignId=True)
            return all_confs_mol
        else:
            # Generate conformers for the molecule without stereochemistry enumeration
            return self._generate_confs_for_mol(mol)

    def _generate_confs_for_mol(self, mol):
        """Generate conformers for a single molecule."""
        # Set up ETKDG parameters
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        params.enforceChirality = True

        # Generate conformers
        confs = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=self.max_confs,
            params=params
        )

        if confs == -1:  # Failed to generate any conformers
            return None

        # Energy minimize conformers
        for conf_id in range(mol.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(mol, confId=conf_id)

        return mol

    def _enumerate_stereoisomers(self, mol):
        """Enumerate stereoisomers of a molecule."""
        opts = EnumerateStereoisomers.StereoEnumerationOptions(
            maxIsomers=self.max_stereo,
            onlyUnassigned=False,
            unique=True
        )
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
        return isomers
```

### 6.2 RDKit-based Shape Similarity Calculation

The core RDKit shape similarity calculation is implemented in the `RDKitShapeSimilarity` class:

```python
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT
from rdkit.Chem import ChemicalFeatures

class RDKitShapeSimilarity(BaseROCSComponent):
    """RDKit-based shape similarity calculator to replace OpenEye's ROCS."""

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        # Extract parameters
        self.shape_weight = self.parameters.specific_parameters.get("shape_weight", 0.5)
        self.color_weight = self.parameters.specific_parameters.get("color_weight", 0.5)
        self.method = self.parameters.specific_parameters.get("method", "usrcat")

        # Load reference molecule
        ref_file = self.parameters.specific_parameters.get("reference_file")
        self.reference_mol = self._load_reference(ref_file)

        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=self.parameters.specific_parameters.get("max_confs", 200),
            energy_window=self.parameters.specific_parameters.get("ewindow", 10),
            max_stereo=self.parameters.specific_parameters.get("max_stereo", 0)
        )

    def _calculate_omega_score(self, smiles, step=-1) -> np.array:
        """Calculate shape similarity scores for a list of SMILES."""
        scores = []

        for smile in smiles:
            # Generate conformers for query molecule
            query_mol = self.conformer_generator.generate_conformers(smile)
            if query_mol is None or query_mol.GetNumConformers() == 0:
                scores.append(0.0)
                continue

            # Calculate similarity based on method
            if self.method == "usrcat":
                best_score = self._calculate_usrcat_similarity(query_mol)
            elif self.method == "o3a":
                best_score = self._calculate_o3a_similarity(query_mol)
            else:
                best_score = 0.0

            scores.append(best_score)

        return np.array(scores)

    def _calculate_usrcat_similarity(self, query_mol):
        """Calculate USRCAT similarity between query and reference molecules."""
        best_score = 0.0

        # Calculate for all conformer pairs
        for q_conf_id in range(query_mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                # Get USRCAT descriptors
                query_descriptor = GetUSRCAT(query_mol, confId=q_conf_id)
                ref_descriptor = GetUSRCAT(self.reference_mol, confId=r_conf_id)

                # Calculate shape similarity
                shape_sim = GetUSRScore(query_descriptor, ref_descriptor)

                # For color/pharmacophore similarity, we could use feature-based USRCAT components
                # or implement a custom pharmacophore matching algorithm
                color_sim = self._calculate_feature_similarity(query_mol, self.reference_mol,
                                                              q_conf_id, r_conf_id)

                # Combine scores
                combined_score = ((self.shape_weight * shape_sim) +
                                 (self.color_weight * color_sim)) / (self.shape_weight + self.color_weight)

                if combined_score > best_score:
                    best_score = combined_score

        return best_score

    def _calculate_o3a_similarity(self, query_mol):
        """Calculate O3A-based similarity between query and reference molecules."""
        best_score = 0.0

        # Calculate for all conformer pairs
        for q_conf_id in range(query_mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                # Create O3A alignment
                pyO3A = AllChem.GetO3A(query_mol, self.reference_mol,
                                      confId1=q_conf_id, confId2=r_conf_id)

                # Get alignment score (shape similarity)
                shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range

                # Align molecules
                pyO3A.Align()

                # Calculate feature similarity after alignment
                color_sim = self._calculate_feature_similarity(query_mol, self.reference_mol,
                                                              q_conf_id, r_conf_id)

                # Combine scores
                combined_score = ((self.shape_weight * shape_sim) +
                                 (self.color_weight * color_sim)) / (self.shape_weight + self.color_weight)

                if combined_score > best_score:
                    best_score = combined_score

        return best_score

    def _calculate_feature_similarity(self, mol1, mol2, conf_id1=0, conf_id2=0):
        """Calculate pharmacophore feature similarity between aligned molecules."""
        # This is a simplified implementation
        # A more sophisticated version would use RDKit's pharmacophore features

        # Get feature factories
        factory = ChemicalFeatures.BuildFeatureFactory()

        # Get features for both molecules
        feats1 = factory.GetFeaturesForMol(mol1, confId=conf_id1)
        feats2 = factory.GetFeaturesForMol(mol2, confId=conf_id2)

        # Calculate feature overlap
        # (This is a simplified approach - a real implementation would be more complex)
        similarity = calculate_feature_overlap(feats1, feats2)

        return similarity

    def _load_reference(self, file_path):
        """Load reference molecule from file."""
        # Implementation depends on file format (SDF, MOL, etc.)
        # Here's a simple example for SDF
        suppl = Chem.SDMolSupplier(file_path)
        ref_mol = next(suppl)

        if ref_mol is None:
            raise ValueError(f"Could not load reference molecule from {file_path}")

        # Generate conformers if needed
        if ref_mol.GetNumConformers() == 0:
            ref_mol = self.conformer_generator.generate_conformers(
                Chem.MolToSmiles(ref_mol)
            )

        return ref_mol
```

### 6.3 Registering New Components

To use the new RDKit-based components, you need to register them in the factory:

```python
# In score_component_factory.py
def _deafult_scoring_component_registry(self) -> dict:
    enum = ScoringFunctionComponentNameEnum()
    component_map = {
        # Existing components...
        enum.ROCS_SIMILARITY: RocsSimilarity,
        enum.PARALLEL_ROCS_SIMILARITY: ParallelRocsSimilarity,

        # Add new RDKit-based components
        enum.RDKIT_SHAPE_SIMILARITY: RDKitShapeSimilarity,
        enum.PARALLEL_RDKIT_SHAPE_SIMILARITY: ParallelRDKitShapeSimilarity,

        # Other components...
    }
    return component_map
```

### 6.4 Updating the Component Enum

```python
# In scoring_function_component_enum.py
@dataclass(frozen=True)
class ScoringFunctionComponentNameEnum:
    # Existing components...
    ROCS_SIMILARITY = "rocs_similarity"
    PARALLEL_ROCS_SIMILARITY = "parallel_rocs_similarity"

    # New RDKit-based components
    RDKIT_SHAPE_SIMILARITY = "rdkit_shape_similarity"
    PARALLEL_RDKIT_SHAPE_SIMILARITY = "parallel_rdkit_shape_similarity"

    # Other components...
```

### 6.5 Using RDKit Components in Configuration

```json
{
    "component_type": "rdkit_shape_similarity",
    "name": "RDKit Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "method": "usrcat",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "max_confs": 50
    }
}
```

By implementing these RDKit-based alternatives, Link-INVENT can operate without requiring OpenEye licenses while still providing similar functionality for 3D shape-based molecular design.


