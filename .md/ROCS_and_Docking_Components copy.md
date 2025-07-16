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
- Open-source docking engines supported by DockStream (AutoDock Vina, etc.)

This document provides a simplified explanation of the following scoring components used in Reinvent's reinforcement learning:

- `rocs_similarity`
- `parallel_rocs_similarity`
- `dockstream`
- `docked_parallel_rocs_similarity`

## 1. What does ROCS/Omega do in Protac-INVENT and Link-INVENT?

In Protac-INVENT and Link-INVENT, ROCS and OMEGA serve critical functions for 3D molecular shape analysis and conformer generation:

### OMEGA's Role:
- Generates 3D conformers of molecules designed by the generative models
- Creates multiple low-energy conformations for each molecule
- Handles stereochemistry enumeration when needed
- Provides realistic 3D structures for shape comparison and docking

### ROCS's Role:
- Evaluates 3D shape similarity between generated molecules and reference structures
- Aligns molecules to reference compounds based on shape and chemical features
- Calculates similarity scores that guide the reinforcement learning process
- Helps identify molecules with similar binding modes to known active compounds
- Enables scaffold hopping by finding molecules with similar 3D properties but different 2D structures

These tools are particularly important in Link-INVENT, where the 3D arrangement of the linker between two binding warheads is critical for proper protein-protein interaction.

## 2. How does it work? (Implementation Details)

### OMEGA Conformer Generation

The code below shows how OMEGA is used to generate conformers in the `ParallelRocsSimilarity` component:

```python
@classmethod
def setup_omega(cls, erange, max_confs):
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetStrictStereo(False)
    omegaOpts.SetEnergyWindow(erange)
    omegaOpts.SetMaxConfs(max_confs)
    cls.omega = oeomega.OEOmega(omegaOpts)
    return cls.omega
```

And here's how conformers are generated for each molecule:

```python
def get_omega_confs(imol, omega, enum_stereo, max_stereo):
    stereo = False
    no_stereo = False
    if enum_stereo:
        enantiomers = list(oeomega.OEFlipper(imol.GetActive(), max_stereo, False, True))
        for k, enantiomer in enumerate(enantiomers):
            enantiomer = oechem.OEMol(enantiomer)
            ret_code = omega.Build(enantiomer)
            if ret_code == oeomega.OEOmegaReturnCode_Success:
                if k == 0:
                    imol = oechem.OEMol(enantiomer.SCMol())
                    imol.DeleteConfs()
                stereo = True
                for x in enantiomer.GetConfs():
                    imol.NewConf(x)
    else:
        no_stereo = omega(imol)
    return no_stereo or stereo, imol
```

### ROCS Shape Similarity Calculation

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

The parallel version distributes these calculations across multiple CPU cores:

```python
def _calculate_omega_score(self, smiles, step) -> np.array:
    inputs = []
    if len(smiles) == 0:
        return np.array(())
    self._prepare_overlay()
    ind = str(step).zfill(4)

    for smile in smiles:
        input = {"smile": smile, "shape_weight": self.shape_weight, "color_weight": self.color_weight,
                 "sim_func_name_set": self.sim_func_name_set, "batch_id": ind,
                 "enum_stereo": self.enum_stereo, "max_stereo": self.max_stereo, "save_overlays": self.save_overlays,
                 "neg_vol_file": self.protein_file, "neg_vol_lig": self.ligand_file
                 }
        inputs.append(input)
    with Pool(processes=min(self.num_cpus, len(inputs))) as pool:
        results = pool.map(self._unfold, inputs)

    scores = []
    # Process results...
    return np.array(scores)
```

## 3. What are the open-source alternatives?

Several open-source alternatives to OpenEye's ROCS and OMEGA are available through RDKit:

### Alternatives to OMEGA (Conformer Generation):

1. **RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry)**
   - Modern conformer generation algorithm with knowledge-based torsion angles
   - Provides high-quality 3D conformers comparable to OMEGA
   - Supports stereochemistry and macrocycle handling
   - Implementation example:
   ```python
   from rdkit import Chem
   from rdkit.Chem import AllChem

   def generate_conformers(mol, num_confs=200, energy_window=10):
       mol = Chem.AddHs(mol)  # Add hydrogens
       params = AllChem.ETKDGv3()
       params.randomSeed = 42
       params.useSmallRingTorsions = True
       params.useMacrocycleTorsions = True
       params.enforceChirality = True
       confs = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
       return mol
   ```

2. **RDKit's Stereochemistry Enumeration**
   - For handling stereochemistry similar to OMEGA:
   ```python
   from rdkit.Chem import EnumerateStereoisomers

   def enumerate_stereo(mol, max_isomers=10):
       opts = EnumerateStereoisomers.StereoEnumerationOptions(
           maxIsomers=max_isomers,
           onlyUnassigned=False,
           unique=True
       )
       isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
       return isomers
   ```

### Alternatives to ROCS (Shape Similarity):

1. **USR and USRCAT (Ultrafast Shape Recognition)**
   - Fast method for shape-based similarity without explicit alignment
   - USRCAT extends USR with chemical feature information
   - Implementation example:
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

2. **RDKit's Shape Alignment**
   - For more ROCS-like explicit alignment and scoring:
   ```python
   from rdkit.Chem import AllChem

   def align_molecules(query_mol, ref_mol):
       # Align query molecule to reference
       pyO3A = AllChem.GetO3A(query_mol, ref_mol)
       score = pyO3A.Score()
       pyO3A.Align()
       return score, query_mol
   ```

## 4. How can we implement these alternatives on top of Link-INVENT?

To implement RDKit-based alternatives in Link-INVENT, we need to create new scoring components that replace the OpenEye functionality:

### Step 1: Create RDKit-based conformer generator class

```python
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

### Step 2: Create RDKit-based shape similarity component

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

### Step 3: Register the new components in the factory

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

### Step 4: Update the component enum

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

### Step 5: Use in Link-INVENT configuration

```json
{
    "component_type": "rdkit_shape_similarity",
    "name": "RDKit Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "rocs_input": "path/to/reference.sdf",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "method": "usrcat",
        "max_confs": 200,
        "max_stereo": 1
    }
}
```

By implementing these RDKit-based alternatives, Link-INVENT can operate without requiring OpenEye licenses while still providing similar functionality for 3D shape-based molecular design.
