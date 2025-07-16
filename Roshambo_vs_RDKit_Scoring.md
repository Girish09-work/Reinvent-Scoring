# Roshambo vs RDKit Shape Similarity Scoring Components

This document explains the flow of scoring linker SMILES in both the Roshambo and RDKit shape similarity scoring components, highlighting their differences and implementation details.

## 1. Input Processing Flow

Both components take SMILES strings as input, typically representing linker molecules in the context of Protac-Invent. The general flow is:

1. SMILES strings are passed to the scoring component
2. Conformers are generated for each molecule
3. Shape similarity is calculated against reference molecules
4. Scores are returned for each input molecule

## 2. Conformer Generation

### 2.1 RoshamboConformerGenerator

The `RoshamboConformerGenerator` class handles conformer generation for the Roshambo component:

```python
def create_temp_sdf_from_smiles(self, smiles_list: List[str], temp_dir: Optional[str] = None) -> str:
    """
    Create a temporary SDF file from SMILES with conformers.
    
    Args:
        smiles_list: List of SMILES strings
        temp_dir: Directory for temporary file
        
    Returns:
        Path to temporary SDF file
    """
    if temp_dir is None:
        temp_dir = tempfile.gettempdir()
    
    Path(temp_dir).mkdir(parents=True, exist_ok=True)
    temp_file = os.path.join(temp_dir, "temp_conformers.sdf")
    
    writer = Chem.SDWriter(temp_file)
    
    for i, smiles in enumerate(smiles_list):
        if not smiles:
            continue
            
        mol = self.generate_conformers_from_smiles(smiles)
        if mol and mol.GetNumConformers() > 0:
            # Set molecule name
            mol.SetProp("_Name", f"mol_{i}")
            
            # Write all conformers
            for conf in mol.GetConformers():
                writer.write(mol, confId=conf.GetId())
    
    writer.close()
    return temp_file
```

The conformer generation process:
1. Converts SMILES to RDKit molecules
2. Adds hydrogens if specified
3. Generates 3D conformers using methods like ETKDGv3
4. Optionally optimizes conformers with force fields (UFF, MMFF94s)
5. Writes the conformers to an SDF file

### 2.2 RDKit Shape Similarity Conformer Generation

The RDKit shape component generates conformers directly in memory:

```python
def _calculate_shape_scores(self, mols):
    """Calculate shape similarity scores for a list of RDKit molecules."""
    scores = []
    
    for mol in mols:
        if mol is None:
            scores.append(0.0)
            continue
            
        # Generate conformers if needed
        if mol.GetNumConformers() == 0:
            mol = self._generate_conformers(mol)
            if mol is None:
                scores.append(0.0)
                continue
                
        # Calculate shape similarity
        score = self._calculate_shape_similarity(mol)
        scores.append(score)
```

## 3. Score Calculation

### 3.1 Roshambo Shape Similarity Scoring

Roshambo uses GPU-accelerated Gaussian molecular shape comparison:

```python
def _run_roshambo_direct(self, ref_file: str, dataset_file: str, n_confs: int, working_dir: str) -> dict:
    """Run roshambo directly and extract scores from CSV."""
    try:
        # Execute roshambo API (always returns None due to bug, but writes CSV)
        self.get_similarity_scores(
            ref_file=ref_file,
            dataset_file=dataset_file,
            ignore_hs=self.ignore_hs,
            n_confs=n_confs,
            use_carbon_radii=self.use_carbon_radii,
            color=self.color_weight > 0,
            sort_by="ComboTanimoto",
            write_to_file=True,
            gpu_id=self.gpu_id,
            working_dir=working_dir
        )

        # Read scores from CSV file
        return self._read_scores_from_csv(working_dir)
    except Exception as e:
        if self.debug:
            print(f"Error running roshambo directly: {e}")
        return {}
```

The scoring process:
1. Roshambo compares the shape of each molecule to the reference using Gaussian functions
2. It calculates both shape and color (pharmacophore) similarity
3. The final score is a weighted combination of shape and color scores
4. Results are written to a CSV file and then read back

### 3.2 RDKit Shape Similarity Scoring

RDKit uses the Open3DAlign (O3A) method for shape comparison:

```python
def _calculate_shape_similarity(self, mol):
    """Calculate shape similarity between a molecule and reference molecules."""
    best_score = 0.0
    best_conf_id = -1
    best_ref_conf_id = -1

    # Check main reference
    if self.reference_mol is not None:
        for q_conf_id in range(mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                try:
                    # Create O3A alignment
                    pyO3A = AllChem.GetO3A(mol, self.reference_mol,
                                          confId1=q_conf_id, confId2=r_conf_id)

                    # Get alignment score (shape similarity)
                    shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range

                    if shape_sim > best_score:
                        best_score = shape_sim
                        best_conf_id = q_conf_id
                        best_ref_conf_id = r_conf_id
                except Exception:
                    continue
```

## 4. Key Differences

| Feature | Roshambo | RDKit Shape |
|---------|----------|-------------|
| **Performance** | GPU-accelerated, much faster | CPU-based, slower |
| **Method** | Gaussian shape comparison | Open3DAlign (O3A) |
| **Output Format** | Writes to SDF/CSV files | In-memory calculation |
| **Pharmacophore** | Supports color (pharmacophore) features | Shape-only by default |
| **Dependencies** | Requires Roshambo package | Built into RDKit |
| **Parallelization** | GPU parallelization | CPU threading |

## 5. Implementation Details

### 5.1 Roshambo-Specific Features

Roshambo has several unique features:
- Dynamic environment management for conda environments
- Subprocess execution for isolation
- File format conversion for compatibility
- GPU device selection

```python
def _get_similarity_scores_with_env(self, **kwargs):
    """
    Execute Roshambo using subprocess with conda environment.
    Simplified version that focuses on CSV output.
    """
    import json
    import tempfile
    import subprocess

    # Create temporary files for input/output
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as input_file:
        json.dump(kwargs, input_file)
        input_path = input_file.name
```

### 5.2 RDKit Shape Similarity Features

RDKit shape similarity is simpler but more integrated:
- Direct in-memory calculation
- No external dependencies
- Integrated with RDKit's conformer generation
- Optional overlay visualization

## 6. Conclusion

Both components serve the same purpose of calculating shape similarity between molecules, but they use different approaches:

- **Roshambo** is optimized for performance using GPU acceleration and is ideal for large-scale calculations, but requires additional setup and dependencies.

- **RDKit Shape** is simpler to use and requires no additional dependencies, but is slower and less feature-rich.

The choice between them depends on performance requirements, available hardware, and the specific needs of the project.