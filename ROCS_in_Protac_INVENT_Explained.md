# ROCS in Protac-INVENT: Reference Structures and Comparison Process

## Understanding Reference Structures in ROCS Comparison

In Protac-INVENT, ROCS (Rapid Overlay of Chemical Structures) is used to compare generated PROTAC conformers to reference structures. However, there's an important distinction to clarify about what constitutes a "reference structure" in this context.

### Reference Structure vs. PTS

- **PTS (PROTAC Ternary Structure)**: This is the complete complex consisting of:
  - The target protein (POI - Protein of Interest)
  - The E3 ligase protein
  - The PROTAC molecule binding both proteins

- **Reference Structure for ROCS**: This is typically just the PROTAC molecule (or a specific part of it) extracted from a known PTS, not the entire protein complex.

ROCS compares the 3D shape and pharmacophore features of molecules, not entire protein complexes. The protein structure information from the PTS is used in other parts of the workflow (particularly docking), but not directly in the ROCS comparison.

## How Reference Structures Are Used

In Protac-INVENT, the reference structure for ROCS comparison is typically:

1. A known effective PROTAC molecule extracted from a crystal structure
2. Sometimes just the linker portion of a known PROTAC
3. Occasionally a designed ideal linker conformation

The workflow described in the paper follows this process:

1. Extract a reference PROTAC from a PTS crystal structure
2. Generate new candidate PROTACs
3. Compare the 3D shape of candidates to the reference using ROCS
4. Score based on similarity to guide the generative model

## Example from Protac-INVENT

The paper specifically mentions using BTK PROTAC (PDB: 6W8I) as a test case:

> "In current study, Bruton's tyrosine kinase (BTK) PROTAC was chosen as the test case. Published PTS crystal structure 6W8I was chosen as the reference structure for structure generation."

In this case:
- The PTS is the entire crystal structure (6W8I)
- The reference for ROCS is the PROTAC molecule extracted from this structure

## Relationship Between Generated PROTACs and Reference Structure

The relationship between generated PROTACs and the reference structure can be:

1. **Same warheads, different linker**: Most common in Protac-INVENT, where the E3 ligand and POI warheads are kept the same, but the linker is varied.
2. **Similar warheads, different linker**: Where warheads are modified slightly but maintain key binding features.
3. **Complete novel design**: Less common, where the entire PROTAC is designed de novo.

## Practical Implementation

In the code, the reference structure is specified in the scoring component configuration:

```json
{
  "component_type": "parallel_rocs_similarity",
  "name": "ROCS Similarity",
  "weight": 1.0,
  "specific_parameters": {
    "rocs_input": "path/to/reference.sdf",  // The reference PROTAC structure
    "input_type": "sdf",
    "shape_weight": 0.5,
    "color_weight": 0.5,
    "sim_measure": "tanimoto"
  }
}
```

The `rocs_input` parameter points to an SDF file containing the reference molecule extracted from the PTS.

## Visualization of the Process

The process can be visualized as follows:

1. **PTS Crystal Structure**: Contains proteins + PROTAC
   ```
   [POI Protein]----[PROTAC]----[E3 Ligase]
   ```

2. **Extract Reference PROTAC**:
   ```
   [PROTAC] = [POI Warhead]--[Linker]--[E3 Ligand]
   ```

3. **Generate New PROTACs**:
   ```
   [New PROTAC] = [POI Warhead]--[New Linker]--[E3 Ligand]
   ```

4. **ROCS Comparison**:
   ```
   Compare: [PROTAC] vs [New PROTAC]
   ```

## Example from BTK System (6W8I)

For the BTK system mentioned in the paper:

1. The PTS is the crystal structure 6W8I containing:
   - BTK protein (the POI)
   - CRBN E3 ligase
   - A PROTAC molecule binding both

2. The reference structure for ROCS would be the PROTAC molecule extracted from this crystal structure.

3. Protac-INVENT would generate new linkers while keeping the BTK inhibitor and CRBN ligand portions the same.

4. ROCS would compare the shape and pharmacophore features of the generated PROTACs to the reference PROTAC.

## References

- PDB entry for 6W8I: [https://www.rcsb.org/structure/6W8I](https://www.rcsb.org/structure/6W8I)
- OpenEye ROCS documentation: [https://www.eyesopen.com/rocs](https://www.eyesopen.com/rocs)
- Protac-INVENT GitHub: [https://github.com/jidushanbojue/Protac-invent](https://github.com/jidushanbojue/Protac-invent)

# PTS in Protac-INVENT: DockStream Integration and Workflow

## Real Use Case of PTS in Protac-INVENT

The Protein Ternary Structure (PTS) plays a crucial role in Protac-INVENT beyond just providing a reference PROTAC structure for ROCS comparison. The PTS is directly used in the DockStream component of the workflow, which handles the molecular docking aspects of PROTAC design.

### DockStream Integration with PTS

DockStream is integrated into Protac-INVENT to evaluate how well generated PROTACs fit into the binding site of the PTS. This is a critical step because a PROTAC must not only have the right shape in isolation but must also fit properly into the binding pocket formed by the protein complex.

As stated in the paper:

> "To further evaluate generated PROTAC, its OPCs were docked into the binding pocket of PTS, which is composed by E3 ligase and POI, to obtain their docking scores."

## Inputs and Outputs for DockStream in Protac-INVENT

### Inputs to DockStream:

1. **Protein Complex Structure**: The PTS without the original PROTAC (just the protein components)
2. **Generated PROTAC Molecules**: The candidate PROTACs in SMILES format
3. **Initial Conformations**: Often derived from OMEGA or other conformer generators
4. **Docking Configuration**: Parameters for the docking process

Here's an example of DockStream configuration from the Protac-INVENT codebase:

```json
{
    "component_type": "dockstream",
    "name": "ADV",
    "specific_parameters": {
        "configuration_path": "/protac-invent/DockStream-master/result/BTK/BTK_scoring.json",
        "docker_script_path": "/protac-invent/DockStream-master/docker.py",
        "environment_path": "D:/conda_envs/DockStream",
        "transformation": {
            "high": 0,
            "k": 0.25,
            "low": -10,
            "transformation_type": "reverse_sigmoid"
        }
    },
    "weight": 1
}
```

### Outputs from DockStream:

1. **Docking Scores**: Numerical values indicating binding affinity
2. **Docked Poses**: 3D conformations of PROTACs in the binding site
3. **Post-Docking Analysis**: Additional metrics like shape similarity after docking

The paper describes a specific component for post-docking analysis:

> "Another component, *Ps*, with a value in the range of [0,1], was set to measure the shape similarity between warheads (including the ligands for E3 ligase and POI) of docking conformation and that of the reference ligand."

## Separation of Reference Structure from PTS

The separation of the reference PROTAC structure from the PTS is a critical step in the Protac-INVENT workflow. This process involves:

### 1. PTS Preparation

1. **Obtain PTS**: Start with a crystal structure of a PROTAC ternary complex (e.g., PDB: 6W8I)
2. **Extract Components**: Separate the structure into:
   - Protein components (POI and E3 ligase)
   - PROTAC molecule

### 2. Reference Structure Extraction

The reference structure extraction depends on what's being designed:

1. **Full PROTAC Reference**: Extract the entire PROTAC molecule from the PTS
2. **Linker-Only Reference**: Extract just the linker portion, removing the warheads
3. **Warhead References**: Extract the warhead portions separately

The paper describes this process:

> "Once the alignment was done, the warheads of the generated PROTAC were removed and only the linker atoms were remained. The warheads of reference ligand were then copied and merged with the linker fragment."

### 3. Practical Implementation

In practice, this separation is typically done using molecular modeling software like PyMOL, Chimera, or programmatically with libraries like RDKit or OpenEye. The extracted components are saved as separate files:

- Protein complex for docking (.pdb)
- Reference PROTAC for ROCS (.sdf)
- Warheads for merging with generated linkers (.sdf)

## Example Workflow with BTK (6W8I)

For the BTK system (PDB: 6W8I) mentioned in the paper:

1. **PTS Preparation**:
   - Download 6W8I PDB file
   - Extract protein components (BTK and CRBN)
   - Extract PROTAC molecule

2. **Reference Structure Extraction**:
   - Save complete PROTAC as reference.sdf
   - Identify attachment points for linker
   - Save warheads separately as warhead_poi.sdf and warhead_e3.sdf

3. **DockStream Configuration**:
   - Prepare protein complex for docking
   - Configure docking parameters for AutoDock Vina
   - Set up scoring function to evaluate docked poses

4. **Integration in Scoring Function**:
   - Use ROCS to compare generated molecules to reference
   - Use DockStream to dock molecules into PTS
   - Evaluate post-docking shape similarity

This workflow is reflected in the configuration files found in the Protac-INVENT repository, such as:

```
/protac-invent/DockStream-master/result/BTK/BTK_scoring.json
```

## Visual Representation of the Process

```
PTS Crystal Structure (6W8I)
        |
        v
+-------------------+
| Extract Components |
+-------------------+
        |
        v
+---------------------+     +----------------------+
| Protein Complex     |     | Reference PROTAC     |
| (for docking)       |     | (for ROCS)           |
+---------------------+     +----------------------+
        |                           |
        v                           v
+---------------------+     +----------------------+
| DockStream          |     | ROCS Similarity      |
| (evaluate binding)  |     | (evaluate shape)     |
+---------------------+     +----------------------+
        |                           |
        v                           v
+---------------------------------------------+
| Combined Score for Reinforcement Learning   |
+---------------------------------------------+
```

This integrated approach allows Protac-INVENT to generate PROTACs that not only have appropriate 3D shapes but also fit well into the binding site of the protein complex, increasing the likelihood of successful ternary complex formation.

## Code for PTS and Reference Structure Separation

After examining the repository, there isn't a dedicated, comprehensive script specifically for separating PTS and reference structures in the Protac-INVENT codebase. This separation appears to be handled through a combination of standard molecular modeling tools and manual preparation steps before the main workflow begins.

However, there are several relevant code snippets that show how the separated components are used:

### 1. Shape Similarity Calculation

The `compute_shape_between_docking_and_reference.py` script shows how the reference structure is used for shape comparison after docking:

```python
def shape_scoring(ref_sdf, warheads_sdf):
    reffs = oechem.oemolistream(ref_sdf)
    fitfs = oechem.oemolistream(warheads_sdf)

    refmol = oechem.OEGraphMol()
    oechem.OEReadMolecule(reffs, refmol)

    prep = oeshape.OEOverlapPrep()
    prep.Prep(refmol)

    shapeFunc = oeshape.OEExactShapeFunc()
    shapeFunc.SetupRef(refmol)

    res = oeshape.OEOverlapResults()

    for idx, fitmol in enumerate(fitfs.GetOEGraphMols()):
        if idx != 0:
            continue
        prep.Prep((fitmol))
        shapeFunc.Overlap(fitmol, res)

        shape_score = res.GetTanimoto()
        print(shape_score)
        return shape_score
```

This function takes a reference SDF file and compares it to warheads from docked structures, calculating a shape Tanimoto score.

### 2. Target Preparation

The DockStream codebase includes target preparation functionality that can work with reference ligands:

```python
def _cavity_by_reference(self, receptor_mol2_path: str, prm_path: str, folder: str):
    # load the reference file (PDB or SDF)
    ref_format = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_FORMAT].upper()
    if ref_format == self._TP.CAVITY_REFERENCE_FORMAT_PDB:
        ref_mol = Chem.MolFromPDBFile(self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH],
                                      sanitize=True)
    elif ref_format == self._TP.CAVITY_REFERENCE_FORMAT_SDF:
        mol_supplier = Chem.SDMolSupplier(self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH])
        for mol in mol_supplier:
            ref_mol = mol
    else:
        raise TargetPreparationFailed("Specified format not supported.")
```

This code from `rDock_target_preparator.py` shows how a reference ligand (which could be the extracted PROTAC) is loaded to define the binding cavity for docking.

### 3. OpenEye Target Preparation

Similarly, the OpenEye target preparator uses reference structures:

```python
ref_istream = oechem.oemolistream()
if not ref_istream.open(self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH]):
    raise TargetPreparationFailed("Could not open reference file.")
ref = oechem.OEGraphMol()
oechem.OEReadMolecule(ref_istream, ref)
oedocking.OEMakeReceptor(target, self._protein, ref)
```

This code from `OpenEye_target_preparator.py` shows how a reference molecule is used to define the binding site.

### 4. Setting References for Ligand Preparation

The ligand preparator has functionality to set reference structures:

```python
def set_references(self, references):
    # usually, references are loaded from files; but this function allows setting them as a list of molecules
    if references is not None:
        if not isinstance(references, list):
            references = [references]
    self._references = references
```

This method from `ligand_preparator.py` allows setting reference molecules for alignment or comparison.

### Typical Workflow for Separation

Based on the code and documentation, the typical workflow for separating PTS and reference structures would involve:

1. Using external tools like PyMOL, Chimera, or RDKit to extract the PROTAC from the PTS
2. Saving the extracted PROTAC as an SDF file
3. Preparing the protein complex separately for docking
4. Configuring the workflow to use these separated files

For example, a user might:
1. Use PyMOL to open the PTS PDB file
2. Select and save the PROTAC molecule as `reference.sdf`
3. Save the protein complex without the PROTAC as `protein_complex.pdb`
4. Configure DockStream to use `protein_complex.pdb` as the receptor and `reference.sdf` as the reference

While there isn't a dedicated script for this separation process in the repository, it's a standard preprocessing step in structure-based drug design workflows.

