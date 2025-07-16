# Detailed Code Explanation of Scoring Components

This document provides a line-by-line explanation of the scoring components you're interested in, focusing on what each component generates and the parameters required.

## 1. `rocs_similarity`

The `RocsSimilarity` class evaluates 3D shape and pharmacophore similarity between generated molecules and reference structures.

### Initialization (Lines 13-27)
```python
def __init__(self, parameters: ComponentParameters):
    super().__init__(parameters)
    self.sim_measure_enum = ROCSSimilarityMeasuresEnum()
    self.input_types_enum = ROCSInputFileTypesEnum()
    self.param_names_enum = ROCSSpecificParametersEnum()
    self.shape_weight = self.parameters.specific_parameters[self.param_names_enum.SHAPE_WEIGHT]
    self.color_weight = self.parameters.specific_parameters[self.param_names_enum.COLOR_WEIGHT]
    self.sim_func_name_set = self.__get_similarity_name_set()
    cff_path = self.parameters.specific_parameters.get(self.param_names_enum.CUSTOM_CFF, None)
    self.prep = self.__set_prep(cff_path)
    self.overlay = self.__prepare_overlay(self.parameters.specific_parameters[self.param_names_enum.ROCS_INPUT],
                                          self.parameters.specific_parameters[self.param_names_enum.INPUT_TYPE])
    self.omega = self.__setup_omega()
    oechem.OEThrow.SetLevel(10000)
```

- **What it does**: Initializes the ROCS similarity component with the necessary parameters
- **Required parameters**:
  - `SHAPE_WEIGHT`: Weight given to shape similarity (typically 0.5)
  - `COLOR_WEIGHT`: Weight given to pharmacophore similarity (typically 0.5)
  - `ROCS_INPUT`: Path to reference structure file (SDF or shape query)
  - `INPUT_TYPE`: Type of reference file ("sdf" or "shape_query")
  - `CUSTOM_CFF` (optional): Path to custom color force field file

### Score Calculation (Lines 29-46)
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

- **What it generates**: An array of similarity scores (0-1) for each input molecule
- **Process**:
  1. Converts SMILES to OpenEye molecule
  2. Generates 3D conformers using OMEGA
  3. Prepares molecule for shape comparison
  4. Aligns molecule to reference structure
  5. Calculates shape and color (pharmacophore) similarity
  6. Combines scores using weighted average

### Similarity Measure Selection (Lines 75-89)
```python
def __get_similarity_name_set(self):
    similarity_collection_name = self.parameters.specific_parameters.get(self.param_names_enum.SIM_MEASURE,
                                                                    self.sim_measure_enum.TANIMOTO)
    similarity_collection = self.__similarity_collection(similarity_collection_name)
    return similarity_collection

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

- **What it does**: Selects the similarity measure to use (Tanimoto or Tversky)
- **Required parameters**:
  - `SIM_MEASURE` (optional): Similarity measure to use (default: "Tanimoto")
  - Options: "Tanimoto", "RefTversky", "FitTversky"

## 2. `parallel_rocs_similarity`

The `ParallelRocsSimilarity` class is a parallel version of ROCS similarity that distributes calculations across multiple CPU cores.

### Initialization (Lines 16-31)
```python
def __init__(self, parameters: ComponentParameters):
    super().__init__(parameters)
    avail_cpus = multiprocessing.cpu_count()
    oechem.OEThrow.SetLevel(10000)
    self.sim_measure_enum = ROCSSimilarityMeasuresEnum()
    self.input_types_enum = ROCSInputFileTypesEnum()
    self.param_names_enum = ROCSSpecificParametersEnum()

    self.num_cpus = min(avail_cpus, self._specific_param("MAX_CPUS"))
    self._set_omega_parameters()
    self._set_rocs_parameters()
    self.shape_weight = self._specific_param("SHAPE_WEIGHT")
    self.color_weight = self._specific_param("COLOR_WEIGHT")
    self.sim_func_name_set = oefuncs.get_similarity_name_set(parameters, self.param_names_enum,
                                                             self.sim_measure_enum)
```

- **What it does**: Initializes the parallel ROCS similarity component
- **Required parameters**:
  - `MAX_CPUS` (optional): Maximum number of CPU cores to use (default: 4)
  - `SHAPE_WEIGHT`: Weight given to shape similarity
  - `COLOR_WEIGHT`: Weight given to pharmacophore similarity
  - Plus parameters for OMEGA and ROCS (see below)

### OMEGA Parameters (Lines 33-40)
```python
def _set_omega_parameters(self):
    self.max_confs = self._specific_param("MAX_CONFS")
    self.erange = self._specific_param("EWINDOW")
    self.enum_stereo = self._specific_param("ENUM_STEREO")
    self.max_stereo = self._specific_param("MAX_STEREO")
    if self.max_stereo == 0:
        self.enum_stereo = False
    self.setup_omega(self.erange, self.max_confs)
```

- **What it does**: Sets up OMEGA conformer generation parameters
- **Required parameters**:
  - `MAX_CONFS` (optional): Maximum number of conformers to generate (default: 200)
  - `EWINDOW` (optional): Energy window for conformer selection in kcal/mol (default: 10)
  - `ENUM_STEREO` (optional): Whether to enumerate stereoisomers (default: False)
  - `MAX_STEREO` (optional): Maximum number of stereoisomers to generate (default: 0)

### ROCS Parameters (Lines 42-57)
```python
def _set_rocs_parameters(self):
    self.file_path = self._specific_param("ROCS_INPUT")
    self.file_type = self._specific_param("INPUT_TYPE")
    self.cff_path = self._specific_param("CUSTOM_CFF")
    self.save_overlays = self._specific_param("SAVE_ROCS_OVERLAYS")
    if self.save_overlays:
        self.dir_name = self._specific_param("ROCS_OVERLAYS_DIR")
        self.overlay_prefix = self._specific_param("ROCS_OVERLAYS_PREFIX")
        Path(self.dir_name).mkdir(parents=True, exist_ok=True)

    self.protein_file = ""
    self.ligand_file = ""
    self.neg_vol = self._specific_param("NEGATIVE_VOLUME")
    if self.neg_vol:
        self.protein_file = self._specific_param("PROTEIN_NEG_VOL_FILE")
        self.ligand_file = self._specific_param("LIGAND_NEG_VOL_FILE")
```

- **What it does**: Sets up ROCS shape comparison parameters
- **Required parameters**:
  - `ROCS_INPUT`: Path to reference structure file
  - `INPUT_TYPE`: Type of reference file ("sdf" or "shape_query")
  - `CUSTOM_CFF` (optional): Path to custom color force field
  - `SAVE_ROCS_OVERLAYS` (optional): Whether to save overlay files (default: False)
  - `ROCS_OVERLAYS_DIR` (optional): Directory to save overlays (required if SAVE_ROCS_OVERLAYS is True)
  - `ROCS_OVERLAYS_PREFIX` (optional): Prefix for overlay filenames (default: "")
  - `NEGATIVE_VOLUME` (optional): Whether to use negative volume constraints (default: False)
  - `PROTEIN_NEG_VOL_FILE` (optional): Protein file for negative volume (required if NEGATIVE_VOLUME is True)
  - `LIGAND_NEG_VOL_FILE` (optional): Ligand file for negative volume (required if NEGATIVE_VOLUME is True)

### Parallel Score Calculation (Lines 59-86)
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
    if self.save_overlays:
        overlay_filename = self.overlay_prefix + ind + ".sdf"
        overlay_file_path = os.path.join(self.dir_name, overlay_filename)
        outfs = oechem.oemolostream(overlay_file_path)
    for result in results:
        score, outmol = result
        scores.append(score)
        if self.save_overlays:
            oechem.OEWriteMolecule(outfs, outmol)
    return np.array(scores)
```

- **What it generates**: An array of similarity scores (0-1) for each input molecule
- **Process**:
  1. Prepares the reference structure for overlay
  2. Creates input dictionaries for each SMILES
  3. Distributes calculations across multiple CPU cores
  4. Collects and processes results
  5. Optionally saves overlay structures to SDF file

## 3. `dockstream`

The `DockStream` class evaluates molecules based on their docking scores against a protein target.

### Initialization (Lines 10-15)
```python
def __init__(self, parameters: ComponentParameters):
    super().__init__(parameters)
    self._configuration_path = self.parameters.specific_parameters[self.component_specific_parameters.DOCKSTREAM_CONFPATH]
    self._docker_script_path = self.parameters.specific_parameters[self.component_specific_parameters.DOCKSTREAM_DOCKERSCRIPTPATH]
    self._environment_path = self.parameters.specific_parameters[self.component_specific_parameters.DOCKSTREAM_ENVPATH]
```

- **What it does**: Initializes the DockStream component
- **Required parameters**:
  - `DOCKSTREAM_CONFPATH`: Path to DockStream configuration file
  - `DOCKSTREAM_DOCKERSCRIPTPATH`: Path to DockStream docker script
  - `DOCKSTREAM_ENVPATH`: Path to Python environment for DockStream

### Command Creation (Lines 23-35)
```python
def _create_command(self, smiles: List[str], step):
    concat_smiles = '"' + ';'.join(smiles) + '"'
    command = ' '.join([self._environment_path,
                        self._docker_script_path,
                        "-conf", self._configuration_path,
                        "-output_prefix", self._get_step_string(step),
                        "-smiles", concat_smiles,
                        "-print_scores"])

    # check, if DockStream is to be executed in debug mode, which will cause its loggers to print out
    # much more detailed information
    command = self._add_debug_mode_if_selected(command)
    return command
```

- **What it does**: Creates the command to run DockStream
- **Required parameters**:
  - `DOCKSTREAM_DEBUG` (optional): Whether to run in debug mode (default: False)

### Score Calculation (Lines 37-59)
```python
def _calculate_score(self, smiles: List[str], step) -> np.array:
    # create the external command
    command = self._create_command(smiles, step)

    # send the batch smiles and retrieve the result as a list of strings
    results = self._send_request_with_stepwize_read(command, len(smiles))

    # note: some ligands might have failed in DockStream (embedding or docking) although they are valid
    #       RDkit molecules -> "docker.py" will return "NA"'s for failed molecules, as '0' could be a perfectly
    #       normal value; anything that cannot be cast to a floating point number will result in '0'
    scores = []
    for score in results:
        try:
            score = float(score)
        except ValueError:
            score = 0
        scores.append(score)
    transform_params = self.parameters.specific_parameters.get(
        self.component_specific_parameters.TRANSFORMATION, {}
    )
    transformed_scores = self._transformation_function(scores, transform_params)

    return np.array(transformed_scores), np.array(scores)
```

- **What it generates**: Two arrays - transformed docking scores and raw docking scores
- **Process**:
  1. Creates and executes the DockStream command
  2. Reads the results (docking scores)
  3. Handles any failed dockings (sets score to 0)
  4. Applies score transformation if specified
  5. Returns both transformed and raw scores

- **Optional parameters**:
  - `TRANSFORMATION`: Score transformation parameters (e.g., sigmoid, reverse_sigmoid)
    - `transformation_type`: Type of transformation
    - `low`: Lower bound for transformation
    - `high`: Upper bound for transformation
    - `k`: Steepness parameter for sigmoid transformations

## 4. `docked_parallel_rocs_similarity`

This component doesn't exist as a separate class in the codebase. It likely refers to a workflow that combines docking and ROCS similarity:

1. First dock molecules using `dockstream`
2. Then use the docked poses as input for `parallel_rocs_similarity`

To implement this workflow, you would need to:
1. Run DockStream and save the docked poses
2. Use those poses as input for ParallelRocsSimilarity
3. Combine the scores from both components

This would require custom scripting to connect the two components, as there's no built-in class for this specific combination.
