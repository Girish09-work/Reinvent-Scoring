from dataclasses import dataclass


@dataclass(frozen=True)
class RoshamboSpecificParametersEnum:
    """Enum for Roshambo-specific parameters."""

    # Basic shape similarity parameters
    REFERENCE_FILE = "reference_file"
    SHAPE_WEIGHT = "shape_weight"
    COLOR_WEIGHT = "color_weight"

    # Multiple reference files for PROTAC design
    WARHEAD1_REFERENCE = "warhead1_reference"
    WARHEAD2_REFERENCE = "warhead2_reference"

    # Conformer generation parameters
    N_CONFS = "n_confs"
    IGNORE_HS = "ignore_hs"
    USE_CARBON_RADII = "use_carbon_radii"

    # GPU and parallelization parameters
    GPU_ID = "gpu_id"
    MAX_GPUS = "max_gpus"
    GPU_IDS = "gpu_ids"

    # Overlay saving parameters
    SAVE_OVERLAYS = "save_overlays"
    OVERLAYS_DIR = "overlays_dir"
    OVERLAY_PREFIX = "overlay_prefix"

    # Environment and execution parameters
    ENVIRONMENT_PATH = "environment_path"
    DEBUG = "debug"

    # Dynamic environment configuration parameters (REQUIRED)
    RDBASE_PATH = "rdbase_path"  # REQUIRED: Path to RDKit build directory
    CONDA_ENV_NAME = "conda_env_name"  # REQUIRED: Name of conda environment with Roshambo
    CONDA_BASE_PATH = "conda_base_path"  # OPTIONAL: Path to conda installation
    CUDA_HOME_PATH = "cuda_home_path"  # OPTIONAL: Path to CUDA installation
    AUTO_SETUP_ENV = "auto_setup_env"  # OPTIONAL: Auto-setup environment variables

    # Advanced Roshambo parameters
    VOLUME_TYPE = "volume_type"
    EPSILON = "epsilon"
    RES = "res"
    MARGIN = "margin"
    PROXY_CUTOFF = "proxy_cutoff"
    N = "n"

    # Feature definition file for color similarity
    FDEF_PATH = "fdef_path"

    # Sorting and output parameters
    SORT_BY = "sort_by"
    WRITE_TO_FILE = "write_to_file"
    MAX_CONFORMERS = "max_conformers"
    FILENAME = "filename"

    # RDKit conformer generation parameters
    RDKIT_METHOD = "rdkit_method"
    RDKIT_RANDOM_SEED = "rdkit_random_seed"
    RDKIT_FF = "rdkit_ff"
    RDKIT_ADD_HS = "rdkit_add_hs"
    RDKIT_OPT_CONFS = "rdkit_opt_confs"
    RDKIT_CALC_ENERGY = "rdkit_calc_energy"
    RDKIT_ENERGY_ITERS = "rdkit_energy_iters"
    RDKIT_ENERGY_CUTOFF = "rdkit_energy_cutoff"
    RDKIT_ALIGN_CONFS = "rdkit_align_confs"
    RDKIT_RMS_CUTOFF = "rdkit_rms_cutoff"
    RDKIT_NUM_THREADS = "rdkit_num_threads"
