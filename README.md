# Reinvent Scoring with Enhanced Shape Similarity Components

This package is a fork of the original [reinvent-scoring](https://github.com/Girish09-work/reinvent-scoring) package with additional shape similarity components including RDKit-based and enhanced Roshambo implementations.

## Added Features

### RDKit Shape Components
- `rdkit_shape_similarity`: An open-source alternative to ROCS similarity using RDKit
- `parallel_rdkit_shape_similarity`: A parallel version for faster processing

### Enhanced Roshambo Shape Components
- `roshambo_shape_similarity`: GPU-accelerated shape similarity using Roshambo with RDKit conformer generation
- `parallel_roshambo_shape_similarity`: Multi-GPU parallel version for high-throughput processing
- `RoshamboConformerGenerator`: Dedicated RDKit-based conformer generator for Roshambo
- Environment path support for conda environments
- Enhanced parameter management with dedicated enums
- Multiple reference file support for PROTAC design

## Why This Fork?

The original package relies on OpenEye's ROCS for 3D shape similarity, which requires a commercial license. This fork adds both RDKit-based alternatives and enhanced Roshambo implementations that provide similar or superior functionality:

- **RDKit components**: Completely open-source alternative to ROCS
- **Enhanced Roshambo**: GPU-accelerated with RDKit conformer generation for optimal performance
- **Environment flexibility**: Support for conda environments and different installation setups
- **PROTAC support**: Multiple reference files for complex molecular design scenarios

## Installation

```bash
pip install reinvent-scoring-rdkit
```

## Usage

### RDKit Shape Similarity

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

### Enhanced Roshambo Shape Similarity

```json
{
    "component_type": "roshambo_shape_similarity",
    "name": "Roshambo Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "shape_weight": 0.6,
        "color_weight": 0.4,
        "n_confs": 20,
        "environment_path": "conda activate roshambo_env &&",
        "gpu_id": 0,
        "save_overlays": true,
        "overlays_dir": "shape_overlays"
    }
}
```

### Parallel Roshambo for Multi-GPU

```json
{
    "component_type": "parallel_roshambo_shape_similarity",
    "name": "Parallel Roshambo",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "path/to/reference.sdf",
        "shape_weight": 0.6,
        "color_weight": 0.4,
        "n_confs": 20,
        "max_gpus": 2,
        "gpu_ids": [0, 1],
        "debug": true
    }
}
```

### PROTAC Design with Multiple References

```json
{
    "component_type": "roshambo_shape_similarity",
    "name": "PROTAC Shape Similarity",
    "weight": 1.0,
    "specific_parameters": {
        "reference_file": "protac_reference.sdf",
        "warhead1_reference": "warhead1.sdf",
        "warhead2_reference": "warhead2.sdf",
        "shape_weight": 0.6,
        "color_weight": 0.4,
        "n_confs": 30,
        "save_overlays": true
    }
}
```

## Original Package

This package is based on the original reinvent-scoring package. All credit for the original code goes to the original authors.

# Developing
## Setup environment
You can use Conda to create an environment with all the necessary packages installed.

```
$ conda env create -f reinvent_scoring
[...]
$ conda activate reinvent_scoring
```

## Run tests
The tests use the `unittest` package testing framework.  Before you can run the tests make sure that you have created a
`config.json`file in the `reinvent_scoring/configs` directory.  There is an example config in the same directory, which
you can base your own config off of.  Make sure that you set `MAIN_TEST_PATH` to a non-existent directory; it is where
temporary files will be written during the tests; if it is set to an existing directory, that directory will be removed
once the tests have finished.

Some tests require a proprietary OpenEye license; you have to set up a few things to make the tests read your
license.  The simple way is to just set the `OE_LICENSE` environment variable to the path of the file containing the
license.  If you just want to set the license in the `reinvent_scoring` Conda environment, it is a bit more complicated,
but you only have to do it once.

```
(reinvent-scoring) $ cd $CONDA_PREFIX
$ mkdir -p etc/conda/activate.d
$ mkdir -p etc/conda/deactivate.d
```

Put the following in `etc/conda/activate.d/env_vars.sh`.

```
#!/bin/sh
export OE_LICENSE='</path/to/your/oe_license/file>'
```

And put the following in `etc/conda/deactivate.d/env_vars.sh`.

```
#!/bin/sh
unset OE_LICENSE
```

Once you have created the files, deactivate and re-activate the environment, and `echo $OE_LICENSE` should output the
path to the license file.

Once you have created and configured your environment, you can run unittests by running

```bash
python main_test.py --unittests
```

If you have a valid Open eye license and other dependencie configured, like Icolos and AZDOCK -
you can also run integration tests, by running command (remember to submit this configuration, since the default one is test):

```bash
python main_test.py --integration --base_config <path to your configuration>
```

# Building
- Building: `python setup.py sdist bdist_wheel`
- Upload build to test: `python -m twine upload --repository testpypi dist/*`
- Upload build: `python -m twine upload dist/*`

