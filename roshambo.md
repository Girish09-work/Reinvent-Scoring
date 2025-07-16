# ROSHAMBO Integration Guide (v2.0)

## Overview

This guide covers the new Roshambo integration with Reinvent Scoring, featuring a Flask API architecture for improved reliability and easier deployment.

## Architecture

The new Roshambo integration consists of:

1. **Simplified Scoring Component** (`roshambo_shape_similarity.py`)
2. **Flask API Server** (`roshambo/app.py`)
3. **Epoch-based File Management**
4. **Comprehensive Configuration Options**

## Installation Instructions

### Step 1: Create Roshambo Environment

```bash
conda create --name roshambo python=3.9.6
conda activate roshambo
```

### Step 2: Install Dependencies

```bash
# Install basic packages
conda install notebook flask requests pandas numpy

# Install RDKit
conda install -c conda-forge rdkit

# Install Roshambo
git clone https://github.com/rashatwi/roshambo.git
cd roshambo
pip install .
```

### Step 3: Set Up Environment Variables

```bash
export RDBASE=/path/to/rdkit
export RDKIT_LIB_DIR=$RDBASE/lib
export RDKIT_INCLUDE_DIR=$RDBASE/Code
export RDKIT_DATA_DIR=$RDBASE/Data
export PYTHONPATH=$PYTHONPATH:$RDBASE
export CUDA_HOME=/usr/local/cuda
```

### Step 4: Start Roshambo Flask API

```bash
cd roshambo
python start_api.py
```

The API will be available at `http://localhost:5000`

## New Roshambo Scoring Component Usage

### Basic Configuration

```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Shape Similarity",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "data/reference.sdf",
    "shape_weight": 0.5,
    "color_weight": 0.5,
    "roshambo_api_url": "http://localhost:5000",
    "save_overlays": true,
    "overlays_dir": "roshambo_overlays"
  }
}
```

### With RDKit Conformer Generation

```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Shape Similarity with Conformers",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "data/reference.sdf",
    "n_confs": 20,
    "rdkit_method": "ETKDGv3",
    "rdkit_add_hs": true,
    "rdkit_opt_confs": true,
    "roshambo_api_url": "http://localhost:5000"
  }
}
```

### Epoch-based File Organization

The component automatically creates epoch folders:

```
roshambo_overlays/
├── epoch_0/
│   ├── reference.sdf
│   ├── dataset_0.sdf
│   ├── mols.sdf
│   ├── hits.sdf
│   └── roshambo.csv
├── epoch_1/
│   └── ...
```

### Flask API Usage

```python
import requests

# Call the API directly
response = requests.post("http://localhost:5000/similarity", json={
    "reference_file": "/path/to/reference.sdf",
    "dataset_file": "/path/to/dataset.sdf",
    "ignore_hs": True,
    "n_confs": 0,
    "use_carbon_radii": True,
    "color": True,
    "write_to_file": True,
    "gpu_id": 0
})

result = response.json()
```

---

## Legacy ROSHAMBO API Usage

### Basic ROSHAMBO Run

```python
from roshambo.api import get_similarity_scores

get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    ignore_hs=True,
    n_confs=0,
    use_carbon_radii=True,
    color=True,
    sort_by="ComboTanimoto",
    write_to_file=True,
    gpu_id=0,
    working_dir="inpdata",
)
```

**Output files:**

* `mols.sdf`
* `hits.sdf`
* `roshambo.csv`

```python
import pandas as pd
df_default = pd.read_csv("inpdata/roshambo.csv", delimiter="\t")
df_default.head(10)
```

---

## Analytic Overlap (Higher Order)

```python
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    ignore_hs=True,
    n_confs=0,
    use_carbon_radii=True,
    color=True,
    sort_by="ComboTanimoto",
    write_to_file=True,
    gpu_id=0,
    volume_type="analytic",
    n=6,
    epsilon=0.5,
    working_dir="data/analytic",
)
```

---

## Grid Method

```python
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    ignore_hs=True,
    n_confs=0,
    use_carbon_radii=True,
    color=True,
    sort_by="ComboTanimoto",
    write_to_file=True,
    gpu_id=0,
    volume_type="gaussian",
    res=0.4,
    margin=0.4,
    working_dir="data/grid",
)
```

---

## Visualizing Overlap Methods

```python
merged_df = pd.merge(df_default, df_sixth, on='Molecule', suffixes=('_Default', '_Sixth'))
merged_df = pd.merge(merged_df, df_grid, on='Molecule')

fig, axs = plt.subplots(1, 2, figsize=(8, 4))
axs[0].scatter(merged_df["ShapeTanimoto_Default"], merged_df["ShapeTanimoto_Sixth"], color="#80B9F9")
axs[1].scatter(merged_df["ShapeTanimoto_Default"], merged_df["ShapeTanimoto"], color="#80B9F9")
```

---

## Custom Color Force Field

```python
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    ignore_hs=True,
    n_confs=0,
    use_carbon_radii=True,
    color=True,
    sort_by="ComboTanimoto",
    write_to_file=True,
    gpu_id=0,
    fdef_path="../data/features.json",
    working_dir="data/color",
)
```

---

## Visualizing Pharmacophores

```python
from rdkit import Chem
from roshambo.pharmacophore import draw_pharm, calc_custom_pharm, load_smarts_from_json, calc_rdkit_pharm

compiled_smarts = load_smarts_from_json("../data/features.json")
rdkit_mol = Chem.MolFromMolFile("data/basic_run/query.sdf")
custom_features = calc_custom_pharm(rdkit_mol, compiled_smarts)
rdkit_features = calc_rdkit_pharm(rdkit_mol)
draw_pharm(rdkit_mol, custom_features, working_dir="data/basic_run")
draw_pharm(rdkit_mol, rdkit_features, working_dir="data/basic_run")
```

---

## Conformer Generation

```python
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    ignore_hs=True,
    n_confs=10,
    keep_mol=True,
    random_seed=109838974,
    opt_confs=True,
    calc_energy=True,
    energy_iters=300,
    energy_cutoff=30,
    align_confs=True,
    rms_cutoff=0.1,
    num_threads=48,
    method="ETKDGv3",
    volume_type="analytic",
    n=2,
    epsilon=0.5,
    use_carbon_radii=True,
    color=True,
    max_conformers=3,
    sort_by="ComboTanimoto",
    write_to_file=True,
    gpu_id=0,
    fdef_path="../data/features.json",
    working_dir="data/conformers",
)
```

---

## Using SMILES Input

```python
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.smi",
    ignore_hs=True,
    n_confs=100,
    keep_mol=True,
    random_seed=109838974,
    opt_confs=False,
    calc_energy=False,
    energy_iters=300,
    energy_cutoff=30,
    align_confs=True,
    rms_cutoff=0.1,
    num_threads=48,
    method="ETKDGv3",
    volume_type="analytic",
    n=2,
    epsilon=0.5,
    use_carbon_radii=True,
    color=True,
    max_conformers=1,
    sort_by="ComboTanimoto",
    write_to_file=True,
    gpu_id=0,
    fdef_path="../data/features.json",
    working_dir="data/smiles",
    # smiles_kwargs={"delimiter": "\t"},
)
```

---

## ROC Analysis

```python
from roshambo.analysis import plot_scores_dist, calc_roc_auc, plot_mult_roc, plot_mult_auc, plot_mult_enrichment
```

### Score Distributions

```python
df = pd.read_csv("data/analysis/roshambo_ligands_CSF1R.csv", delimiter="\t")
plot_scores_dist(
    df,
    columns=["ComboTanimoto", "ShapeTanimoto", "ColorTanimoto"],
    title="CSF1R Score Distributions",
    working_dir="data/analysis",
)
```

### ROC Curves

```python
auc, roce, fig = calc_roc_auc(
    "data/analysis/roshambo_ligands_CXCR4.csv",
    "data/analysis/roshambo_decoys_CXCR4.csv",
    score="ComboTanimoto",
    n_bootstraps=1000,
    interpolation=True,
    eevs=[0.005, 0.01, 0.02, 0.05],
    plot=True,
    log=False,
    working_dir="data/analysis",
)
```

### Multiple ROC Curves

```python
fig = plot_mult_roc(
    rates_dict={
        "CXCR4": "data/analysis/roshambo_roc_CXCR4.csv",
        "CSF1R": "data/analysis/roshambo_roc_CSF1R.csv"
    },
    analysis_dict={
        "CXCR4": "data/analysis/roshambo_analysis_CXCR4.csv",
        "CSF1R": "data/analysis/roshambo_analysis_CSF1R.csv"
    },
    colors_dict={"CXCR4": "#80B9F9", "CSF1R": "#6DAD46"},
    title="CXCR4 vs. CSF1R",
    log=False,
    filename="CXCR4_CSF1R.jpg",
    working_dir="data/analysis",
)
```

### AUC Plot

```python
fig = plot_mult_auc(
    auc_dict={
        "ROSHAMBO": [
            "data/analysis/roshambo_analysis_CSF1R.csv",
            "data/analysis/roshambo_analysis_CXCR4.csv",
        ],
    },
    colors_dict={"ROSHAMBO": "#80B9F9"},
    group_labels=["CSF1R", "CXCR4"],
    working_dir="data/analysis",
)
```

### Enrichment Factors Plot

```python
fig = plot_mult_enrichment(
    enrich_dict={
        "ROSHAMBO": [
            "data/analysis/roshambo_analysis_CSF1R.csv",
            "data/analysis/roshambo_analysis_CXCR4.csv",
        ],
    },
    colors_dict={0: "#80B9F9", 1: "#6DAD46", 2: "gray", 3: "black"},
    hatch_patterns=[None, "+"],
    group_labels=["CSF1R", "CXCR4"],
    working_dir="data/analysis",
)
```
