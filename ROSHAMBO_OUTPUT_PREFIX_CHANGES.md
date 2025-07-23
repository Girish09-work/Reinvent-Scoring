# Roshambo Output Prefix Parameter Implementation

## Overview
This document outlines the changes made to add the `--output_prefix` parameter to Roshambo, allowing users to specify custom prefixes for output files.

## Files Modified

### 1. `roshambo/api.py`
**Function**: `get_similarity_scores()`

**Changes Made**:
```python
def get_similarity_scores(
    ref_file,
    dataset_files_pattern,
    ignore_hs=True,
    n_confs=0,
    use_carbon_radii=True,
    color=True,
    sort_by="ComboTanimoto",
    write_to_file=True,
    working_dir=".",
    volume_type="analytic",
    n=2,
    epsilon=0.1,
    output_prefix="roshambo"  # NEW PARAMETER ADDED
):
```

**Implementation**:
```python
# Inside the function, replace hardcoded "roshambo" with output_prefix
if write_to_file:
    output_file = os.path.join(working_dir, f"{output_prefix}.csv")
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Results written to {output_file}")
```

### 2. `roshambo/cli.py` (if exists)
**Changes Made**:
```python
import argparse

def main():
    parser = argparse.ArgumentParser(description='Roshambo molecular similarity scoring')
    
    # Existing arguments...
    parser.add_argument('--ref_file', required=True, help='Reference molecule file')
    parser.add_argument('--dataset_files_pattern', required=True, help='Dataset files pattern')
    
    # NEW ARGUMENT ADDED
    parser.add_argument('--output_prefix', default='roshambo', 
                       help='Prefix for output files (default: roshambo)')
    
    args = parser.parse_args()
    
    # Call API with new parameter
    get_similarity_scores(
        ref_file=args.ref_file,
        dataset_files_pattern=args.dataset_files_pattern,
        output_prefix=args.output_prefix,  # NEW PARAMETER PASSED
        # ... other parameters
    )
```

### 3. `roshambo/core/similarity.py`
**Function**: `calculate_similarities()` or similar core function

**Changes Made**:
```python
def calculate_similarities(ref_mol, dataset_mols, output_prefix="roshambo", **kwargs):
    """
    Calculate molecular similarities with custom output prefix
    
    Parameters:
    -----------
    ref_mol : rdkit.Chem.Mol
        Reference molecule
    dataset_mols : list
        List of dataset molecules
    output_prefix : str
        Prefix for output files (default: "roshambo")
    """
    
    # Process similarities...
    results = []
    
    # When saving intermediate files or logs
    log_file = f"{output_prefix}_processing.log"
    temp_file = f"{output_prefix}_temp.sdf"
    
    return results
```

### 4. `roshambo/utils/file_handler.py` (if exists)
**Changes Made**:
```python
def save_results(results_df, working_dir=".", output_prefix="roshambo", file_format="csv"):
    """
    Save results with custom prefix
    
    Parameters:
    -----------
    results_df : pandas.DataFrame
        Results dataframe
    working_dir : str
        Working directory
    output_prefix : str
        Prefix for output files
    file_format : str
        Output file format (csv, tsv, etc.)
    """
    
    if file_format.lower() == "csv":
        output_file = os.path.join(working_dir, f"{output_prefix}.csv")
        results_df.to_csv(output_file, sep='\t', index=False)
    elif file_format.lower() == "json":
        output_file = os.path.join(working_dir, f"{output_prefix}.json")
        results_df.to_json(output_file, orient='records', indent=2)
    
    return output_file
```

## Usage Examples

### 1. Python API Usage
```python
from roshambo.api import get_similarity_scores

# Default usage (creates roshambo.csv)
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    working_dir="results/"
)

# Custom prefix usage (creates custom_analysis.csv)
get_similarity_scores(
    ref_file="query.sdf",
    dataset_files_pattern="dataset.sdf",
    working_dir="results/",
    output_prefix="custom_analysis"
)

# Multiple analyses with different prefixes
prefixes = ["benzene_analysis", "phenol_analysis", "aniline_analysis"]
for i, prefix in enumerate(prefixes):
    get_similarity_scores(
        ref_file=f"query_{i}.sdf",
        dataset_files_pattern="dataset.sdf",
        working_dir="results/",
        output_prefix=prefix
    )
```

### 2. Command Line Usage
```bash
# Default usage
python -m roshambo --ref_file query.sdf --dataset_files_pattern dataset.sdf

# Custom prefix usage
python -m roshambo --ref_file query.sdf --dataset_files_pattern dataset.sdf --output_prefix my_analysis

# Batch processing with different prefixes
python -m roshambo --ref_file benzene.sdf --dataset_files_pattern dataset.sdf --output_prefix benzene_similarity
python -m roshambo --ref_file phenol.sdf --dataset_files_pattern dataset.sdf --output_prefix phenol_similarity
```

## Configuration File Support

### `roshambo_config.json`
```json
{
    "reference_file": "query.sdf",
    "dataset_pattern": "dataset.sdf",
    "output_settings": {
        "prefix": "custom_analysis",
        "working_dir": "./results",
        "format": "csv"
    },
    "similarity_settings": {
        "ignore_hs": true,
        "use_carbon_radii": true,
        "color": true,
        "sort_by": "ComboTanimoto"
    }
}
```

## Backward Compatibility
- Default value for `output_prefix` is "roshambo" to maintain backward compatibility
- Existing scripts will continue to work without modification
- New parameter is optional in all function calls

## Testing

### Test Cases to Add
```python
def test_output_prefix_default():
    """Test that default prefix creates roshambo.csv"""
    result_file = get_similarity_scores(...)
    assert os.path.exists("roshambo.csv")

def test_output_prefix_custom():
    """Test that custom prefix creates correctly named file"""
    result_file = get_similarity_scores(..., output_prefix="test_analysis")
    assert os.path.exists("test_analysis.csv")

def test_output_prefix_with_path():
    """Test prefix with directory path"""
    result_file = get_similarity_scores(..., output_prefix="results/analysis")
    assert os.path.exists("results/analysis.csv")
```

## Implementation Checklist

- [ ] Add `output_prefix` parameter to `get_similarity_scores()` function
- [ ] Update all file output operations to use the prefix
- [ ] Add command line argument support
- [ ] Update configuration file parsing
- [ ] Add unit tests for new functionality
- [ ] Update documentation and examples
- [ ] Ensure backward compatibility
- [ ] Test with various file formats (CSV, JSON, etc.)

## Benefits
1. **Organized Output**: Users can organize results by experiment or molecule type
2. **Batch Processing**: Easy to run multiple analyses without file conflicts
3. **Integration**: Better integration with automated pipelines
4. **Flexibility**: Supports both simple names and directory paths as prefixes
