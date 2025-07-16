# Clean Roshambo Shape Similarity Component

## Overview

This is a clean, optimized implementation of the Roshambo shape similarity scoring component for Reinvent Scoring. The component uses the Roshambo Flask API for GPU-accelerated molecular shape comparison with proper epoch-based organization.

## Key Features

### 🧹 Clean Implementation
- Removed unnecessary commented code blocks
- Streamlined API calls using direct Flask API approach
- Clear debug output with emojis for easy identification
- Proper error handling without excessive try-catch blocks

### 📁 Epoch-Based Organization
- Creates `epoch_{step}` folders for each scoring iteration
- Copies reference files to each epoch folder
- Generates dataset SDF files with RDKit molecules
- Organizes all roshambo outputs (CSV, SDF files) in epoch folders

### ⚙️ Configuration
- **Default n_confs**: 50 (as requested)
- **GPU usage**: Explicitly uses GPU 0 by default
- **API endpoint**: Uses localhost (127.0.0.1) for security
- **Debug mode**: Comprehensive debug output with bottleneck identification

## API Data Structure

The component sends the following data to the Roshambo Flask API:

```python
api_data = {
    "reference_file": reference_file,      # Path to reference SDF file
    "dataset_file": dataset_file,          # Path to dataset SDF file  
    "ignore_hs": self.ignore_hs,           # Ignore hydrogens
    "n_confs": self.n_confs,               # Number of conformers (default 50)
    "use_carbon_radii": self.use_carbon_radii,  # Use carbon radii
    "color": self.color_weight > 0,        # Enable color similarity
    "sort_by": "ComboTanimoto",            # Sort by combined score
    "write_to_file": True,                 # Generate CSV output
    "gpu_id": self.gpu_id,                 # GPU ID (default 0)
    "working_dir": epoch_folder            # Working directory
}
```

## File Organization

```
roshambo_overlays/
├── epoch_0/
│   ├── reference.sdf          # Copy of reference file
│   ├── dataset_0.sdf          # Generated dataset SDF
│   ├── roshambo.csv           # Roshambo output CSV
│   ├── mols.sdf               # Roshambo molecules output
│   └── hits.sdf               # Roshambo hits output
├── epoch_1/
│   └── ...
```

## CSV Processing

The component reads the Roshambo CSV output with tab delimiter:

```
Molecule	OriginalName	ComboTanimoto	ShapeTanimoto	ColorTanimoto	...
mol_58_0	mol_58	0.363	0.307	0.056	...
mol_24_0	mol_24	0.354	0.256	0.098	...
```

### Score Calculation
1. **Primary**: Uses `ComboTanimoto` directly from CSV
2. **Fallback**: Calculates weighted combination: `(shape_weight * ShapeTanimoto + color_weight * ColorTanimoto) / (shape_weight + color_weight)`
3. **Multiple conformers**: Takes the maximum score for each molecule

## Debug Output

When `debug=True`, the component provides detailed logging:

```
🔧 Roshambo initialized: reference=ref.sdf, n_confs=50, gpu_id=0
🧮 Roshambo calculate_score called with 113 molecules, step=0
📝 Converted 113 molecules to SMILES
⚙️ Processing 113 SMILES for step 0
📁 Created epoch folder: roshambo_overlays/epoch_0
🧪 Creating SDF file with 113 molecules
📄 Created dataset file: dataset_0.sdf
📄 Copied reference file: reference.sdf
🌐 Calling Roshambo API at http://127.0.0.1:5000/similarity
📡 API response status: 200
✅ API call successful, execution time: 45.67s
📊 Extracting scores from CSV: roshambo.csv
📈 CSV loaded: 113 rows, 11 columns
✅ Extracted 113 scores, mean: 0.245
📊 Calculated 113 scores, mean: 0.245
```

## Usage Example

```python
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity

# Configure parameters
parameters = ComponentParameters(
    component_type="roshambo_shape_similarity",
    name="roshambo_scorer",
    weight=1.0,
    specific_parameters={
        "reference_file": "path/to/reference.sdf",
        "shape_weight": 0.5,
        "color_weight": 0.5,
        "n_confs": 50,
        "ignore_hs": True,
        "use_carbon_radii": True,
        "gpu_id": 0,
        "roshambo_api_url": "http://127.0.0.1:5000",
        "save_overlays": True,
        "overlays_dir": "roshambo_overlays",
        "debug": True
    }
)

# Initialize component
component = RoshamboShapeSimilarity(parameters)

# Calculate scores
molecules = ["c1ccccc1", "CCCc1ccccc1", "CCCCCCCC"]
result = component.calculate_score(molecules, step=0)
print(f"Scores: {result.total_score}")
```

## Requirements

1. **Roshambo Flask API** running on `http://127.0.0.1:5000`
2. **RDKit** for molecule handling and SDF generation
3. **pandas** for CSV processing
4. **numpy** for array operations
5. **requests** for API communication

## Testing

Use the provided test script:

```bash
python test_roshambo_component.py
```

This will test the component with sample molecules and verify the complete workflow.
