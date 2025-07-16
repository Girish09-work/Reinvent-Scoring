# Roshambo Zero Scores Issue - Diagnosis and Fix

## Problem Description

The Roshambo shape similarity component is returning all zero scores during reinforcement learning, even though other scoring components (like Linker Num HBD and Linker Num Rings) are working correctly.

## Root Cause

The issue is caused by **missing required configuration parameters** for the Roshambo component. The component requires two mandatory parameters that are not present in the current configuration:

1. `rdbase_path` - Path to the RDKit build directory
2. `conda_env_name` - Name of the conda environment with Roshambo installed

Without these parameters, the Roshambo component fails to initialize properly and returns zero scores for all molecules.

## Solution

### Step 1: Update Your Configuration

Add the required parameters to your scoring configuration. Here's an example:

```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Roshambo Shape Similarity",
  "weight": 1.0,
  "specific_parameters": {
    "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
    "conda_env_name": "roshambo",
    "conda_base_path": "/home/protacinvent/.conda",
    "reference_file": "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BTK/BTK_sel.sdf",
    "warhead2_reference": "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BTK/BTK_sel.sdf",
    "shape_weight": 0.6,
    "color_weight": 0.4,
    "n_confs": 0,
    "ignore_hs": true,
    "use_carbon_radii": true,
    "gpu_id": 0,
    "debug": true
  }
}
```

### Step 2: Verify Your Paths

Make sure the following paths are correct for your system:

1. **rdbase_path**: Should point to your RDKit build directory
   - Example: `/home/protacinvent/Desktop/roshambo/rdkit`
   - This directory should contain subdirectories like `lib`, `Code`, `Data`

2. **conda_env_name**: Should be the name of your conda environment with Roshambo
   - Example: `roshambo`
   - Verify with: `conda env list`

3. **reference_file**: Should point to your reference SDF file
   - Example: `/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BTK/BTK_sel.sdf`

### Step 3: Test Your Configuration

Use the provided test script to verify your configuration:

```bash
python test_roshambo_configuration.py
```

This script will:
- Test component initialization
- Verify environment setup
- Run a simple scoring test
- Provide detailed error messages if something is wrong

### Step 4: Enable Debug Mode

Add `"debug": true` to your configuration to get detailed output about what's happening:

```json
"specific_parameters": {
  "debug": true,
  // ... other parameters
}
```

## What Was Fixed

1. **Parameter Access**: Fixed inconsistent parameter access in the component code
2. **Error Messages**: Improved error messages to clearly indicate missing required parameters
3. **Debug Output**: Added comprehensive debug output to help diagnose configuration issues
4. **Validation**: Enhanced parameter validation with helpful suggestions

## Files Modified

- `reinvent_scoring/scoring/score_components/roshambo/roshambo_shape_similarity.py`
  - Fixed parameter enum usage
  - Added better error messages
  - Enhanced debug output
  - Improved error handling

## Additional Resources

- `roshambo_config_template.json` - Complete configuration template
- `test_roshambo_configuration.py` - Test script for diagnosis
- `examples/roshambo_configuration_examples.json` - More configuration examples
- `ROSHAMBO_MANDATORY_PARAMETERS.md` - Detailed parameter documentation

## Expected Behavior After Fix

Once the configuration is corrected, you should see:

1. **Initialization**: Component initializes without errors
2. **Debug Output**: Detailed information about environment setup
3. **Non-Zero Scores**: Roshambo returns meaningful similarity scores
4. **CSV Generation**: Diversity filters work correctly with Roshambo scores

## Troubleshooting

If you still get zero scores after updating the configuration:

1. **Check Environment**: Verify conda environment exists and has Roshambo installed
2. **Check Paths**: Ensure all file paths are correct and accessible
3. **Check Permissions**: Make sure the process has read access to reference files
4. **Check Dependencies**: Verify RDKit and CUDA are properly installed
5. **Run Test Script**: Use `test_roshambo_configuration.py` for detailed diagnosis

## Contact

If you continue to experience issues, please provide:
- Your complete configuration file
- Output from the test script
- Debug output from the component
- Your system environment details (conda env list, paths, etc.)
