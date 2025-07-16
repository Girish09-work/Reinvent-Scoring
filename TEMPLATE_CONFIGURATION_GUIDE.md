# Template.json Configuration Guide for Roshambo

## 📋 **Updated Template.json Configuration**

The template.json has been updated to include the enhanced Roshambo shape similarity component with **mandatory parameters**. Here's what you need to know:

## 🔧 **Mandatory Parameters (MUST BE CONFIGURED)**

### **1. `rdbase_path` (REQUIRED)**
```json
"rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit"
```
- **Purpose**: Path to your RDKit build directory
- **What to change**: Update this path to match your actual RDKit build location
- **Example paths**:
  - `/home/protacinvent/Desktop/roshambo/rdkit`
  - `/home/user/rdkit-build`
  - `/opt/rdkit`

### **2. `conda_env_name` (REQUIRED)**
```json
"conda_env_name": "roshambo"
```
- **Purpose**: Name of your conda environment with Roshambo installed
- **What to change**: Update this to match your actual conda environment name
- **Example names**:
  - `"roshambo"`
  - `"my_roshambo_env"`
  - `"protac_env"`

## ⚙️ **Optional Parameters (Can Be Customized)**

### **Environment Paths**
```json
"conda_base_path": "/home/protacinvent/.conda",
"cuda_home_path": "/usr/local/cuda"
```
- **conda_base_path**: Path to your conda installation (auto-discovered if not specified)
- **cuda_home_path**: Path to CUDA installation

### **Reference Files**
```json
"reference_file": "D:\\protac-invent\\Protac-invent\\data\\protac\\BRD9\\protac.sdf",
"warhead1_reference": "D:\\protac-invent\\Protac-invent\\data\\protac\\BRD9\\warhead1.sdf",
"warhead2_reference": "D:\\protac-invent\\Protac-invent\\data\\protac\\BTK\\BTK_sel.sdf"
```
- **reference_file**: Main reference structure for shape similarity
- **warhead1_reference**: First warhead reference (for PROTAC design)
- **warhead2_reference**: Second warhead reference (for PROTAC design)

### **Output Settings**
```json
"save_overlays": true,
"overlays_dir": "D:\\protac-invent\\Protac-invent\\data\\protac\\BRD9\\roshambo_overlays",
"overlay_prefix": "protac_overlay"
```
- **save_overlays**: Whether to save molecular overlays for visualization
- **overlays_dir**: Directory to save overlay files
- **overlay_prefix**: Prefix for overlay filenames

## 📝 **How to Customize Your Template.json**

### **Step 1: Update Mandatory Parameters**
1. **Find your RDKit build path**:
   ```bash
   # Look for your RDKit build directory
   find /home -name "rdkit" -type d 2>/dev/null
   ```

2. **Check your conda environment**:
   ```bash
   conda env list
   ```

3. **Update template.json**:
   ```json
   "rdbase_path": "/your/actual/path/to/rdkit",
   "conda_env_name": "your_actual_env_name"
   ```

### **Step 2: Update File Paths**
1. **Update reference files** to point to your actual SDF files:
   ```json
   "reference_file": "/path/to/your/reference.sdf",
   "warhead1_reference": "/path/to/your/warhead1.sdf",
   "warhead2_reference": "/path/to/your/warhead2.sdf"
   ```

2. **Update output directory**:
   ```json
   "overlays_dir": "/path/to/your/output/directory"
   ```

### **Step 3: Adjust Performance Settings**
```json
"n_confs": 20,          // Number of conformers to generate
"gpu_id": 0,            // GPU to use (0, 1, 2, etc.)
"rdkit_num_threads": 2, // Number of threads for RDKit
"debug": true           // Enable debug output
```

## 🚨 **Common Configuration Errors**

### **Error 1: Missing Mandatory Parameters**
```
ValueError: 'rdbase_path' is required - specify the path to your RDKit build directory
```
**Solution**: Add the missing parameter to your configuration.

### **Error 2: Invalid RDBASE Path**
```
Environment setup failed: Specified RDBASE path not valid: /wrong/path
```
**Solution**: Verify the path exists and contains `lib/`, `Code/`, and `Data/` directories.

### **Error 3: Conda Environment Not Found**
```
Environment setup failed: Conda environment 'wrong_env' not found
```
**Solution**: Check `conda env list` and use the correct environment name.

## ✅ **Validation Checklist**

Before running, ensure:

- [ ] `rdbase_path` points to a valid RDKit build directory
- [ ] `conda_env_name` matches an existing conda environment
- [ ] Reference SDF files exist at specified paths
- [ ] Output directory is writable
- [ ] GPU ID is valid (if using GPU acceleration)

## 🔧 **Environment Setup Commands**

If you need to set up the environment:

```bash
# Create conda environment
conda create -n roshambo python=3.8

# Activate environment
conda activate roshambo

# Install required packages
pip install rdkit
pip install git+https://github.com/molecularinformatics/roshambo.git

# Verify installation
python -c "from rdkit import Chem; import roshambo; print('Setup complete!')"
```

## 📊 **Template.json Structure**

The current template.json includes:

```json
{
  "version": 3,
  "run_type": "reinforcement_learning",
  "parameters": {
    "scoring_function": {
      "name": "custom_product",
      "parallel": true,
      "parameters": [
        {
          "component_type": "linker_num_hbd",
          // ... other components
        },
        {
          "component_type": "roshambo_shape_similarity",
          "name": "Roshambo Shape Similarity",
          "weight": 2.0,
          "specific_parameters": {
            // Mandatory parameters
            "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
            "conda_env_name": "roshambo",
            
            // Reference files
            "reference_file": "path/to/protac.sdf",
            "warhead1_reference": "path/to/warhead1.sdf",
            "warhead2_reference": "path/to/warhead2.sdf",
            
            // Performance settings
            "shape_weight": 0.6,
            "color_weight": 0.4,
            "n_confs": 20,
            "gpu_id": 0,
            
            // Output settings
            "save_overlays": true,
            "overlays_dir": "path/to/output",
            
            // Debug
            "debug": true
          }
        }
      ]
    }
  }
}
```

## 🎯 **Key Benefits of Updated Template**

1. **🔒 Mandatory Parameters**: Ensures reliable configuration
2. **📝 Clear Comments**: Explains each parameter section
3. **🎛️ Organized Structure**: Logical grouping of parameters
4. **🐛 Debug Ready**: Includes debug settings for troubleshooting
5. **🚀 Production Ready**: Suitable for actual PROTAC design workflows

## 📞 **Support**

If you encounter issues:

1. **Check the error message** - it usually indicates what's wrong
2. **Verify file paths** - ensure all files exist
3. **Test environment** - run the test script to validate setup
4. **Enable debug mode** - set `"debug": true` for detailed output

The updated template.json is now configured with mandatory parameters and ready for production use with your Roshambo environment! 🎉
