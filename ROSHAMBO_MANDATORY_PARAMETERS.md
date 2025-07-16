# Roshambo Mandatory Parameters Configuration

## ⚠️ Important: Required Parameters

You are absolutely correct! The Roshambo shape similarity component **requires** users to specify the RDKit build path and conda environment. Auto-discovery configurations are not practical for production use.

## 🔧 **Mandatory Parameters**

### **1. `rdbase_path` (REQUIRED)**
- **Purpose**: Path to your RDKit build directory
- **Example**: `"/home/protacinvent/Desktop/roshambo/rdkit"`
- **Why Required**: Roshambo needs to know where RDKit is built to access libraries and Python modules

### **2. `conda_env_name` (REQUIRED)**
- **Purpose**: Name of the conda environment where Roshambo is installed
- **Example**: `"roshambo"`
- **Why Required**: Roshambo must be run from the correct conda environment with all dependencies

## ✅ **Correct Configuration Examples**

### **Minimal Required Configuration**
```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Roshambo Shape Similarity",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
    "conda_env_name": "roshambo",
    "shape_weight": 0.6,
    "color_weight": 0.4,
    "n_confs": 10,
    "gpu_id": 0
  }
}
```

### **Complete Configuration with Optional Parameters**
```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "Roshambo Complete Configuration",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
    "conda_env_name": "roshambo",
    "conda_base_path": "/opt/conda",
    "cuda_home_path": "/usr/local/cuda",
    "auto_setup_env": true,
    "shape_weight": 0.6,
    "color_weight": 0.4,
    "n_confs": 20,
    "gpu_id": 0,
    "save_overlays": true,
    "debug": true
  }
}
```

### **PROTAC Design Configuration**
```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "PROTAC Roshambo",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "D:\\protac-invent\\Protac-invent\\data\\protac\\BRD9\\protac.sdf",
    "warhead1_reference": "D:\\protac-invent\\Protac-invent\\data\\protac\\BRD9\\warhead1.sdf",
    "warhead2_reference": "D:\\protac-invent\\Protac-invent\\data\\protac\\BTK\\BTK_sel.sdf",
    "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
    "conda_env_name": "roshambo",
    "conda_base_path": "/opt/conda",
    "auto_setup_env": true,
    "shape_weight": 0.6,
    "color_weight": 0.4,
    "n_confs": 20,
    "save_overlays": true,
    "overlays_dir": "D:\\protac-invent\\Protac-invent\\data\\protac\\BRD9\\roshambo_overlays",
    "debug": true
  }
}
```

## ❌ **Invalid Configurations (Removed)**

These configurations are **NOT SUPPORTED** because they lack mandatory parameters:

```json
// ❌ INVALID - Missing rdbase_path and conda_env_name
{
  "component_type": "roshambo_shape_similarity",
  "name": "Roshambo Auto-Discovery",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "auto_setup_env": true,
    "debug": true
  }
}
```

## 🔍 **Parameter Validation**

The component now validates required parameters at initialization:

```python
# Validation in RoshamboShapeSimilarity.__init__()
if not self.rdbase_path:
    raise ValueError("'rdbase_path' is required - specify the path to your RDKit build directory")
if not self.conda_env_name:
    raise ValueError("'conda_env_name' is required - specify the name of your conda environment with Roshambo installed")
```

## 📋 **Setup Requirements**

Before using Roshambo components, ensure:

### **1. RDKit Build Directory**
- RDKit must be built and available at the specified `rdbase_path`
- Directory must contain: `lib/`, `Code/`, `Data/` subdirectories
- Python modules must be accessible

### **2. Conda Environment**
- Create dedicated conda environment: `conda create -n roshambo python=3.8`
- Install required packages:
  ```bash
  conda activate roshambo
  pip install rdkit
  pip install git+https://github.com/molecularinformatics/roshambo.git
  ```

### **3. Environment Variables**
The component automatically sets these when `auto_setup_env: true`:
```bash
export RDBASE="/home/protacinvent/Desktop/roshambo/rdkit"
export RDKIT_LIB_DIR="$RDBASE/lib"
export RDKIT_INCLUDE_DIR="$RDBASE/Code"
export RDKIT_DATA_DIR="$RDBASE/Data"
export PYTHONPATH="$RDBASE:$CONDA_PREFIX/lib/python*/site-packages"
export LD_LIBRARY_PATH="$RDBASE/lib:$CONDA_PREFIX/lib:/usr/local/cuda/lib64"
export CUDA_HOME="/usr/local/cuda"
```

## 🎯 **Benefits of Mandatory Parameters**

1. **Reliability**: No guessing about environment setup
2. **Clarity**: Users know exactly what to configure
3. **Debugging**: Clear error messages when paths are wrong
4. **Consistency**: Same configuration works across different systems
5. **Production Ready**: Suitable for deployment environments

## 🚨 **Error Messages**

If mandatory parameters are missing, you'll see clear error messages:

```
ValueError: 'rdbase_path' is required - specify the path to your RDKit build directory
ValueError: 'conda_env_name' is required - specify the name of your conda environment with Roshambo installed
```

## 📝 **Migration Guide**

If you were using auto-discovery configurations, update them:

### **Before (Invalid)**
```json
{
  "conda_env_name": "roshambo",
  "auto_setup_env": true
}
```

### **After (Valid)**
```json
{
  "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
  "conda_env_name": "roshambo",
  "auto_setup_env": true
}
```

## ✅ **Summary**

- **`rdbase_path`**: REQUIRED - Path to RDKit build directory
- **`conda_env_name`**: REQUIRED - Conda environment name
- **`conda_base_path`**: OPTIONAL - Auto-discovered if not specified
- **`cuda_home_path`**: OPTIONAL - Defaults to "/usr/local/cuda"
- **`auto_setup_env`**: OPTIONAL - Defaults to true

This approach ensures reliable, production-ready Roshambo configurations that work consistently across different user environments.
