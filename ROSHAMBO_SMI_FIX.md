# Roshambo SMI Format Fix

## 🔍 **Problem Identified**

From your debug output:
```
⚠️ Failed to embed molecule 0: CNC(C)C(=O)NC(C(=O)N1Cc2cc([O:1][CH2:1]COCCNc3ccc(...
⚠️ Failed to embed molecule 1: CNC(C)C(=O)NC(C(=O)N1Cc2cc([O:1][C:1](=O)C(C)NC(=O...
⚠️ Failed to embed molecule 2: CNC(C)C(=O)NC(C(=O)N1Cc2cc([O:1][c:1]3ccc(C(=O)NCC...
⚠️ Failed to embed molecule 3: CNC(C)C(=O)NC(C(=O)N1Cc2cc([O:1][C:1](=O)COCC(=O)N...
⚠️ Failed to embed molecule 4: CNC(C)C(=O)NC(C(=O)N1Cc2cc([O:1][CH2:1]c3ccc([NH:0...
```

**All molecules are failing to embed**, which means:
1. SDF file is empty or has very few molecules
2. Roshambo receives empty/invalid input
3. Returns all zero scores
4. Process eventually gets killed

## 🛠️ **Solution: SMI Format**

Added a new parameter `file_format` with two options:

### **Option 1: SMI Format (Recommended)**
```json
{
  "file_format": "smi",
  "n_confs": 50
}
```

**How it works:**
- Creates `.smi` file with SMILES strings
- Roshambo handles all conformer generation
- No memory-intensive RDKit 3D embedding
- Much more reliable for complex PROTAC molecules

### **Option 2: SDF Format (Original)**
```json
{
  "file_format": "sdf",
  "n_confs": 0
}
```

**How it works:**
- Pre-generates 3D coordinates with RDKit
- Creates `.sdf` file with embedded molecules
- More memory-intensive
- Can fail with complex molecules

## 📋 **Updated Configuration**

Use this configuration to fix the issue:

```json
{
  "component_type": "roshambo_shape_similarity",
  "name": "roshambo_scorer",
  "weight": 1.0,
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "file_format": "smi",           // ← KEY FIX: Use SMI format
    "n_confs": 50,                  // ← Roshambo will generate conformers
    "shape_weight": 0.5,
    "color_weight": 0.5,
    "ignore_hs": true,
    "use_carbon_radii": true,
    "gpu_id": 0,
    "roshambo_api_url": "http://127.0.0.1:5000",
    "save_overlays": true,
    "overlays_dir": "roshambo_overlays",
    "debug": true
  }
}
```

## 🔄 **What Changed**

### **Before (Failing):**
1. Convert SMILES → RDKit molecules
2. Add hydrogens
3. Generate 3D coordinates (FAILS for complex molecules)
4. Write SDF file (empty/few molecules)
5. Roshambo gets invalid input → all zeros

### **After (Working):**
1. Validate SMILES strings
2. Write SMI file directly
3. Roshambo reads SMILES
4. Roshambo generates conformers itself
5. Roshambo calculates similarity scores

## 🎯 **Expected Behavior**

With `file_format: "smi"`:

```
🔧 Roshambo initialized: reference=ref.sdf, n_confs=50, gpu_id=0
🔧 File format: SMI, API URL: http://127.0.0.1:5000
📝 Creating SMILES file with 113 molecules
🔧 Roshambo will generate 50 conformers per molecule
✅ Created SMILES file: dataset_2.smi
📊 Valid: 113, Invalid: 0, Total: 113
📋 Using filenames for API: reference='ref.sdf', dataset='dataset_2.smi'
```

## 🚀 **Benefits of SMI Format**

1. **No embedding failures**: Skips problematic RDKit 3D generation
2. **Lower memory usage**: No pre-generation of 3D coordinates
3. **Faster processing**: Direct SMILES validation only
4. **More reliable**: Roshambo handles conformer generation optimally
5. **No "Killed" errors**: Avoids memory exhaustion

## 🔍 **Debug Output to Watch**

Look for:
```
📝 Creating SMILES file with 113 molecules
📊 Valid: 113, Invalid: 0, Total: 113
```

If you see high invalid counts, there might be SMILES format issues.

## ✅ **Quick Fix**

**Just add this one parameter to your configuration:**
```json
"file_format": "smi"
```

This should resolve both the embedding failures and the "Killed" errors!
