# Roshambo Memory Management Fix

## 🔍 **Problem Identified**

The "Killed" message indicates the process was terminated by the system, most likely due to:

1. **Memory exhaustion**: RDKit 3D embedding of large PROTAC molecules consumes excessive RAM
2. **Complex molecules**: Your molecules appear to be large PROTACs with long linkers
3. **Accumulating memory**: Multiple epochs without proper cleanup

## 🛠️ **Fixes Applied**

### 1. **Memory Monitoring & Cleanup**
- Added memory usage monitoring during SDF creation
- Implemented garbage collection every 50 molecules
- Added proper cleanup of RDKit molecule objects
- Limited error message output to prevent memory bloat

### 2. **Robust Embedding**
- Added error handling for embedding failures
- Limited embedding attempts to prevent infinite loops
- Added timeout protection for problematic molecules

### 3. **Fallback Mode**
- Added `skip_3d_embedding` option to create SMILES files instead of SDF
- Roshambo can handle SMILES files and generate conformers itself

## 📋 **Configuration Options**

Add these parameters to your scoring configuration:

```json
{
  "component_type": "roshambo_shape_similarity",
  "specific_parameters": {
    "reference_file": "path/to/reference.sdf",
    "debug": true,
    
    // Memory management options
    "skip_3d_embedding": false,        // Set to true if memory issues persist
    "max_embed_attempts": 1,           // Limit embedding attempts
    "n_confs": 25,                     // Reduce from 50 if needed
    
    // Standard options
    "shape_weight": 0.5,
    "color_weight": 0.5,
    "gpu_id": 0
  }
}
```

## 🚀 **Quick Solutions**

### **Option 1: Try with current fixes**
The memory management improvements should resolve the issue. Run again and monitor the debug output for memory usage.

### **Option 2: Enable fallback mode**
If still getting killed, add this to your config:
```json
"skip_3d_embedding": true
```

This will create SMILES files instead of SDF, letting Roshambo handle conformer generation.

### **Option 3: Reduce batch size**
If using reinforcement learning, reduce the number of molecules per step in your RL configuration.

### **Option 4: Reduce conformers**
Lower the `n_confs` parameter:
```json
"n_confs": 10  // Instead of 50
```

## 🔍 **Debug Output to Watch**

The enhanced debug output will now show:
```
🧪 Creating SDF file with 118 molecules
💾 Current memory usage: 245.3 MB
📊 Success: 85, Failed: 33, Total: 118
💾 Final memory usage: 267.1 MB
```

Watch for:
- **Memory spikes**: If memory jumps dramatically
- **High failure rate**: If many molecules fail to embed
- **Memory not decreasing**: If garbage collection isn't working

## 🆘 **If Still Getting Killed**

1. **Check system memory**: `free -h` or `htop`
2. **Enable fallback mode**: Set `skip_3d_embedding: true`
3. **Reduce batch size**: Process fewer molecules at once
4. **Use smaller molecules**: Test with simpler SMILES first

## 📊 **Expected Behavior**

With these fixes:
- ✅ Memory usage should be monitored and controlled
- ✅ Failed embeddings won't crash the process
- ✅ Garbage collection should prevent memory leaks
- ✅ Fallback mode available if needed

The component should now handle large PROTAC molecules more robustly without getting killed by the system.
