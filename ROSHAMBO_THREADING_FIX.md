# Roshambo Threading/Signal Issue Fix

## 🛠️ Problem Identified

You were absolutely right! The issue was with the `signal` module being used in Flask threads. Here's what was happening:

### The Error
```
ValueError: signal only works in main thread
```

### Root Cause
1. **Flask Development Server**: Runs requests in separate threads (not the main thread)
2. **Signal Module Limitation**: Python's `signal` module can only be used in the main thread
3. **Our Original Code**: Used `signal.SIGALRM` for timeout handling in the Flask request thread
4. **Result**: Python throws the error and aborts execution with HTTP 500

## ✅ Solution Implemented

### Option Chosen: Subprocess Approach
Instead of using threads or signal-based timeouts, we implemented a **subprocess-based solution**:

1. **Separate Subprocess Script**: `run_roshambo_subprocess.py`
2. **Process Isolation**: Roshambo runs in a completely separate process
3. **Built-in Timeout**: Uses `subprocess.run()` with timeout parameter
4. **No Signal Module**: Completely avoids signal-related issues

### Files Modified/Created

#### 1. `roshambo_api/app.py` - Fixed Threading Issues
**Before:**
```python
# Used signal module (BROKEN in Flask threads)
signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(timeout_seconds)
get_similarity_scores(**kwargs)  # Direct call
```

**After:**
```python
# Uses subprocess with built-in timeout (WORKS in Flask threads)
result = subprocess.run(
    [sys.executable, subprocess_script, params_json],
    capture_output=True,
    text=True,
    timeout=timeout_seconds,  # Built-in timeout, no signal needed
    cwd=working_dir
)
```

#### 2. `roshambo_api/run_roshambo_subprocess.py` - New Subprocess Script
- Runs roshambo in isolated process
- Handles JSON parameter passing
- Provides error handling and result reporting
- No threading/signal issues

### Key Improvements

#### ✅ **Threading Safety**
- **Before**: Used `signal` module in Flask threads (BROKEN)
- **After**: Uses subprocess with built-in timeout (WORKS)

#### ✅ **Process Isolation**
- **Before**: Roshambo ran in same process as Flask
- **After**: Roshambo runs in separate process
- **Benefits**: Better error isolation, memory management, timeout handling

#### ✅ **Timeout Handling**
- **Before**: `signal.SIGALRM` (thread-unsafe)
- **After**: `subprocess.run(timeout=...)` (thread-safe)

#### ✅ **Error Handling**
- **Before**: Signal-based timeout errors
- **After**: Subprocess timeout exceptions with proper error messages

## 🚀 How It Works Now

### 1. API Request Flow
```
Flask Request (Thread) 
    ↓
call_roshambo_api() 
    ↓
subprocess.run() with timeout
    ↓
run_roshambo_subprocess.py (Separate Process)
    ↓
roshambo.api.get_similarity_scores()
    ↓
Creates: roshambo.csv, mols.sdf, hits.sdf
    ↓
Returns JSON result to Flask
```

### 2. Timeout Handling
```python
# Thread-safe subprocess timeout
result = subprocess.run(
    [sys.executable, subprocess_script, params_json],
    timeout=timeout_seconds  # No signal module needed!
)
```

### 3. Error Isolation
- **Subprocess crashes**: Flask continues running
- **Timeout occurs**: Clean subprocess termination
- **Import errors**: Isolated to subprocess

## 🧪 Testing

### Updated Test Script
- `test_rebuilt_api.py` now checks for subprocess script
- Tests subprocess-based execution
- Validates file creation and GPU usage

### Manual Testing
```bash
# Start server
cd roshambo_api
python start_api.py

# Test health endpoint
curl http://127.0.0.1:5000/health

# Test similarity endpoint
curl -X POST http://127.0.0.1:5000/similarity \
  -H "Content-Type: application/json" \
  -d '{
    "reference_file": "reference.sdf",
    "dataset_file": "dataset.smi",
    "gpu_id": 1,
    "working_dir": "test_data"
  }'
```

## 📋 Benefits of Subprocess Approach

### ✅ **Thread Safety**
- No signal module usage
- Works perfectly in Flask threads
- No threading-related errors

### ✅ **Process Isolation**
- Roshambo crashes don't affect Flask
- Better memory management
- Clean resource cleanup

### ✅ **Timeout Handling**
- Built-in subprocess timeout
- No signal complications
- Clean timeout error messages

### ✅ **Scalability**
- Multiple requests can run roshambo simultaneously
- Each in its own process
- No shared state issues

### ✅ **Debugging**
- Easier to debug subprocess issues
- Clear separation of concerns
- Better error reporting

## 🎯 GPU Management Maintained

The subprocess approach maintains all GPU preferences:
- **Still uses GPU 1 or 2 by default (not GPU 0)**
- GPU selection passed via JSON parameters
- GPU status reported in responses

## 🔧 Usage

### Starting the Server
```bash
cd roshambo_api
python start_api.py
```

### API Calls (Same as Before)
```python
import requests

response = requests.post("http://127.0.0.1:5000/similarity", json={
    "reference_file": "reference.sdf",
    "dataset_file": "dataset.smi",
    "gpu_id": 1,  # Still prefers GPU 1/2
    "working_dir": "working_dir"
})
```

## ✅ Problem Solved

The **HTTP 500 error** caused by `signal only works in main thread` is now **completely fixed**:

1. ✅ **No more signal module usage**
2. ✅ **Thread-safe subprocess execution**
3. ✅ **Proper timeout handling**
4. ✅ **Process isolation for stability**
5. ✅ **Maintained GPU preferences (1/2, not 0)**
6. ✅ **Comprehensive error handling**

The API now works reliably in Flask's threaded environment!
