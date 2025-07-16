#!/usr/bin/env python3
"""
Debug script to help diagnose the roshambo subprocess issue.
This script tests the roshambo environment setup and execution.
"""

import os
import subprocess
import tempfile
import json

def test_conda_environment(conda_env_name="roshambo"):
    """Test if the conda environment exists and is accessible."""
    print("🔍 Testing Conda Environment")
    print("=" * 50)
    
    # Test conda command
    try:
        result = subprocess.run("conda --version", shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ Conda is available: {result.stdout.strip()}")
        else:
            print(f"❌ Conda command failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"❌ Error running conda command: {e}")
        return False
    
    # Test conda environment
    try:
        result = subprocess.run(f"conda env list", shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            envs = result.stdout
            if conda_env_name in envs:
                print(f"✅ Conda environment '{conda_env_name}' exists")
            else:
                print(f"❌ Conda environment '{conda_env_name}' not found")
                print("Available environments:")
                print(envs)
                return False
        else:
            print(f"❌ Failed to list conda environments: {result.stderr}")
            return False
    except Exception as e:
        print(f"❌ Error listing conda environments: {e}")
        return False
    
    # Test conda run command
    try:
        result = subprocess.run(f"conda run -n {conda_env_name} python --version", 
                              shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ Can execute Python in '{conda_env_name}': {result.stdout.strip()}")
        else:
            print(f"❌ Failed to run Python in '{conda_env_name}': {result.stderr}")
            return False
    except Exception as e:
        print(f"❌ Error testing conda run: {e}")
        return False
    
    return True

def test_roshambo_import(conda_env_name="roshambo"):
    """Test if roshambo can be imported in the conda environment."""
    print("\n🧪 Testing Roshambo Import")
    print("=" * 50)
    
    # Create a simple test script
    test_script = '''
import sys
print(f"Python version: {sys.version}")
print(f"Python path: {sys.executable}")

try:
    import roshambo
    print("✅ Roshambo imported successfully")
    print(f"Roshambo location: {roshambo.__file__}")
    
    # Test the API function
    from roshambo.api import get_similarity_scores
    print("✅ get_similarity_scores function imported successfully")
    
except ImportError as e:
    print(f"❌ Failed to import roshambo: {e}")
    sys.exit(1)
except Exception as e:
    print(f"❌ Error testing roshambo: {e}")
    sys.exit(1)
'''
    
    # Write test script to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
        f.write(test_script)
        script_path = f.name
    
    try:
        # Execute test script in conda environment
        result = subprocess.run(f"conda run -n {conda_env_name} python {script_path}", 
                              shell=True, capture_output=True, text=True)
        
        print("Command output:")
        print(result.stdout)
        
        if result.returncode != 0:
            print("Command errors:")
            print(result.stderr)
            return False
        else:
            print("✅ Roshambo test passed")
            return True
            
    except Exception as e:
        print(f"❌ Error running roshambo test: {e}")
        return False
    finally:
        # Clean up
        try:
            os.unlink(script_path)
        except:
            pass

def test_roshambo_execution(conda_env_name="roshambo", rdbase_path="/home/protacinvent/Desktop/roshambo/rdkit"):
    """Test a simple roshambo execution."""
    print("\n🚀 Testing Roshambo Execution")
    print("=" * 50)
    
    # Create a test script that mimics the actual roshambo call
    test_script = f'''
import os
import sys
import json

# Set up environment variables
os.environ["RDBASE"] = "{rdbase_path}"
os.environ["RDKIT_LIB_DIR"] = "{rdbase_path}/lib"
os.environ["RDKIT_INCLUDE_DIR"] = "{rdbase_path}/Code"
os.environ["RDKIT_DATA_DIR"] = "{rdbase_path}/Data"
os.environ["CUDA_HOME"] = "/usr/local/cuda"

# Set up PYTHONPATH
current_pythonpath = os.environ.get("PYTHONPATH", "")
new_pythonpath = "{rdbase_path}"
if current_pythonpath:
    new_pythonpath += ":" + current_pythonpath
os.environ["PYTHONPATH"] = new_pythonpath

# Set up LD_LIBRARY_PATH
current_ld_path = os.environ.get("LD_LIBRARY_PATH", "")
new_ld_path = "{rdbase_path}/lib"
if current_ld_path:
    new_ld_path += ":" + current_ld_path
os.environ["LD_LIBRARY_PATH"] = new_ld_path

# Add RDKit to Python path
sys.path.insert(0, "{rdbase_path}")

print("Environment setup:")
print(f"RDBASE: {{os.environ.get('RDBASE')}}")
print(f"PYTHONPATH: {{os.environ.get('PYTHONPATH')}}")
print(f"LD_LIBRARY_PATH: {{os.environ.get('LD_LIBRARY_PATH')}}")

try:
    from roshambo.api import get_similarity_scores
    print("✅ Successfully imported get_similarity_scores")
    
    # Test with minimal parameters (this will fail but should give us better error info)
    print("Testing roshambo API call...")
    result = get_similarity_scores(
        ref_file="dummy.sdf",
        dataset_file="dummy.sdf",
        write_to_file=True
    )
    print(f"Result: {{result}}")
    
except Exception as e:
    import traceback
    print(f"❌ Error: {{e}}")
    print(f"Traceback: {{traceback.format_exc()}}")
    sys.exit(1)
'''
    
    # Write test script to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
        f.write(test_script)
        script_path = f.name
    
    try:
        # Execute test script in conda environment
        result = subprocess.run(f"conda run -n {conda_env_name} python {script_path}", 
                              shell=True, capture_output=True, text=True)
        
        print("Command output:")
        print(result.stdout)
        
        if result.returncode != 0:
            print("Command errors:")
            print(result.stderr)
            return False
        else:
            print("✅ Basic roshambo execution test passed")
            return True
            
    except Exception as e:
        print(f"❌ Error running roshambo execution test: {e}")
        return False
    finally:
        # Clean up
        try:
            os.unlink(script_path)
        except:
            pass

def main():
    """Run all diagnostic tests."""
    print("🔧 Roshambo Environment Diagnostic Tool")
    print("=" * 60)
    
    # Test parameters (adjust these to match your setup)
    conda_env_name = "roshambo"
    rdbase_path = "/home/protacinvent/Desktop/roshambo/rdkit"
    
    print(f"Testing conda environment: {conda_env_name}")
    print(f"Testing RDBASE path: {rdbase_path}")
    print()
    
    # Run tests
    tests_passed = 0
    total_tests = 3
    
    if test_conda_environment(conda_env_name):
        tests_passed += 1
    
    if test_roshambo_import(conda_env_name):
        tests_passed += 1
    
    if test_roshambo_execution(conda_env_name, rdbase_path):
        tests_passed += 1
    
    # Summary
    print("\n" + "=" * 60)
    print(f"🏁 Test Results: {tests_passed}/{total_tests} tests passed")
    
    if tests_passed == total_tests:
        print("✅ All tests passed! Roshambo environment should work.")
        print("\n💡 If you're still having issues, the problem might be:")
        print("   - File paths in the roshambo call")
        print("   - Working directory permissions")
        print("   - SDF file format compatibility")
    else:
        print("❌ Some tests failed. Please fix the issues above.")
        print("\n💡 Common solutions:")
        print("   - Install roshambo: conda install -n roshambo roshambo")
        print("   - Check conda environment name")
        print("   - Verify RDBASE path points to valid RDKit installation")
        print("   - Ensure all required dependencies are installed")
    
    return tests_passed == total_tests

if __name__ == "__main__":
    main()
