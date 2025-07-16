#!/usr/bin/env python3
"""
DockStream Debugging Script

This script helps debug issues with the DockStream scoring component by:
1. Testing the DockStream configuration and environment
2. Running a simple docking test
3. Checking the output format
4. Providing detailed diagnostics

Hardcoded paths for protacinvent setup with spaces in directory names.
"""

import subprocess
import os
import json
import tempfile
import sys
from pathlib import Path

def test_environment(env_path):
    """Test if the environment path is valid"""
    print(f"🔍 Testing environment path: {env_path}")
    
    if not os.path.exists(env_path):
        print(f"❌ Environment path does not exist: {env_path}")
        return False
    
    # Test if it's executable
    if not os.access(env_path, os.X_OK):
        print(f"❌ Environment path is not executable: {env_path}")
        return False
    
    # Test basic Python execution
    try:
        result = subprocess.run([env_path, "--version"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print(f"✅ Environment is working: {result.stdout.strip()}")
            return True
        else:
            print(f"❌ Environment test failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"❌ Environment test error: {e}")
        return False

def test_docker_script(docker_path, env_path):
    """Test if the docker script is accessible"""
    print(f"🔍 Testing docker script: {docker_path}")
    
    if not os.path.exists(docker_path):
        print(f"❌ Docker script does not exist: {docker_path}")
        return False
    
    # Test if we can import the script
    try:
        result = subprocess.run([env_path, docker_path, "--help"], 
                              capture_output=True, text=True, timeout=30)
        if result.returncode == 0 or "usage:" in result.stdout.lower() or "usage:" in result.stderr.lower():
            print(f"✅ Docker script is accessible")
            return True
        else:
            print(f"❌ Docker script test failed:")
            print(f"   STDOUT: {result.stdout[:200]}")
            print(f"   STDERR: {result.stderr[:200]}")
            return False
    except Exception as e:
        print(f"❌ Docker script test error: {e}")
        return False

def test_configuration(config_path):
    """Test if the configuration file is valid"""
    print(f"🔍 Testing configuration: {config_path}")
    
    if not os.path.exists(config_path):
        print(f"❌ Configuration file does not exist: {config_path}")
        return False, None
    
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        print(f"✅ Configuration file is valid JSON")
        
        # Check for required fields
        required_fields = ["docking_runs"]
        for field in required_fields:
            if field not in config:
                print(f"⚠️  Missing required field: {field}")
        
        # Check docking runs
        if "docking_runs" in config:
            for i, run in enumerate(config["docking_runs"]):
                print(f"   Docking run {i}: backend = {run.get('backend', 'unknown')}")
                
                # Check for receptor file
                if "parameters" in run and "receptor_pdbqt_path" in run["parameters"]:
                    receptor_paths = run["parameters"]["receptor_pdbqt_path"]
                    if isinstance(receptor_paths, list):
                        for receptor_path in receptor_paths:
                            if os.path.exists(receptor_path):
                                print(f"   ✅ Receptor file found: {receptor_path}")
                            else:
                                print(f"   ❌ Receptor file missing: {receptor_path}")
                    else:
                        if os.path.exists(receptor_paths):
                            print(f"   ✅ Receptor file found: {receptor_paths}")
                        else:
                            print(f"   ❌ Receptor file missing: {receptor_paths}")
        
        return True, config
    except json.JSONDecodeError as e:
        print(f"❌ Configuration file is not valid JSON: {e}")
        return False, None
    except Exception as e:
        print(f"❌ Configuration test error: {e}")
        return False, None

def test_simple_docking(env_path, docker_path, config_path):
    """Test simple docking with a basic molecule"""
    print(f"🔍 Testing simple docking...")
    
    # Simple test molecules
    test_smiles = ["CCO", "CC(=O)O"]  # Ethanol and acetic acid
    smiles_string = ";".join(test_smiles)
    
    # Create command
    command = [
        env_path,
        docker_path,
        "-conf", config_path,
        "-smiles", smiles_string,
        "-print_scores"
    ]
    
    print(f"   Command: {' '.join(command)}")
    
    try:
        # Run with timeout
        result = subprocess.run(command, capture_output=True, text=True, timeout=300)
        
        print(f"   Return code: {result.returncode}")
        print(f"   STDOUT length: {len(result.stdout)} characters")
        print(f"   STDERR length: {len(result.stderr)} characters")
        
        if result.stdout:
            print(f"   STDOUT preview: {result.stdout[:300]}...")
        
        if result.stderr:
            print(f"   STDERR preview: {result.stderr[:300]}...")
        
        # Try to parse scores from output
        lines = result.stdout.strip().split('\n')
        scores = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('INFO'):
                try:
                    score = float(line)
                    scores.append(score)
                except ValueError:
                    continue
        
        print(f"   Parsed scores: {scores}")
        
        if len(scores) >= len(test_smiles):
            print(f"✅ Docking test successful - got {len(scores)} scores for {len(test_smiles)} molecules")
            return True
        else:
            print(f"⚠️  Docking test incomplete - got {len(scores)} scores for {len(test_smiles)} molecules")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"❌ Docking test timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"❌ Docking test error: {e}")
        return False

def main():
    # Hardcoded paths to avoid issues with spaces in command line arguments
    config_path = "/home/protacinvent/Desktop/Getting Started/protac-invent/DockStream-master/example.json"
    docker_path = "/home/protacinvent/Desktop/Getting Started/protac-invent/DockStream-master/docker.py"
    env_path = "/home/protacinvent/.conda/envs/DockStream-master/bin/python"
    
    print("🚀 DockStream Debugging Script")
    print("=" * 50)
    print("Using hardcoded paths:")
    print(f"  Config: {config_path}")
    print(f"  Docker: {docker_path}")
    print(f"  Environment: {env_path}")
    print()
    
    # Test 1: Environment
    env_ok = test_environment(env_path)
    print()
    
    # Test 2: Docker script
    docker_ok = test_docker_script(docker_path, env_path)
    print()
    
    # Test 3: Configuration
    config_ok, config = test_configuration(config_path)
    print()
    
    # Test 4: Simple docking (only if previous tests pass)
    if env_ok and docker_ok and config_ok:
        docking_ok = test_simple_docking(env_path, docker_path, config_path)
        print()
    else:
        print("⏭️  Skipping docking test due to previous failures")
        docking_ok = False
    
    # Summary
    print("📋 Summary")
    print("=" * 50)
    print(f"Environment: {'✅' if env_ok else '❌'}")
    print(f"Docker script: {'✅' if docker_ok else '❌'}")
    print(f"Configuration: {'✅' if config_ok else '❌'}")
    print(f"Simple docking: {'✅' if docking_ok else '❌'}")
    
    if all([env_ok, docker_ok, config_ok, docking_ok]):
        print("\n🎉 All tests passed! DockStream should work correctly.")
        print("\n💡 Next steps:")
        print("   1. Run your REINVENT scoring with debug mode enabled")
        print("   2. Check the detailed output for any issues")
        print("   3. For scaffold.csv generation, use reinforcement_learning mode (not scoring mode)")
    else:
        print("\n❌ Some tests failed. Please fix the issues above.")
        print("\n💡 Common solutions:")
        print("   1. Check file paths are correct and accessible")
        print("   2. Ensure DockStream environment is properly set up")
        print("   3. Verify receptor files exist and are in correct format")
        print("   4. Test DockStream manually outside of REINVENT first")

if __name__ == "__main__":
    main()
