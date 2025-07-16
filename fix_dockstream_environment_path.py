#!/usr/bin/env python3
"""
Fix DockStream Environment Path

This script helps fix the environment_path issue in your REINVENT DockStream configuration.
The issue is that REINVENT is trying to execute the conda environment directory instead of the Python executable.
"""

import json
import os
import sys

def find_and_fix_config_files():
    """Find and fix DockStream configuration files"""
    
    print("🔧 DockStream Environment Path Fixer")
    print("=" * 50)
    
    # The correct environment path
    correct_env_path = "/home/protacinvent/.conda/envs/DockStream-master/bin/python"
    incorrect_env_path = "/home/protacinvent/.conda/envs/DockStream-master"
    
    print(f"❌ Incorrect path: {incorrect_env_path}")
    print(f"✅ Correct path: {correct_env_path}")
    print()
    
    # Check if the correct path exists
    if not os.path.exists(correct_env_path):
        print(f"⚠️  Warning: The correct Python executable doesn't exist at: {correct_env_path}")
        print("   Please check your conda environment installation.")
        return False
    else:
        print(f"✅ Confirmed: Python executable exists at {correct_env_path}")
    
    # Look for configuration files that might need fixing
    config_files_to_check = [
        "dockstream.json",
        "reinvent_scoring/configs/example.config.json",
        "template.json"
    ]
    
    fixed_files = []
    
    for config_file in config_files_to_check:
        if os.path.exists(config_file):
            print(f"\n🔍 Checking {config_file}...")
            
            try:
                with open(config_file, 'r') as f:
                    content = f.read()
                
                # Check if the file contains the incorrect path
                if incorrect_env_path in content:
                    print(f"   ❌ Found incorrect environment path in {config_file}")
                    
                    # Create backup
                    backup_file = f"{config_file}.backup"
                    with open(backup_file, 'w') as f:
                        f.write(content)
                    print(f"   💾 Created backup: {backup_file}")
                    
                    # Fix the path
                    fixed_content = content.replace(incorrect_env_path, correct_env_path)
                    
                    with open(config_file, 'w') as f:
                        f.write(fixed_content)
                    
                    print(f"   ✅ Fixed environment path in {config_file}")
                    fixed_files.append(config_file)
                    
                elif correct_env_path in content:
                    print(f"   ✅ {config_file} already has correct environment path")
                else:
                    print(f"   ℹ️  {config_file} doesn't contain DockStream environment path")
                    
            except Exception as e:
                print(f"   ❌ Error processing {config_file}: {e}")
        else:
            print(f"   ℹ️  {config_file} not found")
    
    return fixed_files

def create_test_config():
    """Create a test configuration with correct paths"""
    
    test_config = {
        "component_type": "dockstream",
        "name": "ADV",
        "weight": 1.0,
        "specific_parameters": {
            "configuration_path": "/home/protacinvent/Desktop/Getting Started/protac-invent/DockStream-master/example.json",
            "docker_script_path": "/home/protacinvent/Desktop/Getting Started/protac-invent/DockStream-master/docker.py",
            "environment_path": "/home/protacinvent/.conda/envs/DockStream-master/bin/python",
            "debug": True,
            "transformation": {
                "transformation_type": "reverse_sigmoid",
                "low": -12,
                "high": -6,
                "k": 0.5
            }
        }
    }
    
    with open("dockstream_component_config_fixed.json", 'w') as f:
        json.dump(test_config, f, indent=4)
    
    print(f"\n📝 Created test configuration: dockstream_component_config_fixed.json")
    print("   Use this configuration in your REINVENT scoring function.")

def main():
    """Main function"""
    
    fixed_files = find_and_fix_config_files()
    
    create_test_config()
    
    print("\n" + "=" * 50)
    print("📋 Summary")
    
    if fixed_files:
        print(f"✅ Fixed {len(fixed_files)} configuration files:")
        for file in fixed_files:
            print(f"   - {file}")
        print("\n💡 Next steps:")
        print("   1. Use the fixed configuration files")
        print("   2. Test your REINVENT scoring again")
        print("   3. You should now get proper docking scores instead of zeros")
    else:
        print("ℹ️  No files needed fixing, but created a test configuration.")
        print("\n💡 Next steps:")
        print("   1. Check your REINVENT configuration manually")
        print("   2. Ensure environment_path points to the Python executable:")
        print("      /home/protacinvent/.conda/envs/DockStream-master/bin/python")
        print("   3. NOT the environment directory:")
        print("      /home/protacinvent/.conda/envs/DockStream-master")
    
    print("\n🎯 Key Point:")
    print("   The environment_path must point to the Python EXECUTABLE, not the environment DIRECTORY")
    print("   ❌ Wrong: /path/to/conda/envs/DockStream-master")
    print("   ✅ Right: /path/to/conda/envs/DockStream-master/bin/python")

if __name__ == "__main__":
    main()
