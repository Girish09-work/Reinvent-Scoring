#!/usr/bin/env python3
"""
Debug script to test the Roshambo API data format.
This script shows what data is being sent to the API.
"""

import os
import tempfile
import shutil
from pathlib import Path

def test_api_data_format():
    """Test the API data format that should be sent to Roshambo."""
    
    # Create a temporary working directory
    temp_dir = tempfile.mkdtemp()
    epoch_folder = os.path.join(temp_dir, "epoch_0")
    Path(epoch_folder).mkdir(parents=True, exist_ok=True)
    
    print(f"📁 Created test directory: {epoch_folder}")
    
    try:
        # Create dummy files
        reference_filename = "reference.sdf"
        dataset_filename = "dataset_0.sdf"
        
        reference_path = os.path.join(epoch_folder, reference_filename)
        dataset_path = os.path.join(epoch_folder, dataset_filename)
        
        # Create dummy SDF files
        with open(reference_path, "w") as f:
            f.write("""
  Mrv2014 01010101

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
""")
        
        with open(dataset_path, "w") as f:
            f.write("""
  Mrv2014 01010101

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
""")
        
        print(f"📄 Created reference file: {reference_path}")
        print(f"📄 Created dataset file: {dataset_path}")
        
        # This is the CORRECT API data format
        correct_api_data = {
            "reference_file": reference_filename,  # FILENAME ONLY
            "dataset_file": dataset_filename,      # FILENAME ONLY
            "ignore_hs": True,
            "n_confs": 50,
            "use_carbon_radii": True,
            "color": True,
            "sort_by": "ComboTanimoto",
            "write_to_file": True,
            "gpu_id": 0,
            "working_dir": epoch_folder  # FULL PATH TO WORKING DIRECTORY
        }
        
        print("\n✅ CORRECT API Data Format:")
        for key, value in correct_api_data.items():
            print(f"   {key}: {value}")
        
        # This is the INCORRECT format (what was happening before)
        incorrect_api_data = {
            "reference_file": reference_path,  # FULL PATH - WRONG!
            "dataset_file": dataset_path,      # FULL PATH - WRONG!
            "ignore_hs": True,
            "n_confs": 50,
            "use_carbon_radii": True,
            "color": True,
            "sort_by": "ComboTanimoto",
            "write_to_file": True,
            "gpu_id": 0,
            "working_dir": epoch_folder
        }
        
        print("\n❌ INCORRECT API Data Format (what was happening before):")
        for key, value in incorrect_api_data.items():
            print(f"   {key}: {value}")
        
        print("\n📋 Key Points:")
        print("   - reference_file and dataset_file should be FILENAMES ONLY")
        print("   - working_dir should be the FULL PATH to the directory containing the files")
        print("   - Roshambo API will look for the files in the working_dir")
        
    finally:
        # Clean up
        shutil.rmtree(temp_dir)
        print(f"\n🧹 Cleaned up: {temp_dir}")

if __name__ == "__main__":
    test_api_data_format()
