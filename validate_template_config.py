#!/usr/bin/env python3
"""
Template.json Configuration Validator for Roshambo Components
This script validates that the template.json has the required mandatory parameters
for Roshambo shape similarity components.
"""

import json
import os
import sys
from pathlib import Path


def validate_roshambo_config(config_data):
    """
    Validate Roshambo configuration in template.json
    
    Args:
        config_data: Parsed JSON configuration
        
    Returns:
        Tuple of (is_valid, issues_list)
    """
    issues = []
    is_valid = True
    
    # Find Roshambo components
    roshambo_components = []
    
    try:
        scoring_params = config_data.get("parameters", {}).get("scoring_function", {}).get("parameters", [])
        
        for component in scoring_params:
            component_type = component.get("component_type", "")
            if "roshambo" in component_type.lower():
                roshambo_components.append(component)
    
    except Exception as e:
        issues.append(f"Error parsing configuration structure: {e}")
        return False, issues
    
    if not roshambo_components:
        issues.append("No Roshambo components found in configuration")
        return True, issues  # Not an error if no Roshambo components
    
    print(f"Found {len(roshambo_components)} Roshambo component(s)")
    
    # Validate each Roshambo component
    for i, component in enumerate(roshambo_components):
        component_name = component.get("name", f"Component {i+1}")
        specific_params = component.get("specific_parameters", {})
        
        print(f"\nValidating: {component_name}")
        
        # Check mandatory parameters
        mandatory_params = {
            "rdbase_path": "Path to RDKit build directory",
            "conda_env_name": "Name of conda environment with Roshambo"
        }
        
        for param, description in mandatory_params.items():
            if param not in specific_params or not specific_params[param]:
                issues.append(f"{component_name}: Missing mandatory parameter '{param}' ({description})")
                is_valid = False
            else:
                print(f"  ✓ {param}: {specific_params[param]}")
        
        # Validate rdbase_path if provided
        rdbase_path = specific_params.get("rdbase_path")
        if rdbase_path:
            if not os.path.exists(rdbase_path):
                issues.append(f"{component_name}: RDBASE path does not exist: {rdbase_path}")
                is_valid = False
            else:
                # Check for required subdirectories
                required_dirs = ["lib", "Code", "Data"]
                for req_dir in required_dirs:
                    dir_path = os.path.join(rdbase_path, req_dir)
                    if not os.path.exists(dir_path):
                        issues.append(f"{component_name}: Missing required directory in RDBASE: {dir_path}")
                        is_valid = False
                
                if all(os.path.exists(os.path.join(rdbase_path, d)) for d in required_dirs):
                    print(f"  ✓ RDBASE path validated: {rdbase_path}")
        
        # Check reference files
        reference_files = [
            ("reference_file", "Main reference file"),
            ("warhead1_reference", "Warhead 1 reference"),
            ("warhead2_reference", "Warhead 2 reference")
        ]
        
        for ref_param, description in reference_files:
            ref_file = specific_params.get(ref_param)
            if ref_file:
                if not os.path.exists(ref_file):
                    issues.append(f"{component_name}: Reference file does not exist: {ref_file}")
                    is_valid = False
                else:
                    print(f"  ✓ {description}: {ref_file}")
        
        # Check output directory
        overlays_dir = specific_params.get("overlays_dir")
        if overlays_dir:
            parent_dir = os.path.dirname(overlays_dir)
            if parent_dir and not os.path.exists(parent_dir):
                issues.append(f"{component_name}: Parent directory for overlays does not exist: {parent_dir}")
                is_valid = False
            else:
                print(f"  ✓ Overlays directory path: {overlays_dir}")
        
        # Check optional but important parameters
        important_params = {
            "shape_weight": "Shape similarity weight",
            "color_weight": "Color similarity weight",
            "n_confs": "Number of conformers",
            "gpu_id": "GPU ID"
        }
        
        for param, description in important_params.items():
            if param in specific_params:
                print(f"  ✓ {description}: {specific_params[param]}")
    
    return is_valid, issues


def main():
    """Main validation function"""
    print("=" * 60)
    print("Template.json Roshambo Configuration Validator")
    print("=" * 60)
    
    # Check if template.json exists
    template_file = "template.json"
    if not os.path.exists(template_file):
        print(f"❌ Error: {template_file} not found in current directory")
        print("Please run this script from the directory containing template.json")
        sys.exit(1)
    
    # Load and parse template.json
    try:
        with open(template_file, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        print(f"✓ Successfully loaded {template_file}")
    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in {template_file}: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading {template_file}: {e}")
        sys.exit(1)
    
    # Validate configuration
    is_valid, issues = validate_roshambo_config(config_data)
    
    print("\n" + "=" * 60)
    print("Validation Results")
    print("=" * 60)
    
    if is_valid:
        print("🎉 Configuration is valid!")
        print("✓ All mandatory parameters are present")
        print("✓ File paths are accessible")
        print("✓ Ready for Roshambo execution")
    else:
        print("❌ Configuration has issues:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nPlease fix these issues before running Reinvent with Roshambo.")
    
    # Additional recommendations
    print("\n" + "=" * 60)
    print("Recommendations")
    print("=" * 60)
    
    print("1. Test your environment setup:")
    print("   conda activate <your_conda_env>")
    print("   python -c \"from rdkit import Chem; import roshambo; print('Environment OK')\"")
    
    print("\n2. Verify GPU availability (if using GPU):")
    print("   python -c \"import torch; print('CUDA available:', torch.cuda.is_available())\"")
    
    print("\n3. Enable debug mode for detailed output:")
    print("   Set \"debug\": true in your Roshambo component configuration")
    
    print("\n4. Test with a small dataset first before full runs")
    
    return 0 if is_valid else 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
