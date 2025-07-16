#!/usr/bin/env python3
"""
Check if the roshambo fix is properly applied.
"""

import sys
import os
import shutil

def clear_cache():
    """Clear Python cache files."""
    print("🧹 Clearing Python cache...")
    
    # Remove __pycache__ directories
    for root, dirs, files in os.walk('.'):
        for dir_name in dirs:
            if dir_name == '__pycache__':
                cache_path = os.path.join(root, dir_name)
                print(f"   Removing: {cache_path}")
                shutil.rmtree(cache_path, ignore_errors=True)
    
    # Clear sys.modules
    modules_to_remove = []
    for module_name in sys.modules:
        if 'roshambo' in module_name:
            modules_to_remove.append(module_name)
    
    for module_name in modules_to_remove:
        print(f"   Removing from sys.modules: {module_name}")
        del sys.modules[module_name]


def check_method_signature():
    """Check the method signature after clearing cache."""
    print("\n🔍 Checking method signature...")
    
    try:
        from reinvent_scoring.scoring.component_parameters import ComponentParameters
        from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity
        
        params = ComponentParameters(
            component_type="roshambo_shape_similarity",
            name="Test",
            weight=1.0,
            specific_parameters={
                "reference_file": "test.sdf",
                "roshambo_api_url": "http://127.0.0.1:5000",
                "overlays_dir": "./test",
                "shape_weight": 0.5,
                "color_weight": 0.5
            }
        )
        
        component = RoshamboShapeSimilarity(params)
        
        import inspect
        sig = inspect.signature(component._extract_scores_from_csv)
        params_list = list(sig.parameters.keys())
        
        print(f"   Method signature: {sig}")
        print(f"   Parameters: {params_list}")
        
        if 'expected_count' in params_list:
            print("   ✅ Fix is applied correctly!")
            return True
        else:
            print("   ❌ Fix is NOT applied - expected_count parameter missing")
            return False
            
    except Exception as e:
        print(f"   ❌ Error checking signature: {e}")
        return False


def test_method_call():
    """Test calling the method with the new signature."""
    print("\n🧪 Testing method call...")
    
    try:
        from reinvent_scoring.scoring.component_parameters import ComponentParameters
        from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity
        
        params = ComponentParameters(
            component_type="roshambo_shape_similarity",
            name="Test",
            weight=1.0,
            specific_parameters={
                "reference_file": "test.sdf",
                "roshambo_api_url": "http://127.0.0.1:5000",
                "overlays_dir": "./test",
                "shape_weight": 0.5,
                "color_weight": 0.5,
                "debug": True
            }
        )
        
        component = RoshamboShapeSimilarity(params)
        
        # Create a temporary CSV file
        import tempfile
        csv_content = "Molecule\tShapeTanimoto\tColorTanimoto\nmol_0_0\t0.8\t0.9\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            csv_file = f.name
        
        try:
            # Test the method call
            scores = component._extract_scores_from_csv(csv_file, 3)
            print(f"   ✅ Method call successful!")
            print(f"   Returned {len(scores)} scores: {scores}")
            
            if len(scores) == 3:
                print("   ✅ Correct number of scores returned!")
                return True
            else:
                print(f"   ❌ Expected 3 scores, got {len(scores)}")
                return False
                
        finally:
            os.unlink(csv_file)
            
    except Exception as e:
        print(f"   ❌ Error testing method: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("🔧 Checking Roshambo Fix Application")
    print("=" * 50)
    
    clear_cache()
    signature_ok = check_method_signature()
    
    if signature_ok:
        test_ok = test_method_call()
        if test_ok:
            print("\n🎉 SUCCESS: The fix is properly applied and working!")
        else:
            print("\n❌ FAILURE: The fix is applied but not working correctly.")
    else:
        print("\n❌ FAILURE: The fix is not applied. You may need to:")
        print("   1. Restart your Python session")
        print("   2. Check if the file was saved correctly")
        print("   3. Verify you're running from the correct directory")


if __name__ == "__main__":
    main()
