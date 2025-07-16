#!/usr/bin/env python3
"""
Verify that the roshambo component has the correct method signature.
"""

import inspect
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity


def verify_method_signatures():
    """Verify that the method signatures are correct."""
    print("🔍 Verifying Roshambo Component Method Signatures")
    print("=" * 60)
    
    # Create test parameters
    params = ComponentParameters(
        component_type="roshambo_shape_similarity",
        name="RoshamboTest",
        weight=1.0,
        specific_parameters={
            "reference_file": "test_ref.sdf",
            "roshambo_api_url": "http://127.0.0.1:5000",
            "overlays_dir": "./test_overlays",
            "shape_weight": 0.5,
            "color_weight": 0.5,
            "debug": True
        }
    )
    
    # Create component
    component = RoshamboShapeSimilarity(params)
    
    # Check _extract_scores_from_csv method signature
    print("📋 Checking _extract_scores_from_csv method:")
    method = getattr(component, '_extract_scores_from_csv')
    sig = inspect.signature(method)
    print(f"   Signature: {sig}")
    print(f"   Parameters: {list(sig.parameters.keys())}")
    
    expected_params = ['csv_file', 'expected_count']
    actual_params = [p for p in sig.parameters.keys() if p != 'self']
    
    if actual_params == expected_params:
        print("   ✅ Method signature is correct!")
    else:
        print(f"   ❌ Method signature is wrong!")
        print(f"      Expected: {expected_params}")
        print(f"      Actual:   {actual_params}")
        return False
    
    # Check _call_roshambo_api method signature
    print("\n📋 Checking _call_roshambo_api method:")
    method = getattr(component, '_call_roshambo_api')
    sig = inspect.signature(method)
    print(f"   Signature: {sig}")
    print(f"   Parameters: {list(sig.parameters.keys())}")
    
    expected_params = ['reference_file', 'dataset_file', 'epoch_folder', 'expected_count']
    actual_params = [p for p in sig.parameters.keys() if p != 'self']
    
    if actual_params == expected_params:
        print("   ✅ Method signature is correct!")
    else:
        print(f"   ❌ Method signature is wrong!")
        print(f"      Expected: {expected_params}")
        print(f"      Actual:   {actual_params}")
        return False
    
    return True


def test_method_calls():
    """Test that the methods can be called with correct parameters."""
    print("\n🧪 Testing Method Calls")
    print("=" * 60)
    
    # Create test parameters
    params = ComponentParameters(
        component_type="roshambo_shape_similarity",
        name="RoshamboTest",
        weight=1.0,
        specific_parameters={
            "reference_file": "test_ref.sdf",
            "roshambo_api_url": "http://127.0.0.1:5000",
            "overlays_dir": "./test_overlays",
            "shape_weight": 0.5,
            "color_weight": 0.5,
            "debug": True
        }
    )
    
    # Create component
    component = RoshamboShapeSimilarity(params)
    
    # Test _call_roshambo_api (should fail gracefully)
    print("📞 Testing _call_roshambo_api with expected_count:")
    try:
        scores = component._call_roshambo_api("dummy_ref", "dummy_dataset", "dummy_folder", 5)
        print(f"   ✅ Method called successfully, returned {len(scores)} scores")
        if len(scores) == 5:
            print("   ✅ Returned correct number of scores!")
        else:
            print(f"   ❌ Expected 5 scores, got {len(scores)}")
            return False
    except Exception as e:
        print(f"   ❌ Method call failed: {e}")
        return False
    
    # Test _extract_scores_from_csv with a dummy file
    print("\n📄 Testing _extract_scores_from_csv with expected_count:")
    import tempfile
    import os
    
    # Create a simple CSV file
    csv_content = "Molecule\tShapeTanimoto\tColorTanimoto\nmol_0_0\t0.8\t0.9\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        csv_file = f.name
    
    try:
        scores = component._extract_scores_from_csv(csv_file, 3)
        print(f"   ✅ Method called successfully, returned {len(scores)} scores")
        if len(scores) == 3:
            print("   ✅ Returned correct number of scores!")
            print(f"   📊 Scores: {scores}")
        else:
            print(f"   ❌ Expected 3 scores, got {len(scores)}")
            return False
    except Exception as e:
        print(f"   ❌ Method call failed: {e}")
        return False
    finally:
        os.unlink(csv_file)
    
    return True


if __name__ == "__main__":
    print("🔧 Verifying Roshambo Component Fix")
    print("=" * 70)
    
    success1 = verify_method_signatures()
    success2 = test_method_calls()
    
    print("\n" + "=" * 70)
    if success1 and success2:
        print("🎉 All verifications passed! The component is correctly fixed.")
        print("\n📋 The fix includes:")
        print("   • _call_roshambo_api now accepts expected_count parameter")
        print("   • _extract_scores_from_csv now accepts expected_count parameter")
        print("   • Both methods return exactly the expected number of scores")
        print("   • This should resolve the AssertionError in base_scoring_function.py")
    else:
        print("❌ Verification failed. There may be an issue with the fix.")
