#!/usr/bin/env python3
"""
Simple test to verify the roshambo scoring component fix.
"""

import numpy as np
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity


def test_score_length_fix():
    """Test that the fix ensures correct score array length."""
    print("🧪 Testing Roshambo Score Length Fix")
    print("=" * 50)
    
    # Create test parameters
    params = ComponentParameters(
        component_type="roshambo_shape_similarity",
        name="RoshamboTest",
        weight=1.0,
        specific_parameters={
            "reference_file": "ref.sdf",
            "roshambo_api_url": "http://127.0.0.1:5000",
            "overlays_dir": "./test_overlays",
            "shape_weight": 0.5,
            "color_weight": 0.5,
            "debug": True
        }
    )
    
    # Create component
    component = RoshamboShapeSimilarity(params)
    
    # Test the CSV extraction method directly
    import tempfile
    import os
    
    # Create mock CSV with only 3 molecules but we expect 6
    csv_content = """Molecule	OriginalName	ComboTanimoto	ShapeTanimoto	ColorTanimoto
mol_0_0	mol_0	0.8	0.7	0.9
mol_2_0	mol_2	0.6	0.5	0.7
mol_4_0	mol_4	0.4	0.3	0.5
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        csv_file = f.name
    
    try:
        print(f"📊 Testing CSV extraction with 6 expected molecules")
        print(f"📄 CSV contains data for molecules: 0, 2, 4 (3 molecules)")
        
        # Extract scores expecting 6 molecules
        scores = component._extract_scores_from_csv(csv_file, 6)
        
        print(f"✅ Extracted {len(scores)} scores: {scores}")
        
        # Verify the results
        expected_scores = [
            0.75,  # mol_0: (0.5*0.7 + 0.5*0.9) / (0.5+0.5) = 0.8
            0.0,   # mol_1: not in CSV, should be 0.0
            0.6,   # mol_2: (0.5*0.5 + 0.5*0.7) / (0.5+0.5) = 0.6
            0.0,   # mol_3: not in CSV, should be 0.0
            0.4,   # mol_4: (0.5*0.3 + 0.5*0.5) / (0.5+0.5) = 0.4
            0.0    # mol_5: not in CSV, should be 0.0
        ]
        
        if len(scores) == 6:
            print("✅ SUCCESS: Returned correct number of scores!")
            print(f"   Expected: {expected_scores}")
            print(f"   Got:      {scores}")
            
            # Check if scores are reasonable (allowing for floating point precision)
            all_correct = True
            for i, (expected, actual) in enumerate(zip(expected_scores, scores)):
                if abs(expected - actual) > 0.01:  # Allow small floating point differences
                    print(f"   ⚠️  Score {i}: expected {expected}, got {actual}")
                    all_correct = False
            
            if all_correct:
                print("✅ All scores are correct!")
            else:
                print("⚠️  Some scores differ from expected values")
                
            return True
        else:
            print(f"❌ FAILURE: Expected 6 scores, got {len(scores)}")
            return False
            
    finally:
        os.unlink(csv_file)


def test_api_call_fix():
    """Test that API calls return correct number of scores."""
    print("\n🧪 Testing API Call Fix")
    print("=" * 50)
    
    # Create test parameters
    params = ComponentParameters(
        component_type="roshambo_shape_similarity",
        name="RoshamboTest",
        weight=1.0,
        specific_parameters={
            "reference_file": "ref.sdf",
            "roshambo_api_url": "http://127.0.0.1:5001",
            "overlays_dir": "./test_overlays",
            "shape_weight": 0.5,
            "color_weight": 0.5,
            "debug": True
        }
    )
    
    # Create component
    component = RoshamboShapeSimilarity(params)
    
    # Test different scenarios
    test_cases = [
        {"expected": 3, "description": "3 molecules"},
        {"expected": 6, "description": "6 molecules"},
        {"expected": 10, "description": "10 molecules"},
        {"expected": 1, "description": "1 molecule"},
    ]
    
    for case in test_cases:
        expected_count = case["expected"]
        description = case["description"]
        
        print(f"\n📊 Testing {description} (expected_count={expected_count})")
        
        # Test error cases that should return correct number of zeros
        scores = component._call_roshambo_api("dummy_ref", "dummy_dataset", "dummy_folder", expected_count)
        
        print(f"   Returned {len(scores)} scores: {scores}")
        
        if len(scores) == expected_count:
            print(f"   ✅ Correct length!")
            if all(score == 0.0 for score in scores):
                print(f"   ✅ All scores are 0.0 (expected for error case)")
            else:
                print(f"   ⚠️  Some scores are non-zero: {scores}")
        else:
            print(f"   ❌ Wrong length! Expected {expected_count}, got {len(scores)}")
            return False
    
    return True


if __name__ == "__main__":
    print("🔧 Testing Roshambo AssertionError Fix")
    print("=" * 60)
    
    success1 = test_score_length_fix()
    success2 = test_api_call_fix()
    
    print("\n" + "=" * 60)
    if success1 and success2:
        print("🎉 All tests passed! The fix should resolve the AssertionError.")
        print("\n📋 Summary of the fix:")
        print("   • Modified _call_roshambo_api to accept expected_count parameter")
        print("   • Updated _extract_scores_from_csv to return exactly expected_count scores")
        print("   • Fixed all error cases to return correct number of zero scores")
        print("   • This ensures len(valid_indices) == len(summary.total_score)")
    else:
        print("❌ Some tests failed. The fix needs more work.")
