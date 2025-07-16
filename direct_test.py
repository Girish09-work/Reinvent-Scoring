#!/usr/bin/env python3
"""
Direct test of the roshambo component to verify the fix.
"""

import sys
import os

# Add the current directory to Python path
sys.path.insert(0, os.getcwd())

# Force reload of the module
if 'reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity' in sys.modules:
    del sys.modules['reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity']

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity


def main():
    print("🔧 Direct Test of Roshambo Component Fix")
    print("=" * 50)
    
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
    print("📦 Creating RoshamboShapeSimilarity component...")
    component = RoshamboShapeSimilarity(params)
    print("✅ Component created successfully")
    
    # Test the _extract_scores_from_csv method
    print("\n📄 Testing _extract_scores_from_csv method...")
    
    # Create a temporary CSV file
    import tempfile
    csv_content = """Molecule	OriginalName	ComboTanimoto	ShapeTanimoto	ColorTanimoto
mol_0_0	mol_0	0.8	0.7	0.9
mol_2_0	mol_2	0.6	0.5	0.7
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        csv_file = f.name
    
    try:
        print(f"📊 Calling _extract_scores_from_csv with expected_count=5...")
        scores = component._extract_scores_from_csv(csv_file, 5)
        print(f"✅ Method call successful!")
        print(f"   Returned {len(scores)} scores: {scores}")
        
        if len(scores) == 5:
            print("✅ SUCCESS: Correct number of scores returned!")
        else:
            print(f"❌ FAILURE: Expected 5 scores, got {len(scores)}")
            
    except Exception as e:
        print(f"❌ ERROR: {e}")
        print(f"   Error type: {type(e).__name__}")
        import traceback
        traceback.print_exc()
    finally:
        os.unlink(csv_file)
    
    # Test the _call_roshambo_api method
    print("\n📞 Testing _call_roshambo_api method...")
    try:
        print(f"📊 Calling _call_roshambo_api with expected_count=3...")
        scores = component._call_roshambo_api("dummy_ref", "dummy_dataset", "dummy_folder", 3)
        print(f"✅ Method call successful!")
        print(f"   Returned {len(scores)} scores: {scores}")
        
        if len(scores) == 3:
            print("✅ SUCCESS: Correct number of scores returned!")
        else:
            print(f"❌ FAILURE: Expected 3 scores, got {len(scores)}")
            
    except Exception as e:
        print(f"❌ ERROR: {e}")
        print(f"   Error type: {type(e).__name__}")
        import traceback
        traceback.print_exc()
    
    # Test the main calculate_score method
    print("\n🧮 Testing calculate_score method...")
    test_molecules = ["CCO", "CC(=O)O", "c1ccccc1"]
    
    try:
        print(f"📊 Calling calculate_score with {len(test_molecules)} molecules...")
        score_summary = component.calculate_score(test_molecules, step=0)
        print(f"✅ Method call successful!")
        print(f"   Input molecules: {len(test_molecules)}")
        print(f"   Returned scores: {len(score_summary.total_score)}")
        print(f"   Score values: {score_summary.total_score}")
        
        if len(score_summary.total_score) == len(test_molecules):
            print("✅ SUCCESS: Score length matches input length!")
            print("🎉 The AssertionError fix is working correctly!")
        else:
            print(f"❌ FAILURE: Expected {len(test_molecules)} scores, got {len(score_summary.total_score)}")
            
    except Exception as e:
        print(f"❌ ERROR: {e}")
        print(f"   Error type: {type(e).__name__}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
