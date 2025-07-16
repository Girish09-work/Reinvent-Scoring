#!/usr/bin/env python3
"""
Force reload the module and test the fix.
"""

import sys
import os
import importlib

def force_reload_and_test():
    """Force reload the module and test the method signature."""
    print("🔄 Force reloading roshambo module...")

    # Remove all related modules from cache
    modules_to_remove = []
    for module_name in list(sys.modules.keys()):
        if 'roshambo' in module_name or 'reinvent_scoring' in module_name:
            modules_to_remove.append(module_name)

    for module_name in modules_to_remove:
        if module_name in sys.modules:
            print(f"   Removing {module_name} from cache")
            del sys.modules[module_name]

    # Clear __pycache__ directories
    import shutil
    for root, dirs, files in os.walk('.'):
        for dir_name in dirs:
            if dir_name == '__pycache__':
                cache_path = os.path.join(root, dir_name)
                try:
                    shutil.rmtree(cache_path)
                    print(f"   Removed cache: {cache_path}")
                except:
                    pass

    # Add the reinvent-scoring directory to Python path if needed
    reinvent_scoring_path = None
    possible_paths = [".", "../reinvent-scoring", "../../reinvent-scoring"]
    for path in possible_paths:
        if os.path.exists(os.path.join(path, "reinvent_scoring")):
            reinvent_scoring_path = os.path.abspath(path)
            break

    if reinvent_scoring_path and reinvent_scoring_path not in sys.path:
        sys.path.insert(0, reinvent_scoring_path)
        print(f"   Added to Python path: {reinvent_scoring_path}")

    # Now import fresh
    print("\n📦 Importing fresh modules...")
    try:
        from reinvent_scoring.scoring.component_parameters import ComponentParameters
        from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity

        print("✅ Import successful")

        # Create component
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

        # Check method signature
        import inspect
        sig = inspect.signature(component._extract_scores_from_csv)
        params_list = list(sig.parameters.keys())

        print(f"\n🔍 Method signature: {sig}")
        print(f"   Parameters: {params_list}")

        if 'expected_count' in params_list:
            print("✅ SUCCESS: Fix is properly applied!")

            # Test the method
            print("\n🧪 Testing method call...")
            import tempfile
            csv_content = "Molecule\tShapeTanimoto\tColorTanimoto\nmol_0_0\t0.8\t0.9\n"

            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                f.write(csv_content)
                csv_file = f.name

            try:
                scores = component._extract_scores_from_csv(csv_file, 5)
                print(f"   ✅ Method call successful!")
                print(f"   Returned {len(scores)} scores: {scores}")

                if len(scores) == 5:
                    print("   ✅ Correct number of scores!")
                    return True
                else:
                    print(f"   ❌ Expected 5 scores, got {len(scores)}")
                    return False
            finally:
                os.unlink(csv_file)

        else:
            print("❌ FAILURE: Fix is NOT applied")
            print("   The method still has the old signature")
            return False

    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def check_file_content():
    """Check the actual file content to verify the fix is there."""
    print("\n📄 Checking file content...")

    # Try different possible paths
    possible_paths = [
        "reinvent_scoring/scoring/score_components/roshambo/roshambo_shape_similarity.py",
        "../reinvent-scoring/reinvent_scoring/scoring/score_components/roshambo/roshambo_shape_similarity.py",
        "../../reinvent-scoring/reinvent_scoring/scoring/score_components/roshambo/roshambo_shape_similarity.py"
    ]

    file_path = None
    for path in possible_paths:
        if os.path.exists(path):
            file_path = path
            print(f"   Found file at: {path}")
            break

    if not file_path:
        print("❌ File not found in any of the expected locations:")
        for path in possible_paths:
            print(f"   Tried: {path}")
        return False

    try:
        with open(file_path, 'r') as f:
            content = f.read()

        # Look for the method definition
        if "def _extract_scores_from_csv(self, csv_file: str, expected_count: int)" in content:
            print("✅ File contains the correct method signature")
            return True
        elif "def _extract_scores_from_csv(self, csv_file: str)" in content:
            print("❌ File contains the OLD method signature")
            return False
        else:
            print("❌ Method definition not found in file")
            return False

    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return False


def main():
    print("🔧 Force Reload and Test Roshambo Fix")
    print("=" * 50)

    # First check the file content
    file_ok = check_file_content()

    if file_ok:
        print("📁 File content is correct, testing module loading...")
        success = force_reload_and_test()

        if success:
            print("\n🎉 SUCCESS: The fix is working correctly!")
        else:
            print("\n❌ FAILURE: There's still an issue with the fix.")
    else:
        print("\n❌ FAILURE: The file doesn't contain the correct fix.")
        print("   You may need to manually edit the file or reapply the changes.")


if __name__ == "__main__":
    main()
