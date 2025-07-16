#!/usr/bin/env python3
"""
Demonstration of how to properly extract scores from roshambo CSV output.
This shows the correct way to get scores after discovering the API bug.
"""

import os
import pandas as pd

def demo_csv_score_extraction():
    """Demonstrate how to extract scores from roshambo CSV output."""
    
    print("🔍 ROSHAMBO CSV SCORE EXTRACTION DEMO")
    print("=" * 50)
    
    # This is the CSV content from your test output
    csv_content = """Molecule	OriginalName	ComboTanimoto	ShapeTanimoto	ColorTanimoto	FitTverskyCombo	FitTversky	FitColorTversky	RefTverskyCombo	RefTversky	RefColorTversky	Overlap
mol_0_0	mol_0	2.0	1.0	1.0	2.0	1.0	1.0	2.0	1.0	1.0	237.352
mol_1_0	mol_1	0.976	0.976	0.0	0.905	0.905	0.0	1.088	1.088	0.0	261.155"""
    
    # Save to temporary file
    with open("demo_roshambo.csv", "w") as f:
        f.write(csv_content)
    
    print("📄 Sample CSV content:")
    print(csv_content)
    print()
    
    # Read the CSV file
    print("📊 Reading CSV with pandas:")
    df = pd.read_csv("demo_roshambo.csv", sep='\t')
    
    print(f"   Shape: {df.shape}")
    print(f"   Columns: {df.columns.tolist()}")
    print()
    
    print("🎯 Extracting individual molecule scores:")
    
    # Extract scores for each molecule
    molecule_scores = {}
    
    for _, row in df.iterrows():
        molecule_name = row['Molecule']
        original_name = row['OriginalName']
        
        # Extract the molecule index from the name
        # "mol_0_0" -> index 0, "mol_1_0" -> index 1
        parts = molecule_name.split('_')
        if len(parts) >= 2:
            mol_index = int(parts[1])
        else:
            continue
            
        shape_score = float(row['ShapeTanimoto'])
        color_score = float(row['ColorTanimoto'])
        combo_score = float(row['ComboTanimoto'])
        
        molecule_scores[mol_index] = {
            'shape': shape_score,
            'color': color_score,
            'combo': combo_score,
            'original_name': original_name
        }
        
        print(f"   Molecule {mol_index} ({original_name}):")
        print(f"      Shape Tanimoto: {shape_score:.3f}")
        print(f"      Color Tanimoto: {color_score:.3f}")
        print(f"      Combo Tanimoto: {combo_score:.3f}")
        print()
    
    print("📈 Final scores array (for reinvent-scoring):")
    
    # Create scores array like reinvent-scoring expects
    max_index = max(molecule_scores.keys()) if molecule_scores else -1
    scores_array = []
    
    for i in range(max_index + 1):
        if i in molecule_scores:
            # Use combo score as final score
            score = molecule_scores[i]['combo']
        else:
            score = 0.0
        scores_array.append(score)
    
    print(f"   Scores: {scores_array}")
    print(f"   Length: {len(scores_array)}")
    
    # Cleanup
    os.remove("demo_roshambo.csv")
    
    return scores_array

def demo_weighted_scoring():
    """Demonstrate custom weighted scoring like in reinvent-scoring."""
    
    print("\n🔧 CUSTOM WEIGHTED SCORING DEMO")
    print("=" * 50)
    
    # Example: 70% shape weight, 30% color weight
    shape_weight = 0.7
    color_weight = 0.3
    
    print(f"Shape weight: {shape_weight}")
    print(f"Color weight: {color_weight}")
    print()
    
    # Sample scores from CSV
    molecules = [
        {'name': 'mol_0', 'shape': 1.0, 'color': 1.0},
        {'name': 'mol_1', 'shape': 0.976, 'color': 0.0}
    ]
    
    print("🧮 Calculating weighted scores:")
    
    for mol in molecules:
        shape_score = mol['shape']
        color_score = mol['color']
        
        # Calculate weighted combination
        weighted_score = (shape_weight * shape_score + color_weight * color_score) / (shape_weight + color_weight)
        
        print(f"   {mol['name']}:")
        print(f"      Shape: {shape_score:.3f} × {shape_weight} = {shape_score * shape_weight:.3f}")
        print(f"      Color: {color_score:.3f} × {color_weight} = {color_score * color_weight:.3f}")
        print(f"      Weighted: {weighted_score:.3f}")
        print()

if __name__ == "__main__":
    # Run the demonstrations
    scores = demo_csv_score_extraction()
    demo_weighted_scoring()
    
    print("✅ SUMMARY:")
    print("   • Roshambo API returns None but writes correct CSV files")
    print("   • CSV files are tab-separated with molecule names and scores")
    print("   • Extract molecule index from 'Molecule' column (e.g., 'mol_1_0' -> index 1)")
    print("   • Use ShapeTanimoto, ColorTanimoto, or ComboTanimoto as needed")
    print("   • Apply custom weighting if required")
    print("   • Convert to numpy array for reinvent-scoring compatibility")
