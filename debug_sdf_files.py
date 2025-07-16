#!/usr/bin/env python3
"""
Debug SDF files to understand why roshambo can't read them.
"""

import os
from rdkit import Chem

def test_sdf_file(filepath, name):
    """Test if an SDF file can be read by RDKit and roshambo."""
    print(f"\n🔍 Testing {name}: {filepath}")
    
    if not os.path.exists(filepath):
        print(f"❌ File does not exist")
        return False
    
    print(f"📁 File size: {os.path.getsize(filepath)} bytes")
    
    # Test with RDKit
    try:
        suppl = Chem.SDMolSupplier(filepath)
        mol_count = 0
        valid_mols = 0
        
        for i, mol in enumerate(suppl):
            mol_count += 1
            if mol is not None:
                valid_mols += 1
                print(f"✅ Molecule {i}: {mol.GetNumAtoms()} atoms, SMILES: {Chem.MolToSmiles(mol)}")
                if i >= 2:  # Show first 3 molecules
                    break
            else:
                print(f"❌ Molecule {i}: Invalid")
                if i >= 2:
                    break
        
        print(f"📊 Total molecules: {mol_count}, Valid: {valid_mols}")
        
        if valid_mols == 0:
            print("❌ No valid molecules found by RDKit")
            return False
        
    except Exception as e:
        print(f"❌ RDKit error: {e}")
        return False
    
    # Test with roshambo utilities
    try:
        from roshambo.utilities import prepare_mols
        
        print("🧪 Testing with roshambo prepare_mols...")
        ref_mol, _, _ = prepare_mols(filepath, ignore_hs=True, n_confs=0)
        print(f"✅ Roshambo can read the file: {ref_mol}")
        return True
        
    except Exception as e:
        print(f"❌ Roshambo error: {e}")
        
        # Try to understand what roshambo is doing
        print("🔍 Debugging roshambo behavior...")
        try:
            from roshambo.utilities import smiles_to_rdmol
            
            # Check if roshambo is trying to read the file as SMILES
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                print(f"First line of file: '{first_line}'")
                
            # Test if roshambo thinks this is a SMILES file
            try:
                result = smiles_to_rdmol([first_line])
                print(f"Roshambo interpreted first line as SMILES: {result}")
            except Exception as e2:
                print(f"Roshambo failed to parse first line as SMILES: {e2}")
                
        except Exception as e2:
            print(f"Error debugging roshambo: {e2}")
        
        return False

def create_simple_sdf():
    """Create a very simple SDF file that should work."""
    print("\n🛠️ Creating simple test SDF...")
    
    # Create a simple molecule (ethanol)
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Write to SDF
    test_sdf = "/tmp/simple_test.sdf"
    writer = Chem.SDWriter(test_sdf)
    mol.SetProp("_Name", "ethanol")
    writer.write(mol)
    writer.close()
    
    print(f"✅ Created simple SDF: {test_sdf}")
    return test_sdf

def test_roshambo_with_simple_sdf():
    """Test roshambo with a simple SDF we create."""
    
    # Create simple SDF files
    ref_sdf = create_simple_sdf()
    
    # Create dataset SDF (same molecule)
    dataset_sdf = "/tmp/simple_dataset.sdf"
    import shutil
    shutil.copy2(ref_sdf, dataset_sdf)
    
    # Test with roshambo
    print("\n🧪 Testing roshambo with simple SDF files...")
    
    try:
        from roshambo.api import get_similarity_scores
        
        # Change to /tmp directory
        original_cwd = os.getcwd()
        os.chdir("/tmp")
        
        results = get_similarity_scores(
            ref_file="simple_test.sdf",
            dataset_files_pattern="simple_dataset.sdf",
            ignore_hs=True,
            n_confs=0,
            use_carbon_radii=True,
            color=True,
            sort_by="ComboTanimoto",
            write_to_file=False,  # Don't write files for this test
            gpu_id=0,
            working_dir=None
        )
        
        print("✅ Roshambo worked with simple SDF!")
        print(f"Results: {results}")
        return True
        
    except Exception as e:
        print(f"❌ Roshambo failed with simple SDF: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        os.chdir(original_cwd)

def main():
    """Run all tests."""
    print("🔍 SDF File Debugging")
    print("=" * 50)
    
    # Test original files
    files_to_test = [
        ("/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BRD9/protac.sdf", "BRD9 PROTAC"),
        ("/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BRD9/warhead.sdf", "BRD9 Warhead"),
        ("/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BTK/BTK_sel.sdf", "BTK Selection"),
    ]
    
    results = []
    for filepath, name in files_to_test:
        success = test_sdf_file(filepath, name)
        results.append((name, success))
    
    # Test with simple SDF
    print("\n" + "=" * 50)
    simple_success = test_roshambo_with_simple_sdf()
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 SUMMARY:")
    for name, success in results:
        print(f"{name}: {'✅ SUCCESS' if success else '❌ FAILED'}")
    print(f"Simple SDF test: {'✅ SUCCESS' if simple_success else '❌ FAILED'}")
    
    if simple_success and not any(success for _, success in results):
        print("\n💡 DIAGNOSIS: Issue with original SDF files")
        print("   - Roshambo works with simple SDF files")
        print("   - Original SDF files have format issues")
        print("   - Need to convert or fix original SDF files")
    elif not simple_success:
        print("\n💡 DIAGNOSIS: Roshambo installation issue")
        print("   - Roshambo cannot read any SDF files")
        print("   - Check roshambo installation")
        print("   - Check RDKit compatibility")

if __name__ == "__main__":
    main()
