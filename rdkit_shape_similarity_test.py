import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Set environment variables if needed
#os.environ['RDBASE'] = r"D:\conda_envs\roshambo\Lib\site-packages\rdkit"

class RDKitShapeSimilarity:
    def __init__(self, reference_file, method="usrcat", shape_weight=0.5, color_weight=0.5, max_confs=50):
        self.reference_file = reference_file
        self.method = method
        self.shape_weight = shape_weight
        self.color_weight = color_weight
        self.max_confs = max_confs
        
        # Load reference molecule
        self.reference_mol = self._load_reference_molecule()
        print(f"Loaded reference molecule from {reference_file}")
        
    def _load_reference_molecule(self):
        """Load reference molecule from SDF file."""
        try:
            # Load the reference molecule from SDF
            reference_mol = Chem.SDMolSupplier(self.reference_file)[0]
            if reference_mol is None:
                print(f"Error: Could not load reference molecule from {self.reference_file}")
                return None
                
            print(f"Reference molecule has {reference_mol.GetNumConformers()} conformers")
            return reference_mol
        except Exception as e:
            print(f"Error loading reference molecule: {e}")
            return None
    
    def generate_conformers(self, mol, max_confs=None):
        """Generate conformers for a molecule."""
        if max_confs is None:
            max_confs = self.max_confs
            
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate conformers
        confs = AllChem.EmbedMultipleConfs(
            mol, 
            numConfs=max_confs,
            pruneRmsThresh=0.5,
            randomSeed=42,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True
        )
        
        # Minimize conformers
        for conf_id in confs:
            AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
            
        return mol
    
    def calculate_similarity(self, smiles):
        """Calculate shape similarity between a SMILES string and the reference molecule."""
        start_time = time.time()
        
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error: Could not parse SMILES: {smiles}")
            return 0.0
            
        # Generate conformers
        print(f"Generating conformers for {smiles}")
        mol_with_confs = self.generate_conformers(mol)
        print(f"Generated {mol_with_confs.GetNumConformers()} conformers")
        
        # Calculate similarity based on method
        if self.method == "usrcat":
            score = self._calculate_usrcat_similarity(mol_with_confs)
        elif self.method == "o3a":
            score = self._calculate_o3a_similarity(mol_with_confs)
        else:
            print(f"Unknown method: {self.method}")
            return 0.0
            
        end_time = time.time()
        print(f"Calculation took {end_time - start_time:.2f} seconds")
        return score
    
    def _calculate_usrcat_similarity(self, mol):
        """Calculate USRCAT similarity."""
        from rdkit.Chem import rdMolDescriptors
        
        best_score = 0.0
        
        for q_conf_id in range(mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                # Get USRCAT descriptors
                query_descriptor = rdMolDescriptors.GetUSRCAT(mol, confId=q_conf_id)
                ref_descriptor = rdMolDescriptors.GetUSRCAT(self.reference_mol, confId=r_conf_id)
                
                # Calculate similarity
                similarity = rdMolDescriptors.GetUSRScore(query_descriptor, ref_descriptor)
                
                if similarity > best_score:
                    best_score = similarity
        
        return best_score
    
    def _calculate_o3a_similarity(self, mol):
        """Calculate O3A similarity."""
        best_score = 0.0
        
        for q_conf_id in range(mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                try:
                    # Create O3A alignment
                    pyO3A = AllChem.GetO3A(mol, self.reference_mol,
                                          confId1=q_conf_id, confId2=r_conf_id)
                    
                    # Get alignment score (shape similarity)
                    shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range
                    
                    if shape_sim > best_score:
                        best_score = shape_sim
                except Exception as e:
                    print(f"Error in O3A calculation: {e}")
                    continue
        
        return best_score

# Example usage
if __name__ == "__main__":
    # Path to your reference SDF file
    reference_file = r"D:\reinvent-scoring\sample_data\query.sdf"
    
    # Create shape similarity calculator
    shape_calculator = RDKitShapeSimilarity(
        reference_file=reference_file,
        method="usrcat",  # or "o3a"
        shape_weight=0.5,
        color_weight=0.5,
        max_confs=10  # Using fewer conformers for faster testing
    )
    
    # Test molecules (SMILES strings)
    test_molecules = [
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",  # Imatinib
        "CC(C)C1=C(C(=CC=C1)C(=O)NC2=CC(=C(C=C2)C(=O)NC3=CC=CC=C3)OC)NC(=O)C4=CC=CC=C4",  # Example molecule
        "C1=CC=C(C=C1)C(=O)NCCOCCOCCN"  # Simple molecule
    ]
    
    # Calculate and print similarity scores
    print("\nCalculating shape similarity scores:")
    print("-" * 50)
    
    for i, smiles in enumerate(test_molecules):
        print(f"\nMolecule {i+1}:")
        print(f"SMILES: {smiles}")
        
        start = time.time()
        score = shape_calculator.calculate_similarity(smiles)
        end = time.time()
        
        print(f"Shape similarity score: {score:.4f}")
        print(f"Total time: {end - start:.2f} seconds")
        print("-" * 50)