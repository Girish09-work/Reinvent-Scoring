# RDKit Alternatives to OpenEye's ROCS and OMEGA in LinkInvent

This document outlines RDKit-based alternatives to OpenEye's proprietary ROCS and OMEGA tools for use in the LinkInvent module of Reinvent 3.2. The goal is to provide open-source replacements that can achieve similar functionality without requiring commercial licenses.

## 1. Overview of OpenEye Tools Used in LinkInvent

Based on the codebase analysis, LinkInvent primarily uses two OpenEye tools:

1. **OMEGA** - For conformer generation
2. **ROCS** - For shape-based similarity calculations

These tools are used in the scoring components, particularly in the `RocsSimilarity` and `ParallelRocsSimilarity` classes, to evaluate the 3D similarity between generated molecules and reference structures.

## 2. RDKit Alternative for OMEGA (Conformer Generation)

### 2.1 ETKDG (Experimental Torsion Knowledge Distance Geometry)

RDKit's ETKDG algorithm is a robust alternative to OpenEye's OMEGA for conformer generation.

#### Implementation Details:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_conformers(mol, num_confs=200, energy_window=10, max_attempts=10):
    """
    Generate conformers using RDKit's ETKDG algorithm.
    
    Args:
        mol: RDKit molecule
        num_confs: Maximum number of conformers to generate
        energy_window: Energy window in kcal/mol for pruning conformers
        max_attempts: Maximum number of attempts for conformer generation
        
    Returns:
        RDKit molecule with conformers
    """
    # Prepare molecule
    mol = Chem.AddHs(mol)  # Add hydrogens
    
    # Set parameters for ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # For reproducibility
    params.numThreads = 0  # Use all available CPUs
    params.useSmallRingTorsions = True  # Improved small ring conformations
    params.useMacrocycleTorsions = True  # Improved macrocycle conformations
    params.useBasicKnowledge = True  # Use basic knowledge
    params.enforceChirality = True  # Enforce chirality
    params.maxIterations = 1000  # Maximum iterations
    
    # Generate conformers
    confs = AllChem.EmbedMultipleConfs(
        mol, 
        numConfs=num_confs,
        params=params
    )
    
    if len(confs) == 0:
        # If no conformers were generated, try again with different parameters
        for attempt in range(max_attempts):
            params.randomSeed += 1
            confs = AllChem.EmbedMultipleConfs(
                mol, 
                numConfs=num_confs,
                params=params
            )
            if len(confs) > 0:
                break
    
    # Energy minimize conformers
    for conf_id in confs:
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
    
    # Prune conformers based on energy window
    if len(confs) > 1:
        # Calculate energies
        energies = []
        for conf_id in confs:
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, 
                AllChem.MMFFGetMoleculeProperties(mol), 
                confId=conf_id
            )
            energies.append((conf_id, ff.CalcEnergy()))
        
        # Sort by energy
        energies.sort(key=lambda x: x[1])
        
        # Keep only conformers within energy window
        min_energy = energies[0][1]
        keep_conf_ids = [conf_id for conf_id, energy in energies if energy - min_energy <= energy_window]
        
        # Remove high energy conformers
        remove_conf_ids = [conf_id for conf_id in confs if conf_id not in keep_conf_ids]
        for conf_id in sorted(remove_conf_ids, reverse=True):
            mol.RemoveConformer(conf_id)
    
    return mol
```

#### Advantages over OMEGA:

1. **Open-source and free** - No license required
2. **Integrated with RDKit** - Seamless integration with other RDKit functionality
3. **Comparable quality** - Recent improvements (ETKDGv3) provide high-quality conformers
4. **Customizable** - Extensive parameters for fine-tuning

#### Limitations compared to OMEGA:

1. **Speed** - Generally slower than OMEGA for large-scale conformer generation
2. **Macrocycles** - Less specialized handling of macrocycles (though improving with recent versions)

### 2.2 Handling Stereochemistry

To replicate OMEGA's stereochemistry enumeration functionality:

```python
from rdkit.Chem import EnumerateStereoisomers

def enumerate_stereoisomers(mol, max_isomers=32):
    """
    Enumerate stereoisomers of a molecule.
    
    Args:
        mol: RDKit molecule
        max_isomers: Maximum number of stereoisomers to generate
        
    Returns:
        List of RDKit molecules representing stereoisomers
    """
    # Get all possible stereoisomers
    isomers = tuple(EnumerateStereoisomers.EnumerateStereoisomers(mol))
    
    # Limit the number of isomers
    if len(isomers) > max_isomers:
        isomers = isomers[:max_isomers]
    
    return isomers
```

## 3. RDKit Alternative for ROCS (Shape-Based Similarity)

### 3.1 USR and USRCAT

Ultrafast Shape Recognition (USR) and its extension USRCAT (USR with Chemical Annotation Types) are efficient alternatives to ROCS for shape-based similarity calculations.

#### Implementation Details:

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT

def calculate_shape_similarity(query_mol, ref_mol, method="usrcat", query_conf_id=-1, ref_conf_id=-1):
    """
    Calculate shape similarity between two molecules using USR or USRCAT.
    
    Args:
        query_mol: Query RDKit molecule with conformers
        ref_mol: Reference RDKit molecule with conformers
        method: Similarity method ('usr' or 'usrcat')
        query_conf_id: Conformer ID for query molecule (-1 for all)
        ref_conf_id: Conformer ID for reference molecule (-1 for all)
        
    Returns:
        Similarity score (0-1, where 1 is identical)
    """
    if method.lower() == "usr":
        # Calculate USR descriptors
        query_descriptor = AllChem.GetUSR(query_mol, confId=query_conf_id)
        ref_descriptor = AllChem.GetUSR(ref_mol, confId=ref_conf_id)
        
        # Calculate similarity
        similarity = GetUSRScore(query_descriptor, ref_descriptor)
    
    elif method.lower() == "usrcat":
        # Calculate USRCAT descriptors
        query_descriptor = AllChem.GetUSRCAT(query_mol, confId=query_conf_id)
        ref_descriptor = AllChem.GetUSRCAT(ref_mol, confId=ref_conf_id)
        
        # Calculate similarity
        similarity = GetUSRScore(query_descriptor, ref_descriptor)
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return similarity
```

#### Multi-Conformer Shape Similarity:

```python
def calculate_best_shape_similarity(query_mol, ref_mol, method="usrcat"):
    """
    Calculate the best shape similarity across all conformer pairs.
    
    Args:
        query_mol: Query RDKit molecule with multiple conformers
        ref_mol: Reference RDKit molecule with multiple conformers
        method: Similarity method ('usr' or 'usrcat')
        
    Returns:
        Best similarity score and conformer IDs
    """
    best_similarity = 0.0
    best_query_conf_id = -1
    best_ref_conf_id = -1
    
    # Iterate through all conformer pairs
    for query_conf_id in range(query_mol.GetNumConformers()):
        for ref_conf_id in range(ref_mol.GetNumConformers()):
            similarity = calculate_shape_similarity(
                query_mol, 
                ref_mol, 
                method=method,
                query_conf_id=query_conf_id,
                ref_conf_id=ref_conf_id
            )
            
            if similarity > best_similarity:
                best_similarity = similarity
                best_query_conf_id = query_conf_id
                best_ref_conf_id = ref_conf_id
    
    return best_similarity, best_query_conf_id, best_ref_conf_id
```

### 3.2 RDKit Shape Alignment

For more accurate shape alignment (though slower than USR/USRCAT):

```python
def align_molecules_by_shape(query_mol, ref_mol, query_conf_id=0, ref_conf_id=0):
    """
    Align a query molecule to a reference molecule based on shape.
    
    Args:
        query_mol: Query RDKit molecule with conformers
        ref_mol: Reference RDKit molecule with conformers
        query_conf_id: Conformer ID for query molecule
        ref_conf_id: Conformer ID for reference molecule
        
    Returns:
        Aligned query molecule and RMSD
    """
    # Create a copy of the query molecule
    aligned_mol = Chem.Mol(query_mol)
    
    # Perform the alignment
    rmsd = AllChem.AlignMol(
        aligned_mol,
        ref_mol,
        prbCid=query_conf_id,
        refCid=ref_conf_id
    )
    
    return aligned_mol, rmsd
```

### 3.3 Shape and Feature Similarity (ROCS-like)

To more closely mimic ROCS functionality, combining shape and pharmacophoric features:

```python
def calculate_shape_and_feature_similarity(query_mol, ref_mol, shape_weight=0.5, feature_weight=0.5):
    """
    Calculate combined shape and pharmacophoric feature similarity.
    
    Args:
        query_mol: Query RDKit molecule with conformers
        ref_mol: Reference RDKit molecule with conformers
        shape_weight: Weight for shape similarity (0-1)
        feature_weight: Weight for feature similarity (0-1)
        
    Returns:
        Combined similarity score (0-1)
    """
    # Calculate shape similarity using USRCAT
    shape_similarity, query_conf_id, ref_conf_id = calculate_best_shape_similarity(
        query_mol, 
        ref_mol, 
        method="usrcat"
    )
    
    # Align molecules based on best conformers
    aligned_mol, _ = align_molecules_by_shape(
        query_mol, 
        ref_mol, 
        query_conf_id=query_conf_id,
        ref_conf_id=ref_conf_id
    )
    
    # Define pharmacophoric features
    feature_factory = AllChem.BuildFeatureFactory()
    query_features = feature_factory.GetFeaturesForMol(aligned_mol, confId=query_conf_id)
    ref_features = feature_factory.GetFeaturesForMol(ref_mol, confId=ref_conf_id)
    
    # Calculate feature similarity based on spatial overlap
    feature_similarity = calculate_feature_similarity(query_features, ref_features)
    
    # Combine similarities
    combined_similarity = (
        (shape_weight * shape_similarity) + 
        (feature_weight * feature_similarity)
    ) / (shape_weight + feature_weight)
    
    return combined_similarity
```

## 4. Implementation Strategy for LinkInvent

### 4.1 Creating RDKit-Based Replacements

To replace OpenEye tools in LinkInvent, we need to create RDKit-based alternatives for:

1. **RDKitConformerGenerator** - To replace OMEGA
2. **RDKitShapeSimilarity** - To replace ROCS

### 4.2 RDKitConformerGenerator Implementation

```python
class RDKitConformerGenerator:
    """RDKit-based conformer generator to replace OpenEye's OMEGA."""
    
    def __init__(self, max_confs=200, energy_window=10, use_stereo_enumeration=False, max_stereo=4):
        """
        Initialize the conformer generator.
        
        Args:
            max_confs: Maximum number of conformers to generate
            energy_window: Energy window in kcal/mol
            use_stereo_enumeration: Whether to enumerate stereoisomers
            max_stereo: Maximum number of stereoisomers to generate
        """
        self.max_confs = max_confs
        self.energy_window = energy_window
        self.use_stereo_enumeration = use_stereo_enumeration
        self.max_stereo = max_stereo
    
    def generate_conformers(self, smiles):
        """
        Generate conformers for a molecule from SMILES.
        
        Args:
            smiles: SMILES string
            
        Returns:
            RDKit molecule with conformers or None if failed
        """
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Handle stereochemistry
        if self.use_stereo_enumeration:
            isomers = enumerate_stereoisomers(mol, max_isomers=self.max_stereo)
            
            # Generate conformers for each isomer
            all_confs_mol = None
            for isomer in isomers:
                confs_mol = generate_conformers(
                    isomer,
                    num_confs=self.max_confs // len(isomers),
                    energy_window=self.energy_window
                )
                
                if confs_mol.GetNumConformers() > 0:
                    if all_confs_mol is None:
                        all_confs_mol = confs_mol
                    else:
                        # Add conformers from this isomer to the main molecule
                        for conf in confs_mol.GetConformers():
                            all_confs_mol.AddConformer(conf, assignId=True)
            
            return all_confs_mol
        else:
            # Generate conformers without stereochemistry enumeration
            return generate_conformers(
                mol,
                num_confs=self.max_confs,
                energy_window=self.energy_window
            )
```

### 4.3 RDKitShapeSimilarity Implementation

```python
class RDKitShapeSimilarity:
    """RDKit-based shape similarity calculator to replace OpenEye's ROCS."""
    
    def __init__(self, reference_file, shape_weight=0.5, feature_weight=0.5, method="usrcat"):
        """
        Initialize the shape similarity calculator.
        
        Args:
            reference_file: Path to reference molecule file (SDF)
            shape_weight: Weight for shape similarity
            feature_weight: Weight for feature similarity
            method: Similarity method ('usr', 'usrcat', or 'alignment')
        """
        self.shape_weight = shape_weight
        self.feature_weight = feature_weight
        self.method = method
        
        # Load reference molecule
        self.reference_mol = self._load_reference(reference_file)
        
        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=50,  # Fewer conformers for faster screening
            energy_window=10
        )
    
    def _load_reference(self, reference_file):
        """Load reference molecule from file."""
        # Check file extension
        if reference_file.endswith('.sdf'):
            # Load SDF file
            suppl = Chem.SDMolSupplier(reference_file)
            if len(suppl) > 0:
                ref_mol = suppl[0]
                # Generate conformers if none exist
                if ref_mol.GetNumConformers() == 0:
                    ref_mol = generate_conformers(ref_mol, num_confs=10)
                return ref_mol
        else:
            # Assume it's a SMILES file
            with open(reference_file, 'r') as f:
                smiles = f.readline().strip()
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return generate_conformers(mol, num_confs=10)
        
        raise ValueError(f"Could not load reference from {reference_file}")
    
    def calculate_similarity(self, smiles):
        """
        Calculate shape similarity between a query molecule and the reference.
        
        Args:
            smiles: SMILES string of query molecule
            
        Returns:
            Similarity score (0-1)
        """
        # Generate conformers for query molecule
        query_mol = self.conformer_generator.generate_conformers(smiles)
        if query_mol is None or query_mol.GetNumConformers() == 0:
            return 0.0
        
        # Calculate similarity based on method
        if self.method == "usrcat":
            similarity, _, _ = calculate_best_shape_similarity(
                query_mol, 
                self.reference_mol, 
                method="usrcat"
            )
            return similarity
        
        elif self.method == "usr":
            similarity, _, _ = calculate_best_shape_similarity(
                query_mol, 
                self.reference_mol, 
                method="usr"
            )
            return similarity
        
        elif self.method == "combined":
            return calculate_shape_and_feature_similarity(
                query_mol,
                self.reference_mol,
                shape_weight=self.shape_weight,
                feature_weight=self.feature_weight
            )
        
        else:
            raise ValueError(f"Unknown method: {self.method}")
```

## 5. Comparison with OpenEye Tools

### 5.1 Advantages of RDKit Alternatives

1. **Open-source and free** - No license costs or restrictions
2. **Integration** - Seamless integration with existing RDKit functionality
3. **Customizability** - Full control over algorithms and parameters
4. **Active development** - Continuous improvements in the RDKit community

### 5.2 Limitations of RDKit Alternatives

1. **Performance** - Generally slower than OpenEye tools, especially for large-scale screening
2. **Conformer quality** - May require more fine-tuning to achieve comparable quality
3. **Feature richness** - Fewer built-in features compared to specialized commercial tools
4. **Documentation** - Less comprehensive documentation for advanced use cases

### 5.3 Benchmark Comparison

| Feature | OpenEye OMEGA/ROCS | RDKit Alternative | Notes |
|---------|-------------------|-------------------|-------|
| Conformer generation speed | ★★★★★ | ★★★☆☆ | OMEGA is typically 2-3x faster |
| Conformer quality | ★★★★★ | ★★★★☆ | ETKDGv3 approaches OMEGA quality |
| Shape similarity accuracy | ★★★★★ | ★★★★☆ | USRCAT is very competitive |
| Shape similarity speed | ★★★★☆ | ★★★★★ | USR/USRCAT can be faster than ROCS |
| Feature handling | ★★★★★ | ★★★☆☆ | ROCS has more sophisticated pharmacophore handling |
| Ease of use | ★★★★☆ | ★★★☆☆ | OpenEye tools have better defaults |
| Customizability | ★★★☆☆ | ★★★★★ | RDKit offers more low-level control |
| License cost | ★☆☆☆☆ | ★★★★★ | RDKit is free and open-source |

## 6. Conclusion

RDKit provides viable alternatives to OpenEye's ROCS and OMEGA tools for use in LinkInvent. While there are some trade-offs in terms of performance and ease of use, the RDKit-based implementations offer sufficient functionality for most applications without requiring commercial licenses.

For LinkInvent specifically, the RDKit alternatives should be able to replace the OpenEye tools with minimal impact on overall performance, especially if the scoring components are optimized for the specific use case.

## 7. References

1. RDKit: Open-source cheminformatics - https://www.rdkit.org/
2. ETKDG: Riniker, S.; Landrum, G. A. Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation. J. Chem. Inf. Model. 2015, 55, 2562–2574.
3. USR: Ballester, P. J.; Richards, W. G. Ultrafast Shape Recognition to Search Compound Databases for Similar Molecular Shapes. J. Comput. Chem. 2007, 28, 1711–1723.
4. USRCAT: Schreyer, A. M.; Blundell, T. USRCAT: Real-time Ultrafast Shape Recognition with Pharmacophoric Constraints. J. Cheminform. 2012, 4, 27.
5. OpenEye ROCS: https://www.eyesopen.com/rocs
6. OpenEye OMEGA: https://www.eyesopen.com/omega
