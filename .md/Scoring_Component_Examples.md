# Understanding Scoring Component Construction in Reinvent 3.2

This document explains how scoring components are constructed in Reinvent 3.2, using two examples:
1. QED Score - A simple component without specific parameters
2. Tanimoto Similarity - A component that uses component-specific parameters

## 1. QED Score Component

The QED (Quantitative Estimate of Drug-likeness) score is a simple scoring component that calculates a drug-likeness score for molecules. It doesn't require any component-specific parameters.

### 1.1 Component Structure

```python
from typing import List
import numpy as np
from rdkit.Chem.QED import qed

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary


class QedScore(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

    def calculate_score(self, molecules: List) -> ComponentSummary:
        score = self._calculate_qed(molecules)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters)
        return score_summary

    def _calculate_qed(self, query_mols) -> np.array:
        qed_scores = []
        for mol in query_mols:
            try:
                qed_score = qed(mol)
            except ValueError:
                qed_score = 0.0
            qed_scores.append(qed_score)
        return np.array(qed_scores, dtype=np.float32)
```

### 1.2 Key Components Explained

1. **Class Definition**:
   ```python
   class QedScore(BaseScoreComponent):
   ```
   - Inherits from `BaseScoreComponent`, which provides common functionality for all scoring components

2. **Initialization**:
   ```python
   def __init__(self, parameters: ComponentParameters):
       super().__init__(parameters)
   ```
   - Takes a `ComponentParameters` object containing component configuration
   - Calls the parent class initializer to set up common attributes
   - No additional initialization needed for this simple component

3. **Main Scoring Method**:
   ```python
   def calculate_score(self, molecules: List) -> ComponentSummary:
       score = self._calculate_qed(molecules)
       score_summary = ComponentSummary(total_score=score, parameters=self.parameters)
       return score_summary
   ```
   - Takes a list of RDKit molecules as input
   - Calls the internal method to calculate QED scores
   - Returns a `ComponentSummary` containing the scores and parameters

4. **Internal Calculation Method**:
   ```python
   def _calculate_qed(self, query_mols) -> np.array:
       qed_scores = []
       for mol in query_mols:
           try:
               qed_score = qed(mol)
           except ValueError:
               qed_score = 0.0
           qed_scores.append(qed_score)
       return np.array(qed_scores, dtype=np.float32)
   ```
   - Iterates through each molecule
   - Calculates the QED score using RDKit's `qed` function
   - Handles errors by assigning a score of 0.0
   - Returns scores as a NumPy array

### 1.3 Using the QED Score Component

```python
qed_params = {
    "component_type": "qed_score",
    "name": "Drug-likeness",
    "weight": 1.0,
    "specific_parameters": {
        "transformation": {
            "transformation_type": "no_transformation"
        }
    }
}
```

- No component-specific parameters are needed
- Transformation is optional (here set to "no_transformation")

## 2. Tanimoto Similarity Component

The Tanimoto Similarity component calculates the similarity between molecules and a set of reference molecules. It requires component-specific parameters to define the reference molecules and fingerprint settings.

### 2.1 Component Structure

```python
from typing import List
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary


class TanimotoSimilarity(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        
        # Get component-specific parameters
        specific_params = parameters.specific_parameters
        
        # Get reference SMILES from parameters
        smiles_key = self.component_specific_parameters.SMILES
        self.smiles = specific_params.get(smiles_key, [])
        
        # Get fingerprint parameters
        self._radius = specific_params.get("radius", 2)
        self._use_counts = specific_params.get("use_counts", True)
        self._use_features = specific_params.get("use_features", True)
        
        # Convert reference SMILES to fingerprints
        self._ref_mols = [Chem.MolFromSmiles(smile) for smile in self.smiles]
        self._ref_fingerprints = self._calculate_fingerprints(self._ref_mols)

    def calculate_score(self, molecules: List) -> ComponentSummary:
        query_fps = self._calculate_fingerprints(molecules)
        score = self._calculate_similarity(query_fps)
        score_summary = ComponentSummary(total_score=score, parameters=self.parameters)
        return score_summary
        
    def _calculate_fingerprints(self, molecules):
        fingerprints = []
        for mol in molecules:
            if mol is not None:
                if self._use_features:
                    fp = AllChem.GetMorganFeaturesFingerprint(
                        mol, self._radius, useChirality=True, useCounts=self._use_counts
                    )
                else:
                    fp = AllChem.GetMorganFingerprint(
                        mol, self._radius, useChirality=True, useCounts=self._use_counts
                    )
                fingerprints.append(fp)
            else:
                fingerprints.append(None)
        return fingerprints
        
    def _calculate_similarity(self, query_fingerprints):
        similarities = []
        for query_fp in query_fingerprints:
            if query_fp is None:
                similarities.append(0.0)
                continue
                
            # Calculate similarity to each reference fingerprint
            max_sim = 0.0
            for ref_fp in self._ref_fingerprints:
                if ref_fp is not None:
                    sim = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
                    max_sim = max(max_sim, sim)
            
            similarities.append(max_sim)
            
        return np.array(similarities, dtype=np.float32)
```

### 2.2 Key Components Explained

1. **Class Definition**:
   ```python
   class TanimotoSimilarity(BaseScoreComponent):
   ```
   - Inherits from `BaseScoreComponent` like all scoring components

2. **Initialization with Component-Specific Parameters**:
   ```python
   def __init__(self, parameters: ComponentParameters):
       super().__init__(parameters)
       
       # Get component-specific parameters
       specific_params = parameters.specific_parameters
       
       # Get reference SMILES from parameters
       smiles_key = self.component_specific_parameters.SMILES
       self.smiles = specific_params.get(smiles_key, [])
       
       # Get fingerprint parameters
       self._radius = specific_params.get("radius", 2)
       self._use_counts = specific_params.get("use_counts", True)
       self._use_features = specific_params.get("use_features", True)
       
       # Convert reference SMILES to fingerprints
       self._ref_mols = [Chem.MolFromSmiles(smile) for smile in self.smiles]
       self._ref_fingerprints = self._calculate_fingerprints(self._ref_mols)
   ```
   - Accesses component-specific parameters from `parameters.specific_parameters`
   - Uses `component_specific_parameters` enum to get standardized parameter names
   - Provides default values for parameters that might not be specified
   - Preprocesses reference molecules into fingerprints for efficiency

3. **Main Scoring Method**:
   ```python
   def calculate_score(self, molecules: List) -> ComponentSummary:
       query_fps = self._calculate_fingerprints(molecules)
       score = self._calculate_similarity(query_fps)
       score_summary = ComponentSummary(total_score=score, parameters=self.parameters)
       return score_summary
   ```
   - Converts input molecules to fingerprints
   - Calculates similarity scores
   - Returns a `ComponentSummary`

4. **Fingerprint Calculation Method**:
   ```python
   def _calculate_fingerprints(self, molecules):
       # Implementation details...
   ```
   - Creates Morgan fingerprints based on the specified parameters
   - Handles different fingerprint types (feature-based or not)
   - Handles count-based or binary fingerprints

5. **Similarity Calculation Method**:
   ```python
   def _calculate_similarity(self, query_fingerprints):
       # Implementation details...
   ```
   - Calculates Tanimoto similarity between query and reference fingerprints
   - Takes the maximum similarity across all reference molecules
   - Returns scores as a NumPy array

### 2.3 Using the Tanimoto Similarity Component

```python
tanimoto_params = {
    "component_type": "tanimoto_similarity",
    "name": "Similarity to References",
    "weight": 1.0,
    "specific_parameters": {
        "smiles": ["CC(=O)OC1=CC=CC=C1C(=O)O", "C1CCCCC1"],  # Aspirin and cyclohexane
        "radius": 3,
        "use_counts": True,
        "use_features": True,
        "transformation": {
            "transformation_type": "sigmoid",
            "high": 1.0,
            "low": 0.0,
            "k": 0.5
        }
    }
}
```

- Component-specific parameters:
  - `smiles`: List of reference SMILES strings
  - `radius`: Morgan fingerprint radius
  - `use_counts`: Whether to use count-based fingerprints
  - `use_features`: Whether to use feature-based fingerprints
- Transformation is set to sigmoid to shape the scoring landscape

## 3. Key Differences Between the Components

1. **Complexity**:
   - QED Score is simple with no component-specific parameters
   - Tanimoto Similarity requires multiple component-specific parameters

2. **Initialization**:
   - QED Score has minimal initialization
   - Tanimoto Similarity processes parameters and precomputes fingerprints

3. **Calculation Logic**:
   - QED Score directly calculates a property
   - Tanimoto Similarity compares molecules to references

4. **Parameter Handling**:
   - QED Score doesn't need special parameters
   - Tanimoto Similarity carefully extracts and uses multiple parameters

## 4. Best Practices for Component Construction

1. **Use the Base Class Properly**:
   - Always inherit from `BaseScoreComponent`
   - Call `super().__init__(parameters)` in your initializer

2. **Handle Component-Specific Parameters Carefully**:
   - Use the `component_specific_parameters` enum for standard parameter names
   - Provide sensible defaults for optional parameters
   - Document required parameters

3. **Optimize for Performance**:
   - Precompute what you can in the initializer
   - Use efficient algorithms for scoring

4. **Handle Errors Gracefully**:
   - Catch exceptions during calculation
   - Provide fallback scores (usually 0.0) for invalid molecules

5. **Return Proper Data Structures**:
   - Always return a `ComponentSummary`
   - Ensure scores are a NumPy array with dtype=np.float32

## 5. Conclusion

These two examples demonstrate the flexibility of the scoring component system in Reinvent 3.2:

- Simple components like QED Score can be implemented with minimal code
- Complex components like Tanimoto Similarity can use component-specific parameters for customization

When creating your own scoring components, choose the appropriate level of complexity based on your needs. For simple property calculations, follow the QED Score pattern. For more configurable components, follow the Tanimoto Similarity pattern with proper component-specific parameter handling.
