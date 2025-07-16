# 🧩 Creating a Custom Scoring Component in Reinvent 3.2

This comprehensive guide explains how to create custom scoring components for Reinvent 3.2, which can be used in reinforcement learning and other scoring contexts.

## 📋 Table of Contents
- [Overview](#overview)
- [Component Hierarchy](#component-hierarchy)
- [Basic Implementation Steps](#basic-implementation-steps)
- [Real-World Examples](#real-world-examples)
  - [QED Score Component](#qed-score-component)
  - [Tanimoto Similarity Component](#tanimoto-similarity-component)
- [Practical Example: Fsp3 Component](#practical-example-fsp3-component)
- [Advanced Topics](#advanced-topics)
  - [Transformations](#transformations)
  - [Component-Specific Parameters](#component-specific-parameters)
- [Best Practices](#best-practices)

## 🔍 Overview

Scoring components in Reinvent 3.2 are modular pieces that calculate specific properties or scores for molecules. These components are combined in scoring functions (like `CustomSum` or `CustomProduct`) to provide a final score that guides the reinforcement learning process.

## 🏗️ Component Hierarchy

Before creating a new component, it's important to understand the component hierarchy:

1. `BaseScoreComponent` - The abstract base class for all scoring components
2. Specialized base components for different types of scoring:
   - 📊 `BasePhysChemComponent` - For physicochemical properties
   - 🔗 `BaseLinkInventComponent` - For Link Invent specific properties
   - 🌐 `BaseRESTComponent` - For components that call REST APIs
   - 🧬 `BaseStructuralComponent` - For structural properties
   - 💻 `BaseConsoleInvokedComponent` - For components that invoke console commands

Choose the appropriate base class based on the type of scoring component you want to create.

## 📝 Basic Implementation Steps

### Step 1: Create the Component File Structure

For this example, let's create a custom physicochemical property component called `MyCustomProperty`. Here's the file structure you'll need to create:

```
📁 reinvent_scoring/
└── 📁 scoring/
    └── 📁 score_components/
        └── 📁 physchem/
            ├── 📄 __init__.py (update this)
            └── 📄 my_custom_property.py (create this)
```

### Step 2: Implement the Component Class

Create the `my_custom_property.py` file with the following content:

```python
from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.physchem.base_physchem_component import BasePhysChemComponent


class MyCustomProperty(BasePhysChemComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        # Initialize any specific attributes or resources needed for your component

    def _calculate_phys_chem_property(self, mol):
        """
        Calculate your custom property for a single molecule.

        Args:
            mol: An RDKit molecule object

        Returns:
            float: The calculated property value
        """
        # Implement your custom property calculation here
        # For example:
        try:
            # Your calculation logic here
            # This is just a placeholder example
            custom_value = mol.GetNumAtoms() / 10.0  # Normalize by dividing by 10
            return custom_value
        except:
            return 0.0
```

The key method to implement is `_calculate_phys_chem_property`, which should calculate your custom property for a single molecule and return a numerical value.

### Step 3: Update the `__init__.py` File

Update the `__init__.py` file in the `physchem` directory to include your new component:

```python
# Existing imports...
from reinvent_scoring.scoring.score_components.physchem.mol_weight import MolWeight
from reinvent_scoring.scoring.score_components.physchem.tpsa import PSA
from reinvent_scoring.scoring.score_components.physchem.rot_bonds import RotatableBonds
# Add your new component
from reinvent_scoring.scoring.score_components.physchem.my_custom_property import MyCustomProperty
```

### Step 4: Register the Component in the Enum

Add your component to the `ScoringFunctionComponentNameEnum` in `scoring/enums/scoring_function_component_enum.py`:

```python
@dataclass(frozen=True)
class ScoringFunctionComponentNameEnum:
    # Existing components...
    MOLECULAR_WEIGHT = "molecular_weight"
    NUM_ROTATABLE_BONDS = "num_rotatable_bonds"
    # Add your new component
    MY_CUSTOM_PROPERTY = "my_custom_property"
    # Other components...
```

### Step 5: Register the Component in the Factory

Add your component to the `ScoreComponentFactory` in `scoring/score_components/score_component_factory.py`:

```python
def _deafult_scoring_component_registry(self) -> dict:
    enum = ScoringFunctionComponentNameEnum()
    component_map = {
        # Existing components...
        enum.MOLECULAR_WEIGHT: MolWeight,
        enum.TPSA: PSA,
        # Add your new component
        enum.MY_CUSTOM_PROPERTY: MyCustomProperty,
        # Other components...
    }
    return component_map
```

### Step 6: Using Your Custom Component

Now you can use your custom component in scoring functions. Here's an example of how to configure it:

```python
my_custom_property_params = {
    "component_type": "my_custom_property",
    "name": "My Custom Property",
    "weight": 1.0,
    "specific_parameters": {
        # Optional transformation parameters
        "transformation": {
            "transformation_type": "sigmoid",
            "high": 10,
            "low": 0,
            "k": 0.5
        }
    }
}

# Add to scoring function parameters
scoring_function_params = ScoringFunctionParameters(
    name="custom_sum",
    parameters=[my_custom_property_params, other_component_params],
    parallel=False
)
```

## 🧪 Real-World Examples

Let's look at two real-world examples of scoring components with different complexity levels:

### 📊 QED Score Component

The QED (Quantitative Estimate of Drug-likeness) score is a simple scoring component that calculates a drug-likeness score for molecules. It doesn't require any component-specific parameters.

#### Component Structure

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

#### Using the QED Score Component

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

### 🔍 Tanimoto Similarity Component

The Tanimoto Similarity component calculates the similarity between molecules and a set of reference molecules. It requires component-specific parameters to define the reference molecules and fingerprint settings.

#### Component Structure

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

#### Using the Tanimoto Similarity Component

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

## 🔬 Practical Example: Fsp3 Component

Let's create a complete, practical example of a custom scoring component that calculates the fraction of sp3 hybridized carbon atoms (Fsp3), which is a common metric in medicinal chemistry.

### What is Fsp3?

Fsp3 is the fraction of carbon atoms that are sp3 hybridized (tetrahedral) relative to the total number of carbon atoms in a molecule. Higher Fsp3 values are often associated with increased solubility and success in clinical development.

### Step 1: Create the Component File

**📄 File: `scoring/score_components/physchem/carbon_sp3_fraction.py`**

```python
from typing import List
import numpy as np
from rdkit import Chem

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.physchem.base_physchem_component import BasePhysChemComponent


class CarbonSP3Fraction(BasePhysChemComponent):
    """
    Calculates the fraction of sp3 hybridized carbon atoms (Fsp3) in a molecule.
    Fsp3 = (number of sp3 hybridized carbon atoms) / (total number of carbon atoms)
    """

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

    def _calculate_phys_chem_property(self, mol):
        """
        Calculate the fraction of sp3 hybridized carbon atoms.

        Args:
            mol: An RDKit molecule object

        Returns:
            float: The Fsp3 value (between 0 and 1)
        """
        if mol is None:
            return 0.0

        # Count carbon atoms
        carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
        total_carbons = len(carbon_atoms)

        if total_carbons == 0:
            return 0.0

        # Count sp3 hybridized carbon atoms
        sp3_carbons = sum(1 for atom in carbon_atoms if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)

        # Calculate fraction
        fsp3 = sp3_carbons / total_carbons

        return fsp3
```

### Step 2: Update the `__init__.py` File

**📄 File: `scoring/score_components/physchem/__init__.py`**

```python
from reinvent_scoring.scoring.score_components.physchem.mol_weight import MolWeight
from reinvent_scoring.scoring.score_components.physchem.tpsa import PSA
# Other imports...
from reinvent_scoring.scoring.score_components.physchem.carbon_sp3_fraction import CarbonSP3Fraction
```

### Step 3: Update the Component Enum

**📄 File: `scoring/enums/scoring_function_component_enum.py`**

```python
@dataclass(frozen=True)
class ScoringFunctionComponentNameEnum:
    # Existing components...
    CARBON_SP3_FRACTION = "carbon_sp3_fraction"
    # Other components...
```

### Step 4: Update the Component Factory

**📄 File: `scoring/score_components/score_component_factory.py`**

```python
def _deafult_scoring_component_registry(self) -> dict:
    enum = ScoringFunctionComponentNameEnum()
    component_map = {
        # Existing components...
        enum.CARBON_SP3_FRACTION: CarbonSP3Fraction,
        # Other components...
    }
    return component_map
```

### Step 5: Using the Component

```python
fsp3_params = {
    "component_type": "carbon_sp3_fraction",
    "name": "Fsp3",
    "weight": 1.0,
    "specific_parameters": {
        "transformation": {
            "transformation_type": "sigmoid",
            "high": 1.0,
            "low": 0.0,
            "k": 0.5
        }
    }
}

# Create scoring function parameters
scoring_function_params = ScoringFunctionParameters(
    name="custom_product",
    parameters=[fsp3_params, qed_params],
    parallel=False
)
```

## 🔧 Advanced Topics

### 🔄 Transformations

Transformations convert raw property values into more useful scores, typically in the range of 0 to 1. The base component already handles this through the `_transformation_function` attribute.

Common transformations include:

- 📈 **Sigmoid**: Creates an S-shaped curve that smoothly transitions from low scores to high scores
  ```python
  "transformation": {
      "transformation_type": "sigmoid",
      "high": 1.0,  # Maximum score
      "low": 0.0,   # Minimum score
      "k": 0.5      # Controls steepness
  }
  ```

- 📊 **Double Sigmoid**: Creates a hill or valley shape with two transitions
  ```python
  "transformation": {
      "transformation_type": "double_sigmoid",
      "high": 1.0,
      "low": 0.0,
      "coef_div": 100.0,
      "coef_si": 150.0,  # Lower inflection point
      "coef_se": 150.0   # Upper inflection point
  }
  ```

- 📏 **Step**: Creates a sharp cutoff between low and high scores
  ```python
  "transformation": {
      "transformation_type": "step",
      "high": 1.0,
      "low": 0.0,
      "threshold": 0.5
  }
  ```

- ➖ **No Transformation**: Uses the raw values directly
  ```python
  "transformation": {
      "transformation_type": "no_transformation"
  }
  ```

### 🔧 Component-Specific Parameters

Component-specific parameters allow you to customize the behavior of your component beyond the standard parameters like name and weight.

1. **Define in the Enum**:
   ```python
   @dataclass(frozen=True)
   class ComponentSpecificParametersEnum:
       # Existing parameters...
       MY_CUSTOM_PARAMETER = "my_custom_parameter"
   ```

2. **Access in Your Component**:
   ```python
   def __init__(self, parameters: ComponentParameters):
       super().__init__(parameters)
       self.my_custom_param = self.parameters.specific_parameters.get(
           self.component_specific_parameters.MY_CUSTOM_PARAMETER, default_value
       )
   ```

3. **Use in Configuration**:
   ```python
   component_params = {
       "component_type": "my_component",
       "name": "My Component",
       "weight": 1.0,
       "specific_parameters": {
           "my_custom_parameter": 42
       }
   }
   ```

## ✅ Best Practices

1. **Use the Base Class Properly**:
   - Always inherit from the appropriate base class
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

## 🏁 Summary

Creating a custom scoring component involves:

1. Choosing the appropriate base class
2. Implementing the required calculation method
3. Updating the necessary registration files
4. Configuring the component for use in scoring functions

By following these steps, you can extend Reinvent 3.2 with custom scoring components tailored to your specific needs.

## Adding RDKit Shape Components

If you want to add RDKit-based shape similarity components, follow these steps:

### Step 1: Create the Conformer Generator

Create a file `rdkit_conformer_generator.py` in the `reinvent_scoring/scoring/score_components/rdkit_shape/` directory:

```python
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers

class RDKitConformerGenerator:
    """RDKit-based conformer generator to replace OpenEye's OMEGA."""
    
    def __init__(self, max_confs=200, energy_window=10, max_stereo=0):
        self.max_confs = max_confs
        self.energy_window = energy_window
        self.max_stereo = max_stereo
        
    def generate_conformers(self, smiles):
        """Generate conformers for a molecule from SMILES."""
        # Implementation details...
```

### Step 2: Create the Shape Similarity Component

Create a file `rdkit_shape_similarity.py` in the same directory:

```python
from typing import List
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_conformer_generator import RDKitConformerGenerator

class RDKitShapeSimilarity(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        # Implementation details...
```

### Step 3: Create the Parallel Version (Optional)

Create a file `parallel_rdkit_shape_similarity.py` for parallel processing:

```python
import os
import multiprocessing
from multiprocessing import Pool

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_shape_similarity import RDKitShapeSimilarity

class ParallelRDKitShapeSimilarity(RDKitShapeSimilarity):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        # Implementation details...
```

### Step 4: Update the Module's `__init__.py`

Create or update the `__init__.py` file in the `rdkit_shape` directory:

```python
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_shape_similarity import RDKitShapeSimilarity
from reinvent_scoring.scoring.score_components.rdkit_shape.parallel_rdkit_shape_similarity import ParallelRDKitShapeSimilarity
```

### Step 5: Update the Main `__init__.py`

Update the main `__init__.py` file in the `score_components` directory to include the new module:

```python
# Existing imports...
from reinvent_scoring.scoring.score_components.rdkit_shape import *
```

### Step 6: Update the Scoring Function Component Enum

Add the new components to the enum in `scoring_function_component_enum.py`:

```python
@dataclass(frozen=True)
class ScoringFunctionComponentNameEnum:
    # Existing components...
    RDKIT_SHAPE_SIMILARITY = "rdkit_shape_similarity"
    PARALLEL_RDKIT_SHAPE_SIMILARITY = "parallel_rdkit_shape_similarity"
```

### Step 7: Register the Components in the Factory

Update the component factory to include the new components:

```python
def _deafult_scoring_component_registry(self) -> dict:
    enum = ScoringFunctionComponentNameEnum()
    component_map = {
        # Existing components...
        enum.RDKIT_SHAPE_SIMILARITY: RDKitShapeSimilarity,
        enum.PARALLEL_RDKIT_SHAPE_SIMILARITY: ParallelRDKitShapeSimilarity,
    }
    return component_map
```

