# REINVENT Scoring Package: Detailed Component Guide

## 1. Standard Components

### TanimotoSimilarity
- **Purpose**: Measures chemical similarity between molecules using fingerprints
- **How it works**: 
  - Converts reference SMILES to fingerprints (Morgan/ECFP by default)
  - Converts query molecules to the same fingerprint type
  - Calculates Tanimoto coefficient between query and reference fingerprints
  - Parameters include radius (default 3), use_counts (default True), and use_features (default True)
- **Use case**: Guiding molecule generation toward similar structures to known actives

### JaccardDistance
- **Purpose**: Measures dissimilarity between molecules
- **How it works**: 
  - Similar to TanimotoSimilarity but returns 1-Tanimoto as a distance metric
  - Useful when you want to penalize similarity rather than reward it
- **Use case**: Promoting diversity away from certain structures

### MatchingSubstructure
- **Purpose**: Checks if molecules contain specific substructures
- **How it works**: 
  - Takes a list of SMARTS patterns
  - Returns 1 if the molecule contains any of the specified substructures, 0 otherwise
  - Acts as a binary filter
- **Use case**: Ensuring molecules contain essential pharmacophore elements

### CustomAlerts
- **Purpose**: Penalizes molecules containing unwanted substructures
- **How it works**: 
  - Takes a list of SMARTS patterns representing problematic substructures
  - Returns 0 if the molecule contains any alerts, 1 otherwise
  - Acts as a penalty component in the scoring function
- **Use case**: Avoiding known toxicophores or reactive groups

### QedScore
- **Purpose**: Quantitative Estimate of Drug-likeness
- **How it works**: 
  - Uses RDKit's QED implementation
  - Combines 8 desirability functions for different molecular properties
  - Returns a score between 0 and 1 (higher is more drug-like)
- **Use case**: General drug-likeness assessment

### PredictivePropertyComponent
- **Purpose**: Uses machine learning models to predict molecular properties
- **How it works**: 
  - Loads a pre-trained model
  - Generates descriptors for input molecules
  - Predicts properties using the model
- **Use case**: Predicting complex properties like activity, ADME, etc.

### SelectivityComponent
- **Purpose**: Evaluates selectivity against off-targets
- **How it works**: 
  - Calculates activity against a primary target
  - Calculates activity against off-targets
  - Computes the difference or ratio between them
- **Use case**: Designing selective compounds that hit target but avoid off-targets

## 2. Physicochemical Properties

### MolWeight
- **Purpose**: Calculates molecular weight
- **How it works**: Uses RDKit's descriptor calculation
- **Use case**: Keeping compounds in drug-like weight range (typically 200-500 Da)

### PSA (TPSA)
- **Purpose**: Calculates Topological Polar Surface Area
- **How it works**: Sums up contributions from polar atoms (O, N, etc.)
- **Use case**: Predicting membrane permeability (lower PSA = better permeability)

### RotatableBonds
- **Purpose**: Counts number of rotatable bonds
- **How it works**: Identifies single bonds between non-terminal heavy atoms
- **Use case**: Assessing molecular flexibility and oral bioavailability

### HBD_Lipinski & HBA_Lipinski
- **Purpose**: Count hydrogen bond donors and acceptors
- **How it works**: Identifies atoms that can donate or accept hydrogen bonds
- **Use case**: Lipinski's Rule of 5 compliance for oral bioavailability

### NumRings, NumAromaticRings, NumAliphaticRings
- **Purpose**: Count different types of rings in molecules
- **How it works**: Uses RDKit's ring detection algorithms
- **Use case**: Controlling molecular complexity and flatness

### SlogP
- **Purpose**: Calculates octanol-water partition coefficient
- **How it works**: Uses atom-based contribution method
- **Use case**: Predicting lipophilicity and membrane permeability

### GraphLength
- **Purpose**: Measures the maximum distance between atoms
- **How it works**: Calculates the longest path in the molecular graph
- **Use case**: Controlling molecular shape and size

### NumberOfStereoCenters
- **Purpose**: Counts stereocenters in a molecule
- **How it works**: Identifies atoms with tetrahedral stereochemistry
- **Use case**: Controlling stereochemical complexity

## 3. Synthetic Accessibility

### SASComponent
- **Purpose**: Evaluates synthetic accessibility
- **How it works**: 
  - Implements Ertl & Schuffenhauer's SA score
  - Considers fragment contributions, ring complexity, stereochemistry
  - Returns a score from 1 (easy to synthesize) to 10 (difficult)
- **Use case**: Ensuring generated molecules are synthetically feasible

## 4. 3D Structure-Based Components

### RocsSimilarity
- **Purpose**: Evaluates 3D shape and pharmacophore similarity
- **How it works**: 
  - Uses OpenEye ROCS to generate 3D conformers
  - Aligns molecules to reference structures
  - Calculates shape and color (pharmacophore) similarity
- **Use case**: Scaffold hopping and 3D similarity-based design

### ParallelRocsSimilarity
- **Purpose**: Parallel version of ROCS similarity for better performance
- **How it works**: Same as RocsSimilarity but with parallel processing
- **Use case**: Large-scale 3D similarity calculations

### AZdock & DockStream
- **Purpose**: Molecular docking to protein targets
- **How it works**: 
  - Interfaces with external docking software
  - Generates 3D conformers of molecules
  - Docks them into protein binding sites
  - Returns docking scores
- **Use case**: Structure-based design targeting specific proteins

## 5. PiP (Prediction in Pipeline) Components

### PiPPredictionComponent & PiPLogPredictionComponent
- **Purpose**: Interface with AstraZeneca's internal prediction services
- **How it works**: 
  - Sends SMILES to REST API endpoints
  - Retrieves predictions for various properties
  - PiPLogPredictionComponent applies log transformation to results
- **Use case**: Accessing proprietary predictive models

### RatPKPiP
- **Purpose**: Predicts rat pharmacokinetic parameters
- **How it works**: Specialized PiP component for PK predictions
- **Use case**: Optimizing compounds for good in vivo exposure

### QptunaPiPModelComponent
- **Purpose**: Uses Optuna-optimized models for predictions
- **How it works**: Interfaces with models optimized using Optuna framework
- **Use case**: Accessing ensemble or optimized predictive models

## 6. Linker-Specific Components

These components analyze properties specific to molecular linkers (parts connecting two fragments):

### LinkerGraphLength & LinkerEffectiveLength
- **Purpose**: Measure the size of linkers
- **How it works**: Calculate path lengths through the linker
- **Use case**: Controlling linker size in fragment linking

### LinkerNumRings, LinkerNumAromaticRings, LinkerNumAliphaticRings
- **Purpose**: Count rings in linker portions
- **How it works**: Identify rings specifically in the linker region
- **Use case**: Controlling linker rigidity and properties

### LinkerNumSPAtoms, LinkerNumSP2Atoms, LinkerNumSP3Atoms
- **Purpose**: Count atoms with different hybridizations in linkers
- **How it works**: Identify sp, sp², and sp³ hybridized atoms in linkers
- **Use case**: Controlling linker flexibility and 3D shape

### LinkerNumHBA, LinkerNumHBD
- **Purpose**: Count H-bond acceptors and donors in linkers
- **How it works**: Identify H-bond forming atoms in linker regions
- **Use case**: Controlling linker polarity and solubility

### LinkerMolWeight & LinkerRatioRotatableBonds
- **Purpose**: Measure linker weight and flexibility
- **How it works**: Calculate molecular weight and rotatable bond ratio in linkers
- **Use case**: Controlling linker properties for optimal joining of fragments

## How These Components Work Together

The scoring system allows you to combine these components in two main ways:

1. **CustomSum**: Weighted average of component scores
   ```
   FinalScore = (w₁·Score₁ + w₂·Score₂ + ... + wₙ·Scoreₙ) / (w₁ + w₂ + ... + wₙ)
   ```

2. **CustomProduct**: Weighted product of component scores
   ```
   FinalScore = Score₁^(w₁/W) · Score₂^(w₂/W) · ... · Scoreₙ^(wₙ/W)
   ```
   where W is the sum of all weights

Additionally, penalty components (like CustomAlerts and MatchingSubstructure) are handled separately and multiply the final score:
```
FinalScore = (Regular Components Score) · (Penalty₁) · (Penalty₂) · ...
```

## Example Configuration

Here's a simplified example of how to configure a scoring function:

```python
# Define component parameters
tanimoto_params = {
    "component_type": "tanimoto_similarity",
    "name": "Similarity to target",
    "weight": 1.0,
    "specific_parameters": {
        "smiles": ["CCO", "CCC"],  # Reference molecules
        "radius": 3,
        "use_counts": True,
        "use_features": True
    }
}

qed_params = {
    "component_type": "qed_score",
    "name": "Drug-likeness",
    "weight": 0.5,
    "specific_parameters": {}
}

alerts_params = {
    "component_type": "custom_alerts",
    "name": "Structural alerts",
    "weight": 1.0,
    "specific_parameters": {
        "smiles": ["[N+](=O)[O-]", "C(=O)Cl"]  # Nitro and acid chloride
    }
}

# Create scoring function parameters
from reinvent_scoring.scoring.scoring_function_parameters import ScoringFunctionParameters

scoring_function_params = ScoringFunctionParameters(
    name="custom_product",
    parameters=[tanimoto_params, qed_params, alerts_params],
    parallel=False
)

# Create scoring function
from reinvent_scoring.scoring.scoring_function_factory import ScoringFunctionFactory

scoring_function = ScoringFunctionFactory(scoring_function_params)

# Use scoring function
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]
results = scoring_function.get_final_score(smiles_list)
```

This example creates a scoring function that:
1. Rewards similarity to ethanol and propane
2. Rewards drug-likeness (QED)
3. Penalizes molecules containing nitro groups or acid chlorides
4. Combines the rewards using a weighted product
