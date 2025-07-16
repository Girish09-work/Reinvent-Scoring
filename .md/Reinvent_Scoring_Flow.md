# Reinvent 3.2 Scoring Components Flow in Link Invent Reinforcement Learning

This document explains how scoring components are used in Reinvent 3.2's reinforcement learning process, specifically focusing on the Link Invent module.

## Overview of Scoring Flow

In Reinvent 3.2, the reinforcement learning process for Link Invent follows these steps:

1. **Sampling**: Generate molecular structures
2. **Scoring**: Evaluate the generated structures using scoring components
3. **Updating**: Update the model based on scores
4. **Logging**: Record results

The scoring process is a critical part of the reinforcement learning loop, as it provides the feedback signal that guides the model's learning.

## Scoring Components Flow

### 1. Entry Point: `_scoring` Method

The scoring process begins in the `LinkInventReinforcementLearning` class with the `_scoring` method:

```python
def _scoring(self, sampled_sequences, step: int) -> FinalSummary:
    return self.scoring_strategy.evaluate(sampled_sequences, step)
```

This method delegates to the `evaluate` method of the `scoring_strategy` object.

### 2. Scoring Strategy Evaluation

The `LinkInventScoringStrategy` class implements the `evaluate` method:

```python
def evaluate(self, sampled_sequences: List[SampledSequencesDTO], step) -> FinalSummary:
    score_summary = self._apply_scoring_function(sampled_sequences, step)
    score_summary = self._clean_scored_smiles(score_summary)
    score_summary.total_score = self.diversity_filter.update_score(score_summary, sampled_sequences, step)
    return score_summary
```

This method:
- Applies the scoring function to the sampled sequences
- Cleans the scored SMILES strings
- Updates the score using the diversity filter
- Returns the final score summary

### 3. Applying the Scoring Function

The `_apply_scoring_function` method:

```python
def _apply_scoring_function(self, sampled_sequences: List[SampledSequencesDTO], step) -> FinalSummary:
    molecules = self._join_linker_and_warheads(sampled_sequences, keep_labels=True)
    smiles = []
    for idx, molecule in enumerate(molecules):
        try:
            smiles_str = self._conversion.mol_to_smiles(molecule) if molecule else "INVALID"
        except RuntimeError as exception:
            smiles_str = "INVALID"
            self.logger.log_message(exception.__str__() + f'\n\tinput: {sampled_sequences[idx].input}'
                                    f'\n\toutput: {sampled_sequences[idx].output}\n')
        finally:
            smiles.append(smiles_str)
    final_score: FinalSummary = self.scoring_function.get_final_score_for_step(smiles, step)
    return final_score
```

This method:
- Joins the linker and warheads to create complete molecules
- Converts the molecules to SMILES strings
- Calls the scoring function's `get_final_score_for_step` method

### 4. Scoring Function Factory

The scoring function is created by the `ScoringFunctionFactory`:

```python
def __new__(cls, sf_parameters: ScoringFunctionParameters) -> BaseScoringFunction:
    enum = ScoringFunctionNameEnum()
    scoring_function_registry = {
        enum.CUSTOM_PRODUCT: CustomProduct,
        enum.CUSTOM_SUM: CustomSum
    }
    return cls.create_scoring_function_instance(sf_parameters, scoring_function_registry)
```

This factory creates either a `CustomProduct` or `CustomSum` scoring function based on the parameters.

### 5. Score Component Factory

The scoring function uses the `ScoreComponentFactory` to create the individual scoring components:

```python
def create_score_components(self) -> [BaseScoreComponent]:
    def create_component(component_params):
        if component_params.component_type in self._current_components:
            component = self._current_components[component_params.component_type]
            component_instance = component(component_params)
        else:
            raise KeyError(f'Component: {component_params.component_type} is not implemented.'
                           f' Consider checking your input.')
        return component_instance

    components = [create_component(component) for component in self._parameters]
    return components
```

The factory maintains a registry of available scoring components and creates instances based on the provided parameters.

### 6. Scoring Component Types

Link Invent has specific scoring components designed for evaluating linkers:

- `LinkerEffectiveLength`
- `LinkerGraphLength`
- `LinkerLengthRatio`
- `LinkerNumRings`
- `LinkerNumAliphaticRings`
- `LinkerNumAromaticRings`
- `LinkerNumSPAtoms`
- `LinkerNumSP2Atoms`
- `LinkerNumSP3Atoms`
- `LinkerNumHBA`
- `LinkerNumHBD`
- `LinkerMolWeight`
- `LinkerRatioRotatableBonds`

Additionally, standard scoring components are available:

- `MatchingSubstructure`
- `RocsSimilarity`
- `PredictivePropertyComponent`
- `TanimotoSimilarity`
- `JaccardDistance`
- `CustomAlerts`
- `QedScore`
- `MolWeight`
- `PSA`
- `RotatableBonds`
- `GraphLength`
- `HBD_Lipinski`
- `HBA_Lipinski`
- `NumRings`
- `NumAromaticRings`
- `NumAliphaticRings`
- `SlogP`
- And many more...

### 7. Scoring Function Calculation

The base scoring function calculates scores in the `get_final_score_for_step` method:

```python
def get_final_score_for_step(self, smiles: List[str], step: int) -> FinalSummary:
    molecules, valid_indices = self._chemistry.smiles_to_mols_and_indices(smiles)
    query_size = len(smiles)
    summaries = [_update_total_score(sc.calculate_score_for_step(molecules, step), query_size, valid_indices) for sc
                 in self.scoring_components]
    return self._score_summary(summaries, smiles, valid_indices)
```

This method:
- Converts SMILES strings to molecules
- Calculates scores for each component
- Combines the scores into a final summary

### 8. Score Aggregation

The scores are aggregated differently depending on the scoring function type:

#### CustomSum

```python
def _compute_non_penalty_components(self, summaries: List[ComponentSummary], smiles: List[str]):
    total_sum = np.full(len(smiles), 0, dtype=np.float32)
    all_weights = 0.

    for summary in summaries:
        if not self._component_is_penalty(summary):
            total_sum = total_sum + summary.total_score * summary.parameters.weight
            all_weights += summary.parameters.weight

    if all_weights == 0:
        return np.full(len(smiles), 1, dtype=np.float32)

    return total_sum / all_weights
```

#### CustomProduct

```python
def _compute_non_penalty_components(self, summaries: List[ComponentSummary], smiles: List[str]):
    product = np.full(len(smiles), 1, dtype=np.float32)
    all_weights = self._get_all_weights(summaries)

    for summary in summaries:
        if not self._component_is_penalty(summary):
            comp_pow = self._calculate_pow(summary.total_score, summary.parameters.weight / all_weights)
            product = product * comp_pow

    return product
```

### 9. Penalty Components

Penalty components (like `CustomAlerts` and `MatchingSubstructure`) are handled separately:

```python
def _compute_penalty_components(self, summaries: List[ComponentSummary], smiles: List[str]):
    penalty = np.full(len(smiles), 1, dtype=np.float32)

    for summary in summaries:
        if self._component_is_penalty(summary):
            penalty = penalty * summary.total_score

    return penalty
```

The final score is the product of the non-penalty component score and the penalty component score.

## Configuration Example

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
scoring_function_params = ScoringFunctionParameters(
    name="custom_product",
    parameters=[tanimoto_params, qed_params, alerts_params],
    parallel=False
)
```

## Summary

1. The reinforcement learning process calls the `_scoring` method
2. The scoring strategy evaluates the sampled sequences
3. The scoring function calculates scores for each component
4. The scores are aggregated using either a sum or product approach
5. Penalty components are applied separately
6. The final score is returned and used to update the model

This scoring process is crucial for guiding the model's learning toward generating molecules with desired properties.
