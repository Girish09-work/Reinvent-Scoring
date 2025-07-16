# Scaffold CSV Generation Solution for Protac-Invent

## Problem Analysis

The issue was that `scaffold.csv` files were not being generated when running protac-invent with roshambo and rdkit_shape components during reinforcement learning. After analyzing the codebase, I found the root cause:

### The Real Issue

**The problem was NOT that scaffold CSV files are only generated during training modes.** The real issue was that:

1. **The diversity filter system** is responsible for generating scaffold_memory.csv files during reinforcement learning
2. **The diversity filter collects scores from ALL scoring components** via the `score_summary.scaffold_log`
3. **The roshambo and rdkit_shape components had a bug** where they didn't properly implement `calculate_score_for_step`
4. **During reinforcement learning**, the system calls `calculate_score_for_step` but the base class implementation was dropping the `step` parameter
5. **This caused the components to not be properly integrated** with the diversity filter system

### Technical Details

The base class `BaseScoreComponent` has this flawed implementation:
```python
def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
    return self.calculate_score(molecules)  # BUG: step parameter is lost!
```

But the scoring components expect:
```python
def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
```

## Solution

I've fixed the **root cause** by properly implementing `calculate_score_for_step` in the affected components. This ensures they integrate correctly with the diversity filter system during reinforcement learning.

## Implementation Details

### Modified Files

1. **`reinvent_scoring/scoring/score_components/roshambo/roshambo_shape_similarity.py`**
   - Fixed `calculate_score` method to accept `step` parameter
   - Added proper `calculate_score_for_step` implementation
   - Removed the custom scaffold CSV generation (no longer needed)

2. **`reinvent_scoring/scoring/score_components/rdkit_shape/rdkit_shape_similarity.py`**
   - Fixed `calculate_score` method to accept `step` parameter
   - Added proper `calculate_score_for_step` implementation

3. **`reinvent_scoring/scoring/score_components/rdkit_shape/parallel_rdkit_shape_similarity.py`**
   - Fixed `calculate_score` method to accept `step` parameter
   - Added proper `calculate_score_for_step` implementation

### Key Changes

The fix was simple but critical:

```python
# BEFORE (broken):
def calculate_score(self, molecules: List) -> ComponentSummary:
    # ... implementation

# AFTER (fixed):
def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
    # ... implementation

def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
    """Calculate scores for a specific step (used in reinforcement learning)."""
    return self.calculate_score(molecules, step)
```

### How the Diversity Filter Works

Now that the components are properly integrated, the diversity filter system will:

1. **Collect scores from all components** during each reinforcement learning step
2. **Calculate scaffolds** for molecules that meet the minimum score threshold
3. **Store molecules in memory** with their scaffolds and component scores
4. **Generate scaffold_memory.csv** automatically at regular intervals

### Generated CSV Format

The generated `scaffold_memory.csv` file will now include data from roshambo and rdkit_shape components:

| Column | Description |
|--------|-------------|
| Step | Reinforcement learning step number |
| Scaffold | Murcko scaffold SMILES |
| SMILES | Original molecule SMILES |
| [Component Name] | Score from each scoring component (e.g., "Roshambo Shape Similarity") |
| raw_[Component Name] | Raw score before transformation (if available) |
| total_score | Final combined score |
| ID | Unique identifier |

## Usage

### 1. Configuration File

Use any existing protac-invent configuration with roshambo or rdkit_shape components:

```json
{
    "run_type": "reinforcement_learning",
    "parameters": {
        "scoring_function": {
            "name": "custom_sum",
            "parameters": [
                {
                    "component_type": "roshambo_shape_similarity",
                    "name": "Roshambo Shape Similarity",
                    "weight": 1.0,
                    "specific_parameters": {
                        "reference_file": "path/to/reference.sdf",
                        "rdbase_path": "/path/to/rdkit",
                        "conda_env_name": "roshambo"
                    }
                }
            ]
        },
        "diversity_filter": {
            "name": "IdenticalMurckoScaffold",
            "minscore": 0.4,
            "bucket_size": 25
        }
    }
}
```

### 2. Running Protac-Invent

```bash
python input.py your_config.json
```

### 3. Output

After running, you'll find:
- **scaffold_memory.csv**: Contains molecules with scores above the diversity filter threshold
- **Normal protac-invent output**: Logs, checkpoints, etc.
- **Component scores included**: Roshambo and RDKit shape scores will appear in the CSV

## Benefits

1. **Proper integration**: Components now work correctly with the diversity filter system
2. **Automatic CSV generation**: scaffold_memory.csv is generated automatically during reinforcement learning
3. **Complete score data**: All component scores are included in the CSV
4. **No configuration changes needed**: Existing configurations will now work correctly

## Example Output

```csv
Step,Scaffold,SMILES,Roshambo Shape Similarity,total_score,ID
0,c1ccc2[nH]c3ccccc3c2c1,CN[C@@H](C)C(=O)N[C@@H]...,0.85,0.85,job_0
1,c1ccc2[nH]c3ccccc3c2c1,CNC(C)C(=O)NC(C(=O)N1Cc2cc...,0.78,0.78,job_1
```

## Troubleshooting

1. **No scaffold_memory.csv generated**:
   - Check that `run_type` is set to `reinforcement_learning` (not `scoring`)
   - Ensure diversity filter is configured with appropriate `minscore`
   - Verify molecules are meeting the minimum score threshold

2. **Empty scaffold_memory.csv**:
   - Lower the diversity filter `minscore` threshold
   - Check that scoring components are returning non-zero scores

3. **Missing component scores in CSV**:
   - This should now be fixed with the proper `calculate_score_for_step` implementation
   - Verify the component is properly configured in the scoring function

4. **Component errors**:
   - Check roshambo environment setup (rdbase_path, conda_env_name)
   - Verify reference files exist and are readable

## Testing the Fix

To verify the fix works:

1. **Run a short reinforcement learning job** with roshambo or rdkit_shape components
2. **Check the results directory** for scaffold_memory.csv
3. **Verify the CSV contains** columns for your scoring components
4. **Confirm scores are non-zero** and reasonable

## Root Cause Summary

The issue was a **method signature mismatch** in the scoring components. The diversity filter system calls `calculate_score_for_step(molecules, step)` during reinforcement learning, but the base class implementation was calling `calculate_score(molecules)` without the step parameter. This caused the components to not integrate properly with the diversity filter system, resulting in missing or incomplete scaffold_memory.csv files.

The fix ensures that the step parameter is properly passed through, allowing the diversity filter to correctly collect and store component scores.
