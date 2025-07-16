# Joining Linkers with Warheads in Protac-invent

## 1. Method Used to Join Linker and Warhead Molecules

The Protac-invent codebase uses a method called `_join_linker_and_warheads` to combine the generated linker with warhead molecules. This method is implemented in several strategy classes, including `LinkInventProductionStrategy` and `LinkInventScoringStrategy`.

```python
def _join_linker_and_warheads(self, sampled_sequences: List[SampledSequencesDTO], keep_labels=False):
    molecules = []
    for sample in sampled_sequences:
        linker = self._attachment_points.add_attachment_point_numbers(sample.output, canonicalize=False)
        molecule = self._bond_maker.join_scaffolds_and_decorations(linker, sample.input,
                                                               keep_labels_on_atoms=keep_labels)
        molecules.append(molecule)
    return molecules
```

The actual joining operation is delegated to the `join_scaffolds_and_decorations` method of the `_bond_maker` object, which is likely an instance of a class from the `reinvent_chemistry.library_design.bond_maker` module.

## 2. Bond Type Determination

Based on the code and comments, the system appears to use attachment point numbers to determine where and how to form bonds between the linker and warheads. The bond types are likely determined by the `BondMaker` class, which handles the chemistry of joining molecular fragments.

From the error handling in the code, we can see that there are limitations in the current implementation:

```python
try:
    smiles_str = self._conversion.mol_to_smiles(molecule) if molecule else "INVALID"
except RuntimeError as exception:
    # NOTE: Current implementation of BondMaker (reinvent_chemistry.library_design.bond_maker) results in
    # impossible conversion of mol to smiles if one single atom has two attachment points and labels are
    # kept. As this case is not relevant in the context of link_invent, then can be discarded as invalid.
    smiles_str = "INVALID"
    self._logger.log_message(exception.__str__() + f'\n\tinput: {sampled_sequences[idx].input}'
                                                   f'\n\toutput: {sampled_sequences[idx].output}\n')
```

This suggests that the bond maker creates covalent bonds at the attachment points, and the bond type is likely determined by the chemistry of the atoms at these attachment points.

## 3. Role of Attachment Points

Attachment points play a crucial role in the joining process. The code shows that before joining, attachment point numbers are added to the linker:

```python
linker = self._attachment_points.add_attachment_point_numbers(sample.output, canonicalize=False)
```

These attachment points serve as markers for where the linker should connect to the warheads. After scoring, these attachment point numbers are removed from the final molecules:

```python
def _clean_scored_smiles(self, score_summary: FinalSummary) -> FinalSummary:
    """
    Remove attachment point numbers from scored smiles
    """
    # Note: method AttachmentPoints.remove_attachment_point_numbers does not work in this context, as it searches
    # for attachment point token ('*')
    score_summary.scored_smiles = [self._conversion.mol_to_smiles(
        self._attachment_points.remove_attachment_point_numbers_from_mol(self._conversion.smile_to_mol(smile))
    ) if idx in score_summary.valid_idxs else smile for idx, smile in enumerate(score_summary.scored_smiles)]
    return score_summary
```

## 4. Code Snippets for Joining Operation

The main joining operation is performed in the `_apply_scoring_function` method, which calls `_join_linker_and_warheads` and then processes the resulting molecules:

```python
def _apply_scoring_function(self, sampled_sequences: List[SampledSequencesDTO], step) -> FinalSummary:
    print("DEBUG: Starting _apply_scoring_function")
    molecules = self._join_linker_and_warheads(sampled_sequences, keep_labels=True)
    print(f"DEBUG: Joined molecules, count: {len(molecules)}")
    smiles = []
    for idx, molecule in enumerate(molecules):
        try:
            smiles_str = self._conversion.mol_to_smiles(molecule) if molecule else "INVALID"
            if idx % 20 == 0:  # Print progress every 20 molecules
                print(f"DEBUG: Processed {idx}/{len(molecules)} molecules")
        except RuntimeError as exception:
            # NOTE: Current implementation of BondMaker (reinvent_chemistry.library_design.bond_maker) results in
            # impossible conversion of mol to smiles if one single atom has two attachment points and labels are
            # kept. As this case is not relevant in the context of link_invent, then can be discarded as invalid.
            smiles_str = "INVALID"
        finally:
            smiles.append(smiles_str)
    final_score: FinalSummary = self._scoring_function.get_final_score_for_step(smiles, step)
    return final_score
```

In contrast, the `LibInventScoringStrategy` uses a different method called `_join_scaffolds_and_decorations`:

```python
def _join_scaffolds_and_decorations(self, sampled_sequences: List[SampledSequencesDTO]):
    molecules = []
    for sample in sampled_sequences:
        scaffold = self._attachment_points.add_attachment_point_numbers(sample.input, canonicalize=False)
        molecule = self._bond_maker.join_scaffolds_and_decorations(scaffold, sample.output)
        molecules.append(molecule)
    return molecules
```

Note the difference: in `LinkInventScoringStrategy`, the linker is the output and the warheads are the input, while in `LibInventScoringStrategy`, the scaffold is the input and the decorations are the output.

## 5. Scoring of Complete PROTAC Molecules

Yes, the scoring functions evaluate the complete joined PROTAC molecule (linker + warheads) during the reinforcement learning iterations. This is evident from the code flow:

1. The linker and warheads are joined to form complete molecules
2. These molecules are converted to SMILES strings
3. The SMILES strings are passed to the scoring function
4. The scoring function returns a `FinalSummary` object containing the scores

The scoring components specifically designed for linkers, such as `LinkerNumHBD`, `LinkerNumHBA`, `LinkerMolWeight`, and `LinkerRatioRotatableBonds`, evaluate properties of the linker portion of the complete molecule:

```python
class LinkerRatioRotatableBonds(BaseLinkInventComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

    def _calculate_linker_property(self, labeled_mol):
        return self._linker_descriptor.ratio_rotatable_bonds(labeled_mol)
```

This confirms that the scoring is performed on the complete molecule, with specific components focusing on the linker portion.