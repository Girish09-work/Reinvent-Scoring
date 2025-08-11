#!/usr/bin/env python3
"""
robust_warhead_extractor.py
Extracts POI/E3 warheads from docked PROTAC poses while fixing valence errors.

Copyright (c) 2025
"""

import os
from typing import List, Tuple, Dict

from rdkit import Chem
from rdkit.Chem import rdFMCS


# --------------------------------------------------------------------------- #
# 1.  Reference Warhead Loader (strict mode)                                  #
# --------------------------------------------------------------------------- #
def load_reference_warheads(path: str) -> Tuple[List[Chem.Mol], List[str]]:
    """Load the two reference warheads in normal (strict) RDKit mode."""
    warheads, names = [], []
    for i, mol in enumerate(Chem.SDMolSupplier(path, removeHs=False)):
        if mol is None:
            raise ValueError(f"Reference warhead #{i} failed to parse!")
        warheads.append(mol)
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"warhead_{i}"
        names.append(name)
    return warheads, names


# --------------------------------------------------------------------------- #
# 2.  Valence-Repair Utilities                                                #
# --------------------------------------------------------------------------- #
# def _fix_atom_valence(atom: Chem.Atom) -> None:
#     """
#     Heuristic patch: turn 4-valent neutral N → N+1 (most frequent docking error).
#     Extend as needed for other hypervalent centers.
#     """
#     explicit_valence = atom.GetExplicitValence()
#     charge = atom.GetFormalCharge()
#     Z = atom.GetAtomicNum()

#     # Nitrogen: allowed valence 3 (neutral) or 4 (quaternary cation).
#     if Z == 7 and explicit_valence == 4 and charge == 0:
#         atom.SetFormalCharge(1)                                      # +1 charge
#         return
#     if Z == 6 and explicit_valence == 4 and charge == 0:
#         atom.SetFormalCharge(1)                                      # +1 charge
#         return
#     # Carbon: occasionally comes in as 5-valent (Texas carbon).
#     # Here we simply mark it unsanitized; more complex chemistry would be needed.
#     # Add similar clauses for P, S, ClO2⁻, PF6⁻, etc. if your dataset contains them.

def _fix_atom_valence(atom: Chem.Atom) -> None:
    """
    Heuristic patch:
    - For nitrogen: turn 4-valent neutral N → N+1 (quaternary ammonium).
    - For carbon: turn 5-valent neutral C → C+1 (Texas carbon cation).
    Extend as needed for other hypervalent centers.
    """
    explicit_valence = atom.GetExplicitValence()
    charge = atom.GetFormalCharge()
    Z = atom.GetAtomicNum()

    # Nitrogen: allowed valence 3 (neutral) or 4 (quaternary cation)
    if Z == 7 and explicit_valence == 4 and charge == 0:
        atom.SetFormalCharge(1)  # Quaternary ammonium
        return

    # Carbon: allowed valence 4 (neutral) or 5 (cationic, very rare/Texas carbon)
    if Z == 6 and explicit_valence == 5 and charge == 0:
        atom.SetFormalCharge(1)  # Texas carbon cation
        return

    # Extend here for phosphorus, sulfur, etc.

def attempt_valence_repair(mol: Chem.Mol) -> Chem.Mol:
    """
    1) Detect chemistry problems, 2) run atom-level fixes, 3) re-sanitize.
    If still unsanitizable, return molecule in 'partially fixed' state.
    """
    mol.UpdatePropertyCache(strict=False)                              # preserves coords
    problems = Chem.DetectChemistryProblems(mol)                      # RDKit ≥2019.09[11]

    for p in problems:
        if p.GetType() == "AtomValenceException":
            _fix_atom_valence(mol.GetAtomWithIdx(p.GetAtomIdx()))

    # Try full sanitization—catch without crashing.
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        # Leave molecule in semi-sanitized state; substructure search will still work.
        print(f"   ↳  residual problem ({exc}); proceeding with unsanitized mol")
    return mol


# --------------------------------------------------------------------------- #
# 3.  PROTAC Loader with Tiered Sanitization                                  #
# --------------------------------------------------------------------------- #
def load_and_repair_protacs(path: str) -> List[Tuple[Chem.Mol, str, Dict[str, str]]]:
    """
    Returns list of tuples  (mol, name, property_dict).
    Applies a three-phase loading strategy:
        A) sanitize=False   • always succeeds (geometry kept)          [22]
        B) attempt atomic repairs (e.g., N⁺)                           [21]
        C) if sanitization still fails, keep mol unsanitized but usable.
    """
    repaired: List[Tuple[Chem.Mol, str, Dict[str, str]]] = []
    supplier = Chem.SDMolSupplier(path, sanitize=False, removeHs=False)

    for idx, mol in enumerate(supplier):
        if mol is None:
            print(f"  parse failure at record {idx}; skipping.")
            continue

        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"pose_{idx}"
        props = {p: mol.GetProp(p) for p in mol.GetPropNames()}

        try:
            mol = attempt_valence_repair(mol)
        except Exception as exc:
            print(f"  could not repair {name}: {exc}")

        repaired.append((mol, name, props))
    return repaired


# --------------------------------------------------------------------------- #
# 4.  Substructure & MCS Matching                                             #
# --------------------------------------------------------------------------- #
def find_warhead_matches(whole: Chem.Mol,
                         refs: List[Chem.Mol],
                         mcs_fallback: bool = True,
                         mcs_threshold: float = 0.8) -> List[Tuple[int, Tuple[int]]]:
    """
    Returns list of (ref_index, atom_idx_tuple) found in 'whole'.
    First exact SMARTS; optional MCS fallback increases tolerance for linker edits.
    """
    matches = []
    for i, ref in enumerate(refs):
        if whole.HasSubstructMatch(ref):
            matches.append((i, whole.GetSubstructMatch(ref)))
        elif mcs_fallback:
            res = rdFMCS.FindMCS([whole, ref],
                                 threshold=mcs_threshold,
                                 timeout=30,
                                 atomCompare=rdFMCS.AtomCompare.CompareElements,
                                 bondCompare=rdFMCS.BondCompare.CompareOrder)
            if res.numAtoms:
                patt = Chem.MolFromSmarts(res.smartsString)
                if whole.HasSubstructMatch(patt):
                    matches.append((i, whole.GetSubstructMatch(patt)))
    return matches


# --------------------------------------------------------------------------- #
# 5.  Coordinate-Preserving Extraction                                        #
# --------------------------------------------------------------------------- #
def extract_fragment(parent: Chem.Mol,
                     atom_indices: Tuple[int]) -> Chem.Mol:
    """
    Build a new ROMol consisting solely of the indices in 'atom_indices',
    retain exact 3-D coordinates.
    """
    emol = Chem.EditableMol(Chem.Mol())                                 # start empty
    old2new = {}

    # 5-A: copy atoms
    for idx in atom_indices:
        new_idx = emol.AddAtom(parent.GetAtomWithIdx(idx))
        old2new[idx] = new_idx

    # 5-B: copy bonds among kept atoms
    for i, ai in enumerate(atom_indices):
        for aj in atom_indices[i + 1:]:
            bond = parent.GetBondBetweenAtoms(ai, aj)
            if bond:
                emol.AddBond(old2new[ai], old2new[aj], bond.GetBondType())

    frag = emol.GetMol()

    # 5-C: transfer coordinates verbatim
    if parent.GetNumConformers():
        pconf = parent.GetConformer()
        conf = Chem.Conformer(frag.GetNumAtoms())
        for new_idx, old_idx in enumerate(atom_indices):
            conf.SetAtomPosition(new_idx, pconf.GetAtomPosition(old_idx))
        frag.AddConformer(conf, assignId=True)

    return frag


# --------------------------------------------------------------------------- #
# 6+.  Production-friendly API                                                #
# --------------------------------------------------------------------------- #

def extract_warheads_to_file(
    warheads_sdf: str,
    protac_sdf: str,
    output_sdf: str,
    *,
    mcs_fallback: bool = True,
    mcs_threshold: float = 0.8,
) -> int:
    """
    Extract warheads from a PROTAC docked SDF using provided reference warheads, and
    write them into output_sdf. Returns the number of fragments written.

    Parameters:
        warheads_sdf: Path to reference warheads SDF (strict RDKit loading)
        protac_sdf:   Path to docked PROTAC poses (coordinates preserved)
        output_sdf:   Path to write extracted warheads SDF
        mcs_fallback: Enable MCS-based relaxed matching when exact substructure fails
        mcs_threshold: Similarity threshold for MCS fallback
    """
    if not os.path.isfile(warheads_sdf):
        raise FileNotFoundError(f"Missing reference warheads file: {warheads_sdf}")
    if not os.path.isfile(protac_sdf):
        raise FileNotFoundError(f"Missing docked PROTAC SDF: {protac_sdf}")

    refs, _ = load_reference_warheads(warheads_sdf)
    protacs = load_and_repair_protacs(protac_sdf)

    extracted = []
    for mol, name, props in protacs:
        hits = find_warhead_matches(mol, refs, mcs_fallback=mcs_fallback, mcs_threshold=mcs_threshold)
        for ref_idx, atom_tuple in hits:
            frag = extract_fragment(mol, atom_tuple)
            extracted.append((name, frag, ref_idx, props))

    if extracted:
        write_warheads(extracted, output_sdf)
    return len(extracted)


# --------------------------------------------------------------------------- #
# 6.  Output Writer                                                           #
# --------------------------------------------------------------------------- #
def write_warheads(records, out_path: str) -> None:
    writer = Chem.SDWriter(out_path)
    for parent_name, frag, ref_idx, parent_props in records:
        frag.SetProp("_Name", f"{parent_name}_warhead_{ref_idx}")
        frag.SetProp("extracted_from", parent_name)
        frag.SetProp("warhead_index", str(ref_idx))
        for k, v in parent_props.items():
            frag.SetProp(f"parent_{k}", str(v))
        writer.write(frag)
    writer.close()


# --------------------------------------------------------------------------- #
# 7.  Main Orchestrator                                                       #
# --------------------------------------------------------------------------- #
def main():
    # 7-A: file paths (EDIT)
    warheads_sdf   = "reference_warheads_postDock.sdf"
    protac_sdf = "docked_pose.sdf"
    output_sdf   = "extracted_warheads3.sdf"

    # 7-B: sanity checks
    for f in (protac_sdf, warheads_sdf):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing file: {f}")

    # 7-C: load data
    refs, ref_names = load_reference_warheads(warheads_sdf)
    print(f"Loaded {len(refs)} reference warheads: {ref_names}")

    protacs = load_and_repair_protacs(protac_sdf)
    print(f"Loaded {len(protacs)} PROTAC poses (after repair stage)")

    # 7-D: extraction loop
    extracted = []
    for mol, name, props in protacs:
        hits = find_warhead_matches(mol, refs)
        if not hits:
            print(f"   no match found in {name}")
            continue

        for ref_idx, atom_tuple in hits:
            frag = extract_fragment(mol, atom_tuple)
            extracted.append((name, frag, ref_idx, props))
            print(f"   extracted warhead {ref_idx} from {name}")

    # 7-E: write results
    if extracted:
        write_warheads(extracted, output_sdf)
        print(f"\n  Saved {len(extracted)} warheads →  {output_sdf}")
    else:
        print("  No warheads were extracted — check reference structures.")


if __name__ == "__main__":
    main()
