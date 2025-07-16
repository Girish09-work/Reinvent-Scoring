from typing import List, Dict, Union, Optional
import os
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
# Import rdMolDescriptors for shape similarity calculations
from rdkit.Chem import rdMolDescriptors

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_conformer_generator import RDKitConformerGenerator


class RDKitShapeSimilarity(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        self.shape_weight = self.parameters.specific_parameters.get("shape_weight", 0.5)
        self.color_weight = self.parameters.specific_parameters.get("color_weight", 0.5)
        self.reference_file = self.parameters.specific_parameters.get("reference_file", "")
        self.method = self.parameters.specific_parameters.get("method", "usrcat")
        self.max_confs = self.parameters.specific_parameters.get("max_confs", 50)

        # Overlay saving parameters
        self.save_overlays = self.parameters.specific_parameters.get("save_overlays", False)
        self.overlays_dir = self.parameters.specific_parameters.get("overlays_dir", "rdkit_shape_overlays")
        self.overlay_prefix = self.parameters.specific_parameters.get("overlay_prefix", "mol_")

        # Multiple reference files
        self.warhead1_reference = self.parameters.specific_parameters.get("warhead1_reference", "")
        self.warhead2_reference = self.parameters.specific_parameters.get("warhead2_reference", "")

        # Minimization configuration
        self.minimization_config = self.parameters.specific_parameters.get("minimization_config", "")

        # Create overlays directory if saving is enabled
        if self.save_overlays:
            Path(self.overlays_dir).mkdir(parents=True, exist_ok=True)

        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=self.max_confs,
            energy_window=self.parameters.specific_parameters.get("ewindow", 10),
            max_stereo=self.parameters.specific_parameters.get("max_stereo", 0),
            minimization_config=self.minimization_config
        )

        # Load reference molecules
        self.reference_mol = self._load_reference_molecule()
        self.warhead1_mol = self._load_reference_molecule(self.warhead1_reference)
        self.warhead2_mol = self._load_reference_molecule(self.warhead2_reference)

    def _load_reference_molecule(self, file_path: str = None):
        """Load the reference molecule from file and generate conformers.

        Args:
            file_path: Path to the reference molecule file. If None, uses self.reference_file.

        Returns:
            RDKit molecule with conformers or None if loading fails.
        """
        # Use the provided file path or fall back to the default reference file
        ref_file = file_path if file_path else self.reference_file

        if not ref_file:
            return None

        # Load reference molecule
        ref_mol = None
        if ref_file.endswith('.sdf'):
            try:
                suppl = Chem.SDMolSupplier(ref_file)
                if suppl and len(suppl) > 0:
                    ref_mol = suppl[0]
            except Exception:
                return None
        else:
            return None

        # Generate conformers for reference molecule if needed
        if ref_mol and ref_mol.GetNumConformers() == 0:
            try:
                ref_smiles = Chem.MolToSmiles(ref_mol)
                ref_mol = self.conformer_generator.generate_conformers(ref_smiles)
            except Exception:
                return None

        return ref_mol

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a list of molecules."""
        # Convert input to RDKit molecules, handling both SMILES strings and RDKit Mol objects
        rdkit_mols = []
        for mol in molecules:
            if isinstance(mol, str):
                # Input is a SMILES string
                rdkit_mols.append(Chem.MolFromSmiles(mol))
            else:
                # Input is already a RDKit Mol object
                rdkit_mols.append(mol)

        # Calculate scores
        scores = self._calculate_shape_scores(rdkit_mols)

        # Create and return component summary
        score_summary = ComponentSummary(total_score=scores, parameters=self.parameters)
        return score_summary

    def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a specific step (used in reinforcement learning)."""
        return self.calculate_score(molecules, step)

    def _calculate_shape_scores(self, mols):
        """Calculate shape similarity scores for a list of RDKit molecules."""
        scores = []

        # Open SDF file for saving overlays if enabled
        overlay_writer = None
        if self.save_overlays:
            overlay_file_path = os.path.join(self.overlays_dir, f"{self.overlay_prefix}.sdf")
            try:
                overlay_writer = Chem.SDWriter(overlay_file_path)
            except Exception as e:
                print(f"Warning: Could not create overlay file: {e}")

        for idx, mol in enumerate(mols):
            if mol is None or self.reference_mol is None:
                scores.append(0.0)
                continue

            # Generate conformers using the conformer generator
            # First ensure we have a SMILES string
            try:
                mol_smiles = Chem.MolToSmiles(mol)
                mol_with_confs = self.conformer_generator.generate_conformers(mol_smiles)

                if mol_with_confs is None or mol_with_confs.GetNumConformers() == 0:
                    scores.append(0.0)
                    continue
            except Exception:
                # If we can't convert to SMILES or generate conformers, assign a zero score
                scores.append(0.0)
                continue

            # Calculate best similarity score across all conformer pairs
            best_score, best_conf_mol = self._calculate_best_similarity(mol_with_confs, return_mol=True)
            scores.append(best_score)

            # Save the best conformer if overlay saving is enabled
            if self.save_overlays and overlay_writer and best_conf_mol is not None:
                # Add properties to the molecule
                best_conf_mol.SetProp("SMILES", mol_smiles)
                best_conf_mol.SetProp("Shape_Score", f"{best_score:.4f}")
                best_conf_mol.SetProp("Molecule_Index", str(idx))

                # Write the molecule to the SDF file
                try:
                    overlay_writer.write(best_conf_mol)
                except Exception as e:
                    print(f"Warning: Could not write molecule to overlay file: {e}")

        # Close the SDF writer if it was opened
        if overlay_writer:
            overlay_writer.close()

        return np.array(scores, dtype=np.float32)

    def _calculate_best_similarity(self, mol, return_mol=False):
        """Calculate the best similarity score across all conformer pairs.

        Args:
            mol: RDKit molecule with conformers
            return_mol: If True, returns the best conformer molecule along with the score

        Returns:
            If return_mol is False: best_score (float)
            If return_mol is True: (best_score, best_conf_mol) tuple
        """
        best_score = 0.0
        best_conf_id = -1
        best_ref_conf_id = -1

        if self.method == "usrcat":
            best_score, best_conf_id, best_ref_conf_id = self._calculate_usrcat_similarity(mol, return_conf_ids=True)
        elif self.method == "o3a":
            best_score, best_conf_id, best_ref_conf_id = self._calculate_o3a_similarity(mol, return_conf_ids=True)

        # If we need to return the molecule with the best conformer
        if return_mol and best_conf_id >= 0:
            # Create a new molecule with just the best conformer
            best_conf_mol = Chem.Mol(mol)

            # Keep only the best conformer
            for conf_id in range(best_conf_mol.GetNumConformers()-1, -1, -1):
                if conf_id != best_conf_id:
                    best_conf_mol.RemoveConformer(conf_id)

            # Apply transformation to align with reference if possible
            if self.method == "o3a" and best_ref_conf_id >= 0 and self.reference_mol is not None:
                try:
                    # Create O3A alignment
                    pyO3A = AllChem.GetO3A(best_conf_mol, self.reference_mol,
                                          confId1=0, confId2=best_ref_conf_id)
                    # Apply the transformation
                    pyO3A.Align()
                except Exception:
                    pass

            return best_score, best_conf_mol

        return best_score

    def _calculate_usrcat_similarity(self, mol, return_conf_ids=False):
        """Calculate similarity using USRCAT method.

        Args:
            mol: RDKit molecule with conformers
            return_conf_ids: If True, returns the best conformer IDs along with the score

        Returns:
            If return_conf_ids is False: best_score (float)
            If return_conf_ids is True: (best_score, best_conf_id, best_ref_conf_id) tuple
        """
        best_score = 0.0
        best_conf_id = -1
        best_ref_conf_id = -1

        for q_conf_id in range(mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                # Get USRCAT descriptors
                query_descriptor = rdMolDescriptors.GetUSRCAT(mol, confId=q_conf_id)
                ref_descriptor = rdMolDescriptors.GetUSRCAT(self.reference_mol, confId=r_conf_id)

                # Calculate similarity
                similarity = rdMolDescriptors.GetUSRScore(query_descriptor, ref_descriptor)

                if similarity > best_score:
                    best_score = similarity
                    best_conf_id = q_conf_id
                    best_ref_conf_id = r_conf_id

        # Check warhead references if available
        if self.warhead1_mol is not None:
            for q_conf_id in range(mol.GetNumConformers()):
                for r_conf_id in range(self.warhead1_mol.GetNumConformers()):
                    # Get USRCAT descriptors
                    query_descriptor = rdMolDescriptors.GetUSRCAT(mol, confId=q_conf_id)
                    ref_descriptor = rdMolDescriptors.GetUSRCAT(self.warhead1_mol, confId=r_conf_id)

                    # Calculate similarity
                    similarity = rdMolDescriptors.GetUSRScore(query_descriptor, ref_descriptor)

                    if similarity > best_score:
                        best_score = similarity
                        best_conf_id = q_conf_id
                        best_ref_conf_id = r_conf_id

        # Check second warhead reference if available
        if self.warhead2_mol is not None:
            for q_conf_id in range(mol.GetNumConformers()):
                for r_conf_id in range(self.warhead2_mol.GetNumConformers()):
                    # Get USRCAT descriptors
                    query_descriptor = rdMolDescriptors.GetUSRCAT(mol, confId=q_conf_id)
                    ref_descriptor = rdMolDescriptors.GetUSRCAT(self.warhead2_mol, confId=r_conf_id)

                    # Calculate similarity
                    similarity = rdMolDescriptors.GetUSRScore(query_descriptor, ref_descriptor)

                    if similarity > best_score:
                        best_score = similarity
                        best_conf_id = q_conf_id
                        best_ref_conf_id = r_conf_id

        if return_conf_ids:
            return best_score, best_conf_id, best_ref_conf_id
        return best_score

    def _calculate_o3a_similarity(self, mol, return_conf_ids=False):
        """Calculate similarity using O3A method.

        Args:
            mol: RDKit molecule with conformers
            return_conf_ids: If True, returns the best conformer IDs along with the score

        Returns:
            If return_conf_ids is False: best_score (float)
            If return_conf_ids is True: (best_score, best_conf_id, best_ref_conf_id) tuple
        """
        best_score = 0.0
        best_conf_id = -1
        best_ref_conf_id = -1

        # Check main reference
        if self.reference_mol is not None:
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
                            best_conf_id = q_conf_id
                            best_ref_conf_id = r_conf_id
                    except Exception:
                        continue

        # Check warhead1 reference if available
        if self.warhead1_mol is not None:
            for q_conf_id in range(mol.GetNumConformers()):
                for r_conf_id in range(self.warhead1_mol.GetNumConformers()):
                    try:
                        # Create O3A alignment
                        pyO3A = AllChem.GetO3A(mol, self.warhead1_mol,
                                              confId1=q_conf_id, confId2=r_conf_id)

                        # Get alignment score (shape similarity)
                        shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range

                        if shape_sim > best_score:
                            best_score = shape_sim
                            best_conf_id = q_conf_id
                            best_ref_conf_id = r_conf_id
                    except Exception:
                        continue

        # Check warhead2 reference if available
        if self.warhead2_mol is not None:
            for q_conf_id in range(mol.GetNumConformers()):
                for r_conf_id in range(self.warhead2_mol.GetNumConformers()):
                    try:
                        # Create O3A alignment
                        pyO3A = AllChem.GetO3A(mol, self.warhead2_mol,
                                              confId1=q_conf_id, confId2=r_conf_id)

                        # Get alignment score (shape similarity)
                        shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range

                        if shape_sim > best_score:
                            best_score = shape_sim
                            best_conf_id = q_conf_id
                            best_ref_conf_id = r_conf_id
                    except Exception:
                        continue

        if return_conf_ids:
            return best_score, best_conf_id, best_ref_conf_id
        return best_score

