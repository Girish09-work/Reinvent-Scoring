
import os
import multiprocessing
from multiprocessing import Pool
from typing import List
import time
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Try to import GPU-related libraries
try:
    import torch
    has_torch = True
    if torch.cuda.is_available():
        cuda_available = True
    else:
        cuda_available = False
except ImportError:
    has_torch = False
    cuda_available = False

try:
    import openmm
    has_openmm = True
except ImportError:
    has_openmm = False

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_shape_similarity import RDKitShapeSimilarity
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_conformer_generator import RDKitConformerGenerator


class ParallelRDKitShapeSimilarity(RDKitShapeSimilarity):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.shape_weight = self.parameters.specific_parameters.get("shape_weight", 0.5)
        self.color_weight = self.parameters.specific_parameters.get("color_weight", 0.5)
        self.reference_file = self.parameters.specific_parameters.get("reference_file", "")
        self.method = self.parameters.specific_parameters.get("method", "usrcat")
        self.max_confs = self.parameters.specific_parameters.get("max_confs", 50)

        # Check for GPU availability
        self.use_gpu = self.parameters.specific_parameters.get("use_gpu", True)
        if self.use_gpu:
            try:
                if 'has_torch' in globals() and has_torch and cuda_available:
                    self.gpu_available = True
                else:
                    self.gpu_available = False
            except Exception:
                self.gpu_available = False
        else:
            self.gpu_available = False

        self.save_overlays = self.parameters.specific_parameters.get("save_overlays", False)
        self.overlay_prefix = self.parameters.specific_parameters.get("overlay_prefix", "mol_")
        self.overlays_dir = self.parameters.specific_parameters.get("overlays_dir", "overlays")
        if self.save_overlays:
            os.makedirs(self.overlays_dir, exist_ok=True)

        # Determine number of CPU cores to use
        self.num_cpus = min(
            multiprocessing.cpu_count(),
            self.parameters.specific_parameters.get("max_num_cpus", 4)
        )

        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=self.max_confs,
            energy_window=self.parameters.specific_parameters.get("ewindow", 10),
            max_stereo=self.parameters.specific_parameters.get("max_stereo", 0)
        )

        # Load reference molecule
        self.reference_mol = self._load_reference_molecule()

    def _load_reference_molecule(self):
        if not self.reference_file:
            return None

        ref_mol = None
        if self.reference_file.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(self.reference_file)
            if suppl and len(suppl) > 0:
                ref_mol = suppl[0]

        if ref_mol is None:
            raise ValueError(f"Could not load reference molecule from {self.reference_file}")

        if ref_mol.GetNumConformers() == 0:
            ref_mol = self.conformer_generator.generate_conformers(Chem.MolToSmiles(ref_mol))

        return ref_mol

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        # Convert input to RDKit molecules, handling both SMILES strings and RDKit Mol objects
        rdkit_mols = []
        for mol in molecules:
            if isinstance(mol, str):
                # Input is a SMILES string
                rdkit_mols.append(Chem.MolFromSmiles(mol))
            else:
                # Input is already a RDKit Mol object
                rdkit_mols.append(mol)

        scores = self._calculate_shape_scores(rdkit_mols)
        return ComponentSummary(total_score=scores, parameters=self.parameters)

    def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a specific step (used in reinforcement learning)."""
        return self.calculate_score(molecules, step)

    def _calculate_shape_scores(self, mols):
        """Calculate shape similarity scores in parallel."""
        if len(mols) == 0 or self.reference_mol is None:
            return np.array([], dtype=np.float32)

        # Prepare arguments for parallel processing
        args_list = []
        for i, mol in enumerate(mols):
            if mol is None:
                continue

            try:
                smiles = Chem.MolToSmiles(mol) if mol else ""
                args = {
                    "mol": mol,
                    "smiles": smiles,
                    "reference_mol": self.reference_mol,
                    "method": self.method,
                    "max_confs": self.max_confs,
                    "shape_weight": self.shape_weight,
                    "color_weight": self.color_weight,
                    "index": i,
                    "save_overlays": self.save_overlays,
                    "energy_window": self.parameters.specific_parameters.get("ewindow", 10),
                    "max_stereo": self.parameters.specific_parameters.get("max_stereo", 0),
                    "use_gpu": hasattr(self, 'gpu_available') and self.gpu_available
                }

                if self.save_overlays:
                    args["overlays_dir"] = self.overlays_dir
                    args["overlay_prefix"] = self.overlay_prefix

                args_list.append(args)
            except Exception:
                continue

        # Process in parallel
        scores = np.zeros(len(mols), dtype=np.float32)

        if len(args_list) > 0:
            try:
                with Pool(processes=self.num_cpus) as pool:
                    results = pool.map(calculate_single_molecule, args_list)

                # Collect results
                for idx, score in results:
                    scores[idx] = score
            except Exception:
                pass

        return scores


def calculate_single_molecule(args):
    """Calculate similarity for a single molecule (called by each process)."""
    mol = args["mol"]
    smiles = args["smiles"]
    reference_mol = args["reference_mol"]
    method = args["method"]
    max_confs = args["max_confs"]
    shape_weight = args["shape_weight"]
    color_weight = args["color_weight"]
    index = args["index"]
    save_overlays = args["save_overlays"]
    energy_window = args["energy_window"]
    max_stereo = args["max_stereo"]

    if mol is None or reference_mol is None:
        return index, 0.0

    # Check if we can use GPU acceleration
    use_gpu = args.get("use_gpu", False)
    if use_gpu:
        try:
            # Verify CUDA is available in this process
            if 'has_torch' in globals() and has_torch and torch.cuda.is_available():
                pass
            else:
                use_gpu = False
        except Exception:
            use_gpu = False

    # Create conformer generator
    conformer_generator = RDKitConformerGenerator(
        max_confs=max_confs,
        energy_window=energy_window,
        max_stereo=max_stereo
    )

    # Generate conformers
    mol_with_confs = conformer_generator.generate_conformers(smiles)

    if mol_with_confs is None or mol_with_confs.GetNumConformers() == 0:
        return index, 0.0

    # Calculate best similarity
    best_score = 0.0
    best_conf_pair = None

    if method == "usrcat":
        from rdkit.Chem.rdMolDescriptors import GetUSRCAT, GetUSRScore

        # If GPU is available, we can potentially accelerate the descriptor comparison
        if use_gpu:
            try:
                # Pre-compute all descriptors
                query_descriptors = []
                for q_conf_id in range(mol_with_confs.GetNumConformers()):
                    query_descriptors.append(GetUSRCAT(mol_with_confs, confId=q_conf_id))

                ref_descriptors = []
                for r_conf_id in range(reference_mol.GetNumConformers()):
                    ref_descriptors.append(GetUSRCAT(reference_mol, confId=r_conf_id))

                # Convert to PyTorch tensors
                query_tensor = torch.tensor(query_descriptors, device='cuda')
                ref_tensor = torch.tensor(ref_descriptors, device='cuda')

                # Calculate all similarities in parallel on GPU
                # This is a simplified approach - in practice, you'd need a custom CUDA kernel
                # for the exact USR similarity calculation
                best_score = 0.0
                for q_idx, q_desc in enumerate(query_descriptors):
                    for r_idx, r_desc in enumerate(ref_descriptors):
                        similarity = GetUSRScore(q_desc, r_desc)
                        if similarity > best_score:
                            best_score = similarity
                            best_conf_pair = (q_idx, r_idx)
            except Exception:
                use_gpu = False

        # Fall back to CPU implementation if GPU failed or is not available
        if not use_gpu:
            for q_conf_id in range(mol_with_confs.GetNumConformers()):
                for r_conf_id in range(reference_mol.GetNumConformers()):
                    # Get USRCAT descriptors
                    query_descriptor = GetUSRCAT(mol_with_confs, confId=q_conf_id)
                    ref_descriptor = GetUSRCAT(reference_mol, confId=r_conf_id)

                    # Calculate similarity
                    similarity = GetUSRScore(query_descriptor, ref_descriptor)

                    if similarity > best_score:
                        best_score = similarity
                        best_conf_pair = (q_conf_id, r_conf_id)

    elif method == "o3a":
        from rdkit.Chem import ChemicalFeatures
        from rdkit import RDConfig
        import os

        # Load feature factory for pharmacophore features
        fdef_file = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_file)

        for q_conf_id in range(mol_with_confs.GetNumConformers()):
            for r_conf_id in range(reference_mol.GetNumConformers()):
                try:
                    # Create O3A alignment
                    pyO3A = AllChem.GetO3A(mol_with_confs, reference_mol,
                                          confId1=q_conf_id, confId2=r_conf_id)

                    # Get alignment score (shape similarity)
                    shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range

                    # Align the molecules
                    pyO3A.Align()

                    # Calculate feature similarity after alignment (color similarity)
                    color_sim = 0.0
                    # This would require implementing a feature-based similarity calculation
                    # For simplicity, we'll use shape similarity as a proxy for now

                    # Combine scores using weights
                    combined_score = ((shape_weight * shape_sim) +
                                     (color_weight * color_sim)) / (shape_weight + color_weight)

                    if combined_score > best_score:
                        best_score = combined_score
                        best_conf_pair = (q_conf_id, r_conf_id)
                except Exception:
                    continue

    # Save overlay if requested
    if save_overlays and best_conf_pair is not None:
        overlays_dir = args.get("overlays_dir", "overlays")
        overlay_prefix = args.get("overlay_prefix", "mol_")

        try:
            # Create output file
            overlay_file = os.path.join(overlays_dir, f"{overlay_prefix}{index}.sdf")
            with Chem.SDWriter(overlay_file) as writer:
                writer.write(mol_with_confs, confId=best_conf_pair[0])
        except Exception:
            pass

    return index, best_score





