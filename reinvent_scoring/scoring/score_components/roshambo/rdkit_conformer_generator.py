"""
RDKit-based conformer generator for Roshambo shape similarity calculations.
This module provides conformer generation functionality using RDKit instead of OpenEye's OMEGA.
"""

import os
import tempfile
from typing import List, Optional, Union
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
import numpy as np


class RoshamboConformerGenerator:
    """RDKit-based conformer generator optimized for Roshambo shape similarity calculations."""

    def __init__(self, 
                 n_confs: int = 10,
                 method: str = "ETKDGv3",
                 random_seed: int = 42,
                 ff: str = "MMFF94s",
                 add_hs: bool = True,
                 opt_confs: bool = True,
                 calc_energy: bool = False,
                 energy_iters: int = 200,
                 energy_cutoff: float = np.inf,
                 align_confs: bool = False,
                 rms_cutoff: Optional[float] = None,
                 num_threads: int = 1,
                 max_stereo: int = 0):
        """
        Initialize the conformer generator.
        
        Args:
            n_confs: Number of conformers to generate
            method: Embedding method (ETDG, ETKDG, ETKDGv2, ETKDGv3)
            random_seed: Random seed for reproducibility
            ff: Force field for optimization (UFF, MMFF94s, MMFF94s_noEstat)
            add_hs: Whether to add hydrogens
            opt_confs: Whether to optimize conformers
            calc_energy: Whether to calculate energies
            energy_iters: Maximum iterations for energy minimization
            energy_cutoff: Energy cutoff for conformer filtering
            align_confs: Whether to align conformers
            rms_cutoff: RMS cutoff for conformer clustering
            num_threads: Number of threads for parallel processing
            max_stereo: Maximum number of stereoisomers to enumerate
        """
        self.n_confs = n_confs
        self.method = method
        self.random_seed = random_seed
        self.ff = ff
        self.add_hs = add_hs
        self.opt_confs = opt_confs
        self.calc_energy = calc_energy
        self.energy_iters = energy_iters
        self.energy_cutoff = energy_cutoff
        self.align_confs = align_confs
        self.rms_cutoff = rms_cutoff
        self.num_threads = num_threads
        self.max_stereo = max_stereo
        
        # Validate parameters
        self._validate_parameters()

    def _validate_parameters(self):
        """Validate initialization parameters."""
        valid_methods = ["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3"]
        if self.method not in valid_methods:
            raise ValueError(f"Method {self.method} not supported. Choose from: {valid_methods}")
        
        valid_ffs = ["UFF", "MMFF94s", "MMFF94s_noEstat"]
        if self.ff not in valid_ffs:
            raise ValueError(f"Force field {self.ff} not supported. Choose from: {valid_ffs}")

    def generate_conformers_from_smiles(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Generate conformers from a SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            RDKit molecule with conformers or None if failed
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        return self.generate_conformers(mol)

    def generate_conformers(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Generate conformers for a molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            RDKit molecule with conformers or None if failed
        """
        if mol is None:
            return None
        
        # Make a copy to avoid modifying the original
        mol = Chem.Mol(mol)
        
        # Add hydrogens if requested
        if self.add_hs:
            mol = Chem.AddHs(mol)
        
        # Handle stereochemistry enumeration if requested
        if self.max_stereo > 0:
            return self._generate_with_stereoisomers(mol)
        else:
            return self._generate_conformers_for_mol(mol)

    def _generate_with_stereoisomers(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Generate conformers with stereoisomer enumeration."""
        isomers = self._enumerate_stereoisomers(mol)
        if not isomers:
            return self._generate_conformers_for_mol(mol)
        
        all_confs_mol = None
        for isomer in isomers:
            confs_mol = self._generate_conformers_for_mol(isomer)
            if confs_mol and confs_mol.GetNumConformers() > 0:
                if all_confs_mol is None:
                    all_confs_mol = confs_mol
                else:
                    # Add conformers from this isomer
                    for conf in confs_mol.GetConformers():
                        all_confs_mol.AddConformer(conf, assignId=True)
        
        return all_confs_mol

    def _enumerate_stereoisomers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """Enumerate stereoisomers of a molecule."""
        try:
            opts = EnumerateStereoisomers.StereoEnumerationOptions(
                maxIsomers=self.max_stereo,
                onlyUnassigned=False,
                unique=True
            )
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
            return isomers
        except Exception:
            return [mol]  # Return original molecule if enumeration fails

    def _generate_conformers_for_mol(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Generate conformers for a single molecule."""
        try:
            # Set up embedding parameters
            params = getattr(AllChem, self.method)()
            params.randomSeed = self.random_seed
            params.clearConfs = True
            params.numThreads = self.num_threads
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions = True
            params.enforceChirality = True
            
            # Generate conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=self.n_confs,
                params=params
            )
            
            if len(conf_ids) == 0:
                return None
            
            # Optimize conformers if requested
            if self.opt_confs:
                self._optimize_conformers(mol, conf_ids)
            
            # Calculate energies if requested
            if self.calc_energy:
                self._calculate_energies(mol, conf_ids)
            
            # Align conformers if requested
            if self.align_confs:
                self._align_conformers(mol, conf_ids)
            
            # Filter by RMS if requested
            if self.rms_cutoff is not None:
                mol = self._filter_by_rms(mol, self.rms_cutoff)
            
            return mol
            
        except Exception as e:
            print(f"Warning: Conformer generation failed: {e}")
            return None

    def _optimize_conformers(self, mol: Chem.Mol, conf_ids: List[int]):
        """Optimize conformers using the specified force field."""
        for conf_id in conf_ids:
            try:
                if self.ff == "UFF":
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=self.energy_iters)
                elif self.ff.startswith("MMFF"):
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=self.energy_iters)
            except Exception:
                continue  # Skip failed optimizations

    def _calculate_energies(self, mol: Chem.Mol, conf_ids: List[int]):
        """Calculate energies for conformers."""
        energies = []
        for conf_id in conf_ids:
            try:
                if self.ff == "UFF":
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                elif self.ff.startswith("MMFF"):
                    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
                
                if ff:
                    energy = ff.CalcEnergy()
                    energies.append((conf_id, energy))
                    mol.GetConformer(conf_id).SetDoubleProp("energy", energy)
            except Exception:
                continue
        
        # Set relative energies
        if energies:
            min_energy = min(energy for _, energy in energies)
            for conf_id, energy in energies:
                rel_energy = energy - min_energy
                mol.GetConformer(conf_id).SetDoubleProp("relative_energy", rel_energy)

    def _align_conformers(self, mol: Chem.Mol, conf_ids: List[int]):
        """Align conformers to the first conformer."""
        if len(conf_ids) > 1:
            try:
                AllChem.AlignMolConformers(mol)
            except Exception:
                pass  # Skip if alignment fails

    def _filter_by_rms(self, mol: Chem.Mol, rms_cutoff: float) -> Chem.Mol:
        """Filter conformers by RMS clustering."""
        try:
            # Get RMS matrix
            rms_matrix = AllChem.GetConformerRMSMatrix(mol)
            
            # Simple clustering: keep conformers that are different enough
            keep_confs = [0]  # Always keep the first conformer
            
            for i in range(1, mol.GetNumConformers()):
                keep = True
                for j in keep_confs:
                    if i < len(rms_matrix) and j < len(rms_matrix[i]):
                        if rms_matrix[i][j] < rms_cutoff:
                            keep = False
                            break
                if keep:
                    keep_confs.append(i)
            
            # Create new molecule with filtered conformers
            new_mol = Chem.Mol(mol)
            new_mol.RemoveAllConformers()
            
            for conf_id in keep_confs:
                conf = mol.GetConformer(conf_id)
                new_mol.AddConformer(conf, assignId=True)
            
            return new_mol
            
        except Exception:
            return mol  # Return original if filtering fails

    def save_conformers_to_sdf(self, mol: Chem.Mol, filename: str):
        """Save conformers to an SDF file."""
        if mol and mol.GetNumConformers() > 0:
            writer = Chem.SDWriter(filename)
            for conf in mol.GetConformers():
                writer.write(mol, confId=conf.GetId())
            writer.close()

    def create_temp_sdf_from_smiles(self, smiles_list: List[str], temp_dir: Optional[str] = None) -> str:
        """
        Create a temporary SDF file from SMILES with conformers.
        
        Args:
            smiles_list: List of SMILES strings
            temp_dir: Directory for temporary file
            
        Returns:
            Path to temporary SDF file
        """
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()
        
        Path(temp_dir).mkdir(parents=True, exist_ok=True)
        temp_file = os.path.join(temp_dir, "temp_conformers.sdf")
        
        writer = Chem.SDWriter(temp_file)
        
        for i, smiles in enumerate(smiles_list):
            if not smiles:
                continue
                
            mol = self.generate_conformers_from_smiles(smiles)
            if mol and mol.GetNumConformers() > 0:
                # Set molecule name
                mol.SetProp("_Name", f"mol_{i}")
                
                # Write all conformers
                for conf in mol.GetConformers():
                    writer.write(mol, confId=conf.GetId())
        
        writer.close()
        return temp_file
