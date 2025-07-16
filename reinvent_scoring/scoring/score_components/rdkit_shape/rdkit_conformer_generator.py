from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
import os

class RDKitConformerGenerator:
    """RDKit-based conformer generator to replace OpenEye's OMEGA."""

    def __init__(self, max_confs=200, energy_window=10, max_stereo=0, minimization_config=None):
        self.max_confs = max_confs
        self.energy_window = energy_window
        self.max_stereo = max_stereo
        self.minimization_config = minimization_config
        self.use_gpu = False

        # Check if OpenMM-enabled MMFF is available in this version of RDKit
        self.has_openmm_mmff = hasattr(AllChem, 'MMFFOptimizeMoleculeOpenMM')

        # If minimization config is provided, parse it
        if minimization_config and os.path.exists(minimization_config):
            self._parse_minimization_config()

    def generate_conformers(self, smiles):
        """Generate conformers for a molecule from SMILES."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Handle stereochemistry if needed
        if self.max_stereo > 0:
            isomers = self._enumerate_stereoisomers(mol)
            if not isomers:
                return None

            # Generate conformers for each stereoisomer
            all_confs_mol = None
            for isomer in isomers:
                confs = self._generate_confs_for_mol(isomer)
                if confs and confs.GetNumConformers() > 0:
                    if all_confs_mol is None:
                        all_confs_mol = confs
                    else:
                        # Add conformers from this isomer to the main molecule
                        for conf in confs.GetConformers():
                            all_confs_mol.AddConformer(conf, assignId=True)
            return all_confs_mol
        else:
            # Generate conformers for the molecule without stereochemistry enumeration
            return self._generate_confs_for_mol(mol)

    def _parse_minimization_config(self):
        """Parse the minimization configuration file."""
        try:
            with open(self.minimization_config, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue

                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip().lower()
                        value = value.strip()

                        if key == 'use_gpu' and value.lower() in ('true', 'yes', '1'):
                            self.use_gpu = True
                        elif key == 'force_field' and value.lower() in ('mmff', 'uff'):
                            self.force_field = value.lower()
                        elif key == 'max_iterations':
                            try:
                                self.max_iterations = int(value)
                            except ValueError:
                                pass
        except Exception as e:
            print(f"Warning: Could not parse minimization config file: {e}")

    def _generate_confs_for_mol(self, mol):
        """Generate conformers for a single molecule."""
        # Set up ETKDG parameters
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        params.enforceChirality = True

        # Generate conformers
        confs = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=self.max_confs,
            params=params
        )

        if confs == -1:  # Failed to generate any conformers
            return None

        # Energy minimize conformers
        try:
            # Standard UFF minimization (CPU-only)
            for conf_id in range(mol.GetNumConformers()):
                try:
                    # Try to use OpenMM if available and GPU is enabled
                    if self.has_openmm_mmff and self.use_gpu:
                        try:
                            AllChem.MMFFOptimizeMoleculeOpenMM(mol, confId=conf_id, useGPU=True)
                            continue  # Skip UFF if MMFF with OpenMM succeeded
                        except Exception as e:
                            print(f"Warning: OpenMM GPU minimization failed: {e}")
                            pass  # Fall back to UFF if MMFF fails

                    # Standard UFF minimization
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                except Exception as e:
                    print(f"Warning: Minimization failed for conformer {conf_id}: {e}")
                    # If minimization fails for this conformer, continue with the next one
                    continue
        except Exception as e:
            print(f"Warning: Minimization process failed: {e}")
            # If any error occurs during minimization, return the molecule with unminimized conformers
            pass

        return mol

    def _enumerate_stereoisomers(self, mol):
        """Enumerate stereoisomers of a molecule."""
        opts = EnumerateStereoisomers.StereoEnumerationOptions(
            maxIsomers=self.max_stereo,
            onlyUnassigned=False,
            unique=True
        )
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
        return isomers