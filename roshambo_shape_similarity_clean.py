"""
Clean Roshambo Shape Similarity Component for Reinvent Scoring.
SMI format only - lets Roshambo handle all conformer generation.
"""

import os
import shutil
import pandas as pd
import numpy as np
import requests
from pathlib import Path
from typing import List

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.enums.roshambo_specific_parameters_enum import RoshamboSpecificParametersEnum


class RoshamboShapeSimilarity(BaseScoreComponent):
    """
    Clean Roshambo Shape Similarity Component - SMI Format Only.
    Creates SMILES files and lets Roshambo handle all conformer generation.
    """

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        # Initialize parameter enum
        self.param_enum = RoshamboSpecificParametersEnum()

        # Basic parameters
        self.reference_file = self.parameters.specific_parameters.get(self.param_enum.REFERENCE_FILE, "")
        self.shape_weight = self.parameters.specific_parameters.get(self.param_enum.SHAPE_WEIGHT, 0.5)
        self.color_weight = self.parameters.specific_parameters.get(self.param_enum.COLOR_WEIGHT, 0.5)
        self.n_confs = self.parameters.specific_parameters.get(self.param_enum.N_CONFS, 50)  # Default 50
        self.ignore_hs = self.parameters.specific_parameters.get(self.param_enum.IGNORE_HS, True)
        self.use_carbon_radii = self.parameters.specific_parameters.get(self.param_enum.USE_CARBON_RADII, True)
        self.gpu_id = self.parameters.specific_parameters.get(self.param_enum.GPU_ID, 0)

        # Flask API configuration
        self.roshambo_api_url = self.parameters.specific_parameters.get("roshambo_api_url", "http://127.0.0.1:5000")

        # Overlay saving parameters
        self.save_overlays = self.parameters.specific_parameters.get(self.param_enum.SAVE_OVERLAYS, True)
        self.overlays_dir = self.parameters.specific_parameters.get(self.param_enum.OVERLAYS_DIR, "roshambo_overlays")

        # Debug mode
        self.debug = self.parameters.specific_parameters.get("debug", False)

        # Create overlays directory
        if self.save_overlays:
            Path(self.overlays_dir).mkdir(parents=True, exist_ok=True)

        # Validate reference file
        if not self.reference_file:
            raise ValueError("Reference file must be provided")
        if not os.path.exists(self.reference_file):
            raise FileNotFoundError(f"Reference file not found: {self.reference_file}")

        if self.debug:
            print(f"🔧 Roshambo SMI-only initialized: reference={self.reference_file}, n_confs={self.n_confs}, gpu_id={self.gpu_id}")

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a list of molecules."""
        if self.debug:
            print(f"🧮 Roshambo calculate_score called with {len(molecules)} molecules, step={step}")

        # Convert input to SMILES strings
        smiles_list = []
        for mol in molecules:
            if isinstance(mol, str):
                smiles_list.append(mol)
            else:
                try:
                    from rdkit import Chem
                    smiles = Chem.MolToSmiles(mol)
                    smiles_list.append(smiles)
                except:
                    smiles_list.append("")  # Empty string will result in zero score

        # Calculate scores
        scores = self._calculate_shape_scores(smiles_list, step)

        if self.debug:
            mean_score = np.mean(scores) if len(scores) > 0 else 0.0
            print(f"📊 Calculated {len(scores)} scores, mean: {mean_score:.3f}")

        # Create and return component summary
        score_summary = ComponentSummary(total_score=scores, parameters=self.parameters)
        return score_summary

    def _calculate_shape_scores(self, smiles_list: List[str], step: int) -> np.array:
        """Calculate shape similarity scores using Roshambo Flask API with SMI format."""
        if self.debug:
            print(f"⚙️ Processing {len(smiles_list)} SMILES for step {step}")

        if not smiles_list:
            return np.array([], dtype=np.float32)

        # Create epoch folder
        epoch_folder = os.path.join(self.overlays_dir, f"epoch_{step}")
        Path(epoch_folder).mkdir(parents=True, exist_ok=True)

        # Create SMILES dataset file
        dataset_file = self._create_dataset_smiles(smiles_list, epoch_folder, step)
        dataset_filename = os.path.basename(dataset_file)

        # Copy reference file to epoch folder
        reference_filename = os.path.basename(self.reference_file)
        epoch_reference_file = os.path.join(epoch_folder, reference_filename)
        shutil.copy2(self.reference_file, epoch_reference_file)

        if self.debug:
            print(f"📄 Created SMI dataset: {dataset_filename}")
            print(f"📄 Copied reference: {reference_filename}")

        # Call Roshambo Flask API
        scores = self._call_roshambo_api(reference_filename, dataset_filename, epoch_folder, len(smiles_list))
        
        # Ensure correct number of scores
        if len(scores) != len(smiles_list):
            if self.debug:
                print(f"⚠️ Score count mismatch: got {len(scores)}, expected {len(smiles_list)}. Padding with zeros.")
            if len(scores) < len(smiles_list):
                scores.extend([0.0] * (len(smiles_list) - len(scores)))
            else:
                scores = scores[:len(smiles_list)]
        
        return np.array(scores, dtype=np.float32)

    def _create_dataset_smiles(self, smiles_list: List[str], epoch_folder: str, step: int) -> str:
        """Create SMILES file for Roshambo to handle conformer generation."""
        dataset_file = os.path.join(epoch_folder, f"dataset_{step}.smi")
        
        if self.debug:
            print(f"📝 Creating SMILES file with {len(smiles_list)} molecules")
            print(f"🔧 Roshambo will generate {self.n_confs} conformers per molecule")
        
        successful_molecules = 0
        with open(dataset_file, "w") as f:
            for i, smiles in enumerate(smiles_list):
                if smiles and smiles.strip():  # Skip empty SMILES
                    f.write(f"{smiles.strip()} mol_{i}\n")
                    successful_molecules += 1
        
        if self.debug:
            print(f"✅ Created SMILES file: {dataset_file} with {successful_molecules} molecules")
        
        return dataset_file

    def _call_roshambo_api(self, reference_file: str, dataset_file: str, epoch_folder: str, expected_count: int) -> List[float]:
        """Call Roshambo Flask API and return scores."""
        
        # Prepare API request data
        api_data = {
            "reference_file": reference_file,  # FILENAME ONLY
            "dataset_file": dataset_file,      # FILENAME ONLY
            "ignore_hs": self.ignore_hs,
            "n_confs": self.n_confs,
            "use_carbon_radii": self.use_carbon_radii,
            "color": self.color_weight > 0,
            "sort_by": "ComboTanimoto",
            "write_to_file": True,
            "gpu_id": self.gpu_id,
            "working_dir": epoch_folder
        }

        if self.debug:
            print(f"🌐 Calling Roshambo API at {self.roshambo_api_url}/similarity")
            print(f"📋 API data: {api_data}")

        try:
            response = requests.post(
                f"{self.roshambo_api_url}/similarity",
                json=api_data,
                timeout=300
            )

            if response.status_code == 200:
                result = response.json()
                if result.get("success"):
                    if self.debug:
                        print(f"✅ API call successful")
                    
                    # Read CSV file and extract scores
                    csv_file = os.path.join(epoch_folder, "roshambo.csv")
                    return self._extract_scores_from_csv(csv_file, expected_count)
                else:
                    if self.debug:
                        print(f"❌ API returned error: {result.get('error', 'Unknown error')}")
                    return [0.0] * expected_count
            else:
                if self.debug:
                    print(f"❌ API request failed with status {response.status_code}")
                return [0.0] * expected_count

        except Exception as e:
            if self.debug:
                print(f"❌ Error calling Roshambo API: {e}")
            return [0.0] * expected_count

    def _extract_scores_from_csv(self, csv_file: str, expected_count: int) -> List[float]:
        """Extract scores from Roshambo CSV output using ROCS-like combo score calculation."""
        if self.debug:
            print(f"📊 Extracting scores from CSV: {csv_file}")

        if not os.path.exists(csv_file):
            if self.debug:
                print(f"❌ CSV file not found: {csv_file}")
            return [0.0] * expected_count

        try:
            # Read CSV with tab delimiter
            df = pd.read_csv(csv_file, sep='\t')

            if self.debug:
                print(f"📈 CSV loaded: {df.shape[0]} rows")

            # Create a mapping from molecule names to scores
            mol_scores = {}
            for _, row in df.iterrows():
                name = row.get("Molecule", "")
                if name and name.startswith("mol_"):
                    try:
                        # Extract molecule index from name like "mol_54_0_3" -> 54
                        parts = name.split("_")
                        if len(parts) >= 2:
                            idx = int(parts[1])
                            
                            # Get individual scores
                            shape_score = float(row.get("ShapeTanimoto", 0.0))
                            color_score = float(row.get("ColorTanimoto", 0.0))
                            combo_score = float(row.get("ComboTanimoto", 0.0))

                            # Use ComboTanimoto directly or calculate ROCS-like weighted combination
                            if combo_score > 0:
                                final_score = combo_score
                            else:
                                # ROCS-like weighted combination
                                total_weight = self.shape_weight + self.color_weight
                                if total_weight > 0:
                                    final_score = (self.shape_weight * shape_score + 
                                                 self.color_weight * color_score) / total_weight
                                else:
                                    final_score = shape_score

                            # Keep the best score for each molecule (multiple conformers)
                            mol_scores[idx] = max(mol_scores.get(idx, 0.0), final_score)

                    except Exception as e:
                        if self.debug:
                            print(f"⚠️ Error processing row for {name}: {e}")
                        continue

            # Convert to ordered list ensuring we return exactly expected_count scores
            scores = [mol_scores.get(i, 0.0) for i in range(expected_count)]

            if self.debug:
                mean_score = np.mean(scores) if len(scores) > 0 else 0.0
                print(f"✅ Extracted {len(scores)} scores, mean: {mean_score:.3f}")

            return scores

        except Exception as e:
            if self.debug:
                print(f"❌ Error reading CSV file {csv_file}: {e}")
            return [0.0] * expected_count
