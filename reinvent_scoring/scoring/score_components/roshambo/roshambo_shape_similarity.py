import os
import shutil
import time
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
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        self.param_enum = RoshamboSpecificParametersEnum()
        self.reference_file = self.parameters.specific_parameters.get(self.param_enum.REFERENCE_FILE, "")
        self.shape_weight = self.parameters.specific_parameters.get(self.param_enum.SHAPE_WEIGHT, 0.5)
        self.color_weight = self.parameters.specific_parameters.get(self.param_enum.COLOR_WEIGHT, 0.5)
        self.n_confs = self.parameters.specific_parameters.get(self.param_enum.N_CONFS, 10)
        self.ignore_hs = self.parameters.specific_parameters.get(self.param_enum.IGNORE_HS, True)
        self.use_carbon_radii = self.parameters.specific_parameters.get(self.param_enum.USE_CARBON_RADII, True)
        self.gpu_id = self.parameters.specific_parameters.get(self.param_enum.GPU_ID, 0)
        self.roshambo_api_url = self.parameters.specific_parameters.get("roshambo_api_url", "http://127.0.0.1:5000")
        self.save_overlays = self.parameters.specific_parameters.get(self.param_enum.SAVE_OVERLAYS, True)
        self.overlays_dir = os.path.abspath(
            self.parameters.specific_parameters.get(self.param_enum.OVERLAYS_DIR, "roshambo_overlays")
        )
        self.debug = self.parameters.specific_parameters.get("debug", False)

        if self.save_overlays:
            Path(self.overlays_dir).mkdir(parents=True, exist_ok=True)

        if not self.reference_file or not os.path.exists(self.reference_file):
            raise FileNotFoundError(f"Reference file not found: {self.reference_file}")

    def calculate_score(self, molecules: List, step=0) -> ComponentSummary:
        smiles_list = []
        for mol in molecules:
            try:
                from rdkit import Chem
                smiles = Chem.MolToSmiles(mol) if not isinstance(mol, str) else mol
                smiles_list.append(smiles)
            except:
                smiles_list.append("")

        scores = self._calculate_shape_scores(smiles_list, step)
        return ComponentSummary(total_score=scores, parameters=self.parameters)

    def _calculate_shape_scores(self, smiles_list: List[str], step: int) -> np.array:
        if not smiles_list:
            return np.array([], dtype=np.float32)

        reference_name = os.path.basename(self.reference_file)
        ref_copy_path = os.path.join(self.overlays_dir, reference_name)
        if step == 0 or not os.path.exists(ref_copy_path):
            shutil.copy2(self.reference_file, ref_copy_path)

        run_id = f"{step}" if step >= 0 else "final"
        dataset_file = self._create_dataset_smiles(smiles_list, self.overlays_dir, run_id)
        dataset_name = os.path.basename(dataset_file)

        scores = self._call_roshambo_api(reference_name, dataset_name, self.overlays_dir, len(smiles_list), run_id)
        return np.array(scores, dtype=np.float32)

    def _create_dataset_smiles(self, smiles_list: List[str], overlays_dir: str, run_id: str) -> str:
        file_path = os.path.join(overlays_dir, f"{run_id}dataset.smi")
        with open(file_path, "w") as f:
            for i, smi in enumerate(smiles_list):
                if smi.strip():
                    f.write(f"{smi.strip()} mol_{i}\n")
        return file_path

    def _call_roshambo_api(self, reference_file: str, dataset_file: str, working_dir: str, expected_count: int, run_id: str) -> List[float]:
        api_data = {
            "reference_file": reference_file,
            "dataset_file": dataset_file,
            "ignore_hs": self.ignore_hs,
            "n_confs": self.n_confs,
            "use_carbon_radii": self.use_carbon_radii,
            "color": self.color_weight > 0,
            "sort_by": "ComboTanimoto",
            "write_to_file": True,
            "gpu_id": self.gpu_id,
            "working_dir": working_dir,
            "output_prefix": f"{run_id}"
        }

        retries = 3
        delay_sec = 2
        csv_path = os.path.join(working_dir, f"{run_id}roshambo.csv")

        for attempt in range(retries):
            try:
                response = requests.post(f"{self.roshambo_api_url}/similarity", json=api_data, timeout=300)
                if response.status_code == 200 and response.json().get("success"):
                    if self.debug:
                        print(f"[Roshambo API] Success on attempt {attempt+1} for run_id={run_id}")
                    
                    # Wait for the CSV file to appear (up to 2s)
                    for _ in range(4):
                        if os.path.exists(csv_path):
                            return self._extract_scores_from_csv(csv_path, expected_count)
                        time.sleep(0.5)
                    if self.debug:
                        print(f"[Roshambo API] CSV not found after success response: {csv_path}")
            except Exception as e:
                if self.debug:
                    print(f"[Roshambo API Error] Attempt {attempt+1}/{retries}: {e}")
                time.sleep(delay_sec)

        if self.debug:
            print(f"[Roshambo API] Failed to compute similarity for run_id={run_id}")
        return [0.0] * expected_count

    def _extract_scores_from_csv(self, csv_file: str, expected_count: int) -> List[float]:
        try:
            df = pd.read_csv(csv_file, sep='\t', usecols=["OriginalName", "ShapeTanimoto", "ColorTanimoto"])
            df["score"] = df["ShapeTanimoto"] * self.shape_weight + df["ColorTanimoto"] * self.color_weight
            grouped = df.groupby("OriginalName")["score"].max()
            score_map = {
                int(name.split("_")[1]): score for name, score in grouped.items() if name.startswith("mol_")
            }
            return [score_map.get(i, 0.0) for i in range(expected_count)]
        except Exception as e:
            if self.debug:
                print(f"[Roshambo CSV Error] Failed to parse CSV {csv_file}: {e}")
            return [0.0] * expected_count
