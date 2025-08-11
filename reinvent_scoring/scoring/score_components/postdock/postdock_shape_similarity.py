import os
import re
import time
import shutil
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import requests

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary


class PostDockShapeSimilarity(BaseScoreComponent):
    """Compute shape similarity scores post docking via Roshambo on extracted warheads.

    Flow per epoch:
      1) Ensure overlays dir exists; copy reference warheads SDF there on first epoch only
      2) Build docked pose filename: e{epoch:04d}_ADV_ligands_docked.sdf under docked_poses_path
      3) Extract warheads from docked poses into overlays dir as e{epoch:04d}_e_warheads.sdf
      4) Call Roshambo API with reference_file (filename only) and dataset_file (filename only)
      5) Parse roshambo CSV and map scores back by OriginalName (expects names containing 'mol_{i}')
    """

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        sp = self.parameters.specific_parameters or {}
        tf = (sp.get("transformation") or {})

        # Inputs
        # Accept both 'reference_file' and legacy 'roshambo_input'
        self.reference_file: str = sp.get("reference_file") or sp.get("roshambo_input", "")
        self.docked_poses_path: str = sp.get("docked_poses_path", "")

        # Roshambo/Scoring params
        self.shape_weight: float = float(sp.get("shape_weight", 1.0))
        self.color_weight: float = float(sp.get("color_weight", 0.0))
        self.similarity_measure: str = sp.get("similarity_measure", "Tanimoto")
        self.n_confs: int = int(sp.get("n_confs", 0))  # for SDF datasets, 0 (use existing conformers)
        self.ignore_hs: bool = bool(sp.get("ignore_hs", True))
        self.use_carbon_radii: bool = bool(sp.get("use_carbon_radii", True))
        self.gpu_id: int = int(sp.get("gpu_id", 0))
        self.roshambo_api_url: str = sp.get("roshambo_api_url", "http://127.0.0.1:5000")

        # Overlays/output management
        self.save_overlays: bool = bool(sp.get("save_overlays", tf.get("save_overlays", True)))
        overlays_dir_raw: str = sp.get("overlays_dir", tf.get("overlays_dir", "roshambo_overlays"))
        self.overlays_dir: str = os.path.abspath(overlays_dir_raw)

        # Misc/debug
        self.debug: bool = bool(sp.get("debug", False))

        # Basic validations
        if not self.reference_file or not os.path.exists(self.reference_file):
            raise FileNotFoundError(f"Reference warheads SDF not found: {self.reference_file}")
        if not self.docked_poses_path or not os.path.isdir(self.docked_poses_path):
            raise FileNotFoundError(f"Docked poses path not found or not a directory: {self.docked_poses_path}")

        if self.save_overlays:
            Path(self.overlays_dir).mkdir(parents=True, exist_ok=True)

        # Copy reference once (idempotent)
        self._ref_basename = os.path.basename(self.reference_file)
        self._ref_copy_path = os.path.join(self.overlays_dir, self._ref_basename)
        if not os.path.exists(self._ref_copy_path):
            shutil.copy2(self.reference_file, self._ref_copy_path)

    def calculate_score(self, molecules: List, step: int = 0) -> ComponentSummary:
        # Determine expected count for mapping back to molecules
        expected_count = len(molecules)

        epoch = self._resolve_epoch_index(step)
        run_prefix = f"e{epoch:04d}_"

        # Locate docked pose file for this epoch
        docked_pose = self._get_docked_pose_file(epoch)
        if not docked_pose:
            if self.debug:
                print(f"[PostDock] Docked pose file not found for epoch {epoch} in {self.docked_poses_path}")
            scores = np.array([0.0] * expected_count, dtype=np.float32)
            return ComponentSummary(total_score=scores, parameters=self.parameters)

        # Extract warheads from docked poses -> dataset SDF
        dataset_sdf = os.path.join(self.overlays_dir, f"{run_prefix}e_warheads.sdf")
        try:
            self._extract_warheads(docked_pose, self._ref_copy_path, dataset_sdf)
        except Exception as e:
            if self.debug:
                print(f"[PostDock] Warhead extraction failed for {docked_pose}: {e}")
            scores = np.array([0.0] * expected_count, dtype=np.float32)
            return ComponentSummary(total_score=scores, parameters=self.parameters)

        # Call Roshambo API
        ref_name = self._ref_basename
        dataset_name = os.path.basename(dataset_sdf)
        scores_list = self._call_roshambo_api(ref_name, dataset_name, self.overlays_dir, run_prefix, expected_count)
        scores = np.array(scores_list, dtype=np.float32)
        return ComponentSummary(total_score=scores, parameters=self.parameters)

    # ------------------------ internals ------------------------ #

    def _resolve_epoch_index(self, step: int) -> int:
        # Epochs are 0-based; DockStream files also start at e0000
        if step is None or step < 0:
            latest = self._find_latest_epoch()
            return latest if latest is not None else 0
        return step

    def _find_latest_epoch(self) -> Optional[int]:
        pattern = re.compile(r"^e(\d{4})_ADV_ligands_docked\.sdf$")
        epochs = []
        try:
            for fn in os.listdir(self.docked_poses_path):
                m = pattern.match(fn)
                if m:
                    epochs.append(int(m.group(1)))
        except Exception:
            return None
        return max(epochs) if epochs else None

    def _get_docked_pose_file(self, epoch: int) -> Optional[str]:
        fname = f"e{epoch:04d}_ADV_ligands_docked.sdf"
        path = os.path.join(self.docked_poses_path, fname)
        return path if os.path.exists(path) else None

    def _extract_warheads(self, protac_sdf: str, warheads_sdf: str, output_sdf: str) -> None:
        """Use the repository's robust warhead extractor to write output_sdf."""
        from reinvent_scoring.scoring.score_components.postdock import extraction_warheads as ew

        n = ew.extract_warheads_to_file(
            warheads_sdf=warheads_sdf,
            protac_sdf=protac_sdf,
            output_sdf=output_sdf,
            mcs_fallback=True,
            mcs_threshold=0.8,
        )
        if n <= 0:
            raise RuntimeError("No warheads extracted from docked poses; check references and docking output.")

    def _call_roshambo_api(self, reference_file: str, dataset_file: str, working_dir: str,
                           run_prefix: str, expected_count: int) -> List[float]:
        api_data = {
            "reference_file": reference_file,           # filename only
            "dataset_file": dataset_file,               # filename only (SDF)
            "ignore_hs": self.ignore_hs,
            "n_confs": self.n_confs,
            "use_carbon_radii": self.use_carbon_radii,
            "color": self.color_weight > 0,
            "sort_by": "ComboTanimoto",
            "write_to_file": True,
            "gpu_id": self.gpu_id,
            "working_dir": working_dir,
            "output_prefix": run_prefix
        }

        csv_path = os.path.join(working_dir, f"{run_prefix}roshambo.csv")

        retries = 3
        delay_sec = 2
        for attempt in range(1, retries + 1):
            try:
                resp = requests.post(f"{self.roshambo_api_url}/similarity", json=api_data, timeout=300)
                if resp.status_code == 200 and resp.json().get("success"):
                    # Wait briefly for CSV
                    for _ in range(10):
                        if os.path.exists(csv_path):
                            return self._extract_scores_from_csv(csv_path, expected_count)
                        time.sleep(0.3)
                if self.debug:
                    print(f"[PostDock] Roshambo API attempt {attempt} returned status {resp.status_code}")
            except Exception as e:
                if self.debug:
                    print(f"[PostDock] Roshambo API error on attempt {attempt}: {e}")
            time.sleep(delay_sec)

        if self.debug:
            print(f"[PostDock] Roshambo API failed; returning zeros. Expected {expected_count} scores.")
        return [0.0] * expected_count

    def _extract_scores_from_csv(self, csv_file: str, expected_count: int) -> List[float]:
        """Parse Roshambo CSV and map best score per molecule index (mol_i)."""
        try:
            df = pd.read_csv(csv_file, sep='\t')
            # Prefer explicit columns, fallback if missing
            if {"OriginalName", "ShapeTanimoto", "ColorTanimoto"}.issubset(df.columns):
                df["score"] = df["ShapeTanimoto"] * self.shape_weight + df["ColorTanimoto"] * self.color_weight
            elif {"OriginalName", "ComboTanimoto"}.issubset(df.columns):
                # Use provided combo if shape/color split not present
                df["score"] = df["ComboTanimoto"].astype(float)
            else:
                if self.debug:
                    print(f"[PostDock] Unexpected CSV columns: {list(df.columns)}")
                return [0.0] * expected_count

            # Extract index from names like 'mol_12_...' or 'mol_12'
            def to_index(name: str) -> Optional[int]:
                if isinstance(name, str) and name.startswith("mol_"):
                    parts = name.split("_")
                    if len(parts) >= 2 and parts[1].isdigit():
                        return int(parts[1])
                return None

            df["mol_index"] = df["OriginalName"].map(to_index)
            grouped = df.dropna(subset=["mol_index"]).groupby("mol_index")["score"].max()
            score_map = {int(k): float(v) for k, v in grouped.to_dict().items()}
            return [score_map.get(i, 0.0) for i in range(expected_count)]
        except Exception as e:
            if self.debug:
                print(f"[PostDock] Failed to parse CSV {csv_file}: {e}")
            return [0.0] * expected_count

