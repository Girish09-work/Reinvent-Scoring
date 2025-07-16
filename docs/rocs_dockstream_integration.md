# ROCS Similarity vs DockStream – How They Interact

This note explains how the **ROCS**‐based 3-D similarity scorer and **DockStream** docking scorer are wired inside the *REINVENT-scoring* package.
It clarifies where files such as `protac.sdf` are produced and why they are **not** consumed by DockStream.

---
## 1. Data flow inside a REINVENT scoring cycle
```
┌────────────┐   SMILES list      ┌────────────────┐  SMILES list    ┌────────────────┐
│ REINVENT   │ ───────────────►  │ ROCS component │ ───────────────► │ DockStream     │
│ agent      │                   │ (shape / color │                  │ component      │
│ (SMILES)   │ ◄─────────────────│  similarity)   │◄────────────────│ (docking)      │
└────────────┘    total score    └────────────────┘   total score    └────────────────┘
```
* Both components receive the **same** list of SMILES directly from the agent.
* No intermediate SDF file is sent from ROCS to DockStream.

---
## 2. ROCS component
Key excerpt from `rocs_similarity.py` (single-CPU variant):
```python
scores = []
for smile in smiles:
    imol = oechem.OEMol()
    if oechem.OESmilesToMol(imol, smile):
        if self.omega(imol):              # <- OpenEye OMEGA 3-D generation
            self.prep.Prep(imol)
            score = oeshape.OEBestOverlayScore()
            self.overlay.BestOverlay(score, imol, predicate)
            best_score_shape = getattr(score, self.sim_func_name_set.shape)()
            best_score_color = getattr(score, self.sim_func_name_set.color)()
            best_score = (self.shape_weight * best_score_shape +
                          self.color_weight * best_score_color) / (self.shape_weight + self.color_weight)
    scores.append(best_score)
```
Highlights:
1. **OMEGA** builds 3-D conformers on-the-fly → *not* written to disk.
2. If `SAVE_ROCS_OVERLAYS` is enabled (see `parallel_rocs_similarity.py`) every best-overlay conformer is **optionally** written to an SDF file, e.g. `<prefix><step>.sdf`.
   ```python
   overlay_filename = self.overlay_prefix + ind + ".sdf"
   outfs = oechem.oemolostream(overlay_file_path)
   oechem.OEWriteMolecule(outfs, outmol)        # Saved for visual analysis only
   ```
3. That SDF (your *protac.sdf*) is merely a **side artefact** for human inspection / debugging.
   The ROCS component already returns the numeric similarity score to REINVENT; nothing else is exported.

---
## 3. DockStream component
Snippet from `structural/dockstream.py`:
```python
concat_smiles = '"' + ';'.join(smiles) + '"'
command = ' '.join([
    self._environment_path,            # e.g. "bash /path/to/env.sh"
    self._docker_script_path,          # dockstream/docker.py
    "-conf", self._configuration_path,
    "-output_prefix", self._get_step_string(step),
    "-smiles", concat_smiles,         # <─ the ligand list (no SDF)
    "-print_scores"
])
```
Observations:
1. DockStream is executed as a **child process** via `subprocess.Popen` (see `BaseStructuralComponent._send_request_with_stepwize_read`).
2. Ligands are forwarded as an *in-line* SMILES string. DockStream then performs its **own** ligand embedding/preparation with the backend chosen in the JSON config (AutoDock-Vina, Glide, GOLD, …).
3. No SDF produced by ROCS is referenced – DockStream starts from scratch or from whatever embedding tool the backend requires.

---
## 4. Why `protac.sdf` is **not** passed to DockStream
| Fact | Explanation |
|------|-------------|
| ROCS writes SDF only if `SAVE_ROCS_OVERLAYS = true` | This is meant for later visual inspection of overlay poses. |
| DockStream interface receives ligands via `-smiles` | It rebuilds 3-D coordinates according to the chosen backend; this ensures consistent protonation, tautomer & charge states. |
| Therefore | The SDF on disk is **ignored** by DockStream. |

If you **really** want DockStream to start from the ROCS-generated 3-D coordinates you would have to modify its configuration to accept an input SDF and change `DockStream._create_command` accordingly. The current integration keeps the two scorers independent.

---
## 5. Summary
1. **No file sharing** occurs between the two scorers in the default pipeline.
2. `protac.sdf` is an *optional* overlay snapshot created by ROCS for inspection, not for DockStream.
3. DockStream wraps various docking backends and generates its own 3-D poses from the SMILES every time.

---
### Further reading
* `scoring/score_components/rocs/parallel_rocs_similarity.py`
* `scoring/score_components/structural/dockstream.py`
* DockStream project docs: <https://github.com/MolecularAI/DockStream>
