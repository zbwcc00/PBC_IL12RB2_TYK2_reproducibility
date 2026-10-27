from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "Fig01_study_overview.py",
    "Fig02_genetic_anchoring.py",
    "Fig03_IL12RB2_multisignal.py",
    "Fig04_pQTL_TYK2.py",
    "Fig05_single_cell_localization.py",
    "Fig06_independent_single_cell_validation.py",
    "Fig07_translation.py",
    "FigS01_genetic_screen_coloc_QC.py",
    "FigS02_PBC_local_LD_SuSiE_QC.py",
    "FigS03_competing_genes.py",
    "FigS04_FinnGen_UKB_QC.py",
    "FigS05_splicing_protein_conflict.py",
    "FigS06_S11.py",
    "FigS12_external_donor_pseudobulk.py",
    "FigS13_mechanistic_context.py",
]

for script in SCRIPTS:
    print(f"[run] {script}", flush=True)
    subprocess.run([sys.executable, str(HERE / script)], check=True)

print("All figure scripts completed.")
