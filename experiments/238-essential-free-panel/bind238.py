"""238 — point 233's instrument at 238's panel registry, WITHOUT copying 233's code.

233's `score233.py` / `fold233.py` are the measurement.  236 reused them verbatim and that is
why its Phase 0 reproduced 233 at max |d| 0.0e+00 on 157/157 rows.  238 does the same, so the
essential-free numbers are produced by the identical code path as the numbers they are compared
against, and the only thing that differs is the panel.

The redirection is done by rebinding the module attributes `fold233` / `score233` read at CALL
time on their shared `panels233` module object:

  BUNDLES / read_panel / read_genes   -> 238's five-panel registry (abc, real, pess, efr, efp)
  op_path / cd4_root                  -> 233's artefacts for abc/real/pess, 238's for efr/efp
  WORK                                -> 238's, so NOTHING is written into 233's directory
  load_frozen_operator_233            -> 238's identical axis-asserting loader

`score233.SHARED` is computed from `panels233.WORK` at IMPORT time, so the rebinding must
happen BEFORE `score233` is imported.  Importing this module first is what guarantees that;
every 238 script does `import bind238` as its first project import.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R233 = ROOT / "experiments/233-october-retention"

for p in (str(HERE), str(R233),
          str(ROOT / "experiments/162-retrodicting-validator"),
          str(ROOT / "experiments/188-magnitude-calibration"),
          str(ROOT / "experiments/203-reopen-125"),
          str(ROOT / "experiments/218-independent-donors"),
          str(ROOT / "experiments/228-theta-gain-plolo"),
          str(ROOT / "experiments/236-october-stack")):
    if p not in sys.path:
        sys.path.insert(0, p)

import panels238 as P238  # noqa: E402
import panels233 as P233  # noqa: E402

if "score233" in sys.modules or "fold233" in sys.modules:
    raise SystemExit("bind238 must be imported BEFORE score233/fold233")

P233.BUNDLES = P238.BUNDLES
P233.WORK = P238.WORK
P233.read_panel = P238.read_panel
P233.read_genes = P238.read_genes
P233.axis_sha256 = P238.axis_sha256
P233.op_path = P238.op_path
P233.cd4_root = P238.cd4_root
P233.bind = P238.bind
P233.load_frozen_operator_233 = P238.load_frozen_operator_238
# 238 keeps PLOLO's fourth fold (`kolf2`) addressable; 233 dropped it.  `CANCER_FOLDS` is
# unchanged, so nothing that defaults to it picks kolf2 up.
P233.FOLD_SOURCES = P238.FOLD_SOURCES

PANELS = P238.BUNDLES
CANCER_FOLDS = P238.CANCER_FOLDS
