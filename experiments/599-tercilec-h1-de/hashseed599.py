"""Score 435's saved base_s0 pre-zs stats with evaluate323.score in this process; print the raw members at full precision."""
import json, os, sys
from pathlib import Path
import numpy as np
ROOT = Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
import run400 as H
E = H.E
r, c, emit_c = E.load_reference()
z = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz", allow_pickle=True)
out = E.score(r, c, z["sums"], z["jk"], z["p"], z["q"], z["lfc"])
raw = out if isinstance(out, dict) else out[0]
print(json.dumps({"PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"), "raw": {k: repr(v) for k, v in (raw.items() if isinstance(raw, dict) else [])}}))
