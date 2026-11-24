"""481: normrel cap (-1.5) weights are relative to each donor's MEDIAN row norm over the panel it covers.
Per panel: each production donor's median raw row norm over its official-axis coverage (the builder's `median_norm`),
the between-donor median ratios (what sets the relative weights of two donors on one target), and the share of
(target, donor) rows that the cap down-weights by more than half (ratio > 2^(1/1.5) = 1.587).
usage: uv run python normrel481.py -> work/normrel481.json  (bank rows by mmap, < 1 GB)"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import proxy481 as P
HERE = Path(__file__).resolve().parent
B = P.banks()
out = {}
for name in ["abc", "h1_2025", "efr", "efp"]:
    panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
    rec = {}
    norms = {}
    for d, (delta, pidx, col) in B.items():
        cols = np.flatnonzero(col >= 0)
        v = np.array([float(np.linalg.norm(np.asarray(delta[pidx[t]], dtype=np.float64)[cols])) for t in panel if t in pidx])
        norms[d] = v
        rec[d] = {"n": int(v.size), "median": round(float(np.median(v)), 3),
                  "frac_downweighted_gt_half": round(float((v / np.median(v) > 2 ** (1 / 1.5)).mean()), 3)}
    rec["median_ratio_hct116_over_k562"] = round(rec["orion_hct116"]["median"] / rec["k562_gwps"]["median"], 3)
    rec["median_ratio_hek293t_over_k562"] = round(rec["orion_hek293t"]["median"] / rec["k562_gwps"]["median"], 3)
    out[name] = rec
    print(name, rec, flush=True)
(HERE / "work/normrel481.json").write_text(json.dumps(out, indent=1))
