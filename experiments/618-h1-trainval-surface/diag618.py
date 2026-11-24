#!/usr/bin/env python3
"""618: the 614 low-power diagnostic on the h1-trainval truth, reported BEFORE the gate runs.
For each harness reference: n_real per target (BH q < 0.05 on the DE universe, the target's own gene excluded; ref614's
count), its quantiles, the targets nmae scores (n_real >= 10, 614's reading), and the fraction of DE calls that are
negative (h1-test ~0.50; the old greedy Orion harness 0.908). Rows: h1-trainval (618), H1 h1-test (614 ref_h1ka), K562 (614 ref_k562s)."""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
P614 = HERE.parents[0] / "614-k562-real-emitter/work"
out = {}
for name, d in (("tv618", HERE / "work/ref_tv"), ("k562_614", P614 / "ref_k562s"), ("h1_614", P614 / "ref_h1ka")):
    js = json.loads((d / "reference.json").read_text())
    v = np.array(list(js["n_real_by_target"].values()), float)
    z = np.load(d / "reference.npz", allow_pickle=False)
    q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
    sig = np.isfinite(q) & (q < 0.05)
    out[name] = {"targets": int(len(v)), "n_real_median": float(np.median(v)), "n_real_q90": float(np.percentile(v, 90)),
                 "n_real_max": int(v.max()), "targets_ge10": int((v >= 10).sum()),
                 "frac_targets_ge10": float((v >= 10).mean()), "de_calls_total": int(sig.sum()),
                 "frac_negative": float((lfc[sig] < 0).mean()) if sig.any() else None,
                 "library_median_ctrl": float(np.median(np.asarray(z["libraries"], float)))}
    print(name, json.dumps(out[name]))
(HERE / "work/diag618.json").write_text(json.dumps(out, indent=1))
