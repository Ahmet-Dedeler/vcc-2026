"""467 summary: raw member deltas vs 435's shipped clean4 metrics (same seed) for the gated variants, with 462l's ungated
Arc programme as the reference; plus the gate's composition by n_real tercile."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
F = "metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.json"
base = {s: json.load(open(ROOT / "experiments/435-zs-on-coarse/work/h1" / F.format(s=s)))["raw"] for s in (0, 1, 2)}
spans = {"pds": 0.45, "mse": -0.955, "nmae": -0.60, "fid": 0.30, "reach": 0.90, "jac": 0.37}   # rough H1-like spans (sign: up is better)
out = {}
variants = [("ungated arcfull (462l)", ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull")] + \
           [(d.name[3:], d) for d in sorted((HERE / "work").glob("h1_*")) if d.name != "h1_record"]
for name, d in variants:
    rows = []
    for s in (0, 1, 2):
        p = d / F.format(s=s)
        if p.exists():
            r = json.load(open(p))["raw"]; rows.append({m: r[m] - base[s][m] for m in r})
    if not rows:
        continue
    mean = {m: float(np.mean([x[m] for x in rows])) for m in rows[0]}
    mean["overall_equiv"] = float(sum(mean[m] / spans[m] for m in spans) / 6)
    out[name] = {k: round(v, 5) for k, v in mean.items()} | {"seeds": len(rows)}
    print(name, out[name])
(HERE / "work/summary467.json").write_text(json.dumps(out, indent=1))
