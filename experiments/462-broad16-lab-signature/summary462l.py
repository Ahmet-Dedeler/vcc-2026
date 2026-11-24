"""462l summary: raw member deltas vs 435's shipped clean4 metrics (same seed), for Arc's programmes and (reference) Pan's."""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
F = "metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.json"
base = {s: json.load(open(ROOT / "experiments/435-zs-on-coarse/work/h1" / F.format(s=s)))["raw"] for s in (0, 1, 2)}
out = {}
for name, d, seeds in (("arcfull", "experiments/462-broad16-lab-signature/work/h1_arcfull", (0, 1, 2)),
                       ("arcbank", "experiments/462-broad16-lab-signature/work/h1_arcbank", (0,)),
                       ("pan (450)", "experiments/450-same-line-programme/work/h1_pan", (0, 1, 2))):
    rows = []
    for s in seeds:
        p = ROOT / d / F.format(s=s)
        if p.exists():
            r = json.load(open(p))["raw"]; rows.append({m: r[m] - base[s][m] for m in r})
    if rows:
        out[name] = {m: round(sum(x[m] for x in rows) / len(rows), 5) for m in rows[0]} | {"seeds": len(rows)}
        print(name, out[name])
(ROOT / "experiments/462-broad16-lab-signature/work/summary462l.json").write_text(json.dumps(out, indent=1))
