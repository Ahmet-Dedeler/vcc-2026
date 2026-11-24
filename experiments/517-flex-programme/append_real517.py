"""517: write (or replace) RESULTS.md §9 from the real-emitter runs (work/real517.json from summ_real517.py, the KA lines in
work/drive517.log, and the zs-rider runs against 435's zs runs), beside price517 / mzgrid517's analytic numbers. Run by
drive517.sh at its end, and by hand."""
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R = HERE / "RESULTS.md"
real = json.loads((HERE / "work/real517.json").read_text()) if (HERE / "work/real517.json").exists() else {}
log = (HERE / "work/drive517.log").read_text() if (HERE / "work/drive517.log").exists() else ""
ka = [l.split(" ", 1)[1] for l in log.splitlines() if "KA PASS" in l or "KA FAIL" in l]
peaks = re.findall(r"peakRSS=([0-9.]+) GB", log)
ANALYTIC = {  # H1 analytic means (3 seeds) for the same arms, from price517 / mzg33_517 / hookgrid517 / grid517 / variantblock517
    "c4-mzg0.05": "nmae −0.0183, reach 0", "c4-arcg0.05": "nmae −0.0148, reach −0.0031", "g33-mzg0.05": "nmae −0.0182, reach 0",
    "g33-arcg0.05": "nmae −0.0157, reach −0.0031", "c4-arct0.05": "nmae −0.0167, reach −0.0175", "c4-arcg0.1": "nmae −0.0106, reach 0",
    "c4-arcg0.05-r0.25": "nmae ≈ −0.0059 (interp.)", "g33-arcg0.1": "nmae −0.0112, reach 0", "c4-c4t0.05": "nmae −0.00146 (= 496's real run)",
    "mz-none": "fid +0.0178, nmae −0.0160, reach +0.0171 (vs CLEAN4)", "mzg-none": "fid +0.0146, nmae −0.0124, reach +0.0213",
    "arc-none-q0.7": "vs 458d's construction: fid −0.0040, nmae +0.0036, reach +0.0308"}
lines = ["## 9. Real emitter runs on H1 (auto-written by `append_real517.py` from `drive517.sh`)", "",
         f"Known answers: {'; '.join(ka) if ka else 'not run yet'}. Peak RSS per run: "
         f"{(min(map(float, peaks)) if peaks else float('nan')):.1f}–{(max(map(float, peaks)) if peaks else float('nan')):.1f} GB.", "",
         "Δ raw vs the no-hook base run of the same seed (c4 = 435, g33 = 467's gated33, arc = 462l; mz/mzg block swaps vs 435).",
         "unsupp = mean `call_set_unsuppressible` per target (base → arm): the leak the analytic path cannot see.", "",
         "| mode | seeds | pds / mse | nmae | fid | reach | jac | overall-equiv (seed range) | unsupp base → arm | analytic (H1) |",
         "|---|---:|---|---:|---:|---:|---:|---|---|---|"]
for mode, rows in sorted(real.items()):
    m = lambda k: float(np.mean([r[k] for r in rows]))     # noqa: E731
    ovs = [r["ov"] for r in rows]
    pm = max(max(abs(r["pds"]), abs(r["mse"])) for r in rows)
    ub = np.mean([r["unsupp"][0] for r in rows]); ua = np.mean([r["unsupp"][1] for r in rows])
    lines.append(f"| {mode} | {len(rows)} | {'exact' if pm == 0 else f'max {pm:.1e}'} | {m('nmae'):+.5f} | {m('fid'):+.5f} | "
                 f"{m('reach'):+.5f} | {m('jac'):+.5f} | {np.mean(ovs):+.5f} ({min(ovs):+.5f} … {max(ovs):+.5f}) | "
                 f"{ub:.0f} → {ua:.0f} | {ANALYTIC.get(mode, '—')} |")
# the zs-rider additivity runs: hook + zs vs 435's coarse + zs, same seed
zs_rows = []
for d in sorted(HERE.glob("work/h1_*_zs")):
    t = "base_s0_bps_hookem_mde0.1cap0.1min0.3coarse_dnk250zs"
    b = ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"
    f = d / f"metrics-{t}.json"
    if f.exists() and b.exists():
        a_, b_ = json.loads(f.read_text())["raw"], json.loads(b.read_text())["raw"]
        zs_rows.append(f"| {d.name[3:]} (seed 0) | " + " | ".join(f"{a_[k] - b_[k]:+.5f}" for k in ("nmae", "fid", "reach", "jac", "pds", "mse")) + " |")
if zs_rows:
    lines += ["", "With the incumbent's zs rider (hook + zs vs 435's coarse + zs, seed 0; additive if nmae ≈ the plain run's):", "",
              "| mode | nmae | fid | reach | jac | pds | mse |", "|---|---:|---:|---:|---:|---:|---:|", *zs_rows]
sec = "\n".join(lines) + "\n"
txt = R.read_text()
if "## 9. Real emitter runs" in txt:
    txt = txt[:txt.index("## 9. Real emitter runs")] + sec
else:
    txt = txt.rstrip("\n") + "\n\n" + sec
R.write_text(txt)
print(sec)
