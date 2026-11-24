"""481: does the H1 union's carrying-row dose mean the same thing on a new panel?

The builder sizes the H1 residual so that its median norm over the rows it CARRIES equals g x the whitened base's median
row norm over ALL rows. The share of the residual on a carrying row therefore depends on how the carrying targets' own
base norms compare with the panel median. Measured two ways on A/B/C (453 dump) and efr (465 dump):
  (a) carriers (the dump's design, H1 x2 included) median norm on H1-carrying rows / median over all rows;
  (b) the donor-side proxy (proxy481: agreeing 3-donor count |mean delta| > 0.2 on context D's gate) median on carrying
      rows / median over all rows -- a base-only measure that does not contain the residual.
usage: uv run python h1share481.py
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import proxy481 as P

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
H1 = set()
for f in (ROOT / "experiments/115-gse295214-h1/work/compact-donor/perts.txt", ROOT / "data/donors/h1_2025/perts.txt"):
    H1 |= {ln.strip() for ln in f.read_text().split("\n") if ln.strip()}
gate = P.gate_D()
B = P.banks()
out = {}
for name, dump in (("abc", ROOT / "experiments/453-abc-dress-rehearsal/work/run/dump_h1x2_cd4x1.npz"),
                   ("efr", ROOT / "experiments/465-both-switches-efr/work/run/dump_h1x2_cd4x1.npz"),
                   ("h1_2025", None)):
    panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
    carry = np.array([t in H1 for t in panel])
    rec = {"n_carry": int(carry.sum())}
    if dump is not None:
        z = np.load(dump, allow_pickle=False)
        car = np.asarray(z["carriers"], dtype=np.float64)
        nrm = np.linalg.norm(car, axis=1)
        rec["carrier_norm_carry_over_all"] = round(float(np.median(nrm[carry]) / np.median(nrm)), 3)
        s = (np.abs(np.where(gate[None, :], car, 0.0)) > 0.2).sum(1)
        rec["strength_median_carry_vs_all"] = [float(np.median(s[carry])), float(np.median(s))]
    px = P.proxy(panel, B, gate)[0.2]
    rec["proxy_median_carry_vs_all"] = [float(np.median(px[carry])), float(np.median(px))]
    rec["proxy_ratio_carry_over_all"] = round(float(np.median(px[carry]) / max(np.median(px), 1e-9)), 3)
    out[name] = rec
    print(name, rec, flush=True)
(HERE / "work/h1share481.json").write_text(json.dumps(out, indent=1))
