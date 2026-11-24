"""426 table: every metrics-*_mde*.json in work/h1 vs the base of the same seed (426's own base run if present, else
416's bit-identical base_sS_bps), members at board slopes."""
import glob
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
SL = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}


def base(seed):
    for p in (HERE / f"work/h1/metrics-base_s{seed}_bps.json",
              HERE.parent / f"416-de-carrier-menu/work/h1/metrics-base_s{seed}_bps.json"):
        if p.exists():
            return json.loads(p.read_text())
    raise FileNotFoundError(seed)


rows = []
for f in sorted(glob.glob(str(HERE / "work/h1/metrics-base_s*_bps_mde*.json"))):
    m = json.loads(Path(f).read_text())
    b = base(m["seed"])
    d = {k: m["raw"][k] - b["raw"][k] for k in SL}
    tot = sum(SL[k] * d[k] / 6 for k in SL)
    nonreach = tot - SL["reach"] * d["reach"] / 6
    rows.append((os.path.basename(f)[8:-5], d, tot, nonreach, m))
for name, d, tot, nr, m in rows:
    print(f"{name:42s} " + " ".join(f"{k} {d[k]:+.5f}" for k in ("nmae", "fid", "jac", "reach", "pds", "mse"))
          + f" | ovr {tot:+.5f} non-reach {nr:+.5f} | unsupp {m.get('unsuppressible_mean', float('nan')):.1f}"
          f" clip {m.get('clip_frac_mean', float('nan')):.3f} mde_lam {m.get('mde_lambda_median')}")
