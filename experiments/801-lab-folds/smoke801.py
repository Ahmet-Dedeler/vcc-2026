"""801 -- the light pricing smoke: on every fold, the incumbent (pct+P2g) against `cap3` and against `pct+S_P2`, with
engine801's per-target arrays saved as work/smoke/<fold>-<arm>.npz (arm priced with the incumbent as its base, so the
change-energy arrays are arm - incumbent) and the fold-level r (relative aggregate cosine, incumbent vs the arm) and
literal / board-emulated Delta pds. No aggregation into board units or bins (the lab lead's check800 does that).
Every fold is loaded with its known answers (engine801.load(check=True)); a failed KA aborts.

    uv run python experiments/801-lab-folds/smoke801.py <fold_id>...       -> work/smoke/<fold>-<arm>.npz, smoke801-<tag>.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import engine801 as E  # noqa: E402

ARMS = ("cap3", "pct+S_P2")
OUT = HERE / "work" / "smoke"


def main(fold_ids, tag):
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    for fid in fold_ids:
        t0 = time.time()
        F = E.load(fid, check=True)
        inc = E.named(F, "incumbent")
        rec = dict(n=F.n, n_carried=int(F.carried.sum()), ka_worst=max(v for k, v in F.ka.items()
                                                                      if not k.startswith("_") and isinstance(v, float)))
        for nm in ARMS:
            arm = E.named(F, nm)
            s, arr = E.price(F, arm, base=inc)
            np.savez_compressed(OUT / f"{fid.replace('/', '-')}-{nm}.npz", targets=np.asarray(F.targets), w=arm.w,
                                q=np.broadcast_to(np.asarray(arm.q, float), (F.n,)), w_base=inc.w, q_base=inc.q, **arr)
            vb = s["vs_base"]
            rec[nm] = dict(cos_agg=s["cos_agg"], pds=s["pds"], pds_lam05=s["pds_lam0.5"], K_mean=s["K_mean"],
                           # incumbent vs the arm (positive = the incumbent is better)
                           r_incumbent_vs_arm=1.0 / (1.0 + vb["r"]) - 1.0, dpds_incumbent_minus_arm=-vb["dpds"],
                           dpds_lam05_incumbent_minus_arm=-vb["dpds_lam05"], dK_mean_arm_minus_incumbent=vb["dK_mean"],
                           share_change_unscored=vb["share_unscored"], share_change_cpm=vb["share_cpm"])
        s0 = E.score_rows_of(F, *E.render(F, inc))[0]
        rec["incumbent"] = dict(cos_agg=s0["cos_agg"], pds=s0["pds"], pds_lam05=s0["pds_lam0.5"], K_mean=s0["K_mean"],
                                hi=float(np.max(inc.w)))
        rec["seconds"] = round(time.time() - t0, 1)
        res[fid] = rec
        print(f"[{fid}] n {F.n}; KA worst {rec['ka_worst']:.1e}; incumbent cos {s0['cos_agg']:.6f} pds {s0['pds']:.6f}; "
              + "; ".join(f"vs {nm}: r {100 * rec[nm]['r_incumbent_vs_arm']:+.3f} % dpds {rec[nm]['dpds_incumbent_minus_arm']:+.5f}"
                          for nm in ARMS) + f" ({rec['seconds']} s)", flush=True)
        (OUT / f"smoke801-{tag}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    tag = sys.argv[1]
    main(sys.argv[2:], tag)
