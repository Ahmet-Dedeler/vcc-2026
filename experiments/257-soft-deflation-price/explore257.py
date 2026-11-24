"""257 EXPLORATORY (post hoc, not part of the fixed grid, not used for the verdict): the k = 5
projector with one direction LEFT IN (`k5_drop{j}`), priced over cd4svd5 with price257's recipe.

    uv run python explore257.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import price257 as R

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"


def main() -> None:
    rows, pt, idn = {}, {}, {}
    for held in R.FOLDS:
        for tag in ("sweep257", "sweep257b"):
            d = json.loads((WORK / f"{tag}-{held}.json").read_text())
            for r in d["rows"]:
                rows[(held, r["arm"])] = r
            z = np.load(WORK / f"pertarget{tag[5:]}-{held}.npz")
            for k in z.files:
                pt[(held, k)] = z[k]
            z = np.load(WORK / f"identity{tag[5:]}-{held}.npz")
            for k in z.files:
                idn[(held, k)] = z[k]
    ka = [(h, r["arm"], r["ka"]["ok"]) for (h, a), r in rows.items() if "ka" in r and a in ("STACK", "k5_s1")]
    print(f"257b known-answer rows: {sum(ok for _, _, ok in ka)}/{len(ka)} PASS")
    out = {}
    print("| cell | left in | E removed | CD4 fraction × | Δ`pds` vs `cd4svd5` | folds + | σ | rel `cos_agg` | rel macro | delivered | Δ`pds` vs STACK | agg vs STACK | **overall point** | **overall conservative** | k_cos 2.15 | `kolf2` Δ`pds` / agg |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for j in range(1, 6):
        cell = f"k5_drop{j}"
        c = R.compare(rows, idn, cell, "k5_s1", R.CANCER)
        cs = R.compare(rows, idn, cell, "STACK", R.CANCER)
        ck = R.compare(rows, idn, cell, "k5_s1", ("kolf2",))
        pr = R.price_both(c, R.COS_CD4SVD5, R.PER_CD4SVD5)
        ret = R.retention(pt, cell, "k5_s1", R.CANCER)
        out[cell] = dict(vs_cd4svd5=c, vs_stack=cs, kolf2=ck, price=pr, retention_vs_cd4svd5=ret,
                         shape_only=abs(c["delivered_ratio"] - 1) <= R.AMP_SHAPE_BOUND)
        print(f"| `{cell}` | direction {j} | {c['energy_removed']:.3f} | {c['fraction_ratio']:.4f} | **{c['d_pds']:+.6f}** | {c['pds_pos']}/3 | {c['sigma']:+.1f} | "
              f"{100*c['rel_agg']:+.3f} % | {100*c['rel_macro']:+.3f} % | {100*(c['delivered_ratio']-1):+.3f} % | {cs['d_pds']:+.6f} | {100*cs['rel_agg']:+.3f} % | "
              f"**{pr['point']['total']:+.5f}** | **{pr['conservative']['total']:+.5f}** | {pr['k_cos_2_15']['total']:+.5f} | {ck['d_pds']:+.5f} / {100*ck['rel_agg']:+.2f} % |")
    print()
    for j in range(1, 6):
        c = out[f"k5_drop{j}"]["vs_cd4svd5"]
        print(f"k5_drop{j}: per-fold d pds {['%+.5f' % x for x in c['per_fold_d_pds']]} agg {['%+.3f%%' % (100*x) for x in c['per_fold_rel_agg']]} "
              f"d nmae89 {c['d_nmae']:+.5f} d prec89 {c['d_prec89']:+.5f} d q_w {c['d_qw']:+.5f}; legs point pds {out[f'k5_drop{j}']['price']['point']['pds']:+.5f} "
              f"mse {out[f'k5_drop{j}']['price']['point']['mse']:+.5f} DE3 {out[f'k5_drop{j}']['price']['point']['de3']:+.5f} nmae {out[f'k5_drop{j}']['price']['point']['nmae']:+.5f}")
    (WORK / "explore257.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote {WORK/'explore257.json'}")


if __name__ == "__main__":
    main()
