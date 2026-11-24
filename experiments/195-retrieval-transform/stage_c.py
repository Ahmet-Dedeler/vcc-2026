"""195 Stage C — the permuted-label falsifier.

The Ahlmann-Eltze zero-shot transfer failed here with 20 of 20 deranged maps beating the real
fit, so a fit is not allowed to claim a gain until it has beaten its own permuted twin.  Every
byte of the pipeline is identical between the two arms except the query -> answer map inside
each batch, which is replaced by a derangement (no query keeps its own answer).

Two comparisons are reported, because the fit has two audiences:

  ddev    donor held-out perturbations -- where the real fit does gain
  dpanel  the held LINE                -- the number that would decide a slot

A real learning signal shows up as a real-vs-permuted separation on `ddev`.  Whether it also
separates on `dpanel` is the whole question of the experiment.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from foldio import Fold  # noqa: E402
from harness import LODO  # noqa: E402
from stage_b import run_fold  # noqa: E402

HERE = Path(__file__).resolve().parent


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=256)
    ap.add_argument("--rank", type=int, default=8)
    ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--out", default=str(HERE / "work" / "stage_c.json"))
    args = ap.parse_args()

    folds = {h: Fold(h) for h in LODO}
    base = dict(k=args.k, rank=args.rank, diag=True, tau=0.05, reg=1e-3, lr=0.01,
                gamma=0.65, batch=256, epochs=args.epochs, log_every=10**6)
    res = {"cfg": base, "arms": []}
    for permute in (False, True):
        n = args.seeds if permute else 3
        for s in range(n):
            cfg = dict(base, seed=195 + 1000 * s)
            per = {h: run_fold(folds[h], cfg, permute, log=lambda *_: None) for h in LODO}
            dpanel = {h: per[h]["best_by_dev"]["panel"] - per[h]["base_panel"] for h in LODO}
            ddev = {h: per[h]["best_by_dev"]["dev"] - per[h]["curve"][0]["dev"] for h in LODO}
            rec = {"permute": permute, "seed": cfg["seed"],
                   "macro_dpanel": float(np.mean(list(dpanel.values()))),
                   "macro_ddev": float(np.mean(list(ddev.values()))),
                   "per_fold_dpanel": dpanel, "per_fold_ddev": ddev}
            res["arms"].append(rec)
            print(f"  permute={int(permute)} seed={cfg['seed']:<6d} "
                  f"ddev {rec['macro_ddev']:+.4f}  dpanel {rec['macro_dpanel']:+.4f}  "
                  f"[{'  '.join(f'{h[:9]}={dpanel[h]:+.4f}' for h in LODO)}]", flush=True)
            Path(args.out).write_text(json.dumps(res, indent=2))

    real = [a for a in res["arms"] if not a["permute"]]
    null = [a for a in res["arms"] if a["permute"]]
    for kind in ("macro_ddev", "macro_dpanel"):
        r = np.array([a[kind] for a in real])
        p = np.array([a[kind] for a in null])
        beat = int((p >= r.mean()).sum())
        print(f"\n{kind}: real {r.mean():+.5f} (sd {r.std(ddof=1):.5f}, n={len(r)})   "
              f"permuted {p.mean():+.5f} (sd {p.std(ddof=1):.5f}, n={len(p)})   "
              f"permuted >= real: {beat}/{len(p)}   p95 {np.percentile(p, 95):+.5f}")
        res.setdefault("summary", {})[kind] = {
            "real_mean": float(r.mean()), "real_sd": float(r.std(ddof=1)),
            "null_mean": float(p.mean()), "null_sd": float(p.std(ddof=1)),
            "null_p95": float(np.percentile(p, 95)), "null_beats_real": beat, "n_null": len(p)}
    Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
