"""195 Stage B — fit the correction to the ranking objective, LOLO, and price it on the line.

Protocol.  For each LODO fold the operator, the subspace, the queries and the galleries all
come from the two non-held donor lines only.  Perturbations are split 80/20 into a fit set and
a dev set; every hyper-parameter is chosen on **dev**, which is a held-out perturbation split
of donor data and never the held line.  The held line's 259 panel targets are scored once per
configuration and never used to choose anything.

Three curves are tracked at every epoch:

    train  in-batch pds on the fit perturbations   (does the optimiser work at all)
    dev    in-batch pds on held-out perturbations  (does it generalise across perturbations)
    panel  259-way pds on the HELD LINE            (does it generalise across cell lines)

`--permute` deranges the query -> answer map inside every batch, refits with everything else
identical, and is the falsifier: a real fit must beat its own permuted twin by more than the
fold spread.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fitlib import Adapter, batch_forward  # noqa: E402
from foldio import Fold  # noqa: E402
from harness import LODO  # noqa: E402

HERE = Path(__file__).resolve().parent


def derange(n: int, rng: np.random.Generator) -> np.ndarray:
    """A permutation with no fixed point, so no query keeps its own true answer."""
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            return p


def run_fold(f: Fold, cfg: dict, permute: bool, log=print) -> dict:
    k, rank, tau = cfg["k"], cfg["rank"], cfg["tau"]
    fit_b = f.batches("fit", cfg["gamma"], k, cfg["batch"], seed=11)
    dev_b = f.batches("dev", cfg["gamma"], k, cfg["batch"], seed=12)
    rng = np.random.default_rng(cfg["seed"] + (9999 if permute else 0))
    labels = [derange(bs.B, rng) for bs in fit_b] if permute else [None] * len(fit_b)

    ad = Adapter(k, rank, cfg["diag"], seed=cfg["seed"])
    base_panel, base_cos = f.panel_scores(cfg["gamma"])
    curve = [{"epoch": 0, "loss": None,
              "train": f.batch_pds(fit_b, None, k), "dev": f.batch_pds(dev_b, None, k),
              "panel": base_panel, "panel_cos": base_cos}]
    order = np.arange(len(fit_b))
    for ep in range(1, cfg["epochs"] + 1):
        rng.shuffle(order)
        losses = []
        for i in order:
            M = ad.M()
            loss, _, dM = batch_forward(fit_b[i], M, tau, labels[i])
            losses.append(loss)
            ad.step(dM, cfg["lr"], cfg["reg"])
        M = ad.M()
        p, c = f.panel_scores(cfg["gamma"], M, k)
        curve.append({"epoch": ep, "loss": float(np.mean(losses)),
                      "train": f.batch_pds(fit_b, M, k), "dev": f.batch_pds(dev_b, M, k),
                      "panel": p, "panel_cos": c})
        if ep % cfg.get("log_every", 5) == 0 or ep == cfg["epochs"]:
            e = curve[-1]
            log(f"      ep{ep:<3d} loss {e['loss']:.4f}  train {e['train']:.4f}  "
                f"dev {e['dev']:.4f}  panel {e['panel']:.4f} ({e['panel']-base_panel:+.4f})  "
                f"cos {e['panel_cos']:.4f}")
    best = max(curve, key=lambda e: e["dev"])
    return {"curve": curve, "base_panel": base_panel, "base_cos": base_cos,
            "best_by_dev": best, "n_params": ad.n_params(),
            "n_fit": int(sum(b.B for b in fit_b)), "n_dev": int(sum(b.B for b in dev_b)),
            "M_absmax": float(np.abs(ad.M()).max())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", default="main")
    ap.add_argument("--permute", action="store_true")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.grid == "main":
        cfgs = [dict(k=k, rank=r, diag=d, tau=t, reg=g, lr=lr)
                for k, r, d in [(64, 0, True), (64, 8, True), (64, 64, True),
                                (256, 8, True), (256, 32, True), (256, 256, False),
                                (1024, 16, True)]
                for t in [0.05]
                for g in [1e-3]
                for lr in [0.01]]
    elif args.grid == "tau":
        cfgs = [dict(k=256, rank=32, diag=True, tau=t, reg=g, lr=lr)
                for t in [0.02, 0.05, 0.12, 0.30]
                for g in [1e-4, 1e-2]
                for lr in [0.003, 0.03]]
    elif args.grid == "best":
        cfgs = [dict(k=64, rank=64, diag=True, tau=0.05, reg=1e-3, lr=0.01)]
    else:
        raise SystemExit(f"unknown grid {args.grid}")

    out = {"grid": args.grid, "permute": args.permute, "runs": []}
    folds = {h: Fold(h) for h in LODO}
    for cfg in cfgs:
        cfg.update(gamma=0.65, batch=256, seed=195, epochs=args.epochs)
        t0 = time.time()
        print(f"\ncfg {cfg}", flush=True)
        per = {}
        for held in LODO:
            print(f"    fold {held}", flush=True)
            per[held] = run_fold(folds[held], cfg, args.permute)
        macro_base = float(np.mean([per[h]["base_panel"] for h in LODO]))
        macro_best = float(np.mean([per[h]["best_by_dev"]["panel"] for h in LODO]))
        macro_cos = float(np.mean([per[h]["best_by_dev"]["panel_cos"] for h in LODO]))
        macro_dev = float(np.mean([per[h]["best_by_dev"]["dev"] for h in LODO]))
        macro_dev0 = float(np.mean([per[h]["curve"][0]["dev"] for h in LODO]))
        macro_tr = float(np.mean([per[h]["best_by_dev"]["train"] for h in LODO]))
        macro_tr0 = float(np.mean([per[h]["curve"][0]["train"] for h in LODO]))
        rec = {"cfg": cfg, "per_fold": per, "macro_base_panel": macro_base,
               "macro_panel": macro_best, "macro_dpanel": macro_best - macro_base,
               "macro_panel_cos": macro_cos, "macro_ddev": macro_dev - macro_dev0,
               "macro_dtrain": macro_tr - macro_tr0,
               "n_params": per[LODO[0]]["n_params"], "secs": time.time() - t0}
        out["runs"].append(rec)
        print(f"  => params {rec['n_params']:6d}  dtrain {rec['macro_dtrain']:+.4f}  "
              f"ddev {rec['macro_ddev']:+.4f}  PANEL {macro_base:.4f} -> {macro_best:.4f} "
              f"({rec['macro_dpanel']:+.4f})  cos {macro_cos:.4f}  [{rec['secs']:.0f}s]",
              flush=True)
        Path(args.out).write_text(json.dumps(out, indent=2))

    print("\n== summary (selection on donor dev, price on the held line) ==")
    for r in out["runs"]:
        c = r["cfg"]
        print(f" k={c['k']:<5d} r={c['rank']:<4d} diag={int(c['diag'])} tau={c['tau']:<5g} "
              f"reg={c['reg']:<7g} lr={c['lr']:<6g} p={r['n_params']:<7d} "
              f"dtrain {r['macro_dtrain']:+.4f} ddev {r['macro_ddev']:+.4f} "
              f"dpanel {r['macro_dpanel']:+.4f}")


if __name__ == "__main__":
    main()
