"""195 wrap-up — the fine gamma grid, an error bar on the arm, and the operator hash table.

Three things the write-up needs and the stages above do not produce:

1. Stage A's grid was coarse near the top (0.55 / 0.65 / 0.75).  "Both objectives peak at 0.65"
   is only worth stating if it survives a fine grid, so refine it.
2. The Stage B/C arm's panel delta needs an interval.  259 panel targets x 3 folds is small, so
   resample perturbations (paired, same targets in both arms) rather than quoting a point.
3. Rule 11: every whitening operator on disk gets a path and a sha256 here, including the ones
   this experiment did NOT use, because eight of them exist with near-identical names and the
   only defence is an explicit table.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fitlib import Adapter, batch_forward  # noqa: E402
from foldio import Fold, _pds  # noqa: E402
from harness import LODO  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINE = [0.45, 0.50, 0.55, 0.58, 0.60, 0.62, 0.65, 0.68, 0.70, 0.75, 0.80, 0.90]


def panel_rows(f: Fold, gamma: float, M=None, k=0) -> np.ndarray:
    w = f.evals ** (-gamma / 2.0)
    Z = f.panel_Y * w
    if M is not None and k:
        Z = Z.copy()
        Z[:, -k:] = Z[:, -k:] + Z[:, -k:] @ M
    Zn = Z / np.linalg.norm(Z, axis=1, keepdims=True)
    Tn = f.panel_T / np.linalg.norm(f.panel_T, axis=1, keepdims=True)
    return _pds(Zn @ Tn.T, np.arange(Z.shape[0]))


def main() -> None:
    out = {}

    # ---- 1. fine gamma grid: panel objective vs donor objective ----
    folds = {h: Fold(h) for h in LODO}
    fine = {"gammas": FINE, "panel_pds": [], "dev_pds": [], "panel_cos": []}
    for g in FINE:
        ps, ds, cs = [], [], []
        for h, f in folds.items():
            p, c = f.panel_scores(g)
            ps.append(p)
            cs.append(c)
            ds.append(f.batch_pds(f.batches("dev", g, 0, 256, seed=2), None, 0))
        fine["panel_pds"].append(float(np.mean(ps)))
        fine["dev_pds"].append(float(np.mean(ds)))
        fine["panel_cos"].append(float(np.mean(cs)))
        print(f"  gamma {g:<5g}  panel {fine['panel_pds'][-1]:.5f}  "
              f"donor-dev {fine['dev_pds'][-1]:.5f}  cos {fine['panel_cos'][-1]:.5f}", flush=True)
    fine["argmax_panel"] = FINE[int(np.argmax(fine["panel_pds"]))]
    fine["argmax_dev"] = FINE[int(np.argmax(fine["dev_pds"]))]
    fine["argmax_cos"] = FINE[int(np.argmax(fine["panel_cos"]))]
    print(f"  argmax panel {fine['argmax_panel']}   argmax donor-dev {fine['argmax_dev']}   "
          f"argmax cosine {fine['argmax_cos']}")
    out["fine_gamma"] = fine

    # ---- 2. paired perturbation bootstrap of the selected arm ----
    k, r = 256, 8
    base, arm = [], []
    for h, f in folds.items():
        fit_b = f.batches("fit", 0.65, k, 256, seed=11)
        dev_b = f.batches("dev", 0.65, k, 256, seed=12)
        ad = Adapter(k, r, True, seed=195)
        rng = np.random.default_rng(195)
        order = np.arange(len(fit_b))
        best, bestM = f.batch_pds(dev_b, None, k), None
        for _ in range(25):
            rng.shuffle(order)
            for i in order:
                _, _, dM = batch_forward(fit_b[i], ad.M(), 0.05)
                ad.step(dM, 0.01, 1e-3)
            d = f.batch_pds(dev_b, ad.M(), k)
            if d > best:
                best, bestM = d, ad.M().copy()
        base.append(panel_rows(f, 0.65))
        arm.append(panel_rows(f, 0.65, bestM, k) if bestM is not None else panel_rows(f, 0.65))
        print(f"  {h:14} base {base[-1].mean():.5f}  arm {arm[-1].mean():.5f}  "
              f"delta {arm[-1].mean()-base[-1].mean():+.5f}"
              f"{'  (dev never improved: arm = incumbent)' if bestM is None else ''}", flush=True)

    d = np.stack([a - b for a, b in zip(arm, base)])          # (3 folds, 259 targets)
    rngb = np.random.default_rng(195)
    boot = np.array([d[:, rngb.integers(0, d.shape[1], d.shape[1])].mean(1).mean()
                     for _ in range(20000)])
    out["arm_bootstrap"] = {
        "cfg": {"k": k, "rank": r, "diag": True, "tau": 0.05, "reg": 1e-3, "lr": 0.01,
                "gamma": 0.65, "epochs": 25, "selection": "max donor-dev pds"},
        "macro_base": float(np.mean([b.mean() for b in base])),
        "macro_arm": float(np.mean([a.mean() for a in arm])),
        "macro_delta": float(d.mean(1).mean()),
        "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "p_positive": float((boot > 0).mean()),
        "per_fold_delta": {h: float(x) for h, x in zip(LODO, d.mean(1))},
        "n_boot": 20000, "resampled": "panel perturbations, paired across arms"}
    b = out["arm_bootstrap"]
    print(f"\n  macro {b['macro_base']:.5f} -> {b['macro_arm']:.5f}   delta {b['macro_delta']:+.5f}"
          f"   95% CI [{b['ci95'][0]:+.5f}, {b['ci95'][1]:+.5f}]   P(positive) {b['p_positive']:.4f}")

    # ---- 3. every whitening operator on disk, with a sha ----
    table = []
    for p in sorted(list((ROOT / "data" / "donors").glob("whitener*.npz"))
                    + list(ROOT.glob("experiments/*/out/whitener*.npz"))
                    + list(ROOT.glob("experiments/*/work/**/whitener*.npz"))
                    + list((HERE / "work" / "operators").glob("*.npz"))):
        z = np.load(p)
        ev = z["evals"] if "evals" in z else None
        table.append({"path": str(p.relative_to(ROOT)),
                      "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                      "keys": list(z.keys()),
                      "max_eigenvalue": float(ev.max()) if ev is not None else None,
                      "n_dim": int(ev.shape[0]) if ev is not None else None,
                      "used_by_195": p.parent.name == "operators"})
        t = table[-1]
        print(f"  {t['path']:72} max_ev={t['max_eigenvalue']}  used={int(t['used_by_195'])}"
              f"  sha={t['sha256'][:16]}")
    out["operator_table"] = table
    (HERE / "work" / "wrapup.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
