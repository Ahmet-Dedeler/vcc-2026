"""195 Stage D — the panel-basis artefact, reproduced on purpose as a positive control.

AGENTS.md 4.2 point 1 says the 259 x 7,760 panel basis is rank-deficient and its gamma sweep
"rises monotonically forever", because at the limit the operator just orthogonalises 259
vectors and a ranking over 259 items improves by construction.  That is the failure signature
this experiment has to be able to detect, so it is reproduced here in two forms:

  D2a  fit the whitening covariance on the 259 panel PREDICTIONS instead of the donor bank
       and sweep gamma;
  D2b  fit the learned adapter directly on the held line's own 259-way retrieval problem
       (query = panel prediction, answer = that line's true effect -- pure leakage) and sweep
       capacity.

If either rises without limit while the legitimate donor-fit does not, the harness can tell
the two apart, which is what makes the Stage B negative trustworthy.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fitlib import Adapter, BatchStats, batch_forward  # noqa: E402
from foldio import Fold, _pds  # noqa: E402
from harness import LODO  # noqa: E402

HERE = Path(__file__).resolve().parent
GAMMAS = [0.0, 0.2, 0.35, 0.55, 0.65, 0.9, 1.3, 2.0, 3.0, 5.0, 8.0]


def panel_basis_sweep(f: Fold) -> dict:
    """D2a: covariance from the panel predictions themselves (still in V's coordinates)."""
    P = f.panel_Y                                   # (259, G), prediction in donor-eigen coords
    C = (P.T @ P) / P.shape[0]
    ev, U = np.linalg.eigh(C)
    ev = np.maximum(ev, ev.max() * 1e-12)
    Tn = f.panel_T / np.linalg.norm(f.panel_T, axis=1, keepdims=True)
    out = {}
    for g in GAMMAS:
        Z = ((P @ U) * (ev ** (-g / 2.0))) @ U.T if g else P
        Zn = Z / np.linalg.norm(Z, axis=1, keepdims=True)
        cos = Zn @ Tn.T
        out[str(g)] = float(_pds(cos, np.arange(cos.shape[0])).mean())
    return out


def panel_fit_sweep(f: Fold, ranks: list[int], epochs: int, k: int = 256) -> dict:
    """D2b: the adapter fit ON the held line's own retrieval problem.  Leakage, by design."""
    Z = (f.panel_Y * (f.evals ** (-0.65 / 2.0)))
    T = f.panel_T
    bs = BatchStats.__new__(BatchStats)
    bs.Zp, bs.Tb = Z[:, -k:].copy(), T[:, -k:].copy()
    bs.R = Z @ T.T - bs.Zp @ bs.Tb.T
    bs.c = np.einsum("ij,ij->i", Z, Z) - np.einsum("ij,ij->i", bs.Zp, bs.Zp)
    bs.tn = np.linalg.norm(T, axis=1)
    bs.B = Z.shape[0]
    out = {}
    for r in ranks:
        ad = Adapter(k, r, True, seed=195)
        for _ in range(epochs):
            _, _, dM = batch_forward(bs, ad.M(), 0.05)
            ad.step(dM, 0.01, 1e-3)
        out[str(r)] = f.panel_scores(0.65, ad.M(), k)[0]
    return out


def donor_fit_rank_sweep(f: Fold, ranks: list[int], epochs: int, k: int = 256) -> dict:
    """D1: the legitimate fit, same capacity ladder, scored on the held LINE and on dev."""
    fit_b = f.batches("fit", 0.65, k, 256, seed=11)
    dev_b = f.batches("dev", 0.65, k, 256, seed=12)
    out = {}
    order = np.arange(len(fit_b))
    for r in ranks:
        ad = Adapter(k, r, True, seed=195)
        rng = np.random.default_rng(195)
        for _ in range(epochs):
            rng.shuffle(order)
            for i in order:
                _, _, dM = batch_forward(fit_b[i], ad.M(), 0.05)
                ad.step(dM, 0.01, 1e-3)
        M = ad.M()
        out[str(r)] = {"panel": f.panel_scores(0.65, M, k)[0],
                       "dev": f.batch_pds(dev_b, M, k),
                       "train": f.batch_pds(fit_b, M, k)}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--panel-epochs", type=int, default=400)
    ap.add_argument("--out", default=str(HERE / "work" / "stage_d.json"))
    args = ap.parse_args()
    ranks = [0, 2, 4, 8, 16, 32, 64, 128, 256]
    res = {"gammas": GAMMAS, "ranks": ranks, "folds": {}}
    for held in LODO:
        t0 = time.time()
        f = Fold(held)
        r = {"panel_basis_gamma": panel_basis_sweep(f),
             "panel_fit_rank": panel_fit_sweep(f, ranks, args.panel_epochs),
             "donor_fit_rank": donor_fit_rank_sweep(f, ranks, args.epochs),
             "base_panel": f.panel_scores(0.65)[0]}
        res["folds"][held] = r
        print(f"{held}  ({time.time()-t0:.0f}s)", flush=True)

    def macro(path):
        keys = list(res["folds"][LODO[0]][path].keys())
        return {kk: float(np.mean([res["folds"][h][path][kk] if not isinstance(
            res["folds"][h][path][kk], dict) else res["folds"][h][path][kk]["panel"]
            for h in LODO])) for kk in keys}

    res["macro"] = {"panel_basis_gamma": macro("panel_basis_gamma"),
                    "panel_fit_rank": macro("panel_fit_rank"),
                    "donor_fit_rank_panel": macro("donor_fit_rank"),
                    "donor_fit_rank_dev": {kk: float(np.mean(
                        [res["folds"][h]["donor_fit_rank"][kk]["dev"] for h in LODO]))
                        for kk in [str(x) for x in ranks]},
                    "base_panel": float(np.mean([res["folds"][h]["base_panel"] for h in LODO]))}
    base = res["macro"]["base_panel"]
    print(f"\nD2a  panel-basis whitener, macro panel pds (donor-basis gamma=0.65 = {base:.4f})")
    for g in GAMMAS:
        print(f"   gamma {g:<5g} -> {res['macro']['panel_basis_gamma'][str(g)]:.4f}")
    print("\nD2b/D1  capacity ladder, macro panel pds")
    print("   rank | panel-FIT (leakage) | donor-FIT panel | donor-FIT dev")
    for r in ranks:
        print(f"   {r:<5d}|      {res['macro']['panel_fit_rank'][str(r)]:.4f}         |"
              f"     {res['macro']['donor_fit_rank_panel'][str(r)]:.4f}      |"
              f"   {res['macro']['donor_fit_rank_dev'][str(r)]:.4f}")
    Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
