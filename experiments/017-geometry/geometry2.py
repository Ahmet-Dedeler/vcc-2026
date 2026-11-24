"""017b — spectral flattening, with the harness corrected and the sweep extended.

## The harness correction

`lodo.py` applied the same treatment (`gene_norm` / `eb`) to the held-out line's effects as
to the predictors. That makes the *ground truth* a normalised object, which it is not: the
scorer computes real effects as plain log-pseudobulk differences from the reference counts.
Normalising both sides can flatter any method that likes normalised geometry -- and
whitening is exactly such a method, so this correction has to be made before believing 017.

Here the truth is always `raw` delta unless overridden. Predictions keep their treatment.

## What is swept

1. `gamma` past 1 -- 017 was monotone right up to the edge of its grid, which usually means
   the grid was too small, not that the optimum sits on the boundary.
2. The **basis** the whitening is computed in. Whitening with the panel's own 259x7760
   prediction matrix estimates a 7,760-dimensional covariance from 259 samples. The donors
   carry 9,787 and 17,767 perturbations; a covariance estimated from all of them is far
   better conditioned and is still target-blind.
3. Ridge `lam`, since the panel covariance is rank-deficient by construction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale  # noqa: E402

DONORS = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
}


def pds_cosine(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    P = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    T = true / (np.linalg.norm(true, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - P @ T.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return 1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)


def load_banks() -> dict:
    out = {}
    for name, path in DONORS.items():
        d = Donor.load(path)
        _, s = d.decompose()
        gl = [str(g) for g in d.genes]
        out[name] = dict(s_raw=s, s_gn=(s / _norm_scale(s)).astype(np.float32),
                         gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                         idx={p: i for i, p in enumerate(d.perts)}, b_ctrl=d.b_ctrl)
    return out


_EIG: dict = {}


def eigbasis(key, basis: np.ndarray):
    """Cache the eigendecomposition of C = basis^T basis / n once per (fold, basis).

    Every (gamma, lam) is then a reweighting of the same eigenvalues, which turns an
    8-value gamma sweep from eight SVDs into one.
    """
    if key not in _EIG:
        n = basis.shape[0]
        _, S, Vt = np.linalg.svd(basis, full_matrices=False)
        ev = S ** 2 / n
        keep = ev > ev.max() * 1e-12
        _EIG[key] = (ev[keep], Vt[keep].T)
    return _EIG[key]


def whiten_apply(P: np.ndarray, key, basis: np.ndarray, gamma: float,
                 lam: float) -> np.ndarray:
    """P (C + lam I)^(-gamma/2), applied in the basis' own non-null subspace.

    Directions the basis never observed are left untouched rather than rescaled, so the
    operator is a reweighting of observed structure and not an extrapolation.
    """
    ev, V = eigbasis(key, basis)
    w = (ev + lam) ** (-gamma / 2.0)
    PV = P @ V
    return P - PV @ V.T + (PV * w) @ V.T


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weights", default="4:1")
    ap.add_argument("--treat", default="gene_norm", choices=["raw", "gene_norm"])
    ap.add_argument("--true-treat", default="raw", choices=["raw", "gene_norm"])
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0])
    ap.add_argument("--lams", type=float, nargs="+", default=[0.0])
    ap.add_argument("--basis", nargs="+", default=["panel", "donor"])
    ap.add_argument("--max-basis-perts", type=int, default=6000)
    args = ap.parse_args()

    banks = load_banks()
    vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
    common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
    print(f"panel {len(common)}  treat={args.treat}  TRUTH={args.true_treat}")

    folds = {}
    rng = np.random.default_rng(0)
    for held in banks:
        others = [n for n in banks if n != held]
        shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
        feat = [g for g in shared if g not in set(common)]
        col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
               for n in [held] + others}
        key = "s_gn" if args.treat == "gene_norm" else "s_raw"
        tkey = "s_gn" if args.true_treat == "gene_norm" else "s_raw"
        T = np.vstack([banks[held][tkey][banks[held]["idx"][p]][col[held]] for p in common])
        ws = dict(zip(others, [float(x) for x in args.weights.split(":")]))
        scale = {n: float(np.median(np.linalg.norm(banks[n][key][:, col[n]], axis=1))) or 1.0
                 for n in others}
        ref = float(np.median(list(scale.values())))
        P = np.zeros((len(common), len(feat)), np.float64)
        # Donor-wide basis: every perturbation the predictor lines measured, panel or not.
        Bl = []
        for n in others:
            P += (ws[n] / sum(ws.values())) * (ref / scale[n]) * np.vstack(
                [banks[n][key][banks[n]["idx"][p]][col[n]] for p in common])
            M = banks[n][key][:, col[n]]
            if M.shape[0] > args.max_basis_perts:
                M = M[np.sort(rng.choice(M.shape[0], args.max_basis_perts, False))]
            Bl.append((ref / scale[n]) * M.astype(np.float64))
        folds[held] = dict(P=P, T=T.astype(np.float64), B=np.vstack(Bl))
        print(f"  {held:14} P{P.shape}  donor-basis {folds[held]['B'].shape}")

    results = {}
    for basis_kind in args.basis:
        for lam in args.lams:
            for g in args.gammas:
                per = {}
                for h, f in folds.items():
                    if g == 0.0:
                        Pg = f["P"]
                    else:
                        B = f["P"] if basis_kind == "panel" else f["B"]
                        Pg = whiten_apply(f["P"], (h, basis_kind), B, g, lam)
                    per[h] = float(pds_cosine(Pg, f["T"]).mean())
                macro = float(np.mean(list(per.values())))
                key = f"{basis_kind:5} lam={lam:g} gamma={g:g}"
                results[key] = {**per, "macro": macro}
                fs = "  ".join(f"{k.replace('orion_','')[:9]}={v:.4f}" for k, v in per.items())
                print(f"  {key:30} macro PDS {macro:.4f}   [{fs}]", flush=True)

    best = max(results, key=lambda k: results[k]["macro"])
    print(f"\nbest: {best} -> {results[best]['macro']:.4f}")
    out = Path(f"experiments/017-geometry/results_{args.treat}_true{args.true_treat}.json")
    out.write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
