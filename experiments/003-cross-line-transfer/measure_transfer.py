"""Experiment 003 — how well does a perturbation response transfer between cell lines?

This is the number the whole project turns on. Experiment 002 showed every metric is a
function of directional accuracy `cos θ` between the predicted and true effect vectors.
The simplest possible model — look up the same gene's response in a donor line and copy it
— achieves exactly `cos(δ_donor, δ_target)`. So measuring that cosine on real data tells us
the score a naive donor transfer earns, and how much headroom a real model has to add.

Replogle 2022 gives us the only clean read on it: ~2,000 perturbations measured in **both**
K562 (from the genome-wide screen) and RPE1, same protocol, same lab, same year. Different
tissue of origin, which is the axis the 2026 Challenge tests.

Three quantities, because a raw cosine means nothing without its null:

  transfer   cos(δ_g^K562, δ_g^RPE1)          same gene, different line   <- the signal
  mismatch   cos(δ_g^K562, δ_h^RPE1), g != h  different gene              <- the null
  meanresp   cos(δ_g^K562, mean_h δ_h^RPE1)   the context-mean baseline   <- what 0 means

If `transfer` is not clearly above `mismatch`, donor transfer cannot work at all. If
`transfer` is not clearly above `meanresp`, it earns nothing over the scored baseline.

    uv run python experiments/003-cross-line-transfer/measure_transfer.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def load(donor_dir: Path):
    delta = np.load(donor_dir / "delta.npy")
    perts = (donor_dir / "perts.txt").read_text().split("\n")
    genes = (donor_dir / "genes.txt").read_text().split("\n")
    return delta, np.asarray(perts), np.asarray(genes)


def unit(M: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(M, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return M / n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--a", default="data/donors/k562_gwps")
    ap.add_argument("--b", default="data/donors/rpe1")
    ap.add_argument("--out", default="experiments/003-cross-line-transfer")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    A, pA, gA = load(Path(args.a))
    B, pB, gB = load(Path(args.b))
    print(f"A {Path(args.a).name}: {A.shape}   B {Path(args.b).name}: {B.shape}")

    # Align gene axes. Effects live in log-pseudobulk space on each line's own axis, so the
    # shared axis is the only space where a cosine between them is meaningful.
    shared_g = np.intersect1d(gA, gB)
    ia = {g: i for i, g in enumerate(gA)}
    ib = {g: i for i, g in enumerate(gB)}
    A = A[:, [ia[g] for g in shared_g]]
    B = B[:, [ib[g] for g in shared_g]]

    shared_p = np.intersect1d(pA, pB)
    pa = {p: i for i, p in enumerate(pA)}
    pb = {p: i for i, p in enumerate(pB)}
    Asub = A[[pa[p] for p in shared_p]]
    Bsub = B[[pb[p] for p in shared_p]]
    print(f"shared: {len(shared_g):,} genes, {len(shared_p):,} perturbations\n")

    Au, Bu = unit(Asub), unit(Bsub)

    transfer = np.sum(Au * Bu, axis=1)

    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(len(shared_p))
    perm[perm == np.arange(len(shared_p))] = (perm[perm == np.arange(len(shared_p))] + 1) % len(shared_p)
    mismatch = np.sum(Au * Bu[perm], axis=1)

    meanB = unit(Bsub.mean(axis=0, keepdims=True))
    meanresp = (Au @ meanB.T).ravel()

    def q(x, name):
        return (f"{name:10} mean {x.mean():+.4f}  median {np.median(x):+.4f}  "
                f"p10 {np.percentile(x,10):+.4f}  p90 {np.percentile(x,90):+.4f}  "
                f"frac>0 {(x>0).mean():.2f}")

    print(q(transfer, "transfer"))
    print(q(mismatch, "mismatch"))
    print(q(meanresp, "meanresp"))

    # Effect size relative to the null, and how much of the transfer is just the shared
    # mean response every perturbation carries.
    resid = transfer - meanresp
    print(f"\ntransfer - mismatch : mean {(transfer-mismatch).mean():+.4f}")
    print(f"transfer - meanresp : mean {resid.mean():+.4f}  frac>0 {(resid>0).mean():.2f}")

    # After removing each line's mean response, is there gene-specific signal left?
    Ac = unit(Asub - Asub.mean(axis=0, keepdims=True))
    Bc = unit(Bsub - Bsub.mean(axis=0, keepdims=True))
    centered = np.sum(Ac * Bc, axis=1)
    cperm = np.sum(Ac * Bc[perm], axis=1)
    print(f"\nafter removing each line's mean response (the part the baseline already gets):")
    print(q(centered, "  centered"))
    print(q(cperm, "  null"))

    strong = np.argsort(-transfer)[:15]
    print(f"\nbest-transferring perturbations:")
    print("  " + ", ".join(f"{shared_p[i]}({transfer[i]:.2f})" for i in strong))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "transfer.npz", perts=shared_p, transfer=transfer,
             mismatch=mismatch, meanresp=meanresp, centered=centered)
    (out / "summary.json").write_text(json.dumps({
        "donor_a": args.a, "donor_b": args.b,
        "n_shared_genes": int(len(shared_g)), "n_shared_perts": int(len(shared_p)),
        "transfer_mean": float(transfer.mean()), "transfer_median": float(np.median(transfer)),
        "mismatch_mean": float(mismatch.mean()), "meanresp_mean": float(meanresp.mean()),
        "centered_mean": float(centered.mean()), "centered_null_mean": float(cperm.mean()),
    }, indent=2))
    print(f"\nwrote {out}/summary.json")


if __name__ == "__main__":
    main()
