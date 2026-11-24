"""Experiment 004 — real K562-gwps -> H1 transfer, measured offline.

Every scored metric is a function of `cos θ` between the predicted and true effect vector
(experiment 002). That cosine is computable in seconds from pseudobulk alone, while a full
`cell-eval2` run costs ~7 minutes. So the model-design search happens here, and only the
winner goes to the real scorer for (alpha, dispersion) calibration.

Experiment 003 measured K562 -> RPE1 at cos 0.128 (gene-specific, centered). H1 hESC is a
different and arguably harder target: embryonic stem cells versus a leukemia line. This is the
closest available proxy for K562 -> the 2026 contexts A/B/C.

Reported for each variant:
  raw       cos(pred, true)                            what pds/mse actually see
  centered  cos(pred - mean_pred, true - mean_true)    the part that earns score
  Because the scored baseline IS the mean response, `centered` is the honest number.

    uv run python experiments/004-donor-transfer-h1/offline_eval.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TransferModel, TS  # noqa: E402

CONTROL = "non-targeting"
PERT_COL = "target_gene"


def load_target(path: Path):
    a = ad.read_h5ad(path)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    labels = a.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(a.var_names)

    Pc = np.asarray(X[labels == CONTROL].sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())

    perts, deltas = [], []
    for g in sorted(set(labels) - {CONTROL}):
        P = np.asarray(X[labels == g].sum(0)).ravel().astype(np.float64)
        deltas.append(np.log1p(TS * P / P.sum()) - b_ctrl)
        perts.append(g)
    return np.asarray(perts), np.vstack(deltas), b_ctrl, genes


def cosines(pred: np.ndarray, true: np.ndarray) -> dict:
    def cs(A, B):
        na = np.linalg.norm(A, axis=1); nb = np.linalg.norm(B, axis=1)
        ok = (na > 0) & (nb > 0)
        out = np.zeros(len(A)); out[ok] = np.sum(A[ok] * B[ok], axis=1) / (na[ok] * nb[ok])
        return out
    raw = cs(pred, true)
    cen = cs(pred - pred.mean(0, keepdims=True), true - true.mean(0, keepdims=True))
    return {"raw": raw, "centered": cen}


def report(name: str, c: dict) -> dict:
    r, x = c["raw"], c["centered"]
    print(f"  {name:34} raw {r.mean():+.4f}   centered {x.mean():+.4f}  "
          f"(median {np.median(x):+.4f}, frac>0 {(x > 0).mean():.2f})")
    return {"variant": name, "raw_mean": float(r.mean()),
            "centered_mean": float(x.mean()), "centered_median": float(np.median(x)),
            "centered_frac_pos": float((x > 0).mean())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="data/harness/h1-test.h5ad")
    ap.add_argument("--donor", default="data/donors/k562_gwps")
    ap.add_argument("--out", default="experiments/004-donor-transfer-h1")
    args = ap.parse_args()

    perts, true, b_ctrl, genes = load_target(Path(args.target))
    print(f"target {Path(args.target).name}: {len(perts)} perturbations, {len(genes)} genes")

    donor = Donor.load(Path(args.donor))
    model = TransferModel(donor, genes)
    cov = model.covered(perts)
    print(f"donor {Path(args.donor).name}: {donor.delta.shape}, "
          f"gene-axis coverage {model.coverage:.1%}, "
          f"perturbation coverage {cov.sum()}/{len(perts)}\n")

    # Restrict to covered perturbations: uncovered ones are a separate problem.
    perts_c, true_c = perts[cov], true[cov]
    results = []

    print("VARIANT SWEEP  (alpha=1, beta=0 unless stated; scale is irrelevant to cosine)")
    m_full = model.m_full

    def build(alpha, beta, gate, space="log1p"):
        out = []
        for p in perts_c:
            eff = model.effect(p, alpha, beta)
            if gate:
                eff = np.where(b_ctrl > 0, eff, 0.0)
            if space == "profile":
                pr = model.profile(p, b_ctrl, alpha, beta, gate_unexpressed=gate)
                eff = np.log1p(TS * pr) - b_ctrl
            out.append(eff)
        return np.vstack(out)

    results.append(report("s only (gene-specific)", cosines(build(1, 0, False), true_c)))
    results.append(report("s only + expression gate", cosines(build(1, 0, True), true_c)))
    results.append(report("m only (mean response)", cosines(build(0, 1, True), true_c)))
    results.append(report("m + s", cosines(build(1, 1, True), true_c)))
    results.append(report("full delta via profile()", cosines(build(1, 1, True, "profile"), true_c)))

    # Does the target's own mean response transfer? (we cannot know it at test time,
    # but this bounds what a perfect m-estimate would buy)
    true_m = true_c.mean(0)
    mm = np.zeros_like(m_full); mm[:] = m_full
    denom = np.linalg.norm(mm) * np.linalg.norm(true_m)
    print(f"\n  cos(donor mean response, H1 mean response) = "
          f"{float(mm @ true_m / denom) if denom else float('nan'):+.4f}")

    # Restricting the comparison to genes the donor actually carries isolates how much of
    # the loss is missing genes rather than wrong biology.
    sub = model.tcols
    print("\nON THE DONOR-COVERED GENE SUBSET ONLY "
          f"({len(sub)} of {len(genes)} genes)")
    results.append(report("s only, covered genes",
                          cosines(build(1, 0, True)[:, sub], true_c[:, sub])))

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    (out / "offline_results.json").write_text(json.dumps(results, indent=2))
    best = max(results, key=lambda r: r["centered_mean"])
    print(f"\nbest variant by centered cos: {best['variant']} -> {best['centered_mean']:+.4f}")
    print(f"(experiment 003 measured K562 -> RPE1 at +0.128 for comparison)")


if __name__ == "__main__":
    main()
