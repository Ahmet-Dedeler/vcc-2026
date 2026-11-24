"""160 — is a cell line's PERTURBATION-RESPONSE subspace visible in its own CONTROL cells?

THE QUESTION
------------
Production builds ONE context-independent effect matrix and applies it to A/B/C, modulated only
by each context's CPM gate and control profile. Experiment 159 measured that ~70% of a target's
response energy is line-specific, and that knowing the CONTEXT is worth ~3x knowing the TARGET
(same-line prediction cos 0.2627 against the best cross-line single donor 0.0897). So the axis
we spend no modelling on is the more valuable one.

The only thing the challenge gives us about a context is 18,400 unperturbed cells. The
hypothesis worth testing is geometric rather than per-gene:

    a perturbation pushes the cell along directions the population ALREADY varies in

Unperturbed populations contain cells at different cell-cycle phases, stress levels and
metabolic states, so the programmes a knockdown engages (growth arrest, ribosome biogenesis,
stress) should be spanned by the principal directions of natural variation. If true, the
control cells hand us a CONTEXT-SPECIFIC basis for free, and the transfer problem collapses
from "predict 18,533 gene directions" to "predict ~50 loadings" -- a far easier object, and one
that can be conditioned on the context by construction.

WHAT IS ALREADY CLOSED, AND WHY THIS IS NOT IT
---------------------------------------------
* 155 regressed every gene on the TARGET gene across control cells (split-half r = 0.028, null).
  That is a per-gene conditional-expectation estimator, not a subspace.
* 067 rescaled effect coordinates by per-gene control variability (macro -0.00025).
* 017 modulated per-gene by the basal profile `p^beta` (+0.003 at beta 0.25).
* 111 fitted DepMap context programmes; Heisenberg (board) ADDED context-shared co-expression
  as a direction and pds fell 0.4509 -> 0.3974.
None of them asked whether the response lives in the SPAN of control variation, and none used
that span as a projection. Adding a shared co-expression direction (Heisenberg) is close to the
opposite operation: it injects one hub direction, where this projects out everything outside a
context-specific subspace.

DESIGN, WITH THE CONTROLS THAT MAKE IT FALSIFIABLE
--------------------------------------------------
For each line: control cells -> top-k right singular vectors `V_k` (the control basis).
True effects `E` (n_targets x G) -> shared mean `m` + target-specific residual `R = E - m`.
Only `R` can earn score (the mean response IS the metric's zero), so `R` is what must be
contained, and reporting containment of `E` alone would be self-deceiving.

Captured fraction `||R V^T||^2 / ||R||^2` is compared against FOUR references:
  own      the line's own control basis            <- the hypothesis
  other    a DIFFERENT line's control basis        <- is it THIS line's structure, or any
                                                      biological co-expression structure?
  random   random directions with the same per-gene sd envelope as the control cells
                                                   <- a k-dim subspace captures k/G by chance,
                                                      but effects and controls share a gene-
                                                      expression envelope, so the naive null is
                                                      far too weak
  oracle   top-k PCs of R itself                   <- the ceiling at this k

`other` is the load-bearing control. If own ~= other the subspace is generic biology and buys
no context conditioning, which is the outcome 017/010/155 would predict.
"""
from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

TS = 5.0e4
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
KS = (5, 10, 20, 50, 100, 200)
SEED = 20260902

SURFACES = {
    "hct116": ("data/harness/orion-hct116-thinned.h5ad",
               "experiments/156-delivered-geometry/work/truth-hct116.npz"),
    "hek293t": ("data/harness/orion-hek293t-thinned.h5ad",
                "experiments/156-delivered-geometry/work/truth-hek293t.npz"),
}
CTRL = "non-targeting"


def _cats(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        c = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                        for x in o["categories"][:]])
        return c[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def control_cells(path: str, label: str = CTRL, chunk: int = 4000):
    """Control cells in the scorer's own per-cell space: x = log1p(TS * counts/sum)."""
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        lab = _cats(h, "obs/target_gene")
        rows = np.flatnonzero(lab == label)
        indptr = h["X/indptr"][:]
        out = np.zeros((len(rows), g), dtype=np.float32)
        pos = 0
        rowset = set(rows.tolist())
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            sel = [r for r in range(r0, r1) if r in rowset]
            if not sel:
                continue
            s, e = int(indptr[r0]), int(indptr[r1])
            X = sp.csr_matrix((h["X/data"][s:e].astype(np.float64),
                               h["X/indices"][s:e], indptr[r0:r1 + 1] - indptr[r0]),
                              shape=(r1 - r0, g))
            sub = X[[r - r0 for r in sel]].toarray()
            tot = sub.sum(1, keepdims=True)
            tot[tot == 0] = 1.0
            out[pos:pos + len(sel)] = np.log1p(TS * sub / tot).astype(np.float32)
            pos += len(sel)
    assert pos == len(rows)
    return out


def top_basis(X: np.ndarray, k: int, seed: int = SEED):
    """Top-k right singular vectors of the cell-centred control matrix, as rows."""
    from sklearn.utils.extmath import randomized_svd
    mu = X.mean(0, keepdims=True)
    _, _, Vt = randomized_svd(X - mu, n_components=k, n_iter=5, random_state=seed)
    return Vt                                    # (k, G), orthonormal rows


def captured(R: np.ndarray, V: np.ndarray) -> float:
    """Fraction of R's Frobenius energy inside span(rows of V)."""
    num = float((np.asarray(R @ V.T) ** 2).sum())
    den = float((R * R).sum())
    return num / max(den, 1e-300)


def random_basis(gene_sd: np.ndarray, k: int, rng) -> np.ndarray:
    """k random directions carrying the control's per-gene sd envelope, orthonormalised.

    A uniformly random k-dim subspace captures k/G by chance, which is a vacuous null here:
    both the effects and the control PCs concentrate on expressed genes, so any envelope-matched
    subspace already captures far more than k/G. This is the null that isolates STRUCTURE.
    """
    A = rng.standard_normal((k, len(gene_sd))) * gene_sd[None, :]
    Q, _ = np.linalg.qr(A.T)
    return Q.T[:k]


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(SEED)
    data = {}
    for name, (h5, truth) in SURFACES.items():
        print(f"[{name}] loading control cells ...", flush=True)
        X = control_cells(os.path.join(REPO, h5))
        z = np.load(os.path.join(REPO, truth), allow_pickle=True)
        E = z["effect"].astype(np.float64)
        tg = np.asarray([str(t) for t in z["targets"]])
        genes = np.asarray([str(g) for g in z["genes"]])
        # pds drops every panel target coordinate; do the same so the geometry matches scoring
        gpos = {g: i for i, g in enumerate(genes)}
        drop = {gpos[t] for t in tg if t in gpos}
        mask = np.asarray([i for i in range(len(genes)) if i not in drop])
        data[name] = dict(X=X, E=E[:, mask], mask=mask, targets=tg,
                          sd=X[:, mask].std(0).astype(np.float64))
        print(f"[{name}] {X.shape[0]} control cells, {E.shape[0]} targets, "
              f"{len(mask)} scored genes")

    print(f"\ncomputing bases up to k={max(KS)} ...", flush=True)
    for name, d in data.items():
        d["V"] = top_basis(d["X"][:, d["mask"]], max(KS))
        print(f"[{name}] control basis {d['V'].shape}")

    out = {}
    for name, d in data.items():
        other = [o for o in data if o != name][0]
        E = d["E"]
        R = E - E.mean(0, keepdims=True)
        Ro = top_basis(R, max(KS))               # oracle basis from the residual itself
        rows = []
        print(f"\n{'=' * 96}\n{name}: fraction of TARGET-SPECIFIC residual energy captured")
        print(f"  {'k':>5} {'own ctrl':>10} {'other ctrl':>11} {'random':>9} "
              f"{'oracle':>9} {'own/other':>10} {'own/random':>11}")
        for k in KS:
            c_own = captured(R, d["V"][:k])
            c_oth = captured(R, data[other]["V"][:k])
            c_rnd = float(np.mean([captured(R, random_basis(d["sd"], k, rng))
                                   for _ in range(3)]))
            c_orc = captured(R, Ro[:k])
            rows.append(dict(k=k, own=c_own, other=c_oth, random=c_rnd, oracle=c_orc))
            print(f"  {k:>5} {c_own:10.4f} {c_oth:11.4f} {c_rnd:9.4f} {c_orc:9.4f} "
                  f"{c_own / max(c_oth, 1e-12):10.3f} {c_own / max(c_rnd, 1e-12):11.3f}")
        # the same table for the SHARED mean, which is the part that cannot earn score
        m = E.mean(0)
        mrow = []
        for k in KS:
            mrow.append(dict(k=k,
                             own=captured(m[None, :], d["V"][:k]),
                             other=captured(m[None, :], data[other]["V"][:k])))
        print(f"  shared mean response, own basis: " +
              " ".join(f"k{r['k']}={r['own']:.3f}" for r in mrow))
        out[name] = dict(residual=rows, mean=mrow, other=other,
                         n_targets=int(E.shape[0]), n_genes=int(E.shape[1]),
                         n_cells=int(d["X"].shape[0]))

    with open(os.path.join(WORK, "subspace.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'subspace.json')}")

    print(f"\n{'=' * 96}\nREADING")
    for name, o in out.items():
        r50 = [r for r in o["residual"] if r["k"] == 50][0]
        print(f"  {name}: at k=50 own {r50['own']:.4f} vs other {r50['other']:.4f} "
              f"vs random {r50['random']:.4f} vs oracle {r50['oracle']:.4f}")
        if r50["own"] <= r50["other"] * 1.10:
            print(f"    -> own basis is NOT better than another line's: the subspace is generic "
                  f"co-expression, NOT context conditioning. Hypothesis REJECTED on this line.")
        else:
            print(f"    -> own basis beats the other line's by "
                  f"{100 * (r50['own'] / r50['other'] - 1):.1f}%: context-specific structure "
                  f"is present and worth a transfer test.")


if __name__ == "__main__":
    sys.exit(main())
