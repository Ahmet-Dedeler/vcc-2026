"""167 — does a target's OWN expression in the recipient predict whether its effect is real?

THE OBSERVATION
---------------
Production builds ONE context-independent effect matrix and applies it to A/B/C. But a panel
target's own expression varies enormously between the three contexts: 69% of the 300 targets
swing more than 2x, 30% more than 4x, and the tail is extreme -- VIM is 1626 / 6132 / **2.02**
CPM in A/B/C, S100A11 is 0.12 / 525 / 2353, KIF21B 107 / 0.22 / 2.00.

CRISPRi represses transcription. A gene that is already silent cannot be repressed further, so
its true downstream effect should be ~0 -- yet we predict a full-magnitude effect for it. That is
a context-specific error we make by construction, and unlike a fitted hyperparameter it is a
mechanism, so it should transfer to the unseen D/E/F panel.

WHY NAIVE ||effect|| WOULD LIE
-----------------------------
The obvious test -- correlate the target's expression with `||true effect||` -- is confounded:
with ~150 cells per target the measured effect of a NULL perturbation still has a large norm,
made entirely of sampling noise. A silent target and a noisy one look identical.

So this measures **reliability** instead: split each target's cells into two disjoint halves,
pseudobulk each half against a DISJOINT half of the control (sharing one control between halves
correlates their noise), and take the cosine between the two halves' effect vectors. That is the
reproducible fraction of the effect, and it goes to zero for a target with no real response no
matter how big its noisy norm is.

Endpoint: reliability as a function of the target's own control CPM, on three surfaces with
truth. If low-expression targets are unreliable, the lever is real and buildable -- the target's
own expression in the recipient is available at deployment from the control cells alone.
"""
from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

TS = 5.0e4
SEED = 20260902
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
CTRL = "non-targeting"
MIN_CELLS = 40                      # need >=20 per half for a meaningful split

SURFACES = {
    "orion_hct116": "data/harness/orion-hct116-thinned.h5ad",
    "orion_hek293t": "data/harness/orion-hek293t-thinned.h5ad",
    "h1_test": "data/harness/h1-test.h5ad",
}


def _labels(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        c = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                        for x in o["categories"][:]])
        return c[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def load_counts(path, chunk=4000):
    """Full count matrix as CSR plus labels and the gene axis."""
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        genes = _labels(h, "var/_index")
        lab = _labels(h, "obs/target_gene")
        indptr = h["X/indptr"][:]
        blocks = []
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            s, e = int(indptr[r0]), int(indptr[r1])
            blocks.append(sp.csr_matrix(
                (h["X/data"][s:e].astype(np.float32), h["X/indices"][s:e],
                 indptr[r0:r1 + 1] - indptr[r0]), shape=(r1 - r0, g)))
        X = sp.vstack(blocks).tocsr()
    return X, lab, genes


def profile(counts_sum):
    t = counts_sum.sum()
    return np.log1p(TS * counts_sum / max(t, 1e-12))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(SEED)
    out = {}
    for name, rel in SURFACES.items():
        X, lab, genes = load_counts(os.path.join(REPO, rel))
        gpos = {g: i for i, g in enumerate(genes)}
        ctrl_rows = np.flatnonzero(lab == CTRL)
        ctot = np.asarray(X[ctrl_rows].sum(0)).ravel().astype(np.float64)
        cpm_all = 1e6 * ctot / ctot.sum()

        # disjoint control halves
        cperm = rng.permutation(ctrl_rows)
        ch = len(cperm) // 2
        cA = np.asarray(X[cperm[:ch]].sum(0)).ravel().astype(np.float64)
        cB = np.asarray(X[cperm[ch:]].sum(0)).ravel().astype(np.float64)
        bA, bB = profile(cA), profile(cB)

        targets = sorted({t for t in set(lab) if t != CTRL})
        # pds drops every panel-target coordinate; do the same
        drop = {gpos[t] for t in targets if t in gpos}
        mask = np.asarray([i for i in range(len(genes)) if i not in drop])

        recs = []
        for t in targets:
            rows = np.flatnonzero(lab == t)
            if len(rows) < MIN_CELLS or t not in gpos:
                continue
            perm = rng.permutation(rows)
            h = len(perm) // 2
            e1 = profile(np.asarray(X[perm[:h]].sum(0)).ravel().astype(np.float64))[mask] - bA[mask]
            e2 = profile(np.asarray(X[perm[h:]].sum(0)).ravel().astype(np.float64))[mask] - bB[mask]
            n1, n2 = np.linalg.norm(e1), np.linalg.norm(e2)
            r = float(e1 @ e2 / max(n1 * n2, 1e-300))
            # on-target knockdown, for context: does the guide actually repress the gene?
            full = np.asarray(X[rows].sum(0)).ravel().astype(np.float64)
            fcpm = 1e6 * full / full.sum()
            kd = float(np.log2((fcpm[gpos[t]] + 0.1) / (cpm_all[gpos[t]] + 0.1)))
            recs.append(dict(target=t, cpm=float(cpm_all[gpos[t]]), rel=r,
                             n_cells=int(len(rows)), kd_log2=kd,
                             norm=float(0.5 * (n1 + n2))))
        out[name] = recs
        c = np.array([x["cpm"] for x in recs])
        rl = np.array([x["rel"] for x in recs])
        kd = np.array([x["kd_log2"] for x in recs])
        nm = np.array([x["norm"] for x in recs])
        print(f"\n{'=' * 100}\n{name}: {len(recs)} targets with >= {MIN_CELLS} cells, "
              f"{len(mask)} scored genes, {len(ctrl_rows)} control cells")
        print(f"  target's own CPM: min {c.min():.2f} median {np.median(c):.1f} max {c.max():.0f}")
        print(f"  split-half reliability: mean {rl.mean():+.4f} median {np.median(rl):+.4f}")
        lo = np.log10(c + 0.1)
        print(f"  Spearman(own CPM, reliability) = "
              f"{_spearman(c, rl):+.4f}   Pearson(log10 CPM, rel) = {np.corrcoef(lo, rl)[0,1]:+.4f}")
        print(f"  Spearman(own CPM, ||effect||)  = {_spearman(c, nm):+.4f}   "
              f"<- the CONFOUNDED statistic, for comparison")
        print(f"  Spearman(own CPM, on-target log2FC) = {_spearman(c, kd):+.4f}")
        # binned table -- the readable form
        qs = np.quantile(c, [0, .2, .4, .6, .8, 1.0])
        print(f"\n  {'own-CPM bin':>22s} {'n':>4s} {'med CPM':>9s} {'mean rel':>9s} "
              f"{'med rel':>8s} {'med ||eff||':>11s} {'med on-tgt log2FC':>18s}")
        for i in range(5):
            m = (c >= qs[i]) & (c <= qs[i + 1]) if i == 4 else (c >= qs[i]) & (c < qs[i + 1])
            if not m.any():
                continue
            print(f"  {f'{qs[i]:.1f}-{qs[i+1]:.1f}':>22s} {int(m.sum()):>4d} "
                  f"{np.median(c[m]):>9.1f} {rl[m].mean():>9.4f} {np.median(rl[m]):>8.4f} "
                  f"{np.median(nm[m]):>11.3f} {np.median(kd[m]):>18.3f}")
        for thr in (5.0, 10.0, 20.0):
            m = c < thr
            if m.sum() >= 3:
                print(f"  targets below {thr:>4.0f} CPM: n={int(m.sum()):>3d}  mean rel "
                      f"{rl[m].mean():+.4f}  vs  above: {rl[~m].mean():+.4f}")

    with open(os.path.join(WORK, "reliability.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'reliability.json')}")


def _spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


if __name__ == "__main__":
    sys.exit(main())
