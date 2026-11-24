"""172 — streaming donor access, and the per-donor gene statistics the key is built from.

DONOR PROCESSING MATCHES THE BUILD.  `scripts/make_submission_pm_cd4.py:2336` constructs
`MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)` and a15 ships `beta = 0.0`
(`submissions/a15-cd4f010.vcc.meta.json`), so the transferred effect is

    consensus(t)  =  alpha * global * mean_over_carrying_donors( renorm_d * s_d[t] )
    s_d           =  delta_d - mean_over_all_perts(delta_d)      (`Donor.decompose`)

on the recipient's gene axis, zeros off the donor's own axis.  `alpha * global` is one scalar
per build, so every feature below is either scale-free or explicitly re-centred per fold.

MEMORY.  Orion's `delta.npy` is 2.7 GB and there are two of them.  Everything here is
memory-mapped and accumulated in row chunks; peak RSS stays under ~1 GB (measured), against
the 5.5 GB a pair of `Donor.load()` + `decompose()` calls would cost.  Two sibling agents are
running concurrently and the standing cap is 10 GB.

The per-gene panel statistics (marginal responsiveness, across-panel sd, the fraction of
perturbations in which a gene is a top responder) are properties of the DONOR ALONE, so they
are computed once on the donor's own axis and cached in `work/donorstats.<name>.npz`.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DONORS = os.path.join(REPO, "data", "donors")
WORK = os.path.join(HERE, "work")
CHUNK = 1000
TOP_FRAC = 0.05          # "is this gene a top responder for this perturbation"


def _lines(path: Path) -> np.ndarray:
    """`Donor.load`'s reader, verbatim, so indices line up with production's."""
    return np.asarray(path.read_text().split("\n"))


class DonorHandle:
    """Memory-mapped access to one donor screen, plus its cached gene statistics."""

    def __init__(self, name: str):
        self.name = name
        d = Path(DONORS) / name
        self.perts = _lines(d / "perts.txt")
        self.genes = _lines(d / "genes.txt")
        self.b_ctrl = np.load(d / "b_ctrl.npy").astype(np.float64)
        self.delta = np.load(d / "delta.npy", mmap_mode="r")
        assert self.delta.shape == (len(self.perts), len(self.genes)), self.delta.shape
        self.pindex = {str(p): i for i, p in enumerate(self.perts)}
        e = np.expm1(self.b_ctrl)
        self.cpm = 1e6 * e / e.sum()
        st = self._stats()
        self.m = st["m"]                       # Donor.decompose()'s gene mean
        self.marg_abs = st["marg_abs"]         # mean_p |s[p, g]|
        self.marg_sd = st["marg_sd"]           # sd_p s[p, g]
        self.resp_frac = st["resp_frac"]       # P(g in the top 5% of |s| for a perturbation)
        self.scale = float(st["scale"])        # MultiDonorModel._scale's input, own axis
        self.n_perts = len(self.perts)

    # ---- cached per-gene statistics ---------------------------------------------------

    def _stats(self) -> dict:
        cache = os.path.join(WORK, f"donorstats.{self.name}.npz")
        if os.path.exists(cache):
            return dict(np.load(cache))
        n, G = self.delta.shape
        acc = np.zeros(G, dtype=np.float64)
        for i in range(0, n, CHUNK):
            acc += np.asarray(self.delta[i:i + CHUNK], dtype=np.float64).sum(0)
        m = acc / n
        s_abs = np.zeros(G, dtype=np.float64)
        s_sq = np.zeros(G, dtype=np.float64)
        resp = np.zeros(G, dtype=np.float64)
        norms = []
        k = G - max(1, int(round(TOP_FRAC * G)))
        for i in range(0, n, CHUNK):
            S = np.asarray(self.delta[i:i + CHUNK], dtype=np.float32) - m.astype(np.float32)
            A = np.abs(S)
            s_abs += A.sum(0, dtype=np.float64)
            s_sq += np.square(A, dtype=np.float32).sum(0, dtype=np.float64)
            thr = np.partition(A, k, axis=1)[:, k]
            resp += (A >= thr[:, None]).sum(0)
            norms.append(np.linalg.norm(S, axis=1))
        norms = np.concatenate(norms)
        out = dict(m=m, marg_abs=s_abs / n, marg_sd=np.sqrt(np.maximum(s_sq / n, 0.0)),
                   resp_frac=resp / n, scale=np.float64(np.median(norms)))
        os.makedirs(WORK, exist_ok=True)
        np.savez_compressed(cache, **out)
        return out

    # ---- projection onto a recipient axis ---------------------------------------------

    def project(self, axis_genes: np.ndarray):
        """(tcols, dcols): recipient columns and the donor columns that feed them."""
        gpos = {str(g): j for j, g in enumerate(self.genes)}
        tcols, dcols = [], []
        for t, g in enumerate(axis_genes):
            j = gpos.get(str(g))
            if j is not None:
                tcols.append(t)
                dcols.append(j)
        return np.asarray(tcols, np.int64), np.asarray(dcols, np.int64)

    def rows(self, targets, axis_genes, tcols=None, dcols=None):
        """`s_d[t]` for each target, on the recipient axis; zeros off the donor's own axis."""
        if tcols is None:
            tcols, dcols = self.project(axis_genes)
        S = np.zeros((len(targets), len(axis_genes)), dtype=np.float32)
        carried = np.zeros(len(targets), dtype=bool)
        for r, t in enumerate(targets):
            i = self.pindex.get(str(t))
            if i is None:
                continue
            carried[r] = True
            S[r, tcols] = (np.asarray(self.delta[i], dtype=np.float64) - self.m)[dcols]
        return S, carried

    def raw_rows(self, targets, axis_genes, tcols=None, dcols=None):
        """`delta_d[t]` (NOT mean-centred) -- the analogue of the scorer's `real_lfc`."""
        if tcols is None:
            tcols, dcols = self.project(axis_genes)
        S = np.full((len(targets), len(axis_genes)), np.nan, dtype=np.float32)
        carried = np.zeros(len(targets), dtype=bool)
        for r, t in enumerate(targets):
            i = self.pindex.get(str(t))
            if i is None:
                continue
            carried[r] = True
            S[r, tcols] = np.asarray(self.delta[i], dtype=np.float64)[dcols]
        return S, carried


def overlap_audit(h1test_targets) -> None:
    """AGENTS.md trap #1.  `h1_2025` covers 85/85 h1-test targets and is pure oracle leakage."""
    tset = {str(t) for t in h1test_targets}
    print(f"{'donor':<16}{'perts':>8}{'genes':>8}{'h1-test targets covered':>26}")
    for name in sorted(os.listdir(DONORS)):
        d = Path(DONORS) / name
        if not (d / "perts.txt").exists():
            continue
        perts = {str(p) for p in _lines(d / "perts.txt")}
        genes = _lines(d / "genes.txt")
        cov = len(tset & perts)
        flag = "   <-- BANNED (oracle)" if cov >= len(tset) else ""
        print(f"{name:<16}{len(perts):>8}{len(genes):>8}"
              f"{cov:>18}/{len(tset):<4}{100*cov/len(tset):>5.1f}%{flag}")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, HERE)
    from harness172 import H1Test
    T = H1Test()
    overlap_audit(T.perts)
    print()
    for n in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        h = DonorHandle(n)
        print(f"{n:<16} perts {h.n_perts:>6}  genes {len(h.genes):>6}  "
              f"scale {h.scale:.4f}  marg_abs median {np.median(h.marg_abs):.4f}  "
              f"resp_frac median {np.median(h.resp_frac):.4f}")
    print()
    print("pairwise perturbation overlap (LODO folds use the two donors NOT held out):")
    P = {n: {str(p) for p in _lines(Path(DONORS) / n / "perts.txt")}
         for n in ("k562_gwps", "orion_hct116", "orion_hek293t")}
    for held in P:
        src = [n for n in P if n != held]
        both = P[src[0]] | P[src[1]]
        print(f"  held-out {held:<16} targets also in a source: "
              f"{len(P[held] & both):>6}/{len(P[held]):<6}  "
              f"(in BOTH sources: {len(P[held] & P[src[0]] & P[src[1]])})")
