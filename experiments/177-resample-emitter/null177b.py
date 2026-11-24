"""177 Q0b — how much SYMMETRIC noise the resampling law can carry before its null breaks,
and the calibration table the delivery schedule is built from.

`null177.py` established the two facts this depends on:
  * the resampling law's zero-effect null is CLEAN (1/9929 BH-significant on challenge ctx A,
    against 2,818 for a flat quasi 0.25/phi4 and 4,302 for the shipped mixture), and
  * MULTIPLICATIVE jitter destroys it (`mul` sigma 0.60 -> 72.2% BH, 100% DOWN), because at a
    pinned column total a multiplicative spread lowers the MEDIAN, which is the only thing a
    rank-sum reads.  That is 171 section 12c's mechanism reappearing in a new costume.

Only the ADDITIVE SYMMETRIC PAIR TRANSFER survived (sigma 0.60 -> 0.1% BH), and it broke at
sigma 1.00 (39.2% BH, 100% DOWN).  The suspected bound is the ZERO FLOOR: a cell cannot give
away more counts than it holds, so once `sigma * m` approaches the gene's own counts/cell the
down-moves clip, the transfer stops being symmetric, and the median falls.  If that is the
mechanism the bound must be a function of the gene's COUNTS PER CELL, not of sigma alone, and
high-expression genes -- exactly the ones the expression-compensating schedule needs to noise --
should tolerate far more.

This measures that surface directly and emits `sigma_cap(m)`, which the schedule then respects.

KNOWN-ANSWER ROW: the coordinator's confound, reproduced exactly (donor cells rescaled to
`EMIT_L = 100_000` while the pinned panel total stays at the artifact's ~20k/cell depth, so the
exact-total repair peels ~80% of every column one count per cell at a time).  `work/realize177.log`
recorded 10,776/10,779 genes significant; this must reproduce that, or the null is not measuring
the thing that went wrong.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os
import sys
import time

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, HERE)
from null177 import Ref, _fast_repair_totals, jitter  # noqa: E402
import vcc2026.donor_transfer as _dt  # noqa: E402

_dt._repair_totals = _fast_repair_totals
from vcc2026.donor_transfer import emit_counts  # noqa: E402

N, PHI, MIN_CPM, ZTHR = 400, 4.0, 5.0, 2.8
SURF = os.environ.get("SURF", "h1test")


def load_surface():
    if SURF == "h1test":
        R = ad.read_h5ad(f"{ROOT}/data/harness/h1-test.h5ad")
        lab = np.asarray(R.obs["target_gene"].astype(str).values, dtype="U64")
        X = R.X[np.flatnonzero(lab == "non-targeting")].tocsr().astype(np.float64)
        del R
        return X
    A = ad.read_h5ad(f"{ROOT}/data/challenge/context_{SURF}.h5ad")
    X = A.X.tocsr().astype(np.float64)
    del A
    return X


def antipair(C, sigma, rng):
    """Exact-total, exactly i<->j symmetric variance AMPLIFIER.

    For a random pairing, `v'_i = v_i + t*(v_i - v_j)`, `v'_j = v_j + t*(v_j - v_i)`.  The pair
    SUM is invariant (so every gene's column total is preserved exactly, with no repair), the
    pair DIFFERENCE is multiplied by `(1 + 2t)`, and the construction is invariant under
    swapping the two cells -- so no direction is privileged.  Clipped at zero, which is where
    the symmetry (and hence the null) is expected to fail.
    """
    n, G = C.shape
    out = C.copy()
    perm = rng.permutation(n)
    h = n // 2
    i, j = perm[:h], perm[h:2 * h]
    t = np.asarray(sigma, dtype=np.float64)[None, :]
    a, b = out[i, :], out[j, :]
    na = np.clip(a + t * (a - b), 0.0, None)
    nb = a + b - na                                  # exact pair-sum preservation
    bad = nb < 0.0
    if bad.any():
        na = np.where(bad, a + b, na)
        nb = np.where(bad, 0.0, nb)
    out[i, :], out[j, :] = na, nb
    return out


def fracpair(C, sigma, rng):
    """Pair transfer of a FRACTION of the loser's own count -- clips by construction, never
    forces a negative, so it isolates whether clipping is what breaks the pair kernel."""
    n, G = C.shape
    out = C.copy()
    perm = rng.permutation(n)
    h = n // 2
    i, j = perm[:h], perm[h:2 * h]
    u = rng.uniform(0.0, 1.0, size=(h, G)) * np.clip(np.asarray(sigma)[None, :], 0.0, 1.0)
    d = u * out[j, :]
    out[i, :] += d
    out[j, :] -= d
    return out


def main() -> None:
    t0 = time.time()
    X = load_surface()
    m, G = X.shape
    lib = np.asarray(X.sum(1)).ravel()
    lib[lib == 0] = 1.0
    D = (sp.diags(1e6 / lib) @ X).tocsc()
    cpm = np.asarray(D.mean(0)).ravel()
    gated = np.flatnonzero(cpm > MIN_CPM)
    colsum = np.asarray(X.sum(0)).ravel()
    p = colsum / colsum.sum()
    L = float(np.median(lib))
    ref = Ref(D, gated, m)
    del D
    rng0 = np.random.default_rng(11)
    don = np.sort(rng0.choice(m, N, replace=False))
    Craw = np.asarray(X[don].todense(), dtype=np.float64)
    dlib = Craw.sum(1)
    # counts per cell at the emission depth -- the quantity the zero floor is a function of
    mcell = N * L * p / N
    print(f"[{SURF}] {m} controls x {G} genes, {gated.size} gated, median lib {L:.0f}, "
          f"counts/cell at gated genes: q10 {np.percentile(mcell[gated],10):.2f} "
          f"med {np.median(mcell[gated]):.2f} q90 {np.percentile(mcell[gated],90):.2f} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    e1 = cpm
    e2 = np.asarray((sp.diags(1e6 / lib) @ X).power(2).mean(0)).ravel()
    cv = np.sqrt(np.maximum(e2 - e1 ** 2, 0.0)) / np.maximum(e1, 1e-12)
    del X
    # expression terciles of the gated set, by counts/cell
    q = np.percentile(mcell[gated], [33.3, 66.7])
    tier = np.digitize(mcell[gated], q)
    print(f"[{SURF}] counts/cell tercile edges {q[0]:.2f} / {q[1]:.2f}; "
          f"median control CV per tercile: "
          + " ".join(f"T{k}={np.median(cv[gated][tier==k]):.2f}" for k in range(3)), flush=True)

    pm = np.zeros(G, dtype=bool)
    flat = np.full(G, 1.0)

    def emit(profile, libs, C, cprof, seed):
        rng = np.random.default_rng(seed)
        v, c, ip = emit_counts(profile, libs, 0.1, rng, mode="resample_mix", spread=flat,
                               overdispersion=PHI, control_counts=C, control_profile=cprof,
                               param_mask=pm)
        B = np.zeros((N, G))
        for r in range(N):
            sl = slice(ip[r], ip[r + 1])
            B[r, c[sl]] = v[sl]
        ss = B.sum(1)
        ss[ss == 0] = 1.0
        return B[:, gated] * (1e6 / ss[:, None])

    hdr = (f"{'law':>36}{'|z|med':>8}{'sig':>7}{'rate':>7}{'up':>7}{'BH':>7}{'BHrate':>8}"
           f"{'  rate by counts/cell tercile':>30}")
    print(hdr, flush=True)

    def row(label, E):
        z = ref.z(E)
        sig = np.abs(z) > ZTHR
        bh = false_discovery_control(2 * norm.sf(np.abs(z)), method="bh") < 0.05
        per = " ".join(f"{bh[tier == k].mean():>7.1%}" for k in range(3))
        print(f"{label:>36}{np.median(np.abs(z)):>8.2f}{int(sig.sum()):>7d}{sig.mean():>7.1%}"
              f"{((z[sig] > 0).mean() if sig.any() else 0):>7.1%}{int(bh.sum()):>7d}"
              f"{bh.mean():>8.1%}   {per}", flush=True)
        return z

    # ------------------------------------------------------- KNOWN ANSWER: the confound
    rngz = np.random.default_rng(5)
    C100 = np.floor(Craw * (100_000.0 / dlib)[:, None] + rngz.random(Craw.shape))
    row("KA CONFOUND coord 100k vs 20k pin",
        emit(p, np.full(N, 100_000.0), C100, p, 202))
    del C100

    # ------------------------------------------------------- the clean law
    row("rs raw (no jitter)", emit(p, dlib, Craw, p, 202))
    row("rs raw (no jitter) seed b", emit(p, dlib, Craw, p, 909))

    # ------------------------------------------------------- kernel x sigma grid
    kern = {"pair": lambda C, s, r: jitter(C, s, "pair", r),
            "antipair": antipair, "fracpair": fracpair}
    for kname, fn in kern.items():
        for s in (0.4, 0.8, 1.2, 2.0):
            sv = np.zeros(G)
            sv[gated] = s
            row(f"{kname} sigma {s:.1f} (all gated)",
                emit(p, dlib, fn(Craw, sv, np.random.default_rng(41)), p, 202))

    # ------------------------------------------------------- does the zero floor explain it?
    for s in (0.8, 1.2, 2.0, 3.0):
        for floor in (1.0, 4.0):
            sv = np.zeros(G)
            sel = gated[mcell[gated] >= floor]
            sv[sel] = s
            row(f"pair sigma {s:.1f} only m>={floor:.0f}/cell ({sel.size})",
                emit(p, dlib, jitter(Craw, sv, "pair", np.random.default_rng(41)), p, 202))

    print(f"\n[{SURF}] done [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
