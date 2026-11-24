"""177 Q0 — the ZERO-EFFECT NULL for the resampling emission law, on the challenge's own controls.

WHY THIS RUNS FIRST
-------------------
The whole experiment rests on one property: **a resampling emitter has no median artefact.**
Real control cells carrying a predicted fold change are exchangeable with the reference controls
at zero effect, so `Fbar = 0.5` identically rather than asymptotically, and adding SYMMETRIC noise
on top of them can then only *reduce* |z| -- which is the thing 171 section 12c proved is
impossible through the quasi law, whose two noise channels (even allocator at s->0, zero-inflated
gamma-Poisson at s->1) both shift the MEDIAN and therefore MANUFACTURE significance.

This is exactly 157's `zero_effect.py` surface (`data/challenge/context_{A,B,C}.h5ad`, the only
place the production emission depth and the production control distribution both live), with two
changes: the emission goes through the UNMODIFIED production `emit_counts` rather than a closed
form, and the resampling laws are added.

KNOWN-ANSWER ROWS (AGENTS.md 11.1 rule 10) -- from `experiments/157-reach-order/zero_effect.log`:
  * quasi, flat spread 0.25 / phi 4   ->  ctx A: 19.0% significant, 87.9% of them UP
  * quasi, flat spread 1.00 / phi 4   ->  ctx A: 78.9% significant,  0.0% UP  (the median-collapse)
  * quasi, flat spread 1.00 / phi 0   ->  ctx A:  7.7% significant, 100% UP
  * quasi, flat spread 0.25 / phi 0   ->  ctx A: 44.9% significant, 100% UP
157 used a closed form with the threshold |z| > 2.8; that column is reproduced here alongside BH.

THE COORDINATOR'S CONFOUND is also carried as a known-answer row (`rs_100k`): rescaling each real
control cell from its own ~20k library to the harness's nominal `EMIT_L = 100_000` before
resampling, which produced 10,776/10,779 genes significant in `work/realize177.log`.  It must
reproduce that pathology here, or this null is not measuring what it claims to.

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
import vcc2026.donor_transfer as _dt  # noqa: E402

N, PHI, KEEP_TOP = 400, 4.0, 8000
MIN_CPM = 5.0
ZTHR = 2.8               # 157's threshold: BH-0.05 over ~10k genes needs |z| ~ 2.8
CTXS = os.environ.get("CTX", "A").split(",")


# --------------------------------------------------------------------------- fast exact repair
def _fast_repair_totals(block, totals, rng):
    """Vectorised equivalent of `donor_transfer._repair_totals` (same contract: each gene's column
    total forced to `totals` by moving single counts between cells; removals uniform over the
    cells that still hold a count, exactly as the original's `rng.choice(nz)` loop)."""
    out = np.asarray(block, dtype=np.int64).copy()
    n = out.shape[0]
    resid = np.asarray(totals, dtype=np.int64) - out.sum(0)
    for j in np.flatnonzero(resid > 0):
        np.add.at(out[:, j], rng.integers(0, n, int(resid[j])), 1)
    for j in np.flatnonzero(resid < 0):
        d = int(-resid[j])
        col = out[:, j]
        while d > 0:
            nz = np.flatnonzero(col > 0)
            if nz.size == 0:
                break
            take = rng.choice(nz, size=min(d, nz.size), replace=False)
            col[take] -= 1
            d -= take.size
        out[:, j] = col
    return out


_dt._repair_totals = _fast_repair_totals
from vcc2026.donor_transfer import emit_counts  # noqa: E402


# --------------------------------------------------------------------------- noise families
def jitter(C, sigma, kind, rng):
    """Per-gene symmetric spread applied to the donor cells BEFORE the fold.

    `emit_counts(mode='resample_mix')` renormalises each cell to its own library and then pins
    every gene's column total exactly, so a mean-1 kernel here changes only the cell-to-cell
    spread -- the pseudobulk is restored by the repair either way.  The families differ in what
    they do to the MEDIAN, which is the only thing a rank-sum reads:

      "mul"   M = 1 + sigma*u,  u ~ U(-sqrt3, sqrt3).   mean 1 AND median 1, symmetric in the
              LINEAR domain -- the domain the rank test lives in after its monotone CPM map.
      "logn"  M = exp(sigma*eps - sigma^2/2).  mean 1, median exp(-sigma^2/2) < 1: it BUYS its
              mean-preservation with a downward median shift, i.e. exactly the pathology that
              killed the quasi schedule.  Carried to show the difference is real.
      "pair"  additive symmetric transfer: cells are paired at random and `d` counts move from
              one to the other.  Preserves each gene's column total EXACTLY (no repair needed)
              and is symmetric by construction; `d` is clipped so counts stay >= 0.
    """
    n, G = C.shape
    s = np.asarray(sigma, dtype=np.float64)
    if kind == "none" or not np.any(s > 0):
        return C
    if kind == "mul":
        M = 1.0 + s[None, :] * rng.uniform(-np.sqrt(3.0), np.sqrt(3.0), size=(n, G))
        return C * np.clip(M, 0.0, None)
    if kind == "logn":
        M = np.exp(s[None, :] * rng.standard_normal((n, G)) - 0.5 * s[None, :] ** 2)
        return C * M
    if kind == "pair":
        out = C.copy()
        m = C.mean(0)                                   # per-gene mean count per cell
        perm = rng.permutation(n)
        i, j = perm[: n // 2], perm[n // 2: 2 * (n // 2)]
        d = np.abs(rng.standard_normal((n // 2, G))) * (s * m)[None, :]
        d = np.minimum(d, out[j, :])                    # never take a cell below zero
        out[i, :] += d
        out[j, :] -= d
        return out
    raise ValueError(kind)


# --------------------------------------------------------------------------- Wilcoxon vs controls
class Ref:
    """Per-gene control CPM distribution, for exact mid-rank Mann-Whitney of 400 emitted cells."""

    def __init__(self, D_csc, gated, m):
        self.m = m
        self.nz = [np.sort(D_csc.data[D_csc.indptr[g]:D_csc.indptr[g + 1]]) for g in gated]
        self.n_zero = [m - a.size for a in self.nz]
        self.fac = np.sqrt(12.0 * N * m / (m + N + 1.0))

    def z(self, E):
        """E: (N, n_gated) emitted CPM."""
        out = np.empty(E.shape[1])
        for k in range(E.shape[1]):
            e = E[:, k]
            a, nz0 = self.nz[k], self.n_zero[k]
            lo = np.searchsorted(a, e, "left")
            hi = np.searchsorted(a, e, "right")
            below = np.where(e == 0.0, 0, nz0 + lo)
            eq = np.where(e == 0.0, nz0, hi - lo)
            out[k] = ((below + 0.5 * eq) / self.m).mean()
        return (out - 0.5) * self.fac


def spread_vec(G, gated_idx, profile, cpm, tail, ktop, kspread):
    """`emission_spread_vector`: tail everywhere, `keep_top` genes by |log2 fold| at keep_spread.
    At zero effect the fold is 1 for every gene, so the keep block is the lexsort's first
    `ktop` gated genes -- which is what production would do on a zero prediction."""
    sv = np.full(G, float(tail))
    if ktop <= 0:
        return np.full(G, float(kspread))
    pcpm = 1e6 * profile / profile.sum()
    eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
    score = np.full(G, -np.inf)
    score[gated_idx] = eff[gated_idx]
    order = np.lexsort((np.arange(G), -score))[:min(ktop, gated_idx.size)]
    sv[order] = float(kspread)
    return sv


def main() -> None:
    t0 = time.time()
    print(f"{'ctx':>4}{'law':>34}{'Fbar med':>10}{'|z| med':>9}{'sig|z|>2.8':>11}"
          f"{'rate':>8}{'up':>8}{'BH n_pred':>11}{'BH rate':>9}{'BH up':>8}")
    for ctx in CTXS:
        A = ad.read_h5ad(f"{ROOT}/data/challenge/context_{ctx}.h5ad")
        X = A.X.tocsr().astype(np.float64)
        G = X.shape[1]
        lib = np.asarray(X.sum(1)).ravel()
        lib[lib == 0] = 1.0
        m = X.shape[0]
        D = (sp.diags(1e6 / lib) @ X).tocsc()
        cpm = np.asarray(D.mean(0)).ravel()                 # mean per-cell CPM (the scorer's gate)
        gated = np.flatnonzero(cpm > MIN_CPM)
        colsum = np.asarray(X.sum(0)).ravel()
        p = colsum / colsum.sum()                           # POOLED control profile = zero effect
        mprof = cpm / cpm.sum()                             # the mean-per-cell-CPM profile
        ref = Ref(D, gated, m)
        del D
        rng0 = np.random.default_rng(11)
        don = np.sort(rng0.choice(m, N, replace=False))
        Craw = np.asarray(X[don].todense(), dtype=np.float64)
        dlib = Craw.sum(1)
        print(f"  [ctx {ctx}] {m} controls, {gated.size} gated, median lib {np.median(lib):.0f}, "
              f"donor median lib {np.median(dlib):.0f} [{time.time()-t0:.0f}s]", flush=True)

        # per-gene control CV in linear CPM space -- the schedule's key
        e1 = cpm
        e2 = np.asarray(sp.csc_matrix((sp.diags(1e6 / lib) @ X)).power(2).mean(0)).ravel()
        cv = np.sqrt(np.maximum(e2 - e1 ** 2, 0.0)) / np.maximum(e1, 1e-12)
        mcell = colsum / m * 0 + (p * np.median(lib))        # counts/cell at the emission depth

        def row(label, E):
            z = ref.z(E)
            sig = np.abs(z) > ZTHR
            pv = 2 * norm.sf(np.abs(z))
            q = false_discovery_control(pv, method="bh")
            bh = q < 0.05
            print(f"{ctx:>4}{label:>34}{np.median(z / ref.fac + 0.5):>10.4f}"
                  f"{np.median(np.abs(z)):>9.2f}{int(sig.sum()):>11d}{sig.mean():>8.1%}"
                  f"{((z[sig] > 0).mean() if sig.any() else 0):>8.1%}"
                  f"{int(bh.sum()):>11d}{bh.mean():>9.1%}"
                  f"{((z[bh] > 0).mean() if bh.any() else 0):>8.1%}", flush=True)

        def emit(mode, seed, **kw):
            rng = np.random.default_rng(seed)
            v, c, ip = emit_counts(p, kw.pop("libs"), 0.1, rng, mode=mode, **kw)
            B = np.zeros((N, G))
            for r in range(N):
                sl = slice(ip[r], ip[r + 1])
                B[r, c[sl]] = v[sl]
            ss = B.sum(1)
            ss[ss == 0] = 1.0
            return B[:, gated] * (1e6 / ss[:, None])

        # ------------------------------------------------------------ quasi known-answer rows
        flat = np.full(G, 1.0)
        for s_, phi_ in ((0.25, PHI), (1.00, PHI), (1.00, 0.0), (0.25, 0.0)):
            row(f"KA quasi flat {s_:.2f}/phi{phi_:g}",
                emit("quasi", 101, libs=np.full(N, float(np.median(lib))),
                     spread=np.full(G, s_), overdispersion=phi_))
        sv = spread_vec(G, gated, p, cpm, 1.0, KEEP_TOP, 0.25)
        row("KA quasi SHIPPED tail1/kt8000/ks.25",
            emit("quasi", 101, libs=np.full(N, float(np.median(lib))),
                 spread=sv, overdispersion=PHI))

        # ------------------------------------------------------------ resampling laws
        pm = np.zeros(G, dtype=bool)          # all genes resampled

        def rs(label, C, libs, cprof, seed=202):
            row(label, emit("resample_mix", seed, libs=libs, spread=flat, overdispersion=PHI,
                            control_counts=C, control_profile=cprof, param_mask=pm))

        # the coordinator's confound, as a known-answer row
        rngz = np.random.default_rng(5)
        C100 = np.floor(Craw * (100_000.0 / dlib)[:, None] + rngz.random(Craw.shape))
        rs("KA CONFOUND rs 100k rescale", C100, np.full(N, 100_000.0), mprof)
        del C100
        rs("rs raw, ctrl_prof = mean-CPM", Craw, dlib, mprof)
        rs("rs raw, ctrl_prof = POOLED", Craw, dlib, p)
        rs("rs raw, POOLED, seed b", Craw, dlib, p, seed=909)
        # totals-matched donors: bring the block's total to the emitted panel total by BINOMIAL
        # thinning (a depth change, which is what real shallower sequencing is) instead of
        # letting the exact-total repair peel it off one count per cell at a time
        Ltot = float(np.median(dlib)) * N
        thin = min(1.0, Ltot / dlib.sum())
        Cth = rngz.binomial(Craw.astype(np.int64), thin).astype(np.float64)
        rs(f"rs binom-thinned x{thin:.4f}", Cth, dlib, p)
        del Cth

        # ------------------------------------------------------------ noise families
        sig_flat = np.full(G, 0.0)
        sig_flat[gated] = 0.60
        for kind in ("mul", "logn", "pair"):
            rs(f"rs + {kind} jitter sigma 0.60",
               jitter(Craw, sig_flat, kind, np.random.default_rng(31)), dlib, p)
        sig_hi = np.full(G, 0.0)
        sig_hi[gated] = 1.00
        for kind in ("mul", "pair"):
            rs(f"rs + {kind} jitter sigma 1.00",
               jitter(Craw, sig_hi, kind, np.random.default_rng(32)), dlib, p)

        # the real schedule shape: sigma increasing with expression (INBOX-4 section 3b's
        # direction), floored at 15 CPM where the control's own Poisson noise is unmatchable
        sig_sched = np.zeros(G)
        ok = gated[(cpm[gated] > 15.0)]
        T = float(np.percentile(cv[ok], 60.0))
        sig_sched[ok] = np.clip(np.sqrt(np.maximum(T ** 2 - cv[ok] ** 2, 0.0)), 0.0, 1.2)
        print(f"  [ctx {ctx}] schedule: T(cv q60) {T:.3f}, sigma med {np.median(sig_sched[ok]):.3f}, "
              f"max {sig_sched.max():.3f}, nonzero {(sig_sched > 0).sum()}", flush=True)
        for kind in ("mul", "pair"):
            rs(f"rs + {kind} EXPR-COMP schedule",
               jitter(Craw, sig_sched, kind, np.random.default_rng(33)), dlib, p)

        del A, X, ref, Craw
        print(f"  [ctx {ctx}] done [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
