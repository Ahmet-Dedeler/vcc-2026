"""177 Q0c — the SIGNAL response: does symmetric noise on the resampling law move |z| the way
171 section 12c proved it cannot move it on the quasi law?

171 section 12c closed the spread-schedule family with a mechanism: `|z|` is U-shaped in the
quasi law's `s`, minimised at s ~ 0.25, because BOTH of that law's noise channels shift the
MEDIAN (even allocator at s->0, zero-inflated gamma-Poisson at s->1).  "Suppress unearned
significance by adding noise" was therefore unavailable -- past s ~ 0.25 adding noise
MANUFACTURES significance.

`null177.py` / `null177b.py` established that the resampling law's zero-effect null is clean
(0/10,779 BH-significant on h1-test controls) and that the only noise kernel that KEEPS it clean
is `fracpair` -- a pair transfer of a uniform fraction of the loser's own count, which preserves
every gene's column total exactly, is symmetric under swapping the pair, and can never drive a
count negative.  The remaining question, and the one the whole delivery arm rests on, is whether
that kernel actually has a LEVER: does `|z_g|` fall monotonically in the kernel strength, and
does it fall DIFFERENTIALLY enough by expression to flatten the +0.485 expression dependence?

Measured here on real predicted effects (169's recovered `H`, production mass fix, the exact
profile the realisation will emit), against the real h1-test controls, for a handful of targets.

KNOWN-ANSWER ROW: the coordinator's confound, reproduced exactly -- donor cells inflated to
`EMIT_L = 100_000` while the pinned panel total stays at the artifact's own ~20k/cell depth, so
the exact-total repair peels ~80% of every column.  `work/realize177.log` recorded median
`n_pred` 10,776 of 10,779.  If this row does not reproduce that, the diagnosis is wrong.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import importlib.util
import os
import sys
import time

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import false_discovery_control, norm, rankdata

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = f"{ROOT}/experiments/169-sharpening-six-member/work"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, HERE)
from null177 import Ref, _fast_repair_totals  # noqa: E402
import vcc2026.donor_transfer as _dt  # noqa: E402

_dt._repair_totals = _fast_repair_totals
from vcc2026.donor_transfer import emit_counts  # noqa: E402

N, PHI, MIN_CPM = 400, 4.0, 5.0
NT = int(os.environ.get("NT", "6"))          # targets to average over


def fracpair(C, strength, rng, reps=1):
    """Exact-total, pair-symmetric variance amplifier that cannot clip.

    For a random pairing and per-gene strength `s_g` in [0,1]: `d = u*v_j` with
    `u ~ U(0, s_g)`, then `v_i += d`, `v_j -= d`.  The transferred amount is a FRACTION of what
    the losing cell actually holds, so no count is ever driven negative and the kernel never
    has to clip -- which is what breaks the additive kernels (`null177b.py`: an absolute
    pair transfer at sigma 0.8 lands 32.4% of genes BH-significant at ZERO effect, 100% DOWN,
    because the losing half piles up on the zero floor and drags the median under the
    control's).  `reps` independent passes compound the spread past the `s <= 1` ceiling.
    """
    out = np.asarray(C, dtype=np.float64).copy()
    n = out.shape[0]
    s = np.clip(np.asarray(strength, dtype=np.float64), 0.0, 1.0)[None, :]
    if not np.any(s > 0):
        return out
    h = n // 2
    for _ in range(int(reps)):
        perm = rng.permutation(n)
        i, j = perm[:h], perm[h:2 * h]
        d = rng.uniform(0.0, 1.0, size=(h, out.shape[1])) * s * out[j, :]
        out[i, :] += d
        out[j, :] -= d
    return out


def load_builder():
    path = f"{ROOT}/scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_177s", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_177s"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    t0 = time.time()
    E = np.load(os.path.join(WORK, "h1test-effect.npz"), allow_pickle=True)
    H0, b, a_gate = E["H"].astype(np.float64), E["b"].astype(np.float64), E["gate"]
    esums = E["sums"].astype(np.float64)
    etargets = np.asarray([str(t) for t in E["targets"]])
    prod = load_builder()

    R = ad.read_h5ad(f"{ROOT}/data/harness/h1-test.h5ad")
    lab = np.asarray(R.obs["target_gene"].astype(str).values, dtype="U64")
    X = R.X[np.flatnonzero(lab == "non-targeting")].tocsr().astype(np.float64)
    del R
    m, G = X.shape
    lib = np.asarray(X.sum(1)).ravel()
    lib[lib == 0] = 1.0
    L = float(np.median(lib))
    D = (sp.diags(1e6 / lib) @ X).tocsc()
    cpm = np.asarray(D.mean(0)).ravel()
    gated = np.flatnonzero(cpm > MIN_CPM)
    colsum = np.asarray(X.sum(0)).ravel()
    ctrl_rel = colsum / colsum.sum()                     # POOLED: the space `b` lives in
    ref = Ref(D, gated, m)
    del D
    rng0 = np.random.default_rng(11)
    don = np.sort(rng0.choice(m, N, replace=False))
    Craw = np.asarray(X[don].todense(), dtype=np.float64)
    dlib = Craw.sum(1)
    del X
    mcell = L * ctrl_rel                                  # counts/cell at the emission depth
    tq = np.percentile(mcell[gated], [33.3, 66.7])
    tier = np.digitize(mcell[gated], tq)
    print(f"h1-test controls {m} x {G}; {gated.size} gated; median lib {L:.0f}; "
          f"counts/cell tercile edges {tq[0]:.2f}/{tq[1]:.2f} [{time.time()-t0:.0f}s]", flush=True)

    pm = np.zeros(G, dtype=bool)
    flat = np.full(G, 1.0)

    def emit(profile, libs, C, seed):
        rng = np.random.default_rng(seed)
        v, c, ip = emit_counts(profile, libs, 0.1, rng, mode="resample_mix", spread=flat,
                               overdispersion=PHI, control_counts=C, control_profile=ctrl_rel,
                               param_mask=pm)
        B = np.zeros((N, G))
        for r in range(N):
            sl = slice(ip[r], ip[r + 1])
            B[r, c[sl]] = v[sl]
        ss = B.sum(1)
        ss[ss == 0] = 1.0
        return B[:, gated] * (1e6 / ss[:, None])

    # ------------------------------------------------------------------ KA: the confound
    rngz = np.random.default_rng(5)
    C100 = np.floor(Craw * (100_000.0 / dlib)[:, None] + rngz.random(Craw.shape))
    k0 = 0
    h0 = prod.mass_consistent_effect(H0[k0].copy(), b, a_gate, mode="minnorm")
    lin0 = np.clip(np.expm1(b + h0), 0.0, None)
    p0 = lin0 / lin0.sum()
    # the coordinator's exact geometry: cells at 100k, panel total pinned at the artifact's
    # own ~20k/cell depth, so `profile` is scaled by S/(N*EMIT_L)
    S0 = float(esums[k0].sum())
    zc = ref.z(emit(p0 * (S0 / (N * 100_000.0)), np.full(N, 100_000.0), C100, 202))
    qc = false_discovery_control(2 * norm.sf(np.abs(zc)), method="bh")
    print(f"KA CONFOUND (coordinator's geometry): n_pred {(qc < 0.05).sum()}/{gated.size} "
          f"  expect ~10776/10779 from work/realize177.log", flush=True)
    del C100

    # ------------------------------------------------------------------ signal response
    POLICIES = [("sigma 0", None, 0), ("flat 0.50 x1", 0.50, 1), ("flat 1.00 x1", 1.0, 1),
                ("flat 1.00 x2", 1.0, 2), ("flat 1.00 x4", 1.0, 4), ("flat 1.00 x8", 1.0, 8)]
    absH_all, z_all, q_all = {}, {}, {}
    keys = [p[0] for p in POLICIES]
    for k in range(NT):
        hq = prod.mass_consistent_effect(H0[k].copy(), b, a_gate, mode="minnorm")
        lin = np.clip(np.expm1(b + hq), 0.0, None)
        prof = lin / lin.sum()
        for name, s, reps in POLICIES:
            if s is None:
                C = Craw
            else:
                sv = np.zeros(G)
                sv[gated] = s
                C = fracpair(Craw, sv, np.random.default_rng(700 + k), reps=reps)
            z = ref.z(emit(prof, dlib, C, 202 + k))
            z_all.setdefault(name, []).append(z)
            q_all.setdefault(name, []).append(
                false_discovery_control(2 * norm.sf(np.abs(z)), method="bh"))
        absH_all[k] = np.abs(hq[gated])
        print(f"  target {etargets[k]} done [{time.time()-t0:.0f}s]", flush=True)

    lc = np.log10(np.maximum(cpm[gated], 1e-9))
    print(f"\nSIGNAL RESPONSE, mean over {NT} targets. |z| by counts/cell tercile, and the two "
          f"Spearman\ndiagnostics 171 section 12c used (a15 reads +0.485 expr / -0.391 |H| "
          f"INSIDE its call set).")
    print(f"{'policy':>16}{'n_pred':>8}{'|z|T0':>8}{'|z|T1':>8}{'|z|T2':>8}{'T2/T0':>7}"
          f"{'rho expr':>10}{'rho |H|':>9}{'rho expr in':>12}{'rho |H| in':>11}")
    for name in keys:
        npd, zt, rex, rh, rexi, rhi = [], [[], [], []], [], [], [], []
        for k in range(NT):
            z, q = z_all[name][k], q_all[name][k]
            sig = q < 0.05
            npd.append(int(sig.sum()))
            lp = -np.log10(np.maximum(q, 1e-300))
            for t in range(3):
                zt[t].append(np.median(np.abs(z)[tier == t]))
            aH = absH_all[k]
            rex.append(np.corrcoef(rankdata(lp), rankdata(lc))[0, 1])
            rh.append(np.corrcoef(rankdata(lp), rankdata(aH))[0, 1])
            if sig.sum() > 50:
                rexi.append(np.corrcoef(rankdata(lp[sig]), rankdata(lc[sig]))[0, 1])
                rhi.append(np.corrcoef(rankdata(lp[sig]), rankdata(aH[sig]))[0, 1])
        t0m, t1m, t2m = (float(np.mean(zt[t])) for t in range(3))
        print(f"{name:>16}{int(np.median(npd)):>8}{t0m:>8.2f}{t1m:>8.2f}{t2m:>8.2f}"
              f"{t2m/max(t0m,1e-9):>7.2f}{np.mean(rex):>+10.3f}{np.mean(rh):>+9.3f}"
              f"{(np.mean(rexi) if rexi else np.nan):>+12.3f}"
              f"{(np.mean(rhi) if rhi else np.nan):>+11.3f}", flush=True)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
