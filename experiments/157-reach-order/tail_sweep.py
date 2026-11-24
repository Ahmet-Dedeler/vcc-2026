"""157 Q3 — localise what the SHIPPED tier spread policy does to the confidence ORDER.

`arms.py` found, on the repaired HCT116 harness, that the shipped a12/a14/a15 policy
(`--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`) reads raw
reach 0.066-0.075 while a8's FLAT `spread 0.25 --overdispersion 4.0` reads 0.152, with
`n_pred` going the wrong way (2,867 -> 4,555).  The board agrees in sign: a8's raw reach was
0.1890 and a15's is 0.1811.

Hypothesis: at `spread = 1.0, phi = 4` a low-expression gene's 400 emitted cells are so
zero-inflated that the rank-sum calls a large DOWN shift at zero predicted effect.  Those
calls land at the FRONT of `reach`'s order (significant calls first, ascending p_adj) and
they are drawn from exactly the 2,261 gated genes the policy declared least informative.

This sweeps the tail policy, holding phi = 4 and the profile fixed, and additionally records
what fraction of the first ten pool genes in the emitted order come from the tail block.
Deterministic emission with a re-seeded twin as the noise control.  No cell-eval2, no build,
no submission.

`pds` and `nmae` are structurally invariant across every arm here (exact-total repair, fixed
per-cell library L), so they are not measured.  `mse` is NOT invariant -- the tail policy
exists to raise the emitted dispersion that `mse`'s jackknife credit reads -- and it is not
measured here either.  Any revert must be priced on `mse` first.
"""
from __future__ import annotations
import os
import sys
import time

import anndata as ad
import numpy as np
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts        # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z        # noqa: E402

REAL = f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"
PRED = f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
CTRL, COL, MIN_CPM, ALPHA = "non-targeting", "target_gene", 5.0, 0.05
N_EMIT, EMIT_L = 400, 100_000
NP = int(os.environ.get("N_PERTS", "132"))
B_R, R_R = 0.079022, 0.964145

# name -> (tail_spread, keep_top, keep_spread, gated_only, phi, seed)
# `gated_only` keeps the sub-gate genes at 0.25 instead of at `tail`, which separates the
# per-cell library-variance channel from the gated tail genes' own p-values.
ARMS = [
    ("a8_flat_0.25",            None, 0,     None, False, 4.0, 90_000),
    ("a15_tier8000_tail1.0",     1.0, 8000,  0.25, False, 4.0, 90_000),
    ("a15_tier8000_tail1.0_b",   1.0, 8000,  0.25, False, 4.0, 7_000_000),
    ("tier8000_tail1.0_gatedonly", 1.0, 8000, 0.25, True, 4.0, 90_000),
    ("tier8000_tail0.50",        0.5, 8000,  0.25, False, 4.0, 90_000),
    ("tier9500_tail1.0",         1.0, 9500,  0.25, False, 4.0, 90_000),
]


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    t0 = time.time()
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes_full = np.asarray([str(g) for g in R.var_names])
    G_full = genes_full.size
    ci = np.flatnonzero(rlab == CTRL)
    cpm = np.zeros(G_full)
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1); s[s == 0] = 1.0
        cpm += (m * (1e6 / s[:, None])).sum(0)
    cpm /= ci.size
    keep = cpm > MIN_CPM
    G = int(keep.sum())
    cs = np.empty((G, ci.size))
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1); s[s == 0] = 1.0
        cs[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    ref_lin = cpm[keep]
    b_ctrl = np.log1p(cpm)
    gk = genes_full[keep]
    gpos = {g: i for i, g in enumerate(gk)}
    fr = np.empty(G, np.int64)
    fr[np.argsort(gk)] = np.arange(G)
    print(f"harness: {ci.size} controls, {G}/{G_full} gated [{time.time()-t0:.0f}s]")

    ref = np.load(f"{ROOT}/experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
    assert np.array_equal(ref["genes"].astype(str), gk)
    real_padj, real_lfc = ref["real_padj"], ref["real_lfc"]
    rperts = list(ref["perts"].astype(str))

    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    targets = [t for t in rperts if (plab == t).sum() > 0][:NP]
    print(f"{len(targets)} perturbations, {len(ARMS)} arms")

    reach = {a[0]: [] for a in ARMS}
    npred = {a[0]: [] for a in ARMS}
    headtail = {a[0]: [] for a in ARMS}     # fraction of the top-10 drawn from the tail block
    headpure = {a[0]: [] for a in ARMS}
    for i, t in enumerate(targets):
        ri = rperts.index(t)
        rs = (real_padj[ri] < ALPHA) & np.isfinite(real_padj[ri])
        j = gpos.get(t)
        if j is not None:
            rs[j] = False
        gs = np.flatnonzero(rs)
        nc = gs.size
        if nc == 0:
            continue
        rl = real_lfc[ri]

        totals = dense(P, np.flatnonzero(plab == t)).sum(0).astype(np.int64)
        profile = (totals + 0.5) / (N_EMIT * EMIT_L)
        pcpm = 1e6 * totals / max(totals.sum(), 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        gated = (cpm > MIN_CPM) & (b_ctrl > 0)
        score = np.where(gated, eff, -np.inf)
        rank = np.lexsort((np.arange(G_full), -score))

        for name, tail, ktop, kspread, gonly, phi, seed in ARMS:
            if tail is None:
                sv = np.full(G_full, 0.25)
                in_tail = np.zeros(G_full, bool)
            else:
                sv = np.full(G_full, 0.25 if gonly else float(tail))
                k = min(int(ktop), int(gated.sum()))
                in_tail = np.zeros(G_full, bool)
                in_tail[rank[k:]] = True
                if gonly:
                    sv[in_tail & gated] = float(tail)
                    in_tail &= gated
                else:
                    sv[rank[:k]] = float(kspread)
                    in_tail &= gated
            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                                          mode="quasi", spread=sv, overdispersion=phi)
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            s = B.sum(1); s[s == 0] = 1.0
            nb = np.log1p(B[:, keep] * (1e6 / s[:, None]))
            del B
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
            pp = 2 * norm.sf(np.abs(z))
            pq = false_discovery_control(pp, method="bh")
            pl = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
            del nb
            ps = pq < ALPHA
            if j is not None:
                ps[j] = False
            npred[name].append(int(ps.sum()))
            order = gs[np.lexsort((fr[gs], -np.abs(np.nan_to_num(pl[gs], nan=-np.inf)),
                                  pp[gs], pq[gs], ~ps[gs]))]
            adj = np.isfinite(rl[order]) & (rl[order] != 0)
            order = order[adj]
            if order.size == 0:
                reach[name].append(0.0)
                continue
            ok = np.sign(pl[order]) == np.sign(rl[order])
            pur = np.cumsum(ok) / np.arange(1, ok.size + 1)
            d = np.flatnonzero(pur >= 0.9)
            reach[name].append((int(d[-1]) + 1) / nc if d.size else 0.0)
            h = min(10, order.size)
            headpure[name].append(float(ok[:h].mean()))
            headtail[name].append(float(in_tail[keep][order[:h]].mean()))
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(targets)} [{time.time()-t0:.0f}s]", flush=True)

    base = float(np.mean(reach["a8_flat_0.25"]))
    print(f"\n{'arm':<30}{'n_pred':>8}{'reach':>9}{'d_reach':>10}{'d_overall':>11}"
          f"{'top10 pure':>12}{'top10 from tail':>17}")
    for name, *_ in ARMS:
        r = float(np.mean(reach[name]))
        print(f"{name:<30}{np.median(npred[name]):>8.0f}{r:>9.5f}{r-base:>+10.5f}"
              f"{(r-base)/(R_R-B_R)/6:>+11.5f}{np.mean(headpure[name]):>12.4f}"
              f"{np.mean(headtail[name]):>17.4f}")
    print("\nbaseline = a8_flat_0.25.  `pds` and `nmae` are invariant across all arms; "
          "`mse` is NOT and is unmeasured here.")


if __name__ == "__main__":
    main()
