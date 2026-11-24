"""Does the ANALYTIC predicted lfc reproduce the REALISED (emitted + rank-summed) one?

Rule 3: before believing any analytic sweep, confirm the quantity it computes is the quantity
the scorer sees.  Drives production `emit_counts` with a12's exact tiered-spread policy and
phi=4, twice with different seeds, and compares per-row NMAE against the analytic value.

    N_PERTS=40 uv run python validate.py <h1|orion>
"""
from __future__ import annotations
import os, sys, time
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from msekit import group_stats
from vcc2026.donor_transfer import emit_counts

W_ = sys.argv[1]; NP_ = int(os.environ.get("N_PERTS", "40"))
KEEP_TOP, KEEP_SPREAD, TAIL_SPREAD, PHI = 8000, 0.25, 1.0, 4.0
H = K.Harness(W_); eff, s = K.a12_eff(H)
rp, rs, rc, rj, rg = group_stats(K.REAL[W_])
pp, ps, pc, pj, pg = group_stats(K.PRED[W_])
S = ps.sum(1).astype(np.float64); NC = int(round(float(np.median(pc))))
expression = np.expm1(H.bctrl); cpm_b = 1e6 * expression / max(expression.sum(), 1e-12)
gated_b = (cpm_b > 5.0) & (H.bctrl > 0)
G = H.G
sel = H.rows[:NP_]; selr = H.rrows[:NP_]
Lan = H.lfc_hat(eff[sel])
print(f"# {W_}: {sel.size} rows, {NC} cells/row, cliff s={s:.4f}", flush=True)
res = {0: [], 1: []}; ana = []
t0 = time.time()
for k, (i, j) in enumerate(zip(sel, selr)):
    gs = H.sig[j]; r = H.rlfc[j][gs]; den = np.abs(r).sum()
    ana.append(np.abs(Lan[k][gs] - r).sum() / den)
    e = eff[i]
    prof = np.expm1(H.bctrl + e); np.clip(prof, 0.0, None, out=prof); prof /= prof.sum()
    sc_ = np.where(gated_b, np.abs(e), -np.inf)
    order = np.lexsort((np.arange(G), -sc_))
    sv = np.full(G, TAIL_SPREAD); sv[order[:min(KEEP_TOP, int(gated_b.sum()))]] = KEEP_SPREAD
    for sd in (0, 1):
        v, col, ip = emit_counts(prof, np.full(NC, S[i] / NC), 0.0,
                                 np.random.default_rng(7919 + k + 10_000 * sd),
                                 mode="quasi", spread=sv, overdispersion=PHI)
        blk = np.zeros((NC, G))
        for r_ in range(NC):
            blk[r_, col[ip[r_]:ip[r_+1]]] = v[ip[r_]:ip[r_+1]]
        t = blk.sum(1); t[t == 0] = 1.0
        pl = np.log2((blk[:, H.gate] * (1e6 / t[:, None])).mean(0) + 1e-9) - \
             np.log2(H.ref_lin + 1e-9)
        res[sd].append(np.abs(pl[gs] - r).sum() / den)
    if (k + 1) % 10 == 0: print(f"  {k+1}/{sel.size} [{time.time()-t0:.0f}s]", flush=True)
ana = np.array(ana); r0 = np.array(res[0]); r1 = np.array(res[1])
print(f"\n{W_}  analytic {ana.mean():.5f}   realised seed0 {r0.mean():.5f}  "
      f"seed1 {r1.mean():.5f}   seed twin gap {abs(r0.mean()-r1.mean()):.5f}")
print(f"  realised - analytic:  {r0.mean()-ana.mean():+.5f} / {r1.mean()-ana.mean():+.5f}")
print(f"  per-row corr(analytic, realised) = {np.corrcoef(ana, r0)[0,1]:.5f}; "
      f"max |row diff| {np.abs(r0-ana).max():.4f}, median {np.median(np.abs(r0-ana)):.4f}")
