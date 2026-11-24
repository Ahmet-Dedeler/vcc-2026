"""GATE 2: does `pds_cosine` turn negative past gamma_gene ~0.7, as mse-mechanism reports?

The claim I am asked to confirm or refute independently: raw `mse` is monotone in gamma_gene but
`pds` improves only for gamma_gene <= 0.5 and turns negative past ~0.7, so a gamma chosen to
maximise the mse multiplier alone walks past the binding constraint. Their evidence is a single
surface (h1-test at three donors).

Why this needs no emission, which is what makes a two-surface answer affordable. `pds_cosine` is
computed from the PSEUDOBULK in `bulk_lognorm` space, and the quasi emitter's exact-total repair
preserves each gene's panel total exactly -- so the pseudobulk, and therefore `pds`, is invariant
to (spread, phi). 096 asserts this; 097's `emit_invariance.py` verified it for X/A/B bit-exactly;
a8's board result measured `pds` delta = 0.000000 under phi=4. And at s = 1 the render round-trip
`dh -> profile -> log1p(TS*profile) - b_ctrl` is exact, so the predicted pseudobulk IS
`b_ctrl + dh'`. Every arm is therefore compared at the SAME amplitude, which isolates the
shrinkage -- the thing in dispute -- from the amplitude, which is not.

Scorer's own `discrimination_score` is called, not a reimplementation.
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import rankdata
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor"); sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats, CTRL, TS
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor

CFG = [("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        ["k562_gwps", "orion_hek293t"], "orion_hct116"),
       ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        ["k562_gwps", "orion_hek293t"], None)]
GG = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
GP = [0.0, 1.0, 2.0]
FLOOR = 0.3          # a13's --shrink-pert-floor


def donor_mats(dons, genes, tg):
    P, G = len(tg), genes.size
    M, H = {}, {}
    for nd in dons:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cc = np.array([gi.get(str(g), -1) for g in genes]); ok = cc >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        Mx = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(tg):
            if x in pi: Mx[r, ok] = s[pi[x]][cc[ok]]; h[r] = True
        M[nd] = Mx / (float(np.median(np.linalg.norm(s[:, cc[ok]], axis=1))) or 1.0)
        H[nd] = h; del dn, s
    return M, H


def keys(M, H, dons, P, G):
    """alpha (per-gene midrank of the all-donor-unanimity agreement) and z (per-pert cosine)."""
    ag = np.zeros((P, G)); dc = np.zeros(P); nd = np.zeros(P, int)
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        nd[i] = len(have)
        if len(have) < 2:
            continue
        st = np.vstack([M[n][i] for n in have]); sg = np.sign(st)
        un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
        ag[i] = np.where(un, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)
        tot, cnt = 0.0, 0
        for a in range(len(have)):
            for b in range(a + 1, len(have)):
                x, y = M[have[a]][i], M[have[b]][i]
                tot += float(x @ y) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12); cnt += 1
        dc[i] = tot / cnt
    al = np.empty((P, G))
    for i in range(P):
        al[i] = rankdata(ag[i], method="average") / G
        if np.unique(ag[i]).size < 2: al[i] = 1.0
    z = (rankdata(dc, method="average") - 0.5) / P
    z[nd < 2] = 0.0                     # a13's `--shrink-pert-missing damp`: rank to the bottom
    z = np.maximum(z, FLOOR)            # ...then a13's `--shrink-pert-floor 0.3` lifts it.
    # The floor is what keeps a missing row's effect NON-zero. That matters beyond n_pred:
    # `discrimination.py`'s own docstring records that a ZERO-norm predicted effect ties its
    # whole pds row at distance 1.0 and scores exactly 0.5 under the midrank tie policy. With
    # `damp` and no floor the 6 challenge targets carrying one donor would each cost
    # (0.833 - 0.5)/300 = 0.00111 raw pds, i.e. -0.0025 overall. The floor removes that.
    return al, z, nd


for nm, real, pred, base_d, third in CFG:
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    dh0, w, cols, sD = t["dh"], t["w"], t["cols"], t["sumD"]
    keep = t["keep"]; tg = np.array([str(x) for x in pp[keep]])
    P, G = dh0.shape
    ci = int(np.flatnonzero(rp == CTRL)[0])
    bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
    rmeans = bulk_lognorm_means(rs.astype(np.float64), TS)
    combos = [("2-donor", base_d)] + ([("3-donor", base_d + [third])] if third else [])
    for lbl, dons in combos:
        M, H = donor_mats(dons, rg, tg)
        al, z, ndon = keys(M, H, dons, P, G)
        print(f"\n##### {nm} / {lbl} key: P={P}, rows with <2 donors {int((ndon<2).sum())}",
              flush=True)
        print(f"  {'':<10}" + "".join(f"{'gg='+str(g):>19}" for g in GG))
        for gp in GP:
            row = []
            for gg in GG:
                dh = (z ** gp)[:, None] * (al ** gg) * dh0
                pm = np.tile(bctrl, (P + 1, 1)); pm[:P] += dh
                perts = np.append(tg, CTRL)
                pdsv = discrimination_score(
                    pred_bulk=(perts, pm), real_bulk=(rp, rmeans),
                    pert_col="target_gene", control=CTRL, distance="cosine",
                    exclude_target_gene=True, exclusion_scope="panel",
                    control_source="pred", genes=rg)
                pds = float(np.mean(list(pdsv.values()))) if isinstance(pdsv, dict) else float(pdsv)
                X = float((np.einsum("ij,ij->i", dh, t["d"]) * w).sum())
                B = float(_across_pert_budget(dh, cols, w))
                kap = X / np.sqrt(B * sD)
                row.append((pds, kap))
            print(f"  gp={gp:<7.1f}" + "".join(f"  pds {p:.4f} kap {k:.5f}" for p, k in row),
                  flush=True)
