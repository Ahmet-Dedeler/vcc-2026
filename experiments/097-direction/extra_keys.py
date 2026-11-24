"""Route (d), the remainder: two more per-gene weightings with a mechanism, both measured.

 (d1) CONTEXT-CONTROL keys.  The brief asks about "better use of the target context's own control
      cells".  The only per-gene quantities they give are expression level and cell-to-cell
      dispersion.  096's framework predicts these fail -- they condition on SURFACE-SPECIFIC
      structure (this context's depth) rather than on the model's confidence -- and 096 already
      saw a CPM-keyed emission tail lose to a plain top-K.  Tested here against X directly.
 (d2) A PANEL-LEVEL per-gene donor confidence: how consistently the two donors agree on gene g
      ACROSS all perturbations.  This is a per-gene confidence (a model property, which 096
      measured transfers at 1.17x) rather than a per-gene slope (a cell-line property, 1.02x),
      and it is much cheaper to build than the per-(p,g) matrix.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor"); sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats, CTRL, TS
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor

HARNESS = [("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
            f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
           ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
            f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")]
DON = ("k562_gwps", "orion_hek293t")
for nm, real, pred in HARNESS:
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    dh, d, w = t["dh"], t["d"], t["w"]
    ci = int(np.flatnonzero(rp == CTRL)[0])
    bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
    cpm = np.expm1(bctrl); cpm = 1e6 * cpm / cpm.sum()
    tg = np.array([str(x) for x in pp[t["keep"]]])
    P, G = dh.shape

    def shape(DH):
        X = float((np.einsum("ij,ij->i", DH, d) * w).sum())
        A = float((np.einsum("ij,ij->i", DH, DH) * w).sum())
        B = float(_across_pert_budget(DH, t["cols"], w))
        return X / np.sqrt(B * t["sumD"]), (A - B) / B
    k0, e0 = shape(dh)

    # donor matrices, panel-wide
    Ms, hv = {}, {}
    for nd in DON:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}")
        m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in rg]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        shared = np.array([pi[p] for p in dn.perts if p in pi])  # all
        Sfull = np.zeros((s.shape[0], G), np.float32)
        Sfull[:, ok] = s[:, cols[ok]]
        Ms[nd] = Sfull / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nd] = {p: i for i, p in enumerate(dn.perts)}
        del dn, s
    # (d2) per-gene agreement ACROSS the donors' shared perturbations
    common = [p for p in hv[DON[0]] if p in hv[DON[1]]]
    ia = np.array([hv[DON[0]][p] for p in common]); ib = np.array([hv[DON[1]][p] for p in common])
    Aq, Bq = Ms[DON[0]][ia], Ms[DON[1]][ib]
    na = np.linalg.norm(Aq, axis=0); nb = np.linalg.norm(Bq, axis=0)
    gene_cos = np.einsum("ij,ij->j", Aq, Bq) / np.maximum(na * nb, 1e-12)   # per-GENE cosine
    gene_sign = (np.sign(Aq) == np.sign(Bq)).mean(0)                        # per-GENE sign concord
    del Aq, Bq, Ms

    def r01(v):
        return (np.argsort(np.argsort(v)).astype(float) + 0.5) / v.size
    KEYS = {
        "(d1) control CPM rank": r01(cpm),
        "(d1) control CPM rank, inverted": 1 - r01(cpm),
        "(d1) log CPM, gated >5 CPM": (cpm > 5).astype(float) * 0.9 + 0.1,
        "(d2) per-GENE donor cosine": r01(gene_cos),
        "(d2) per-GENE donor sign concord": r01(gene_sign),
    }
    print(f"\n##### {nm}: base kappa={k0:.5f} eps={e0:+.5f}   "
          f"({len(common)} donor-shared perturbations built the per-gene keys)")
    print(f"  {'per-gene key':<38}" + "".join(f"{'^'+str(g):>9}" for g in (0.5, 1, 2, 3)))
    for lbl, v in KEYS.items():
        row = []
        for g in (0.5, 1.0, 2.0, 3.0):
            k, e = shape((v ** g)[None, :] * dh)
            row.append(k / k0)
        print(f"  {lbl:<38}" + "".join(f"{x:>9.3f}" for x in row))
