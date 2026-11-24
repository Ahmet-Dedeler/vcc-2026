"""164 — the two things `decomp.py` left unsettled.

(A) POOLED vs PER-PERTURBATION.  `fid` is `mean_p( k_p / max(n_pred_p, n_real_p) )`.  Pooling
    all (perturbation, gene) calls and reporting one k/n is a DIFFERENT statistic, and
    AGENTS.md 11.1 rule 8 records that this exact substitution has already produced an
    authoritative-looking wrong answer in this repo (exp 154, pooled 1.0000 vs per-pert 4/70).
    `decomp.py` printed pooled block precisions.  They disagree with the metric: pooled
    keep-block precision is 0.5122 while `fid` after dropping the tail is 0.5077.  This
    resolves which is which and why.

(B) THE h1-TEST SUB-GATE BAND.  `decomp.py`'s modal-ratio detector for zeroed genes returned
    98 genes with control pseudobulk CPM up to 105 -- that is integer-quantisation
    coincidence, not a CPM gate, so the detector FAILED and its output is discarded.  The
    replacement compares `log2(pred_pb / ctrl_pb)` -- both ratio-of-sums, so the
    normalisation artefact cancels exactly and what is left is the builder's own effect --
    binned by control pseudobulk CPM.  A `--min-cpm 10` build must show a hard floor at 10.

    FAILURE MODE, so the check is not vacuous: if the 095 h1-test artifact was built with
    `--min-cpm 5.0`, the 5-10 CPM band and the 10-20 CPM slice show the SAME effect size and
    the h1-test surface simply has no sub-gate band to decompose -- which would mean the
    challenge's 794/958/997-gene block is NOT represented on this surface at all.
"""
from __future__ import annotations

import sys

import anndata as ad
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = f"{ROOT}/experiments/164-callset-artefact"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface  # noqa: E402

ARMS = ["a15", "a15_b", "armA_tailphi0", "armA_tailphi0_b", "armB_tail050",
        "armC_a8_flat", "armC_a8_flat_b"]
REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = (f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
CTRL, COL = "non-targeting", "target_gene"


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def control_stats(genes):
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    gfull = np.asarray([str(g) for g in R.var_names])
    ci = np.flatnonzero(rlab == CTRL)
    csum = np.zeros(gfull.size)
    cellcpm = np.zeros(gfull.size)
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        csum += m.sum(0)
        cellcpm += (m * (1e6 / s[:, None])).sum(0)
    cellcpm /= ci.size
    ctrl_pb = 1e6 * csum / csum.sum()
    del R
    gpos = {g: i for i, g in enumerate(gfull)}
    idx = np.array([gpos[g] for g in genes])
    return cellcpm, ctrl_pb, gfull, idx


def pred_pb_all(gfull):
    """Predicted pseudobulk CPM for a few targets (the gate is context-level, so any target
    exposes it; several are used so a per-target donor gap cannot be mistaken for the gate)."""
    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    pg = np.asarray([str(g) for g in P.var_names])
    tg = sorted(set(plab) - {CTRL})[:4]
    pos = {g: i for i, g in enumerate(pg)}
    cols = np.array([pos.get(g, -1) for g in gfull])
    ok = cols >= 0
    out = []
    for t in tg:
        tot = dense(P, np.flatnonzero(plab == t)).sum(0)
        v = np.zeros(gfull.size)
        v[ok] = 1e6 * tot[cols[ok]] / tot.sum()
        out.append(v)
    del P
    return tg, np.array(out)


def main() -> None:
    S = {a: Surface(f"{HERE}/de164.{a}.npz") for a in ARMS}
    spread = {a: np.load(f"{HERE}/de164.{a}.npz")["spread"] for a in ARMS}
    s0 = S["a15"]
    genes = s0.genes

    # ---------------------------------------------------------------------------------
    print("=" * 100)
    print("(B) does the h1-test prediction carry a sub-gate band?  (both sides ratio-of-sums)")
    print("=" * 100)
    cellcpm, ctrl_pb, gfull, idx = control_stats(genes)
    tg, ppb = pred_pb_all(gfull)
    print(f"targets sampled: {tg}")
    eff = np.log2((ppb + 1e-9) / (ctrl_pb[None, :] + 1e-9))     # the builder's own effect
    scored = np.zeros(gfull.size, bool)
    scored[idx] = True
    print(f"{'control pb-CPM bin':<22}{'n scored':>10}{'median |effect|':>17}"
          f"{'frac |effect|<1e-6':>21}")
    BINS = [(0, 5), (5, 7), (7, 10), (10, 12), (12, 20), (20, 50), (50, 1e9)]
    for lo, hi in BINS:
        m = scored & (ctrl_pb > lo) & (ctrl_pb <= hi)
        if not m.any():
            continue
        e = np.abs(eff[:, m])
        print(f"{f'({lo:g}, {hi:g}]':<22}{int(m.sum()):>10}{np.median(e):>17.5f}"
              f"{float((e < 1e-6).mean()):>21.1%}")
    band = scored & (cellcpm > 5.0) & (ctrl_pb <= 10.0)
    above = scored & (ctrl_pb > 10.0) & (ctrl_pb <= 20.0)
    print(f"\nh1-test sub-gate band (cell-CPM>5 & pb-CPM<=10): {int(band.sum())} genes")
    print(f"   median |effect| in band {np.median(np.abs(eff[:, band])):.5f}   "
          f"in the 10-20 CPM slice {np.median(np.abs(eff[:, above])):.5f}")
    verdict = ("ZEROED -> the band is represented on this surface"
               if np.median(np.abs(eff[:, band])) < 0.1 * np.median(np.abs(eff[:, above]))
               else "NOT ZEROED -> the 095 h1-test artifact was NOT built at --min-cpm 10; "
                    "the challenge's sub-gate block has NO h1-test analogue")
    print(f"   VERDICT: {verdict}")

    bidx = np.flatnonzero(band[idx])
    if bidx.size:
        sv = spread["a15"][:, bidx]
        print(f"   where the band lands in a15's spread vector: "
              f"tail(1.0) {float((sv > 0.5).mean()):.1%}, keep(0.25) {float((sv <= 0.5).mean()):.1%}")

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("(A) POOLED vs PER-PERTURBATION precision  (rule 8)")
    print("=" * 100)
    ok = np.sign(s0.pred_lfc) == np.sign(s0.real_lfc)
    called = s0.pred_sig & s0.adj
    tail = spread["a15"] > 0.5
    nr = np.array([int(s0.real_sig[p].sum()) for p in range(len(s0.perts))])

    def per_pert(mask):
        n = (called & mask).sum(1)
        k = (called & mask & ok).sum(1)
        keep = n > 0
        return n, k, float((k[keep] / n[keep]).mean()), float(k.sum() / max(n.sum(), 1))

    print(f"{'block':<16}{'n calls':>10}{'pooled k/n':>12}{'mean_p(k/n)':>13}"
          f"{'mean n/pert':>12}{'corr(n_p, prec_p)':>19}")
    for nm, mask in (("keep (0.25)", ~tail), ("tail (1.0)", tail),
                     ("ALL", np.ones_like(tail))):
        n, k, mp, pooled = per_pert(mask)
        keep = n > 0
        c = np.corrcoef(n[keep], k[keep] / n[keep])[0, 1] if keep.sum() > 2 else np.nan
        print(f"{nm:<16}{int(n.sum()):>10}{pooled:>12.4f}{mp:>13.4f}"
              f"{n.mean():>12.0f}{c:>19.3f}")

    print(f"\n  a15 fid as scored                     {s0.fid():.5f}")
    n_all = called.sum(1)
    k_all = (called & ok).sum(1)
    print(f"  mean_p( k_p / max(n_p, n_real_p) )    "
          f"{float(np.mean(k_all / np.maximum(n_all, nr))):.5f}   <- identical, as it must be")
    print(f"  pooled  sum(k)/sum(n)                 {float(k_all.sum()/n_all.sum()):.5f}")

    nk = (called & ~tail).sum(1)
    kk = (called & ~tail & ok).sum(1)
    print(f"\n  counterfactual: drop the tail block post hoc (arm D, NOT shippable)")
    print(f"    mean_p( k/max(n,n_real) )           "
          f"{float(np.mean(kk / np.maximum(nk, nr))):.5f}")
    print(f"    pooled                              {float(kk.sum()/nk.sum()):.5f}")
    print(f"    -> the pooled figure overstates the gain by "
          f"{float(kk.sum()/nk.sum()) - float(np.mean(kk/np.maximum(nk, nr))):.4f}, "
          f"because precision rises with n_pred across perturbations")

    # per-perturbation precision of the tail block, distribution
    nt = (called & tail).sum(1)
    kt = (called & tail & ok).sum(1)
    m = nt > 0
    pt = kt[m] / nt[m]
    print(f"\n  tail-block per-perturbation precision: n={int(m.sum())} perts, "
          f"mean {pt.mean():.4f}, median {np.median(pt):.4f}, sd {pt.std(ddof=1):.4f}")
    print(f"    sd of the mean {pt.std(ddof=1)/np.sqrt(m.sum()):.4f}  -> "
          f"{'indistinguishable from 0.5' if abs(pt.mean()-0.5) < 2*pt.std(ddof=1)/np.sqrt(m.sum()) else 'DIFFERENT from 0.5'}")
    pk = kk[nk > 0] / nk[nk > 0]
    print(f"  keep-block per-perturbation precision: mean {pk.mean():.4f}, "
          f"sd of the mean {pk.std(ddof=1)/np.sqrt((nk>0).sum()):.4f}")
    d = pk - pt
    print(f"  paired keep - tail: mean {d.mean():+.4f}, sd of the mean "
          f"{d.std(ddof=1)/np.sqrt(d.size):.4f}, t = {d.mean()/(d.std(ddof=1)/np.sqrt(d.size)):.2f}")

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("(C) what the arms actually do to the call set, per perturbation")
    print("=" * 100)
    print(f"{'arm':<20}{'mean n_pred':>13}{'mean_p prec':>13}{'fid':>10}"
          f"{'jac':>10}{'overlap w/ real':>17}")
    for a in ARMS:
        s = S[a]
        c = s.pred_sig & s.adj
        o = np.sign(s.pred_lfc) == np.sign(s.real_lfc)
        n = c.sum(1)
        k = (c & o).sum(1)
        ov = (s.pred_sig & s.real_sig).sum(1)
        print(f"{a:<20}{n.mean():>13.0f}{float((k[n>0]/n[n>0]).mean()):>13.4f}"
              f"{s.fid():>10.5f}{s.jac():>10.5f}{ov.mean():>17.1f}")
    print(f"{'reference n_real':<20}{nr.mean():>13.0f}")


if __name__ == "__main__":
    main()
