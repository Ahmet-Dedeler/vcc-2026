"""164 — decompose `fid` precision on our own call set, and price the arms on the DE members.

Reads the seven `de164.*.npz` arm tables written by `realize164.py` (h1-test, donors
`k562_gwps` + `orion_hek293t`, no `h1_2025`).  Everything here is deterministic given those
files: no emitter, no cell-eval2, no build, no submission.

WHAT IS BEING TESTED
--------------------
`n_pred` (median 6035 here, ~4,400 on the challenge) far exceeds `n_real` (median 89 here,
~240 there), so `fid = k / max(n_pred, n_real)` collapses to **precision on our own call
set** in BOTH regimes.  The level does not transfer (rule 9) but the DECOMPOSITION does:
which block of our call set supplies the calls, and what is each block's precision.

Blocks, from the per-target `spread` vector `realize164.py` saved:
  * keep  -- spread 0.25, the top-8000 gated genes by |emitted effect|
  * tail  -- spread 1.0, everything else
and, cut across those, the genes whose PREDICTED EFFECT IS EXACTLY ZERO because the
builder's `--min-cpm 10.0` gate zeroed them (`profile_from_effect`, make_submission_pm_cd4.py
:1281).  Those are detected here from the prediction itself, not assumed.

KNOWN-ANSWER GATES (AGENTS.md 11.1 rule 10's corollary).  All three could fail:
  1. `a15`  native-BH reach == 0.18283   (157 `power_h1test.log` f=1.0)
  2. `armC` native-BH reach == 0.20291   (157 `power_h1test.log` f=1.0, a8_flat_s0.25_phi4)
  3. `pred_lfc` invariance across arms -- asserted, NOT assumed (095 saw `nmae` invariant
     only to 4 dp).  If the emitter's per-cell library sums drifted between arms this would
     read a large number and the "spread/phi never changes a scored direction" claim dies.
"""
from __future__ import annotations

import os
import sys

import anndata as ad
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = f"{ROOT}/experiments/164-callset-artefact"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface, ANCHOR  # noqa: E402

ARMS = ["a15", "a15_b", "armA_tailphi0", "armA_tailphi0_b", "armB_tail050",
        "armC_a8_flat", "armC_a8_flat_b"]
GATE_REACH = {"a15": 0.18283, "armC_a8_flat": 0.20291}

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = (f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
CTRL, COL = "non-targeting", "target_gene"


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


# ---------------------------------------------------------------------------------------
# 1. the zero-predicted-effect gene set, measured from the prediction
# ---------------------------------------------------------------------------------------
def zero_effect_mask(genes):
    """Genes whose predicted pseudobulk proportion equals the control's exactly.

    `profile_from_effect` zeroes the effect outside `(pb_cpm > min_cpm) & (b_ctrl > 0)`, so a
    zeroed gene's predicted abundance is `expm1(b_ctrl_g)` -- i.e. proportional to the control
    pseudobulk with ONE normalising constant shared by every zeroed gene.  So the ratio
    `pred_pb / ctrl_pb` has a spike at that constant, and its members are the zeroed set.
    The gate depends only on the CONTEXT's control profile, so it is the same for every
    target; one target suffices.

    This is a measurement with a visible failure mode: if the builder had used `--min-cpm 5.0`
    (or no gate) the spike would be empty and the function returns an all-False mask.
    """
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

    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    pg = np.asarray([str(g) for g in P.var_names])
    tgt = sorted(set(plab) - {CTRL})[0]
    blk = dense(P, np.flatnonzero(plab == tgt))
    ptot = blk.sum(0)
    del blk, P
    pos = {g: i for i, g in enumerate(pg)}
    cols = np.array([pos.get(g, -1) for g in gfull])
    pred_pb = np.zeros(gfull.size)
    ok = cols >= 0
    pred_pb[ok] = 1e6 * ptot[cols[ok]] / ptot.sum()

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where((ctrl_pb > 0) & (pred_pb > 0), pred_pb / ctrl_pb, np.nan)
    fin = ratio[np.isfinite(ratio)]
    # the modal ratio, to 6 significant figures, is the normalising constant
    key = np.round(np.log10(fin), 6)
    vals, cts = np.unique(key, return_counts=True)
    mode = vals[np.argmax(cts)]
    spike = np.isfinite(ratio) & (np.abs(np.round(np.log10(np.where(
        np.isfinite(ratio), ratio, 1.0)), 6) - mode) < 1e-9)

    gpos = {g: i for i, g in enumerate(gfull)}
    idx = np.array([gpos[g] for g in genes])
    print(f"zero-effect detection: modal pred_pb/ctrl_pb spike holds {int(spike.sum())} of "
          f"{int(np.isfinite(ratio).sum())} genes ({cts.max()} at the mode)")
    if spike.sum():
        cp = ctrl_pb[spike]
        print(f"   their control pseudobulk CPM: max {cp.max():.4f}, "
              f"99.9th pct {np.percentile(cp, 99.9):.4f}  -> the builder's --min-cpm")
        nz = ctrl_pb[~spike & (ctrl_pb > 0)]
        print(f"   non-zeroed genes' control pb CPM: min {nz.min():.4f}")
    return spike[idx], cellcpm[idx], ctrl_pb[idx]


# ---------------------------------------------------------------------------------------
# 2. members and decomposition
# ---------------------------------------------------------------------------------------
def members(s):
    return dict(fid=s.fid(), reach=s.reach(), jac=s.jac(), nmae=s.nmae(), pds=s.pds())


def decompose(s, spread, zmask):
    """Per-perturbation precision, split into keep / tail x real-effect / zero-effect."""
    tail = spread > 0.5                      # 1.0 (or 0.5 for armB) vs the 0.25 keep block
    if not tail.any():                       # flat arm -- everything is 'keep'
        tail = np.zeros_like(tail, dtype=bool)
    blocks = {
        "keep":        ~tail & ~zmask[None, :],
        "tail_real":    tail & ~zmask[None, :],
        "tail_zero":    tail &  zmask[None, :],
        "keep_zero":   ~tail &  zmask[None, :],
    }
    ok = np.sign(s.pred_lfc) == np.sign(s.real_lfc)
    called = s.pred_sig & s.adj
    out = {}
    for name, m in blocks.items():
        m = np.broadcast_to(m, called.shape) & called
        n = m.sum(1)
        k = (m & ok).sum(1)
        out[name] = (int(n.sum()), int(k.sum()),
                     float(k.sum() / max(n.sum(), 1)), float(np.median(n)))
    n_all = called.sum(1)
    k_all = (called & ok).sum(1)
    out["ALL"] = (int(n_all.sum()), int(k_all.sum()),
                  float(k_all.sum() / max(n_all.sum(), 1)), float(np.median(n_all)))
    return out


def fid_dropping(s, drop):
    """Arm D's bound: `fid` if the genes in `drop` were removed from the call set post hoc.
    Not shippable (we emit counts, not calls) -- this is the ceiling the arms chase."""
    ok = np.sign(s.pred_lfc) == np.sign(s.real_lfc)
    called = s.pred_sig & s.adj & ~np.broadcast_to(drop, s.pred_sig.shape)
    out = []
    for p in range(len(s.perts)):
        n_real = int(s.real_sig[p].sum())
        n_pred = int(called[p].sum())
        if n_pred == 0 and n_real == 0:
            continue
        out.append(int((called[p] & ok[p]).sum()) / max(n_pred, n_real))
    return float(np.mean(out))


def main() -> None:
    S = {a: Surface(f"{HERE}/de164.{a}.npz") for a in ARMS}
    spread = {a: np.load(f"{HERE}/de164.{a}.npz")["spread"] for a in ARMS}
    s0 = S["a15"]
    nr, npd = s0.counts()
    print(f"h1-test surface: {s0.perts.size} perturbations x {s0.genes.size} gated genes")
    print(f"n_real  median {np.median(nr):.0f}   n_pred(a15) median {np.median(npd):.0f}")
    print(f"perturbations where n_pred > n_real: "
          f"{int((npd > nr).sum())}/{npd.size}  -> fid == precision\n")

    # ---- GATE 1 & 2 : reach ------------------------------------------------------------
    print("--- known-answer gates ---")
    bad = []
    for arm, want in GATE_REACH.items():
        got = S[arm].reach()
        okg = abs(got - want) < 5e-5
        print(f"  {arm:<16} reach {got:.5f}  expect {want:.5f}  "
              f"{'PASS' if okg else 'FAIL'}")
        if not okg:
            bad.append(arm)
    if bad:
        sys.exit(f"ABORT: known-answer reach gate failed for {bad}")

    # ---- GATE 3 : pred_lfc invariance --------------------------------------------------
    base = S["a15"].pred_lfc
    fin = np.isfinite(base)
    print("  pred_lfc invariance vs a15 (max |delta| over finite entries):")
    for a in ARMS[1:]:
        d = np.abs(S[a].pred_lfc[fin] - base[fin])
        print(f"    {a:<18} max {d.max():.3e}   median {np.median(d):.3e}   "
              f"sign flips {int((np.sign(S[a].pred_lfc[fin]) != np.sign(base[fin])).sum())}"
              f" / {int(fin.sum())}")
    print()

    # ---- the zero-effect block ---------------------------------------------------------
    zmask, cellcpm, ctrl_pb = zero_effect_mask(s0.genes)
    print(f"zero-predicted-effect genes inside the scored set: {int(zmask.sum())} "
          f"of {zmask.size} ({zmask.mean():.1%})\n")

    # the normalisation artefact: what lfc does a zero-effect gene get scored with?
    lfc0 = np.log2((ctrl_pb + 1e-9) / (cellcpm + 1e-9))
    print("--- the normalisation artefact, h1-test (mean-of-ratios vs ratio-of-sums) ---")
    for nm, m in (("zero-effect genes", zmask), ("all scored genes", np.ones_like(zmask))):
        v = lfc0[m]
        print(f"  {nm:<20} n {v.size:>6}  median lfc0 {np.median(v):+.4f}  "
              f"DOWN {float((v < 0).mean()):.1%}")
    # does the SCORED pred_lfc of a zero-effect gene equal lfc0?
    pl = s0.pred_lfc[:, zmask]
    print(f"  scored pred_lfc on those genes vs lfc0: max |delta| "
          f"{np.nanmax(np.abs(pl - lfc0[zmask][None, :])):.4f}, "
          f"scored DOWN {float((pl < 0).mean()):.1%}")
    # and how often is the REFERENCE down on them, i.e. what precision does DOWN buy?
    rl = s0.real_lfc[:, zmask]
    rs = s0.real_sig[:, zmask]
    fin2 = np.isfinite(rl) & (rl != 0)
    print(f"  reference DOWN rate on those genes: all {float((rl[fin2] < 0).mean()):.1%}, "
          f"reference-DE only {float((rl[rs & fin2] < 0).mean()):.1%}\n")

    # ---- six members per arm -----------------------------------------------------------
    print("--- DE members, h1-test raw (rule 9: LEVELS not transferable, deltas/ordering are) ---")
    print(f"{'arm':<20}{'n_pred':>8}{'fid':>9}{'reach':>9}{'jac':>9}{'nmae':>9}{'pds':>9}")
    M = {}
    for a in ARMS:
        M[a] = members(S[a])
        n = float(np.median(S[a].counts()[1]))
        print(f"{a:<20}{n:>8.0f}" + "".join(f"{M[a][k]:>9.5f}"
              for k in ("fid", "reach", "jac", "nmae", "pds")))

    print("\n--- deltas vs a15 (seed-paired where a twin exists) ---")
    pairs = {"a15": "a15_b", "armA_tailphi0": "armA_tailphi0_b",
             "armC_a8_flat": "armC_a8_flat_b"}
    def avg(a):
        b = pairs.get(a)
        ks = ("fid", "reach", "jac", "nmae", "pds")
        return {k: (M[a][k] + M[b][k]) / 2 if b else M[a][k] for k in ks}
    A15 = avg("a15")
    print(f"{'arm':<20}{'d fid':>10}{'d reach':>10}{'d jac':>10}{'d nmae':>10}{'d pds':>10}")
    for a in ("armA_tailphi0", "armB_tail050", "armC_a8_flat"):
        v = avg(a)
        print(f"{a:<20}" + "".join(f"{v[k]-A15[k]:>+10.5f}"
              for k in ("fid", "reach", "jac", "nmae", "pds")))
    print(f"{'seed spread (a15)':<20}"
          + "".join(f"{abs(M['a15'][k]-M['a15_b'][k]):>10.5f}"
                    for k in ("fid", "reach", "jac", "nmae", "pds")))

    # ---- precision decomposition -------------------------------------------------------
    print("\n--- PRECISION DECOMPOSITION of the call set (pooled k/n over all perturbations) ---")
    print(f"{'arm':<20}{'block':<12}{'n calls':>10}{'share':>8}{'precision':>11}"
          f"{'med n/pert':>12}")
    D = {}
    for a in ARMS:
        D[a] = decompose(S[a], spread[a], zmask)
        tot = D[a]["ALL"][0]
        for blk in ("keep", "keep_zero", "tail_real", "tail_zero", "ALL"):
            n, k, prec, med = D[a][blk]
            if n == 0 and blk != "ALL":
                continue
            print(f"{a if blk == 'keep' else '':<20}{blk:<12}{n:>10d}"
                  f"{n/max(tot,1):>8.1%}{prec:>11.4f}{med:>12.0f}")

    # ---- arm D bound and the conserved-headroom projection ------------------------------
    print("\n--- arm D (post-hoc drop; NOT shippable, a bound only) ---")
    a15 = S["a15"]
    tail15 = spread["a15"] > 0.5
    for nm, drop in (("drop tail block", tail15),
                     ("drop zero-effect genes", np.broadcast_to(zmask, tail15.shape)),
                     ("drop tail AND zero", tail15 | zmask[None, :])):
        print(f"  {nm:<26} fid {fid_dropping(a15, drop):.5f}   (a15 {M['a15']['fid']:.5f})")

    print("\n--- board projection, 096's conserved-headroom form ---")
    print("  p_after - 0.5 = (p0 - 0.5) * n_before / n_after   [removed genes at EXACTLY 0.5]")
    BOARD_FID, BOARD_NPRED = 0.5022, 4400.0
    for a in ("armA_tailphi0", "armB_tail050", "armC_a8_flat"):
        r = np.median(S[a].counts()[1]) / np.median(S["a15"].counts()[1])
        p = 0.5 + (BOARD_FID - 0.5) / r
        print(f"  {a:<20} n_pred ratio {r:.4f} -> board n_pred {BOARD_NPRED*r:>7.0f}, "
              f"fid {p:.5f} (raw {p-BOARD_FID:+.5f}, scaled "
              f"{(p-BOARD_FID)/(ANCHOR['fid'][1]-ANCHOR['fid'][0]):+.5f}, "
              f"overall {(p-BOARD_FID)/(ANCHOR['fid'][1]-ANCHOR['fid'][0])/6:+.6f})")
    print("  NOTE: this is the OPTIMISTIC bound -- it credits the arm with removing only "
          "at-chance calls.")

    np.savez_compressed(f"{HERE}/decomp.npz", zmask=zmask, lfc0=lfc0,
                        cellcpm=cellcpm, ctrl_pb=ctrl_pb)


if __name__ == "__main__":
    main()
