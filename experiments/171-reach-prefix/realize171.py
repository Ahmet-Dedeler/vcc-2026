"""171 — realize the emission arms the diagnosis points at, on h1-test.

Derived from `experiments/164-callset-artefact/realize164.py` (same code path, same donors,
same gate convention, so 157's board-validated `reach` calibration still applies).  Three
additions:

  * the ranking score `eff` = |log2(pred_cpm / ctrl_cpm)| is SAVED per target, so the block
    boundary is inspectable rather than inferred;
  * a sha256 of the DELIVERED per-gene panel total is recorded per (arm, target).  Every arm
    must produce the same digest.  This is the rule-12 pseudobulk bit-identity check in its
    non-vacuous form: if `emit_counts`'s exact-total repair were absent the Poisson residual
    would land in the totals and the digests would differ, so the check reads differently when
    the thing it tests is missing;
  * the arms are the ones 164 did NOT open -- 164 only moved the tail DOWN in dispersion
    (which is what kills `claim`).  These move `keep_spread` and `keep_top`, i.e. they change
    which block sits at the FRONT of the emitted p-order, and two of them RAISE total emitted
    dispersion rather than lowering it.

WHY THESE ARMS.  `diag171.py` measures that 61-79% of the head positions of the reach walk
are TAIL-block genes at per-perturbation sign accuracy 0.5212, while KEEP-block calls run
0.6629.  The tail wins the head because `spread 1.0 / phi 4` produces median |z| ~ 7 at zero
effect while `spread 0.25 / phi 4` produces ~1.7 (157 section 5b) -- the emitted significance
curve is U-shaped in `spread` and the shipped policy parks the SIGNAL block at its minimum.
Lowering `keep_spread` moves the keep block back to the high-|z| end WITHOUT removing any of
the tail's dispersion, which is the part `mse` is paid for.

KNOWN-ANSWER GATES (rule 10's corollary):
  * `a15` must reproduce median n_pred 6035 and reach 0.18283 (157/164, bit-identical);
  * `armC` (a8 flat) must reproduce median n_pred 3581 and reach 0.20291.
Both re-derived here from a fresh emission on this machine.

DONORS: 095's h1-test artifact -- `k562_gwps` + `orion_hek293t`.  `data/donors/h1_2025`
covers 85/85 h1-test targets and never enters.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import hashlib
import os
import sys
import time

import anndata as ad
import numpy as np
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

REAL = os.environ.get("REAL", f"{ROOT}/data/harness/h1-test.h5ad")
PRED = os.environ.get(
    "PRED",
    f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
    f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
)
OUT = os.environ.get("OUT", f"{ROOT}/experiments/171-reach-prefix/work/de171")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000

# name, tail_spread (None = flat keep_spread everywhere), keep_top, keep_spread,
#                                                          phi_keep, phi_tail, seed
ARMS = [
    # --- known-answer gates -----------------------------------------------------------
    ("a15", 1.0, 8000, 0.25, 4.0, 4.0, 90_000),
    ("armC_a8_flat", None, 0, 0.25, 4.0, 4.0, 90_000),
    # --- the new family: move the SIGNAL block to the high-|z| end of the U ------------
    ("ks000", 1.0, 8000, 0.00, 4.0, 4.0, 90_000),          # keep_spread 0
    ("ks000_b", 1.0, 8000, 0.00, 4.0, 4.0, 7_000_000),     # re-seeded twin
    ("ks010", 1.0, 8000, 0.10, 4.0, 4.0, 90_000),
    ("kt2000_ks000", 1.0, 2000, 0.00, 4.0, 4.0, 90_000),   # smaller protected block
    ("kt2000_ks025", 1.0, 2000, 0.25, 4.0, 4.0, 90_000),
    ("flat100_phi4", None, 0, 1.00, 4.0, 4.0, 90_000),     # everything at the tail law
    # --- ROUND 2: push the whole tail block OUT of the DE universe ---------------------
    # `keep_top` >= the gated count puts every DE-visible gene at `keep_spread` and leaves
    # `tail_spread` acting only on sub-gate genes, which no DE member can see.
    ("kt20000", 1.0, 20000, 0.25, 4.0, 4.0, 90_000),
    ("kt20000_b", 1.0, 20000, 0.25, 4.0, 4.0, 7_000_000),
    ("kt10000", 1.0, 10000, 0.25, 4.0, 4.0, 90_000),
    ("armC_a8_flat_b", None, 0, 0.25, 4.0, 4.0, 7_000_000),
    # --- ROUND 3: the KNEE arm.  claim171b puts the mse knee at challenge keep_top ~8330
    # of 8976 gated = 92.8% of the gated axis promoted; the fraction-matched h1-test value
    # is 0.928 * 10779 = 10000.  Two seeds, plus a bracketing point at 9000 (83.5%).
    ("kt10000_b", 1.0, 10000, 0.25, 4.0, 4.0, 7_000_000),
    ("kt9000", 1.0, 9000, 0.25, 4.0, 4.0, 90_000),
    ("kt9000_b", 1.0, 9000, 0.25, 4.0, 4.0, 7_000_000),
]
if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in _only]


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
        s = m.sum(1)
        s[s == 0] = 1.0
        cpm += (m * (1e6 / s[:, None])).sum(0)
    cpm /= ci.size
    keep = cpm > MIN_CPM

    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    pg = np.asarray([str(g) for g in P.var_names])
    if np.array_equal(pg, genes_full):
        PIDX = None
    else:
        pp = {g: i for i, g in enumerate(pg)}
        PIDX = np.array([pp.get(g, -1) for g in genes_full])
        keep &= PIDX >= 0
    genes = genes_full[keep]
    G = int(keep.sum())
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past the {MIN_CPM}-CPM gate "
          f"[{time.time()-t0:.0f}s]", flush=True)

    cs = np.empty((G, ci.size))
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    ref_lin = cpm[keep]
    b_ctrl = np.log1p(cpm)
    print(f"control block ready [{time.time()-t0:.0f}s]", flush=True)

    targets = [t for t in sorted(set(rlab) - {CTRL}) if (plab == t).sum() > 0]
    T = len(targets)
    print(f"{T} perturbations, arms: {[a[0] for a in ARMS]}", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    svec = {n: np.zeros((T, G), dtype=np.float32) for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}
    digest = {n: [] for n, *_ in ARMS}
    effs = np.zeros((T, G), dtype=np.float32)

    for i, t in enumerate(targets):
        Br = dense(R, np.flatnonzero(rlab == t))
        s = Br.sum(1)
        s[s == 0] = 1.0
        Br = np.log1p(Br[:, keep] * (1e6 / s[:, None]))
        z = wilcoxon_z(cs, tie, np.ascontiguousarray(Br.T))
        real["p"][i] = 2 * norm.sf(np.abs(z))
        real["q"][i] = false_discovery_control(real["p"][i], method="bh")
        real["l"][i] = np.log2((np.expm1(Br).mean(0) + 1e-9) / (ref_lin + 1e-9))
        del Br

        blk = dense(P, np.flatnonzero(plab == t))
        totals_p = blk.sum(0).astype(np.int64)
        del blk
        if PIDX is None:
            totals = totals_p
        else:
            totals = np.zeros(G_full, np.int64)
            m = PIDX >= 0
            totals[m] = totals_p[PIDX[m]]
        profile = (totals + 0.5) / (N_EMIT * EMIT_L)
        pcpm = 1e6 * totals / max(totals.sum(), 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        effs[i] = eff[keep]
        gated = (cpm > MIN_CPM) & (b_ctrl > 0) & (PIDX >= 0 if PIDX is not None else True)
        score = np.where(gated, eff, -np.inf)
        order = np.lexsort((np.arange(G_full), -score))

        for name, tail, ktop, kspread, phi_keep, phi_tail, seed in ARMS:
            if tail is None:
                sv = np.full(G_full, float(kspread))
                is_keep = np.ones(G_full, dtype=bool)
            else:
                sv = np.full(G_full, float(tail))
                is_keep = np.zeros(G_full, dtype=bool)
                k = min(int(ktop), int(gated.sum()))
                if k > 0:
                    sv[order[:k]] = float(kspread)
                    is_keep[order[:k]] = True
            od = float(phi_keep) if phi_keep == phi_tail else np.where(
                is_keep, float(phi_keep), float(phi_tail))
            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(
                profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                mode="quasi", spread=sv, overdispersion=od,
            )
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            delivered = B.sum(0).astype(np.int64)
            # rule-12 pseudobulk identity, in the form that FAILS if the repair is absent
            digest[name].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            assert np.array_equal(delivered, totals), name
            s = B.sum(1)
            s[s == 0] = 1.0
            nb = np.log1p(B[:, keep] * (1e6 / s[:, None]))
            del B
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
            p = 2 * norm.sf(np.abs(z))
            q = false_discovery_control(p, method="bh")
            out[name]["p"][i], out[name]["q"][i] = p, q
            out[name]["l"][i] = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
            svec[name][i] = sv[keep]
            npred[name].append(int((q < 0.05).sum()))
            del nb
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{T} [{time.time()-t0:.0f}s]", flush=True)

    print()
    print("PSEUDOBULK BIT-IDENTITY (sha256 of the delivered per-gene panel total, per target)")
    ref = digest[ARMS[0][0]]
    for name, *_ in ARMS:
        same = sum(a == b for a, b in zip(digest[name], ref))
        print(f"  {name:<16} identical to {ARMS[0][0]} on {same}/{T} targets   "
              f"first digest {digest[name][0]}")
    print()

    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            spread=svec[name], eff=effs,
        )
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}", flush=True)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
