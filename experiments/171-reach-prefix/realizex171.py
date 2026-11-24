"""171x — the EXPRESSION-COMPENSATING spread schedule, realized on h1-test.

The diagnosis (RESULTS.md section 2) says the head of the reach walk is 61-79% zero-effect
genes at chance accuracy.  The coordinator's falsifier (`INBOX-4`, reproduced in
`work/falsifier.log` to the digit) says WHY: our emitted p-order tracks control EXPRESSION at
Spearman +0.485 inside the call set and ANTI-tracks our own predicted effect at -0.391.  So
`keep_top`, which tiers spread by |effect| rank, cannot fix the order -- the disorder is
caused by expression, and 093 / 095 / the a15 tier all tiered on the wrong variable.

MECHANISM.  `_quasi_block` (`donor_transfer.py:274`) gives each gene a per-cell variance of
exactly

    sigma^2_ours = s*m + phi*s^2*m^2          m = per-cell mean count = totals/400

against a control side ~ m*(1 + od_ctrl).  So the rank test's signal-to-noise runs

    z/h  ~  sqrt(m) / sqrt(A + s + phi*s^2*m),        A = 1 + od_ctrl

which RISES with m through the low regime and plateaus at 1/(s*sqrt(phi)) once
phi*s^2*m >> 1.  At the shipped s = 0.25, phi = 4 the crossover is m* ~ 1/(phi s^2) = 4
counts/cell ~ 200 CPM, and essentially the whole 5-CPM-gated universe sits below it.  That is
the +0.485.

THE SCHEDULE.  z/h is monotonically DECREASING in s, and the control's Poisson noise is a
floor we cannot lower, so flattening means pulling the high-m plateau DOWN to the low-m level
by RAISING s on high-m genes.  Solve `sqrt(m)/sqrt(A + s + phi s^2 m) = c` for s:

    (phi c^2 m) s^2  +  (c^2) s  +  (c^2 A - m)  =  0
    s = ( -c^2 + sqrt( c^4 - 4 phi c^2 m (c^2 A - m) ) ) / (2 phi c^2 m)

clipped to [S_LO, S_HI].  `c` is set as a quantile of each gene's own best-achievable
`sqrt(m)/sqrt(A + S_LO)`, so a stated fraction of gated genes is exactly feasible and the
weakest ones simply stay at S_LO.

WHY THIS IS THE ONE ARM `mse` CANNOT VETO.  The schedule only ever RAISES s above S_LO = 0.25.
RESULTS.md section 5 measured that a15's `claim` (0.24581) sits +0.01313 ABOVE the knee
s^2 B = 0.23268, and that every policy which raises dispersion lands at the reference cap
0.24653 with raw `mse` EXACTLY a15's.  Every arm 164 and this experiment killed on `mse` did
so by LOWERING s.  This one cannot.  Still measured, not asserted (`claimx171.py`).

THE CONTROL, and what it reads if the mechanism is absent.  `xcomp10_perm` uses the SAME
multiset of per-gene spreads, randomly permuted across the gated axis, so the emitted
dispersion budget and the spread histogram are identical and only the gene<->spread PAIRING
is destroyed.  If the gain came from "having a spread distribution" rather than from
compensating expression, the permuted arm would score the same.  It is the falsifiable half.

KNOWN-ANSWER GATES: `a15` -> reach 0.18283, median n_pred 6035 (re-emitted here).

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

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = (f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
OUT = f"{ROOT}/experiments/171-reach-prefix/work/de171"
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
PHI, S_LO, S_HI, A_CTRL = 4.0, 0.25, 1.0, 2.0


def xcomp_spread(m, gated, q, rng=None, permute=False):
    """Expression-compensating per-gene spread.  `m` = per-cell mean count."""
    sv = np.full(m.size, S_HI)                       # sub-gate genes keep the tail law
    g = np.flatnonzero(gated)
    mg = np.maximum(m[g], 1e-12)
    c = float(np.quantile(np.sqrt(mg) / np.sqrt(A_CTRL + S_LO), q))
    c2 = c * c
    disc = c2 * c2 - 4.0 * PHI * c2 * mg * (c2 * A_CTRL - mg)
    s = np.full(mg.size, S_LO)
    ok = disc >= 0
    s[ok] = (-c2 + np.sqrt(disc[ok])) / (2.0 * PHI * c2 * mg[ok])
    s = np.clip(s, S_LO, S_HI)
    if permute:
        s = s[rng.permutation(s.size)]               # same multiset, pairing destroyed
    sv[g] = s
    return sv


def rank_spread(m, gated):
    """Crude control: s linear in EXPRESSION RANK from S_LO (lowest m) to S_HI (highest)."""
    sv = np.full(m.size, S_HI)
    g = np.flatnonzero(gated)
    o = np.lexsort((g, m[g]))
    r = np.empty(g.size)
    r[o] = np.arange(g.size) / max(g.size - 1, 1)
    sv[g] = S_LO + (S_HI - S_LO) * r
    return sv


# name -> callable(m, gated, seed_rng) -> spread vector  (None = the a15 gate arm)
ARMS = {
    "a15":          (None, 90_000),
    "xcomp10":      (lambda m, g, r: xcomp_spread(m, g, 0.10), 90_000),
    "xcomp10_b":    (lambda m, g, r: xcomp_spread(m, g, 0.10), 7_000_000),
    "xcomp25":      (lambda m, g, r: xcomp_spread(m, g, 0.25), 90_000),
    "xcomp10_perm": (lambda m, g, r: xcomp_spread(m, g, 0.10, r, True), 90_000),
    "sgrad_expr":   (lambda m, g, r: rank_spread(m, g), 90_000),
}
if os.environ.get("ONLY"):
    ARMS = {k: v for k, v in ARMS.items() if k in os.environ["ONLY"].split(",")}


def dense(A, rows):
    x = A.X[rows]
    return np.asarray(x.todense() if hasattr(x, "todense") else x, dtype=np.float64)


def main() -> None:
    t0 = time.time()
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes_full = np.asarray([str(g) for g in R.var_names])
    G_full = genes_full.size
    ci = np.flatnonzero(rlab == CTRL)

    cpm = np.zeros(G_full)
    for i in range(0, ci.size, 2000):
        x = dense(R, ci[i:i + 2000])
        s = x.sum(1)
        s[s == 0] = 1.0
        cpm += (x * (1e6 / s[:, None])).sum(0)
    cpm /= ci.size
    keep = cpm > MIN_CPM

    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    pg = np.asarray([str(g) for g in P.var_names])
    pp = {g: i for i, g in enumerate(pg)}
    PIDX = np.array([pp.get(g, -1) for g in genes_full])
    keep &= PIDX >= 0
    genes, G = genes_full[keep], int(keep.sum())
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past {MIN_CPM} CPM "
          f"[{time.time()-t0:.0f}s]", flush=True)

    cs = np.empty((G, ci.size))
    for i in range(0, ci.size, 2000):
        x = dense(R, ci[i:i + 2000])
        s = x.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + x.shape[0]] = np.log1p(x[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    ref_lin, b_ctrl = cpm[keep], np.log1p(cpm)
    targets = [t for t in sorted(set(rlab) - {CTRL}) if (plab == t).sum() > 0]
    T = len(targets)
    print(f"{T} perturbations, arms: {list(ARMS)} [{time.time()-t0:.0f}s]", flush=True)

    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n in ARMS}
    svec = {n: np.zeros((T, G), np.float32) for n in ARMS}
    npred = {n: [] for n in ARMS}
    digest = {n: [] for n in ARMS}
    effs = np.zeros((T, G), np.float32)

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
        tp = blk.sum(0).astype(np.int64)
        del blk
        totals = np.zeros(G_full, np.int64)
        mm = PIDX >= 0
        totals[mm] = tp[PIDX[mm]]
        profile = (totals + 0.5) / (N_EMIT * EMIT_L)
        pcpm = 1e6 * totals / max(totals.sum(), 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        effs[i] = eff[keep]
        gated = (cpm > MIN_CPM) & (b_ctrl > 0) & (PIDX >= 0)
        m_cell = np.floor(N_EMIT * EMIT_L * profile) / N_EMIT       # exactly what emit_counts lays down
        order = np.lexsort((np.arange(G_full), -np.where(gated, eff, -np.inf)))

        for name, (fn, seed) in ARMS.items():
            if fn is None:                                          # a15 as shipped
                sv = np.full(G_full, 1.0)
                sv[order[:min(8000, int(gated.sum()))]] = 0.25
            else:
                sv = fn(m_cell, gated, np.random.default_rng(11 + i))
            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                                          mode="quasi", spread=sv, overdispersion=PHI)
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            deliv = B.sum(0).astype(np.int64)
            digest[name].append(hashlib.sha256(deliv.tobytes()).hexdigest()[:16])
            assert np.array_equal(deliv, totals), name
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

    print("\nPSEUDOBULK BIT-IDENTITY (sha256 of the delivered per-gene panel total)")
    ref = digest[list(ARMS)[0]]
    for n in ARMS:
        print(f"  {n:<14} identical on {sum(a == b for a, b in zip(digest[n], ref))}/{T} "
              f"targets   first {digest[n][0]}")
    print()
    for n in ARMS:
        np.savez_compressed(
            f"{OUT}.{n}.npz", genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[n]["p"], pred_padj=out[n]["q"], pred_lfc=out[n]["l"],
            spread=svec[n], eff=effs)
        print(f"wrote de171.{n}.npz   median n_pred {np.median(npred[n]):.0f}   "
              f"median spread {np.median(svec[n]):.3f}", flush=True)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
