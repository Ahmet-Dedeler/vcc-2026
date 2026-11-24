"""164 — realize the incumbent and the three call-set arms on h1-test, saving the BLOCK MASK.

Derived from `experiments/157-reach-order/realize_h1test.py`, same code path, same donors,
same gate convention (so 157's board-validated `reach` calibration still applies).  Three
additions:

  * arm A -- per-gene `overdispersion`: `phi = 4` on the keep block, `phi = 0` on the tail.
    `emit_counts` has accepted a per-gene phi since the 2026-08-27 freeze
    (`donor_transfer.py:403-411`), so this needs no builder change.
  * arm B -- `tail_spread = 0.5`.
  * the per-target spread vector is SAVED, so the call set can be split into the tail block
    and the keep block exactly rather than by reconstructing the `|effect|` ranking.

Known-answer rows, per AGENTS.md 11.1 rule 10's corollary:
  * `a15` must reproduce 157's median `n_pred` 6035 and its `reach` 0.1828;
  * `a8` (arm C) must reproduce 157's median `n_pred` and its `reach` 0.2029.
Both are re-derived here from a fresh emission, so agreement is a real check.

DONORS: the prediction is 095's h1-test artifact, whose donors are `k562_gwps` +
`orion_hek293t`.  `data/donors/h1_2025` covers 85/85 h1-test targets and must never enter.

No cell-eval2, no build, no submission.
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
from vcc2026.donor_transfer import emit_counts  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

REAL = os.environ.get("REAL", f"{ROOT}/data/harness/h1-test.h5ad")
PRED = os.environ.get(
    "PRED",
    f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
    f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
)
OUT = os.environ.get("OUT", f"{ROOT}/experiments/164-callset-artefact/de164")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000

# name, tail_spread (None = flat SPREAD), keep_top, keep_spread, phi_keep, phi_tail, seed
ARMS = [
    ("a15", 1.0, 8000, 0.25, 4.0, 4.0, 90_000),
    ("a15_b", 1.0, 8000, 0.25, 4.0, 4.0, 7_000_000),
    ("armA_tailphi0", 1.0, 8000, 0.25, 4.0, 0.0, 90_000),
    ("armA_tailphi0_b", 1.0, 8000, 0.25, 4.0, 0.0, 7_000_000),
    ("armB_tail050", 0.5, 8000, 0.25, 4.0, 4.0, 90_000),
    ("armC_a8_flat", None, 0, None, 4.0, 4.0, 90_000),
    ("armC_a8_flat_b", None, 0, None, 4.0, 4.0, 7_000_000),
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
    # the emitted spread on the SCORED axis, per target per arm -- this is the block label
    svec = {n: np.zeros((T, G), dtype=np.float32) for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}

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
        gated = (cpm > MIN_CPM) & (b_ctrl > 0) & (PIDX >= 0 if PIDX is not None else True)
        score = np.where(gated, eff, -np.inf)
        order = np.lexsort((np.arange(G_full), -score))

        for name, tail, ktop, kspread, phi_keep, phi_tail, seed in ARMS:
            if tail is None:
                sv = np.full(G_full, 0.25)
                is_keep = np.ones(G_full, dtype=bool)
            else:
                sv = np.full(G_full, float(tail))
                is_keep = np.zeros(G_full, dtype=bool)
                k = min(int(ktop), int(gated.sum()))
                if k > 0:
                    sv[order[:k]] = float(kspread)
                    is_keep[order[:k]] = True
            if phi_keep == phi_tail:
                od = float(phi_keep)
            else:
                od = np.where(is_keep, float(phi_keep), float(phi_tail))
            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(
                profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                mode="quasi", spread=sv, overdispersion=od,
            )
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            assert np.array_equal(B.sum(0).astype(np.int64), totals), name
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

    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            spread=svec[name],
        )
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}", flush=True)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
