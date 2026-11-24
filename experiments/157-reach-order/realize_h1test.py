"""157 — build the reference AND shipped-emission DE tables on the h1-test harness.

WHY THIS SURFACE.  Experiment 154 screened confidence keys on `orion-hct116-thinned`, and
094's A8-BOARD-RESULT records that this surface is **reach-blind**: production reads raw
reach 0.0838 there against 0.1701 on the board, while `h1-test` reads 0.1878 — within 11%
of the board.  088 §5 named h1-test the reach-representative harness for exactly this
reason.  So 154's headline was measured on the one surface whose reach magnitudes are known
not to transfer, and the closure has never been checked where they do.

Builds, on the SAME code path as `realize_shipped.py`:
  * the reference DE table (real_p / real_padj / real_lfc), which does not exist for h1-test;
  * the shipped a15 emission (tail 1.0 / keep_top 8000 / keep_spread 0.25 / phi 4) and its
    re-seeded twin, plus the pre-a8 phi=0 emitter as the known-answer control.

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
from vcc2026.donor_transfer import emit_counts        # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z        # noqa: E402

REAL = os.environ.get("REAL", f"{ROOT}/data/harness/h1-test.h5ad")
PRED = os.environ.get(
    "PRED", f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
            f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
OUT = os.environ.get("OUT", f"{ROOT}/experiments/157-reach-order/de_h1test")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000

ARMS = [
    ("prod_phi0_s0.25", None, 0, None, 0.0, 90_000),
    ("a15_tier8000_phi4", 1.0, 8000, 0.25, 4.0, 90_000),
    ("a15_tier8000_phi4_b", 1.0, 8000, 0.25, 4.0, 7_000_000),
    ("a8_flat_s0.25_phi4", None, 0, None, 4.0, 90_000),
    ("a8_flat_s0.25_phi4_b", None, 0, None, 4.0, 7_000_000),
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
        s = m.sum(1); s[s == 0] = 1.0
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
        print(f"gene-axis align: prediction {pg.size} vs harness {G_full}; "
              f"{int((PIDX < 0).sum())} harness genes absent from the prediction, dropped")
    genes = genes_full[keep]
    G = int(keep.sum())
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past the {MIN_CPM}-CPM gate "
          f"[{time.time()-t0:.0f}s]")

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
    print(f"control block ready [{time.time()-t0:.0f}s]")

    targets = [t for t in sorted(set(rlab) - {CTRL}) if (plab == t).sum() > 0]
    print(f"{len(targets)} perturbations, arms: {[a[0] for a in ARMS]}")
    T = len(targets)
    real = {k: np.full((T, G), np.nan) for k in ("p", "q", "l")}
    out = {n: {k: np.full((T, G), np.nan) for k in ("p", "q", "l")} for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}

    for i, t in enumerate(targets):
        Br = dense(R, np.flatnonzero(rlab == t))
        s = Br.sum(1); s[s == 0] = 1.0
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

        for name, tail, ktop, kspread, phi, seed in ARMS:
            if tail is None:
                sv = np.full(G_full, 0.25)
            else:
                sv = np.full(G_full, float(tail))
                k = min(int(ktop), int(gated.sum()))
                if k > 0:
                    sv[np.lexsort((np.arange(G_full), -score))[:k]] = float(kspread)
            rng = np.random.default_rng(seed + i)
            vals, cols, ptr = emit_counts(profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                                          mode="quasi", spread=sv, overdispersion=phi)
            B = np.zeros((N_EMIT, G_full))
            for r in range(N_EMIT):
                sl = slice(ptr[r], ptr[r + 1])
                B[r, cols[sl]] = vals[sl]
            assert np.array_equal(B.sum(0).astype(np.int64), totals), name
            s = B.sum(1); s[s == 0] = 1.0
            nb = np.log1p(B[:, keep] * (1e6 / s[:, None]))
            del B
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
            p = 2 * norm.sf(np.abs(z))
            q = false_discovery_control(p, method="bh")
            out[name]["p"][i], out[name]["q"][i] = p, q
            out[name]["l"][i] = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
            npred[name].append(int((q < 0.05).sum()))
            del nb
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{T} [{time.time()-t0:.0f}s]", flush=True)

    for name, *_ in ARMS:
        path = f"{OUT}.{name}.npz"
        np.savez_compressed(path, genes=genes, perts=np.array(targets),
                            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
                            pred_p=out[name]["p"], pred_padj=out[name]["q"],
                            pred_lfc=out[name]["l"])
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}")
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
