"""157 Q1/Q3 — realise the SHIPPED (a15) emission law and save its DE table.

`experiments/095-clean-harness/de_thinned.npz` is the only repaired-harness table with a
prediction side, and its prediction was emitted at **spread 0.25, phi 0** (grid.json:
`g0.65_q1.0_a1.08_spread0.25`, built 2026-08-27 20:32, before a8 shipped over-dispersion).
Every ordering statistic read off it therefore describes the PRE-a8 emitter.

a15 emits with `--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`.
This drives the UNMODIFIED production `emit_counts` at that setting, on the same profile and
the same repaired harness, and writes a full (pred_p, pred_padj, pred_lfc) table so the two
emission laws can be compared as a clean single-variable change.

The profile round-trips exactly: `emit_counts(mode="quasi")` preserves each gene's panel
total, so the totals read off the existing pred.h5ad ARE the profile, and
`floor(n*L*(T+0.5)/(n*L)) == T`.

Memory: controls are streamed in chunks; peak is the gated sorted control block
(~10.3k genes x 18.4k cells float64 = 1.5 GB).  No cell-eval2, no build, no submission.
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
OUT = os.environ.get("OUT", f"{ROOT}/experiments/157-reach-order/de_shipped.npz")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000
NP = int(os.environ.get("N_PERTS", "132"))

# (name, tail_spread, keep_top, keep_spread, phi, seed)
# `prod_phi0` reproduces de_thinned's own prediction and is the known-answer row
# (AGENTS.md 11.1 rule 10 corollary): if it does not land near reach 0.077 the
# instrument is wrong and nothing else here is readable.
ARMS = [
    ("prod_phi0_s0.25", None, 0, None, 0.0, 90_000),
    ("a15_tier8000_phi4", 1.0, 8000, 0.25, 4.0, 90_000),
    ("a15_tier8000_phi4_b", 1.0, 8000, 0.25, 4.0, 7_000_000),
    ("a8_flat_s0.25_phi4", None, 0, None, 4.0, 90_000),
]


def chunked_control_stats(A, rows, G, chunk=2000):
    """CPM per gene over control cells, streamed."""
    acc = np.zeros(G)
    for i in range(0, rows.size, chunk):
        m = A.X[rows[i:i + chunk]]
        m = np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)
        s = m.sum(1)
        s[s == 0] = 1.0
        acc += (m * (1e6 / s[:, None])).sum(0)
    return acc / rows.size


def main() -> None:
    t0 = time.time()
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes_full = np.asarray([str(g) for g in R.var_names])
    G_full = genes_full.size
    ci = np.flatnonzero(rlab == CTRL)
    cpm = chunked_control_stats(R, ci, G_full)
    keep = cpm > MIN_CPM
    genes = genes_full[keep]
    G = int(keep.sum())
    print(f"harness: {ci.size} controls, {G}/{G_full} genes past the {MIN_CPM}-CPM gate "
          f"[{time.time()-t0:.0f}s]")

    # gated, log1p-CPM control block, transposed and sorted -- the Wilcoxon reference
    cs = np.empty((G, ci.size), dtype=np.float64)
    for i in range(0, ci.size, 2000):
        m = R.X[ci[i:i + 2000]]
        m = np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)
        s = m.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    ref_lin = cpm[keep]
    b_ctrl = np.log1p(cpm)          # the builder's own control profile, full axis
    print(f"control block ready [{time.time()-t0:.0f}s]")

    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    pg = np.asarray([str(g) for g in P.var_names])
    assert np.array_equal(pg, genes_full), "prediction and harness gene axes differ"
    targets = [t for t in sorted(set(rlab) - {CTRL}) if (plab == t).sum() > 0][:NP]
    print(f"{len(targets)} perturbations, arms: {[a[0] for a in ARMS]}")

    out = {n: {k: np.full((len(targets), G), np.nan) for k in ("p", "q", "l")}
           for n, *_ in ARMS}
    npred = {n: [] for n, *_ in ARMS}
    for i, t in enumerate(targets):
        rows = np.flatnonzero(plab == t)
        blk = P.X[rows]
        blk = np.asarray(blk.todense() if hasattr(blk, "todense") else blk, dtype=np.float64)
        totals = blk.sum(0).astype(np.int64)
        del blk
        profile = (totals + 0.5) / (N_EMIT * EMIT_L)

        # |predicted effect| on the full axis, buildable: our own emitted pseudobulk CPM
        # against the control's. Nothing from the reference enters it.
        pcpm = 1e6 * totals / max(totals.sum(), 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        gated = (cpm > MIN_CPM) & (b_ctrl > 0)
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
            s = B.sum(1)
            s[s == 0] = 1.0
            nb = np.log1p(B[:, keep] * (1e6 / s[:, None]))
            del B
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
            p = 2 * norm.sf(np.abs(z))
            q = false_discovery_control(p, method="bh")
            out[name]["p"][i] = p
            out[name]["q"][i] = q
            out[name]["l"][i] = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
            npred[name].append(int((q < 0.05).sum()))
            del nb
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{len(targets)} [{time.time()-t0:.0f}s]", flush=True)

    ref = np.load(f"{ROOT}/experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
    assert np.array_equal(ref["genes"].astype(str), genes), "gene axis drift vs de_thinned"
    rperts = list(ref["perts"].astype(str))
    assert rperts[:len(targets)] == targets, "pert axis drift vs de_thinned"
    sl = slice(0, len(targets))
    for name, *_ in ARMS:
        path = OUT.replace(".npz", f".{name}.npz")
        np.savez_compressed(
            path, genes=genes, perts=np.array(targets),
            real_p=ref["real_p"][sl], real_padj=ref["real_padj"][sl],
            real_lfc=ref["real_lfc"][sl],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"])
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}")
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
