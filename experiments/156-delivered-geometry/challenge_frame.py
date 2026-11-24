"""156c — the mass ratio on the CHALLENGE frame, and a real verification that a23 shipped it.

TWO JOBS
--------
1. Every number in `RESULTS.md` was measured on the two held Orion artifacts. The mass ratio
   there is 1.00181 / 1.00168 and the prediction in `A23-PREDICTION.md` is scaled off it. It has
   never been measured on contexts A/B/C, which is what the board scores. Measure it.

2. Verify the shipped a23 artifact actually carries the fix. The naive check -- recover
   `delta_hat = log1p(TS*p) - b_ctrl` from the emitted counts and test
   `sum(expm1(b_ctrl + delta_hat)) == TS` -- is a TAUTOLOGY: `expm1(log1p(TS*p)) = TS*p`, which
   sums to TS for any normalised `p` whatsoever. A pre-fix build (v9) reads exactly 1.00000000
   on it too. So the artifact alone cannot distinguish, and the only real check is to rebuild
   the intended effect matrix and compare both renderings against what was emitted.

This script therefore reproduces the a15/a23 effect stage exactly (same donors, same operator,
same PerturbMatch, same CD4 solve) and requires `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10`.
It does NOT emit counts and does NOT write any submission.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[2]
GATE_CPM = 10.0


def load_builder():
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_156c", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_156c"] = mod
    spec.loader.exec_module(mod)
    return mod


def context_bctrl(M, bundle: Path, context: str):
    """The builder's own control profile for one context (run() lines ~1975-1982)."""
    import anndata as ad
    import scipy.sparse as sp
    spec = M.ChallengeSpec.load(bundle)
    controls = ad.read_h5ad(bundle / f"context_{context}.h5ad")
    X = controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
    totals = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    return np.log1p(M.TS * totals / totals.sum()), spec


def emitted_profile(path: Path, context: str, targets, genes_n: int):
    """Per-gene emitted totals -> normalised profile, for named targets in one context."""
    out = {}
    with h5py.File(path, "r") as h:
        ip = h["X/indptr"][:]
        tg = h["obs/target_gene"]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in tg["categories"][:]])
        codes = tg["codes"][:]
        cx = h["obs/context"]
        ccats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                            for c in cx["categories"][:]])
        ccodes = cx["codes"][:]
        want = int(np.flatnonzero(ccats == context)[0])
        inctx = ccodes == want
        for t in targets:
            code = int(np.flatnonzero(cats == t)[0])
            rows = np.flatnonzero(inctx & (codes == code))
            acc = np.zeros(genes_n)
            r0, r1 = int(rows.min()), int(rows.max()) + 1
            assert r1 - r0 == len(rows), "target rows are not contiguous"
            s, e = int(ip[r0]), int(ip[r1])
            np.add.at(acc, h["X/indices"][s:e], h["X/data"][s:e].astype(np.float64))
            out[t] = acc / acc.sum()
    return out


def main() -> None:
    assert os.environ.get("VCC_ALPHA") == "0.982", "set VCC_ALPHA=0.982"
    assert os.environ.get("VCC_CD4_TARGET_F") == "0.10", "set VCC_CD4_TARGET_F=0.10"
    assert not os.environ.get("VCC_MASS_FIX"), "leave VCC_MASS_FIX unset here"
    M = load_builder()
    bundle = REPO / "data/challenge"

    print("rebuilding the a15/a23 effect matrix (this is the 17 GiB stage)...", flush=True)
    spec = M.ChallengeSpec.load(bundle)
    donors = {n: M.Donor.load(REPO / p) for n, p in {
        M.K562: "data/donors/k562_gwps",
        "orion_hct116": "data/donors/orion_hct116",
        "orion_hek293t": "data/donors/orion_hek293t"}.items()}
    model = M.MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)
    weights = M.donor_weights(model)
    pm = M.EXP061.prepare_perturbmatch(
        REPO / "data/external/perturbmatch", tuple(map(str, spec.genes)),
        chunk_rows=20000, max_rss_gib=40.0)
    operator = M.load_production_operator(
        model, tuple(map(str, spec.perts)),
        REPO / "data/donors/whitener_union_pruned.npz",
        REPO / "data/donors/whitener_union_pruned.json",
        expected_sha256=M.WHITENER_REGISTRY["union_pruned"][0])
    panel_positions = np.asarray(
        [{str(g): i for i, g in enumerate(spec.genes)}[str(t)] for t in spec.perts],
        dtype=np.int64)
    _, cd4 = M.prepare_frozen_cd4(
        REPO / "data/external/cd4t/cloud-compact/cd4_compact", tuple(map(str, spec.perts)),
        tuple(map(str, spec.genes)), tuple(map(str, spec.perts)),
        chunk_rows=20000, max_rss_gib=40.0)
    effects, audit = M.build_production_effects(
        model, tuple(map(str, spec.perts)), operator.whitener, operator.columns,
        panel_positions, pm=pm, cd4=cd4, weights=weights)
    print(f"  effect matrix {effects.shape}, post_whiten_sha {audit['post_whiten_sha256'][:16]}, "
          f"after_all_median_norm {audit['whitening']['after_all_median_norm']!r}")
    del donors, model, pm
    perts = [str(p) for p in spec.perts]

    print("\n=== 1. mass ratio on the CHALLENGE frame, after a15's CPM>10 gate ===")
    ratios = {}
    for ctx in ("A", "B", "C"):
        b, _ = context_bctrl(M, bundle, ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / expr.sum()
        gate = (cpm > GATE_CPM) & (b > 0)
        G = np.where(gate[None, :], effects, 0.0)
        lin = np.clip(np.expm1(b[None, :] + G), 0.0, None)
        r = lin.sum(1) / M.TS
        resid = np.median(np.linalg.norm(
            np.log1p(M.TS * lin / lin.sum(1, keepdims=True)) - b[None, :] - G, axis=1)
            / np.linalg.norm(G, axis=1))
        ratios[ctx] = (r, gate, b, G)
        print(f"  context {ctx}: gated genes {int(gate.sum()):5d} | mass ratio mean "
              f"{r.mean():.6f} sd {r.std():.6f} min {r.min():.6f} max {r.max():.6f} | "
              f"render residual {resid * 100:.2f}% of row norm")
    print("  (held Orion artifacts, for comparison: 1.00181 / 1.00168 and 3.96% / 4.53%)")

    print("\n=== 2. does the shipped a23 artifact carry the fix? ===")
    ctx = "A"
    r, gate, b, G = ratios[ctx]
    probe = perts[:8]
    idx = {p: i for i, p in enumerate(perts)}
    emitted = emitted_profile(REPO / "submissions/a23-massfix.h5ad", ctx, probe, effects.shape[1])
    print(f"  {'target':<12} {'||p_emit - p_FIXED||1':>22} {'||p_emit - p_UNFIXED||1':>24} verdict")
    ok = 0
    for t in probe:
        g = G[idx[t]]
        p_unfixed = M.profile_from_effect(effects[idx[t]], b, min_cpm=GATE_CPM)
        h = M.mass_consistent_effect(g, b, gate, mode="minnorm")
        lin = np.clip(np.expm1(b + h), 0.0, None)
        p_fixed = lin / lin.sum()
        pe = emitted[t]
        d_fix = float(np.abs(pe - p_fixed).sum())
        d_un = float(np.abs(pe - p_unfixed).sum())
        good = d_fix < d_un
        ok += good
        print(f"  {t:<12} {d_fix:22.3e} {d_un:24.3e} "
              f"{'FIXED (ratio %.0fx closer)' % (d_un / max(d_fix, 1e-300)) if good else 'UNFIXED'}")
    print(f"\n  {ok}/{len(probe)} probed targets match the FIXED rendering more closely.")
    print("  Note the residual against p_FIXED is integer-emission rounding, not the renderer:")
    print("  the emitter pins per-gene totals to round(400*L*p), so ~1e-6 of L1 mass is expected.")


if __name__ == "__main__":
    main()
