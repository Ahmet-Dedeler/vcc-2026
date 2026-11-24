"""224 — the production question: does the per-gene gain survive spectral whitening?

Production does NOT ship the equal-weight donor mean; it ships a whitened consensus
(AGENTS.md §4.2, `src/vcc2026/whiten.py`, `P -> P (C+lam I)^(-gamma/2)`, gamma 0.65,
lam 0).  A gain measured against the RAW consensus is only a live lever if it is not
something whitening already collects.

The shipped operator's own gene axis is `union(donor coverage) - panel` and is not
recorded in its sidecar json (AGENTS.md §8.3 is precisely about that), so the shipped
operator cannot be restricted to this surface's gene axis without rebuilding it.  What is
done instead is an exact re-fit of the SAME family on THIS surface, from the same kind of
basis (donor effect rows, held line excluded), which is what §4.2 prescribes and what
`build_union_whitener.py` does.

Four arms, all evaluated on held-line PANEL targets:
    A  raw consensus                     P
    B  raw consensus x gain              P * g
    C  whitened consensus                P W
    D  whitened consensus x gain         P W * g_w     (g_w refit on WHITENED inner pairs)

`pds` matters more than cosine here: §4.2 measures whitening as LOWERING transfer cosine
32% while RAISING pds 8%, so a cosine-only read of arm C would be a category error.

    uv run python whiten224.py
"""
from __future__ import annotations

import itertools
import json
import sys

import numpy as np

from core224 import (DIAG_SHRINK, WORK, decoy_cos, diag_fit, donor_mean,  # noqa: E402
                     inner_pairs, macro_cos, pds, unit)
from crosspanel224 import FIT_LINES, PANEL_LINES, genes, meta, perts
from substrate224 import PANEL_SET, ROOT, load_effects

sys.path.insert(0, str(ROOT / "src"))
from vcc2026.whiten import Whitener  # noqa: E402

SEED = 20260904
GAMMAS = (0.0, 0.3, 0.65, 1.0)


def main() -> None:
    rng = np.random.default_rng(SEED)
    panel_lines = PANEL_LINES + ["cd4_rest"]
    allnames = sorted(set(panel_lines) | set(FIT_LINES))
    gaxis = sorted(set.intersection(*[genes(n) for n in allnames]) - PANEL_SET)
    fit_targets = sorted(set.intersection(*[perts(n) for n in FIT_LINES]) - PANEL_SET)
    eval_targets = sorted(set.intersection(*[perts(n) for n in panel_lines]) & PANEL_SET)
    G = len(gaxis)
    gene_bins = np.arange(G)
    print(f"gene axis {G}; fit {len(fit_targets)} non-panel; eval {len(eval_targets)} panel")

    Cf = {n: (lambda v: v - v.mean(0, keepdims=True))(
        load_effects(meta(n)["path"], fit_targets, gaxis)) for n in FIT_LINES}
    Ce = {n: (lambda v: v - v.mean(0, keepdims=True))(
        load_effects(meta(n)["path"], eval_targets, gaxis)) for n in panel_lines}
    Cea = np.stack([Ce[n] for n in panel_lines])

    q = rng.permutation(len(fit_targets))
    tsel, tval = np.sort(q[: int(0.6 * len(q))]), np.sort(q[int(0.6 * len(q)):])
    rows = []
    for A in panel_lines:
        F = [n for n in FIT_LINES if n != A]
        Cfa = np.stack([Cf[n] for n in F])
        # the operator is fitted on the SAME rows the gain is, held line excluded
        basis = np.concatenate([Cfa[i] for i in range(len(F))])
        W = Whitener.fit(basis)
        Pa, Ya, _ = inner_pairs(Cfa, list(range(len(F))), np.arange(len(fit_targets)))
        Ps, Ys, Fs = inner_pairs(Cfa, list(range(len(F))), tsel)
        Pv, Yv, Fv = inner_pairs(Cfa, list(range(len(F))), tval)
        hold = int(rng.integers(len(F)))
        ms, mv = Fs != hold, Fv == hold
        ai = panel_lines.index(A)
        T = Cea[ai]
        pool = [panel_lines.index(n) for n in panel_lines if n != A]
        for gamma in GAMMAS:
            def wh(X):
                return X if gamma == 0.0 else unit(W.apply(X, gamma=gamma, renorm=False))
            sh = max(DIAG_SHRINK, key=lambda s: macro_cos(
                wh(Pv[mv]) * diag_fit(wh(Ps[ms]), Ys[ms], gene_bins, G, s), Yv[mv]))
            gw = diag_fit(wh(Pa), Ya, gene_bins, G, sh)
            for k in range(1, len(pool) + 1):
                for D in itertools.combinations(pool, k):
                    P = donor_mean(Cea, list(D), np.arange(len(eval_targets)))
                    Pw = wh(P)
                    rows.append(dict(
                        held=A, gamma=gamma, k=k,
                        donors=[panel_lines[i] for i in D],
                        cos_plain=macro_cos(Pw, T), cos_gain=macro_cos(Pw * gw, T),
                        pds_plain=pds(Pw, T), pds_gain=pds(Pw * gw, T),
                        decoy_plain=decoy_cos(Pw, T), decoy_gain=decoy_cos(Pw * gw, T),
                        shrink=float(sh)))
        print(f"  held {A:16s} operator fitted on {len(F)} lines x "
              f"{len(fit_targets)} rows; max eval {W.evals.max():.4g}")

    print(f"\n{'gamma':>6}{'k':>3}{'n':>4}{'cos':>10}{'cos+gain':>10}{'rel%':>8}"
          f"{'pds':>9}{'pds+gain':>10}{'dpds':>9}")
    out = {}
    for gamma in GAMMAS:
        for k in sorted({r["k"] for r in rows}):
            v = [r for r in rows if r["gamma"] == gamma and r["k"] == k]
            cp = np.mean([r["cos_plain"] for r in v])
            cg = np.mean([r["cos_gain"] for r in v])
            pp = np.mean([r["pds_plain"] for r in v])
            pg = np.mean([r["pds_gain"] for r in v])
            out[f"{gamma}|{k}"] = dict(n=len(v), cos=cp, cos_gain=cg,
                                       rel=100 * (cg - cp) / cp, pds=pp, pds_gain=pg)
            print(f"{gamma:6.2f}{k:3d}{len(v):4d}{cp:10.5f}{cg:10.5f}"
                  f"{100 * (cg - cp) / cp:+8.2f}{pp:9.4f}{pg:10.4f}{pg - pp:+9.4f}")

    print("\nWHAT §4.2 PREDICTS AND WHAT THIS SURFACE SHOWS "
          "(whitening: cosine DOWN, pds UP)")
    for k in sorted({r["k"] for r in rows}):
        a = out[f"0.0|{k}"]
        b = out[f"0.65|{k}"]
        print(f"  k={k}: cos {a['cos']:.5f} -> {b['cos']:.5f} "
              f"({100 * (b['cos'] - a['cos']) / a['cos']:+.1f}%)   "
              f"pds {a['pds']:.4f} -> {b['pds']:.4f} ({b['pds'] - a['pds']:+.4f})")

    (WORK / "whiten224.json").write_text(json.dumps(
        dict(gene_axis=G, n_fit=len(fit_targets), n_eval=len(eval_targets),
             gammas=list(GAMMAS), summary=out, rows=rows), indent=1))


if __name__ == "__main__":
    main()
