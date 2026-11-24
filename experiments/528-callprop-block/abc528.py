#!/usr/bin/env python3
"""528 — transfer checks on the challenge contexts: is the propensity feature stable across A/B/C (from their own controls),
and do the swapped genes differ materially between contexts? (No A/B/C truth: descriptive only.)

  * excess dropout x and the key lc + beta x across H1 / A / B / C on genes >= 5 CPM in both (517's r's reproduced);
  * per context, the 2,200-gene CPM block (among non-panel genes of the 5-CPM universe; the per-target |H| head is not
    excluded here, so swap counts are approximate) vs the propensity block at beta 22 (lcap) and 11.54 (ltv), with raw x and
    with x quantile-mapped onto H1's distribution (qm: the same x scale as the fit's); the swapped-in sets' size, their
    overlap across contexts and with H1's, and their expression and dropout.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/abc528.py
"""
from __future__ import annotations

import csv
import json
import sys
from itertools import combinations
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
N_CPM = 2200


def jac(a, b):
    return len(a & b) / max(len(a | b), 1)


def qmap(x, ref):
    """map x onto the distribution of ref by rank (ties by order)."""
    r = np.argsort(np.argsort(x, kind="stable"), kind="stable") / max(x.size - 1, 1)
    return np.quantile(ref, r)


def main():
    official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
    F = {c: dict(np.load(HERE / f"work/feat528_{c}.npz", allow_pickle=True)) for c in ("H1", "A", "B", "C")}
    for c in "ABC":
        assert list(F[c]["genes"].astype(str)) == official, f"{c} axis is not the official order"
    # everything on the official axis (H1 mapped by name)
    pos = {g: i for i, g in enumerate(F["H1"]["genes"].astype(str))}
    hcol = np.array([pos.get(g, -1) for g in official])
    V = {}
    for c, d in F.items():
        if c == "H1":
            V[c] = {k: np.where(hcol >= 0, d[k][np.maximum(hcol, 0)], np.nan) for k in ("x", "mean_cpm")}
            V[c]["keep"] = (hcol >= 0) & (np.nan_to_num(V[c]["mean_cpm"]) > 5)
        else:
            V[c] = {k: d[k] for k in ("x", "mean_cpm", "keep")}
    out = {"x_r": {}, "key_r": {}}
    ctx = ("H1", "A", "B", "C")
    for a, b in combinations(ctx, 2):
        m = V[a]["keep"] & V[b]["keep"]
        out["x_r"][f"{a}-{b}"] = float(np.corrcoef(V[a]["x"][m], V[b]["x"][m])[0, 1])
        ka_ = np.log10(V[a]["mean_cpm"][m]) + 22 * V[a]["x"][m]
        kb_ = np.log10(V[b]["mean_cpm"][m]) + 22 * V[b]["x"][m]
        out["key_r"][f"{a}-{b}"] = float(np.corrcoef(ka_, kb_)[0, 1])
    print("x r:", json.dumps({k: round(v, 3) for k, v in out["x_r"].items()}))
    print("key r (lc + 22 x):", json.dumps({k: round(v, 3) for k, v in out["key_r"].items()}))
    notpanel = np.array([g not in panel for g in official])
    h1x = V["H1"]["x"][V["H1"]["keep"]]
    blocks = {}
    for c in ctx:
        e = V[c]["keep"] & notpanel
        lc = np.log10(np.maximum(np.nan_to_num(V[c]["mean_cpm"]), 1e-3))
        x = np.nan_to_num(V[c]["x"])
        xq = x.copy()
        xq[e] = qmap(x[e], h1x)
        res = {}

        def top(score):
            s = np.where(e, score, -np.inf)
            return set(np.argsort(-s, kind="stable")[:N_CPM].tolist())
        b0 = top(lc)
        blocks[(c, "cpm")] = b0
        for nm, beta, xx in (("lcap", 22.0, x), ("lcap_qm", 22.0, xq), ("ltv", 11.54, x), ("ltv_qm", 11.54, xq)):
            b = top(lc + beta * xx)
            sw = b - b0
            blocks[(c, nm)] = sw
            idx = np.array(sorted(sw))
            res[nm] = dict(swapped=len(sw), swapped_median_cpm=float(np.median(V[c]["mean_cpm"][idx])) if idx.size else None,
                           swapped_median_x=float(np.median(x[idx])) if idx.size else None)
        res["block_edge_cpm"] = float(np.sort(V[c]["mean_cpm"][e])[::-1][N_CPM - 1])
        res["x_median_universe"] = float(np.median(x[V[c]["keep"]]))
        out[c] = res
        print(c, json.dumps(res))
    # overlap of swapped-in sets across contexts
    for nm in ("lcap", "lcap_qm", "ltv", "ltv_qm"):
        out[f"swapped_jaccard_{nm}"] = {f"{a}-{b}": jac(blocks[(a, nm)], blocks[(b, nm)]) for a, b in combinations(ctx, 2)}
        print(f"swapped-in Jaccard {nm}:", json.dumps({k: round(v, 3) for k, v in out[f'swapped_jaccard_{nm}'].items()}))
    out["cpm_block_jaccard"] = {f"{a}-{b}": jac(blocks[(a, "cpm")], blocks[(b, "cpm")]) for a, b in combinations(ctx, 2)}
    print("CPM block Jaccard:", json.dumps({k: round(v, 3) for k, v in out["cpm_block_jaccard"].items()}))
    (HERE / "work/abc528.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
