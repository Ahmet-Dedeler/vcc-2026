#!/usr/bin/env python3
"""622f — does |I| (identity strength) predict DE MEMBERSHIP beyond control expression? Pooled AUROC of gate membership
over tested pairs on 636's four surfaces for: control log-CPM, |I|, and a rank-mix; plus |I|'s partial association given
expression (AUROC within control-log-CPM deciles). Flex references (A/B/C [+h1]). Analysis only, < 2 GB."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst  # noqa: E402


def auroc(score, y):
    r = rankdata(score); n1 = int(y.sum()); n0 = y.size - n1
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)) if n1 and n0 else float("nan")


def main() -> None:
    I = Inst()
    prof = {}
    for s, f in SURF.items():
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    z = np.load(HERE / "work/abc_cpm622.npz", allow_pickle=True)
    for k in ("ctxA", "ctxB", "ctxC"):
        prof[k] = dict(zip(z["genes"].astype(str), z[k].astype(float)))
    out = {}
    for s in I.S:
        D = I.S[s]
        refs = [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != s]
        common = sorted(set(prof[s]).intersection(*(set(prof[r]) for r in refs)))
        Iv = I.gene_vec(s, common, np.log1p([prof[s][g] for g in common]) - np.mean([np.log1p([prof[r][g] for g in common]) for r in refs], 0))
        lc = I.gene_vec(s, list(prof[s]), np.log10(np.array(list(prof[s].values())) + 1e-3))
        m = D["tested"] & np.isfinite(Iv)[None, :] & np.isfinite(lc)[None, :]
        y = D["gate"][m]
        LC = np.broadcast_to(lc, m.shape)[m]; AI = np.broadcast_to(np.abs(Iv), m.shape)[m]
        res = {"auroc_lcpm": auroc(LC, y), "auroc_absI": auroc(AI, y), "auroc_mix": auroc(rankdata(LC) + rankdata(AI), y)}
        dec = np.digitize(LC, np.quantile(LC, np.linspace(0.1, 0.9, 9)))
        within = [auroc(AI[dec == k], y[dec == k]) for k in range(10)]
        res["auroc_absI_within_lcpm_deciles"] = [round(x, 4) for x in within]
        res["mean_within"] = float(np.nanmean(within))
        out[s] = res
        print(f"{s:8s} AUROC lcpm {res['auroc_lcpm']:.4f}  |I| {res['auroc_absI']:.4f}  mix {res['auroc_mix']:.4f}  | |I| within expression deciles mean {res['mean_within']:.4f}", flush=True)
    (HERE / "work/member622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
