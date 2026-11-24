#!/usr/bin/env python3
"""622n: identity's TAILS on uncalled gate genes. Is the signal concentrated at the extremes of I, where a delivery
route (630's UP pass; a top-I DOWN pass) could carry it with real magnitude?

H1 (625d's `none` s0, the real production emitter: called = q < 0.05) gave on uncalled gate genes:
  bottom-10% I -> UP   0.61 (f 1) / 0.71 (f 0.35);
  top-10% I    -> DOWN 0.60 / 0.58;
  today's delivered sign accuracy on both tails 0.52-0.55;
  mean |l_real| 0.5-0.76.
Here: the Broad-16 (12 working lines; the 4 reserved lines are never parsed). The "uncalled" region is 622m's
approximate production mask: genes outside the gated-consensus head (400) and the top-2,200-expression block. Gate =
the authors' fdr < 0.05. Per line: the truth-UP fraction on the bottom 5 % / 10 % of I and the truth-DOWN fraction on the
top 5 % / 10 %, against CLEAN4's sign accuracy on the same pairs and the line's base rate of DOWN on uncalled gate genes.
Analysis only, < 4 GB.
usage: uv run python experiments/622-identity-track/tails622.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import broad622 as B  # noqa: E402

HEAD, BLOCK = 400, 2200


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    all_t = set(k for _, k in cl)
    ctrl_cols = {ln: [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t] for ln in lines}
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    cpm = {}
    for ln in lines:
        s = np.nansum(PB[:, ctrl_cols[ln]], axis=1)
        absent = np.all(~np.isfinite(PB[:, ctrl_cols[ln]]), axis=1)
        v = 1e6 * s / np.nansum(s); v[absent] = np.nan
        cpm[ln] = np.where(rowmap >= 0, v[np.maximum(rowmap, 0)], np.nan)
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc))
    C4 = np.array([cmap.get(g, np.nan) for g in gl])
    colline = np.array([ln for ln, _ in cl]); coltgt = np.array([k for _, k in cl])
    out = {}
    for ln in lines:
        others = [o for o in lines if o != ln]
        expressed = cpm[ln] > 5
        I = np.where(expressed, np.log1p(cpm[ln]) - np.nanmean([np.log1p(cpm[o]) for o in others], 0), np.nan)
        jl = np.flatnonzero(colline == ln)
        rank_expr = np.argsort(-np.nan_to_num(cpm[ln], nan=-1.0), kind="stable")
        unc = np.ones((len(gl), len(jl)), bool)
        for c, j in enumerate(jl):
            t = coltgt[j]
            oj = np.flatnonzero((coltgt == t) & (colline != ln))
            notg = gl != t
            head = np.zeros(len(gl), bool)
            if len(oj) >= 2:
                cons = np.abs(np.mean(np.where((FDR[:, oj] < 0.05) & np.isfinite(LFC[:, oj]), LFC[:, oj], 0.0), axis=1))
                cons = np.where(notg & np.isfinite(cons), cons, -np.inf)
                head[np.argsort(-cons, kind="stable")[:HEAD]] = True
            blk = [g for g in rank_expr if not head[g] and notg[g] and np.isfinite(cpm[ln][g])][:BLOCK]
            unc[head, c] = False
            unc[blk, c] = False
        gate = (FDR[:, jl] < 0.05) & np.isfinite(LFC[:, jl]) & (LFC[:, jl] != 0) & (gl[:, None] != coltgt[jl][None, :])
        m0 = gate & unc & np.isfinite(I)[:, None]
        dn = LFC[:, jl] < 0
        res = {"uncalled_gate_pairs": int(m0.sum()), "base_down": float((dn & m0).sum() / max(m0.sum(), 1))}
        for q in (0.05, 0.10):
            lo, hi = np.nanquantile(I, q), np.nanquantile(I, 1 - q)
            mlo, mhi = m0 & (I <= lo)[:, None], m0 & (I >= hi)[:, None]
            c4lo = np.isfinite(C4)[:, None] & mlo; c4hi = np.isfinite(C4)[:, None] & mhi
            res[f"bottom{int(q * 100)}_up"] = float(((~dn) & mlo).sum() / max(mlo.sum(), 1))
            res[f"bottom{int(q * 100)}_n"] = int(mlo.sum())
            res[f"bottom{int(q * 100)}_clean4_acc"] = float((((C4 > 0.5)[:, None] == dn) & c4lo).sum() / max(c4lo.sum(), 1))
            res[f"top{int(q * 100)}_down"] = float((dn & mhi).sum() / max(mhi.sum(), 1))
            res[f"top{int(q * 100)}_n"] = int(mhi.sum())
            res[f"top{int(q * 100)}_clean4_acc"] = float((((C4 > 0.5)[:, None] == dn) & c4hi).sum() / max(c4hi.sum(), 1))
        out[ln] = res
        print(f"{ln:9s} uncalled gate {res['uncalled_gate_pairs']:5d} base DOWN {res['base_down']:.3f} | bottom10% I: UP {res['bottom10_up']:.3f} "
              f"(n {res['bottom10_n']}, CLEAN4 {res['bottom10_clean4_acc']:.3f}) | top10% I: DOWN {res['top10_down']:.3f} "
              f"(n {res['top10_n']}, CLEAN4 {res['top10_clean4_acc']:.3f}) | bottom5 UP {res['bottom5_up']:.3f} top5 DOWN {res['top5_down']:.3f}", flush=True)
    keys = ("base_down", "bottom10_up", "bottom10_clean4_acc", "top10_down", "top10_clean4_acc", "bottom5_up", "top5_down")
    summ = {k: float(np.mean([out[ln][k] for ln in lines])) for k in keys}
    summ["lines_bottom10_up_gt_0.5"] = int(sum(out[ln]["bottom10_up"] > 0.5 for ln in lines))
    summ["lines_top10_down_gt_base"] = int(sum(out[ln]["top10_down"] > out[ln]["base_down"] for ln in lines))
    out["_summary"] = summ
    print("\nmean over 12 lines:", {k: round(v, 3) for k, v in summ.items()})
    (HERE / "work/tails622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
