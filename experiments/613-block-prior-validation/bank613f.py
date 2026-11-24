#!/usr/bin/env python3
"""613f: CLEAN4's BANK COMPOSITION for D/E/F, read on the 12 working Broad/DepMap lines with 613's block readings.
Committed before its first run.

Why: 613c found Nadig-Jurkat alone beating 4-bank CLEAN4 on both Orion folds' null block (+0.012 / +0.033), and the K562
banks alone far worse there (k562_gwps −0.036 / −0.053): a line-specific programme diluting CLEAN4 off-line, like 458's H1
programme. But H1 disagreed (−0.0065), and the two Orion folds share a platform. The Broad-16 screen (another lab;
Cas9 KO; 3′ chemistry) has 12 working lines with per-target truth. The 4 reserved lines are dropped by column before
parsing (broad622's loader; never read).

Per working line L (broad622's definitions):
  gate  = authors' fdr < 0.05, finite lfc != 0, gene != the KO target;
  null  = finite lfc != 0, fdr >= 0.05 (finite), gene != the KO target;
  head  = each target's top 400 genes by the other 11 lines' mean |lfc| for the same KO (the D/E/F-like key);
  block = the top 2,200 of the rest by L's OR-control CPM.
Sources: CLEAN4 (4 banks, equal weight; the reference); each bank alone; NADIG2 = the mean of Jurkat + HepG2;
NOGWPS = CLEAN4 without k562_gwps. DOWN where P > 1/2. Paired cells only (both the source and CLEAN4 signed).
Primary: the pooled paired precision Δ vs CLEAN4 on the NULL block (fid's filler), per line. Reported beside it: the
gate block.
LEAD (fixed now): a composition is a lead only if it beats CLEAN4 on the null block on >= 9 of 12 lines AND on both
Orion folds (613c), AND is not below CLEAN4 by > 0.005 on H1 (613c). Caveat stated in advance: KO gates lean UP on these
lines (622g: all-down −0.24), so only relative differences between banks are read here.
Output: work/bank613f.out.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as BR  # noqa: E402  (loaders; RESERVED dropped by usecols)

BANKS = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")


def main():
    gl, LFC, cl = BR.load_matrix(BR.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = BR.load_matrix(BR.P / "fdr_matrix.csv", 1)
    assert list(gl) == list(gf) and cl == cf
    gp, PB, cp = BR.load_matrix(BR.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert not (set(lines) & BR.RESERVED) and len(lines) == 12, lines
    all_t = set(k for _, k in cl)
    ctrl_cols = {ln: [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t] for ln in lines}
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    c4 = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(c4["bank_names"].astype(str))
    Pb = {b: c4["P_banks"].astype(float)[names.index(b)] for b in BANKS}
    cg = c4["genes"].astype(str)
    cmap = {g: i for i, g in enumerate(cg)}
    col = np.array([cmap.get(g, -1) for g in gl])
    def on_axis(v):
        out = np.full(len(gl), np.nan)
        out[col >= 0] = v[col[col >= 0]]
        return out
    P = {b: on_axis(Pb[b]) for b in BANKS}
    srcs = {"CLEAN4": np.nanmean(np.stack([P[b] for b in BANKS]), 0)}
    for b in BANKS:
        srcs[f"{b} alone"] = P[b]
    srcs["NADIG2 (jurkat+hepg2)"] = np.nanmean(np.stack([P["nadig_jurkat"], P["nadig_hepg2"]]), 0)
    srcs["NOGWPS (3 banks)"] = np.nanmean(np.stack([P[b] for b in BANKS if b != "k562_gwps"]), 0)
    C4 = srcs["CLEAN4"]
    res = {}
    for ln in lines:
        jl = [j for j, (l2, _) in enumerate(cl) if l2 == ln]
        tl = [cl[j][1] for j in jl]
        L, F = LFC[:, jl], FDR[:, jl]
        own = gl[:, None] == np.array(tl)[None, :]
        fin = np.isfinite(L) & (L != 0) & ~own
        gate = fin & (F < 0.05)
        null = fin & np.isfinite(F) & (F >= 0.05)
        s = np.nansum(PB[:, ctrl_cols[ln]], axis=1)
        cpm_p = 1e6 * s / np.nansum(s)
        cpm = np.where(rowmap >= 0, cpm_p[np.maximum(rowmap, 0)], np.nan)
        # head per target: top 400 by the other lines' mean |lfc| for the same KO (where the KO exists elsewhere)
        head = np.zeros(L.shape, bool)
        for a, t in enumerate(tl):
            jo = [j for j, (l2, k) in enumerate(cl) if k == t and l2 != ln]
            if jo:
                k = np.nanmean(np.abs(LFC[:, jo]), 1)
                head[np.argsort(-np.nan_to_num(k, nan=-1.0), kind="stable")[:400], a] = True
        block = np.zeros(L.shape, bool)
        rank_c = np.nan_to_num(cpm, nan=-1.0)
        for a in range(L.shape[1]):
            block[np.argsort(-np.where(head[:, a], -np.inf, rank_c), kind="stable")[:2200], a] = True
        block &= ~head
        dn = L < 0
        res[ln] = {}
        for name, S in srcs.items():
            if name == "CLEAN4":
                continue
            both = np.isfinite(S) & np.isfinite(C4) & (S != 0.5) & (C4 != 0.5)
            r = {}
            for lab, M in (("null", null & block), ("gate", gate & block)):
                m = M & both[:, None]
                ok_s = (((S > 0.5)[:, None] == dn) & m).sum()
                ok_c = (((C4 > 0.5)[:, None] == dn) & m).sum()
                r[lab] = float((ok_s - ok_c) / max(m.sum(), 1))
            res[ln][name] = r
    lines_out = []
    for name in [k for k in srcs if k != "CLEAN4"]:
        nv = [res[ln][name]["null"] for ln in lines]
        gv = [res[ln][name]["gate"] for ln in lines]
        lines_out.append(f"{name:24s} null: wins {sum(v > 0 for v in nv):2d}/12 mean {np.mean(nv):+.4f} (min {min(nv):+.4f}, max {max(nv):+.4f}) | "
                         f"gate: wins {sum(v > 0 for v in gv):2d}/12 mean {np.mean(gv):+.4f}")
    print("Δ vs 4-bank CLEAN4 on the BLOCK, 12 working Broad lines (+ = the composition is right more often):")
    print("\n".join(lines_out))
    print("\nper line, null block:")
    for ln in lines:
        print(f"  {ln:9s} " + "  ".join(f"{k.split(' ')[0]} {res[ln][k]['null']:+.4f}" for k in res[ln]))
    (HERE / "work/bank613f.out").write_text("\n".join(lines_out) + "\n")


if __name__ == "__main__":
    main()
