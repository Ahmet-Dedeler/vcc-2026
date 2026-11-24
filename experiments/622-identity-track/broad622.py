#!/usr/bin/env python3
"""622g — independent confirmation of the identity principle on the Broad/DepMap 16-line Perturb-seq (268's source).

One lab, many cancer lines, 90 Cas9-KO targets each: the same structure as Arc's A-F, but a different lab, modality and
chemistry from every surface 636/639/622 used, and never tuned on by the identity work. ONLY the 12 working lines are
read; the 4 reserved lines (JVE127, MG63, SLR23, UMRC3; 268's reserve.json) are excluded by column BEFORE any value is
parsed (usecols), so their numbers never enter memory.
Per held-out working line L:
  gate   = authors' fdr < 0.05, finite lfc != 0, gene != the KO target (authors' lfc_matrix / fdr_matrix)
  I_L    = log1p(CPM of L's OR-control pseudobulk) - mean over the OTHER 11 working lines' log1p(CPM)  (same-lab refs)
Sources scored on gate genes (sign accuracy, paired vs CLEAN4 on genes both cover, 2,000-draw target bootstrap):
  identity alone (z(I) > 0 -> down), all-down, CLEAN4 (CRISPRi public banks, all four: no bank shares these lines),
  z(C4) + lam z(I) for lam 0.5 / 1, and (reference row, not available for D/E/F) the same-lab LOLO programme =
  fraction of the other 11 lines' gated responses that are down, per gene.
Analysis only, < 4 GB.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
P = ROOT / "data/external/depmap26q1/extracted/processed"
RESERVED = {"JVE127", "MG63", "SLR23", "UMRC3"}
RNG = np.random.default_rng(6227)


def header(f: Path, nrows: int):
    with open(f) as fh:
        r = csv.reader(fh)
        return [next(r) for _ in range(nrows)]


def load_matrix(f: Path, skip_label_rows: int):
    h = header(f, 2)
    lines, kos = h[0], h[1]
    keep = [0] + [i for i in range(1, len(lines)) if lines[i] not in RESERVED]
    assert not any(lines[i] in RESERVED for i in keep[1:])
    df = pd.read_csv(f, header=None, skiprows=2 + skip_label_rows, usecols=keep, low_memory=False)
    genes = df.iloc[:, 0].astype(str).to_numpy()
    X = df.iloc[:, 1:].to_numpy(dtype=np.float64)
    cols = [(lines[i], kos[i]) for i in keep[1:]]
    return genes, X, cols


def main() -> None:
    gl, LFC, cl = load_matrix(P / "lfc_matrix.csv", 1)
    gf, FDR, cf = load_matrix(P / "fdr_matrix.csv", 1)
    assert list(gl) == list(gf) and cl == cf
    gp, PB, cp = load_matrix(P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert not (set(lines) & RESERVED) and len(lines) == 12, lines
    targets_by_line = {ln: [k for l2, k in cl if l2 == ln] for ln in lines}
    # OR controls: KOs in the pseudobulk that are not analysis targets
    all_t = set(k for _, k in cl)
    ctrl_cols = {ln: [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t] for ln in lines}
    print("controls per line:", {ln: [cp[j][1] for j in v] for ln, v in list(ctrl_cols.items())[:2]}, flush=True)
    # control CPM per line on the pseudobulk gene axis (NaN-safe: genes absent in a line are NaN in the authors' file)
    cpm = {}
    for ln in lines:
        s = np.nansum(PB[:, ctrl_cols[ln]], axis=1)
        absent = np.all(~np.isfinite(PB[:, ctrl_cols[ln]]), axis=1)
        v = 1e6 * s / np.nansum(s)
        v[absent] = np.nan
        cpm[ln] = v
    gpos = {g: i for i, g in enumerate(gp)}
    # CLEAN4 (all four CRISPRi banks; no bank shares a Broad-16 line)
    c4 = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    c4map = dict(zip(c4["genes"].astype(str), Pc))
    C4 = np.array([c4map.get(g, np.nan) for g in gl])
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    # identity per line on the lfc gene axis
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    out = {}
    for ln in lines:
        others = [o for o in lines if o != ln]
        Ls = np.log1p(cpm[ln]); Lr = np.nanmean([np.log1p(cpm[o]) for o in others], 0)
        Ip = Ls - Lr
        I = np.where(rowmap >= 0, Ip[np.maximum(rowmap, 0)], np.nan)
        expressed = np.where(rowmap >= 0, cpm[ln][np.maximum(rowmap, 0)], np.nan) > 5
        I = np.where(expressed, I, np.nan)
        zi = (I - np.nanmean(I)) / np.nanstd(I)
        # same-lab LOLO programme (reference row only): fraction down among the other lines' gated responses per gene
        jj = [j for j, (l2, _) in enumerate(cl) if l2 in others]
        gated_o = (FDR[:, jj] < 0.05) & np.isfinite(LFC[:, jj]) & (LFC[:, jj] != 0)
        dn_o = (LFC[:, jj] < 0) & gated_o
        prog = (dn_o.sum(1) + 0.5) / (gated_o.sum(1) + 1.0)
        zp = np.where(gated_o.sum(1) >= 3, (prog - 0.5) / np.nanstd(prog - 0.5), np.nan)
        jl = [j for j, (l2, _) in enumerate(cl) if l2 == ln]
        tl = [cl[j][1] for j in jl]
        gate = (FDR[:, jl] < 0.05) & np.isfinite(LFC[:, jl]) & (LFC[:, jl] != 0)
        gate &= (gl[:, None] != np.array(tl)[None, :])
        truth_dn = LFC[:, jl] < 0
        srcs = {"identity": zi, "all_down": np.ones_like(zi), "clean4": zc, "c4+0.5I": zc + 0.5 * np.nan_to_num(zi, nan=0.0) * np.isfinite(zi) + np.where(np.isfinite(zi), 0, np.nan),
                "c4+1I": zc + 1.0 * np.nan_to_num(zi, nan=0.0) + np.where(np.isfinite(zi), 0, np.nan),
                "samelab_programme(ref)": zp}
        res = {"targets": len(tl), "gate_pairs": int(gate.sum())}
        for name, sc in srcs.items():
            both = np.isfinite(sc) & np.isfinite(zc) & np.isfinite(zi)
            m = gate & both[:, None]
            okS = (((sc > 0)[:, None] == truth_dn) & m).sum(0).astype(float)
            okC = (((zc > 0)[:, None] == truth_dn) & m).sum(0).astype(float)
            n = m.sum(0).astype(float)
            idx = np.flatnonzero(n > 0)
            bs = [(okS[b].sum() - okC[b].sum()) / n[b].sum() for b in (RNG.choice(idx, idx.size) for _ in range(2000))]
            res[name] = dict(acc=float(okS.sum() / max(n.sum(), 1)), d_vs_c4=float((okS.sum() - okC.sum()) / max(n.sum(), 1)),
                             ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], n=int(n.sum()))
        out[ln] = res
        print(f"{ln:9s} targets {len(tl):3d} gate {res['gate_pairs']:6d} | " + "  ".join(
            f"{k}: {res[k]['acc']:.3f} ({res[k]['d_vs_c4']:+.3f})" for k in srcs), flush=True)
    # pooled over lines: mean of per-line deltas
    summ = {k: float(np.mean([out[ln][k]["d_vs_c4"] for ln in lines])) for k in ("identity", "all_down", "c4+0.5I", "c4+1I", "samelab_programme(ref)")}
    wins = {k: int(sum(out[ln][k]["ci"][0] > 0 for ln in lines)) for k in summ}
    losses = {k: int(sum(out[ln][k]["ci"][1] < 0 for ln in lines)) for k in summ}
    out["_summary"] = {"mean_delta_vs_clean4": summ, "lines_CI_above_0": wins, "lines_CI_below_0": losses}
    print("\nmean Δ vs CLEAN4 over 12 lines:", {k: round(v, 4) for k, v in summ.items()})
    print("lines with CI > 0:", wins, " CI < 0:", losses)
    (HERE / "work/broad622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
