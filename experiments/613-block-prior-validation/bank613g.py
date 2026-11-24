#!/usr/bin/env python3
"""613g: the ENSEMBLE hypothesis for the block prior. Committed before its first run.

613f: on 12 independent Broad lines, every single bank, and every 2- or 3-bank subset tried, loses to the equal-weight
4-bank CLEAN4. 458d / 613: a single line's programme does not transfer. Hypothesis: averaging MORE lines' programmes
generalises better. For D/E/F all 7 public CRISPRi banks in 371's consensus file are legal (none is a D/E/F line):
CLEAN4's four, plus h1_2025, orion_hct116, orion_hek293t.

Compositions, each using only banks NOT forbidden on the surface (inst636.BANK_FORBIDDEN; on Broad lines nothing is
forbidden):
  ENS_ALL   = the mean of every allowed bank among the 7;
  ENS_NOH1  = the same without h1_2025 (458d's harmful single-line programme);
against CLEAN4 = the mean of the allowed among its 4 (clean4_lolo's rule).
Surfaces:
  - the 12 working Broad lines (613f's definitions; reserved lines never read);
  - 636's h1 / HCT116 / HEK293T / K562 (613c's definitions: block = decomp639's head-excluded call set; null = q ≥ 0.05).
Readings: the pooled paired precision Δ vs CLEAN4 on the NULL block (primary) and the gate block.
LEAD (fixed now): beats CLEAN4 on the null block on >= 9/12 Broad lines, with the mean gate-block Δ over Broad >= 0, AND
is not below CLEAN4 by > 0.005 on the null block of any 636 surface.
Try count: 613f's 6 compositions + these 2. Output: work/bank613g.out.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as BR  # noqa: E402
import block613 as B  # noqa: E402  (puts 636 / 639 on the path)
import decomp639 as DC  # noqa: E402
from inst636 import BANK_FORBIDDEN, G480, SURF, Inst  # noqa: E402

ALL7 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2", "h1_2025", "orion_hct116", "orion_hek293t")
C4B = ALL7[:4]


def bank_P():
    c = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(c["bank_names"].astype(str))
    return c["genes"].astype(str), {b: c["P_banks"].astype(float)[names.index(b)] for b in ALL7}


def comps(allowed):
    a4 = [b for b in C4B if b in allowed]
    a7 = [b for b in ALL7 if b in allowed]
    return {"CLEAN4": a4, "ENS_ALL": a7, "ENS_NOH1": [b for b in a7 if b != "h1_2025"]}


def mean_of(Pax, banks):
    return np.nanmean(np.stack([Pax[b] for b in banks]), 0)


def paired_delta(S, C4, truth_dn, M):
    both = np.isfinite(S) & np.isfinite(C4) & (S != 0.5) & (C4 != 0.5)
    m = M & both
    ok_s = (((S > 0.5) == truth_dn) & m).sum()
    ok_c = (((C4 > 0.5) == truth_dn) & m).sum()
    return float((ok_s - ok_c) / max(m.sum(), 1))


def main():
    cg, PB7 = bank_P()
    out = []
    # ---- Broad 12 (613f's construction)
    gl, LFC, cl = BR.load_matrix(BR.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = BR.load_matrix(BR.P / "fdr_matrix.csv", 1)
    assert list(gl) == list(gf) and cl == cf
    gp, PBk, cp = BR.load_matrix(BR.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert not (set(lines) & BR.RESERVED) and len(lines) == 12
    all_t = set(k for _, k in cl)
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    cmap = {g: i for i, g in enumerate(cg)}
    col = np.array([cmap.get(g, -1) for g in gl])
    Pax = {}
    for b in ALL7:
        v = np.full(len(gl), np.nan)
        v[col >= 0] = PB7[b][col[col >= 0]]
        Pax[b] = v
    cmp_b = comps(set(ALL7))
    S = {k: mean_of(Pax, v) for k, v in cmp_b.items()}
    bres = {k: {"null": [], "gate": []} for k in ("ENS_ALL", "ENS_NOH1")}
    for ln in lines:
        jl = [j for j, (l2, _) in enumerate(cl) if l2 == ln]
        tl = [cl[j][1] for j in jl]
        L, F = LFC[:, jl], FDR[:, jl]
        own = gl[:, None] == np.array(tl)[None, :]
        fin = np.isfinite(L) & (L != 0) & ~own
        gate, null = fin & (F < 0.05), fin & np.isfinite(F) & (F >= 0.05)
        ctrl = [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t]
        s = np.nansum(PBk[:, ctrl], axis=1)
        cpm = np.where(rowmap >= 0, (1e6 * s / np.nansum(s))[np.maximum(rowmap, 0)], np.nan)
        head = np.zeros(L.shape, bool)
        for a, t in enumerate(tl):
            jo = [j for j, (l2, k) in enumerate(cl) if k == t and l2 != ln]
            if jo:
                head[np.argsort(-np.nan_to_num(np.nanmean(np.abs(LFC[:, jo]), 1), nan=-1.0), kind="stable")[:400], a] = True
        block = np.zeros(L.shape, bool)
        rc = np.nan_to_num(cpm, nan=-1.0)
        for a in range(L.shape[1]):
            block[np.argsort(-np.where(head[:, a], -np.inf, rc), kind="stable")[:2200], a] = True
        block &= ~head
        dn = L < 0
        for k in bres:
            bres[k]["null"].append(paired_delta(S[k][:, None], S["CLEAN4"][:, None], dn, null & block))
            bres[k]["gate"].append(paired_delta(S[k][:, None], S["CLEAN4"][:, None], dn, gate & block))
    for k, r in bres.items():
        out.append(f"Broad-12 {k:9s} null: wins {sum(v > 0 for v in r['null']):2d}/12 mean {np.mean(r['null']):+.4f} "
                   f"(min {min(r['null']):+.4f}) | gate: wins {sum(v > 0 for v in r['gate']):2d}/12 mean {np.mean(r['gate']):+.4f}")
    # ---- 636 surfaces (613c's construction)
    I = Inst()
    sres = {}
    for s in ("h1",) + B.CANCER:
        D = I.S[s]
        lfc, gate, tested = D["lfc"], D["gate"], D["tested"]
        q = np.load(G480 / SURF[s], allow_pickle=True)["q"].astype(float)
        null = tested & np.isfinite(q) & (q >= 0.05) & (lfc != 0)
        _, block = DC.callset_for(I, s)
        allowed = {b for b in ALL7 if s not in BANK_FORBIDDEN.get(b, set())}
        cm = comps(allowed)
        Px = {b: I.gene_vec(s, cg, PB7[b]) for b in allowed}
        Ss = {k: mean_of(Px, v) for k, v in cm.items()}
        dn = lfc < 0
        sres[s] = {}
        for k in ("ENS_ALL", "ENS_NOH1"):
            sres[s][k] = (paired_delta(Ss[k][None, :], Ss["CLEAN4"][None, :], dn, null & block),
                          paired_delta(Ss[k][None, :], Ss["CLEAN4"][None, :], dn, gate & block), cm[k])
        out.append(f"{s:8s} " + " | ".join(f"{k}: null {v[0]:+.4f} gate {v[1]:+.4f} ({len(v[2])} banks)" for k, v in sres[s].items()))
    for k in ("ENS_ALL", "ENS_NOH1"):
        lead = (sum(v > 0 for v in bres[k]["null"]) >= 9 and np.mean(bres[k]["gate"]) >= 0
                and all(sres[s][k][0] >= -0.005 for s in sres))
        out.append(f"LEAD {k}: {'YES' if lead else 'no'}")
    print("\n".join(out))
    (HERE / "work/bank613g.out").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
