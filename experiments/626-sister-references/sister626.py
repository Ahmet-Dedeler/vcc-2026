#!/usr/bin/env python3
"""626 — is the identity rider robust to SISTER references (a reference of the context's lineage)? ([4d0207]'s item)

On Oct 22 the guard drops only an exact same-call reference. A D/E/F line that is a sister of a reference (two T-ALL lines,
two cervical / squamous lines) shrinks the identity contrast for lineage genes. Broad-16 (12 working lines; the reserved 4
never read) has natural sisters: KYSE140/KYSE450 (esophagus), C4I/SKGII (cervix). Per held-out line L:
  - similarity of L to each other line: Spearman of their log1p control-CPM DEVIATIONS from the 12-line mean (truth-free,
    computable on the day among Flex contexts);
  - identity rider gain (gate sign of z(C4 - .5) + 0.5 z(I) vs CLEAN4, authors' fdr < .05 gates) with reference pools:
      all11; no-sister (same-lineage lines dropped); 3 refs WITH the sister (+2 random others) vs 3 refs WITHOUT (20 draws);
      guard "drop refs with similarity > t" for t in a grid;
  - the identity-only accuracy too (the sign of I alone), where a sister's damage shows most directly.
Then the guard's cost on 622's / 636's surfaces: pairwise similarity among {A, B, C, h1} and the 3' surfaces, and which
references a threshold would drop. Analysis only, < 4 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(6260)
LINEAGE = {"C4I": "cervix", "HCT15": "bowel", "HKGZCC": "biliary", "KMRC20": "kidney", "KYSE140": "esophagus",
           "KYSE450": "esophagus", "ONE58": "pleura", "RPMI7951": "skin", "SKGII": "cervix", "SNU8": "ovary",
           "T3M5": "head_neck", "U343": "cns"}
THRESH = (0.3, 0.4, 0.5, 0.6, 0.7)


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl}); assert set(lines) == set(LINEAGE) and not (set(lines) & B.RESERVED)
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc)); C4 = np.array([cmap.get(g, np.nan) for g in gl])
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    all_t = set(k for _, k in cl)
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    Lg = {}
    for ln in lines:
        cc = [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t]
        s = np.nansum(PB[:, cc], 1); v = 1e6 * s / np.nansum(s)
        v[np.all(~np.isfinite(PB[:, cc]), 1)] = np.nan
        Lg[ln] = np.where(rowmap >= 0, np.log1p(v[np.maximum(rowmap, 0)]), np.nan)
    allok = np.all([np.isfinite(Lg[l]) for l in lines], 0)
    mean12 = np.mean([Lg[l] for l in lines], 0)
    dev = {l: Lg[l] - mean12 for l in lines}
    sim = {(a, b): float(spearmanr(dev[a][allok], dev[b][allok]).correlation) for a in lines for b in lines if a != b}
    gated = (FDR < 0.05) & np.isfinite(LFC) & (LFC != 0)
    colline = np.array([ln for ln, _ in cl]); tgt = np.array([k for _, k in cl])

    def gain(L, refs):
        I = Lg[L] - np.nanmean([Lg[r] for r in refs], 0)
        expressed = np.exp(Lg[L]) - 1 > 5
        I = np.where(expressed, I, np.nan)
        zi = (I - np.nanmean(I)) / np.nanstd(I)
        jl = np.flatnonzero(colline == L)
        gate = gated[:, jl] & (gl[:, None] != tgt[jl][None, :]); tdn = LFC[:, jl] < 0
        both = np.isfinite(zc) & np.isfinite(zi)
        m = gate & both[:, None]
        S = zc + 0.5 * zi
        okS = (((S > 0)[:, None] == tdn) & m).sum(); okC = (((zc > 0)[:, None] == tdn) & m).sum()
        okI = (((zi > 0)[:, None] == tdn) & m).sum()
        n = max(m.sum(), 1)
        return dict(d_rider=float((okS - okC) / n), acc_identity=float(okI / n), n=int(m.sum()))

    out = {"similarity": {f"{a}~{b}": v for (a, b), v in sim.items() if a < b}, "lines": {}}
    for L in lines:
        others = [o for o in lines if o != L]
        sisters = [o for o in others if LINEAGE[o] == LINEAGE[L]]
        nonsis = [o for o in others if o not in sisters]
        nearest = max(others, key=lambda o: sim[(L, o)])
        rec = {"lineage": LINEAGE[L], "sisters": sisters, "nearest": nearest, "nearest_sim": sim[(L, nearest)],
               "all11": gain(L, others), "no_sister": gain(L, nonsis) if sisters else None}
        if sisters:
            w, wo = [], []
            for _ in range(20):
                pick = list(RNG.choice(nonsis, 2, replace=False))
                w.append(gain(L, sisters[:1] + pick)); wo.append(gain(L, list(RNG.choice(nonsis, 3, replace=False))))
            rec["three_with_sister"] = {k: float(np.mean([x[k] for x in w])) for k in ("d_rider", "acc_identity")}
            rec["three_without_sister"] = {k: float(np.mean([x[k] for x in wo])) for k in ("d_rider", "acc_identity")}
        rec["guard"] = {}
        for t in THRESH:
            keep = [o for o in others if sim[(L, o)] <= t]
            rec["guard"][str(t)] = dict(n_refs=len(keep), dropped=[o for o in others if o not in keep],
                                        **(gain(L, keep) if len(keep) >= 2 else {}))
        out["lines"][L] = rec
        s = f"{L:9s} [{LINEAGE[L]:9s}] nearest {nearest:9s} sim {sim[(L, nearest)]:+.3f} | all11 rider {rec['all11']['d_rider']:+.4f} I-acc {rec['all11']['acc_identity']:.3f}"
        if sisters:
            s += (f" | no-sister {rec['no_sister']['d_rider']:+.4f} | 3 refs with sister {rec['three_with_sister']['d_rider']:+.4f} "
                  f"(I {rec['three_with_sister']['acc_identity']:.3f}) vs without {rec['three_without_sister']['d_rider']:+.4f} "
                  f"(I {rec['three_without_sister']['acc_identity']:.3f})")
        print(s, flush=True)
    # guard summary: mean rider gain over lines at each threshold vs all11, and how many refs dropped
    print()
    for t in THRESH:
        d = [out["lines"][L]["guard"][str(t)].get("d_rider", np.nan) - out["lines"][L]["all11"]["d_rider"] for L in lines]
        nd = [len(out["lines"][L]["guard"][str(t)]["dropped"]) for L in lines]
        print(f"guard sim > {t}: mean Δ(rider gain) vs all11 {np.nanmean(d):+.4f}; refs dropped per line mean {np.mean(nd):.2f}, lines with any drop {sum(x > 0 for x in nd)}/12")
    sis_pairs = [("KYSE140", "KYSE450"), ("C4I", "SKGII")]
    nonsis_sims = [v for (a, b), v in sim.items() if LINEAGE[a] != LINEAGE[b] and a < b]
    print("sister-pair similarities:", {f"{a}~{b}": round(sim[(a, b)], 3) for a, b in sis_pairs},
          " non-sister: median", round(float(np.median(nonsis_sims)), 3), "max", round(float(np.max(nonsis_sims)), 3))
    out["sister_pairs"] = {f"{a}~{b}": sim[(a, b)] for a, b in sis_pairs}
    out["nonsister_sim_median_max"] = [float(np.median(nonsis_sims)), float(np.max(nonsis_sims))]
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/sister626.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
