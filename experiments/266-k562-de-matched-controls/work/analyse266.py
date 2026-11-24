"""266 — tables from runs266.csv + de_tables266.npz: KA check vs 260, three-scheme table, medians / gate
fractions with seed spread, diagnostics, DE-set overlap on the five best-covered targets."""
import json, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
W = Path("experiments/266-k562-de-matched-controls/work")
R260 = Path("experiments/260-raw-de-diagnostic-feasibility/work")
df = pd.read_csv(W / "runs266.csv"); ref = pd.read_csv(R260 / "candidates_de_power.csv").set_index("target")
T = np.load(W / "de_tables266.npz"); genes = T["genes"].astype(str)
out = []; P = out.append
SCH = ["pooled", "gem", "gemdepth"]; SEEDS = [0, 1, 2]

# ---- 1. KA: pooled seed 0 vs 260 on all 50
ka = df[(df.scheme == "pooled") & (df.seed == 0)].set_index("target").loc[ref.index]
exact = {c: int((ka[c].values == ref[c].values).sum()) for c in ("n_real_tc", "n_real_uc", "gate_tc", "gate_uc", "cells")}
fn = np.nan_to_num(ka.frac_neg_tc.values, nan=-1) ; fr = np.nan_to_num(ref.frac_neg_tc.values, nan=-1)
kd = np.nan_to_num(ka.target_knockdown_lfc.values, nan=-99); kr = np.nan_to_num(ref.target_knockdown_lfc.values, nan=-99)
P("## KA — pooled seed 0 (260's exact control draw) vs `candidates_de_power.csv`\n")
P("| field | targets identical / 50 |\n|---|---:|")
for c, v in exact.items(): P(f"| `{c}` | {v} |")
P(f"| `frac_neg_tc` (|Δ| < 1e-9) | {int((np.abs(fn - fr) < 1e-9).sum())} |")
P(f"| `target_knockdown_lfc` (|Δ| < 1e-4; means in float64 here) | {int((np.abs(kd - kr) < 1e-4).sum())} |")
P(f"| `median_umi` | {int((ka.median_umi_t.values == ref.median_umi.values).sum())} |")
P(f"\nControl median UMI (pooled s0) {ka.median_umi_c.iloc[0]:.0f} (260: 11,406). Genes passing the CPM > 5 gate: "
  f"{sorted(df.genes_gated.unique())} in every run; min control-mean CPM over all 450 runs {df.min_ctrl_mean_cpm.min():.2f}.\n")

# ---- 2. per-target three-scheme table
piv = df.pivot_table(index="target", columns=["scheme", "seed"], values="n_real_tc")
piv = piv.loc[ref.index]
P("## Per-target `n_real` (tie-corrected, target gene excluded) — 50 candidates, three schemes × three seeds\n")
P("| target | cells | pooled s0 (=260) | pooled s1 | pooled s2 | gem s0 | gem s1 | gem s2 | gem+depth s0 | s1 | s2 | gem mean − pooled mean | gem+depth mean − pooled mean |")
P("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for t in ref.index:
    r = piv.loc[t]; pm = r["pooled"].mean(); gm = r["gem"].mean(); dm = r["gemdepth"].mean()
    P(f"| {t} | {int(ref.loc[t].cells)} | {int(r[('pooled',0)])} | {int(r[('pooled',1)])} | {int(r[('pooled',2)])} | "
      f"{int(r[('gem',0)])} | {int(r[('gem',1)])} | {int(r[('gem',2)])} | {int(r[('gemdepth',0)])} | {int(r[('gemdepth',1)])} | {int(r[('gemdepth',2)])} | "
      f"{gm - pm:+.1f} | {dm - pm:+.1f} |")

# ---- 3. summaries per scheme × seed
P("\n## Medians and gate fractions per scheme and seed (50 candidates)\n")
P("| scheme | seed | median `n_real` | IQR | ≥ 10 (`nmae` gate) | ≥ 20 | ≥ 30 | ≥ 89 | zero | median paired Δ vs pooled s0 | targets up / down / tie vs pooled s0 |")
P("|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|")
summ = {}
base = piv[("pooled", 0)]
for s in SCH:
    for sd in SEEDS:
        v = piv[(s, sd)]; d = v - base
        rec = dict(median=float(v.median()), q1=float(v.quantile(.25)), q3=float(v.quantile(.75)), ge10=int((v >= 10).sum()),
                   ge20=int((v >= 20).sum()), ge30=int((v >= 30).sum()), ge89=int((v >= 89).sum()), zero=int((v == 0).sum()),
                   med_paired=float(d.median()), up=int((d > 0).sum()), down=int((d < 0).sum()), tie=int((d == 0).sum()))
        summ[f"{s}_s{sd}"] = rec
        P(f"| {s} | {sd} | {rec['median']:.1f} | [{rec['q1']:.0f}, {rec['q3']:.0f}] | {rec['ge10']} ({rec['ge10']*2} %) | {rec['ge20']} | {rec['ge30']} | {rec['ge89']} | {rec['zero']} | "
          f"{rec['med_paired']:+.1f} | {rec['up']} / {rec['down']} / {rec['tie']} |")
P("\nSeed-averaged per target (mean of the three seeds), then the median over targets and the gate fraction on the seed-mean:\n")
P("| scheme | median of seed-mean `n_real` | IQR | seed-mean ≥ 10 | seed-mean ≥ 20 | median paired Δ vs pooled seed-mean | up / down |\n|---|---:|---|---:|---:|---:|---|")
pm = piv["pooled"].mean(1)
for s in SCH:
    v = piv[s].mean(1); d = v - pm
    P(f"| {s} | {v.median():.1f} | [{v.quantile(.25):.0f}, {v.quantile(.75):.0f}] | {int((v >= 10).sum())} ({(v >= 10).mean()*100:.0f} %) | {int((v >= 20).sum())} | {d.median():+.1f} | {int((d > 0).sum())} / {int((d < 0).sum())} |")
# paired log-ratio (n+1) for a scale-free read
P("\nPaired log2((n+1)/(n_pooled+1)) per target, seed-means: " + "; ".join(
    f"{s}: median {np.log2((piv[s].mean(1)+1)/(pm+1)).median():+.3f}, mean {np.log2((piv[s].mean(1)+1)/(pm+1)).mean():+.3f}" for s in SCH[1:]))
# uncorrected p as sensitivity
pivu = df.pivot_table(index="target", columns=["scheme", "seed"], values="n_real_uc").loc[ref.index]
P("\nUncorrected-variance sensitivity (`n_real_uc`), median over targets per scheme (seeds 0/1/2): " + "; ".join(
    f"{s}: {pivu[(s,0)].median():.0f}/{pivu[(s,1)].median():.0f}/{pivu[(s,2)].median():.0f}" for s in SCH))

# ---- 4. diagnostics
P("\n## Diagnostics per scheme (median over 50 targets × 3 seeds unless stated)\n")
P("| scheme | ctrl pool UMI CV | ctrl UMI mean | depth ratio median [min, max] | outside [0.95, 1.05] / 150 | zero-fraction gap (pert − ctrl) median | gemgroups in ctrl pool | gemgroup TV distance target vs ctrl | n_ctrl |")
P("|---|---:|---:|---|---:|---:|---:|---:|---:|")
diag = {}
for s in SCH:
    d = df[df.scheme == s]
    rec = dict(cv=float(d.ctrl_umi_cv.median()), umi_mean=float(d.ctrl_umi_mean.median()), dr_med=float(d.depth_ratio.median()),
               dr_min=float(d.depth_ratio.min()), dr_max=float(d.depth_ratio.max()), dr_out=int(((d.depth_ratio < .95) | (d.depth_ratio > 1.05)).sum()),
               zero_gap=float(d.zero_gap.median()), zero_gap_mean_abs=float(d.zero_gap.abs().mean()), n_gem_c=float(d.n_gem_c.median()), tv=float(d.gem_tv.median()),
               n_ctrl_min=int(d.n_ctrl.min()), n_ctrl_max=int(d.n_ctrl.max()))
    diag[s] = rec
    P(f"| {s} | {rec['cv']:.3f} | {rec['umi_mean']:.0f} | {rec['dr_med']:.4f} [{rec['dr_min']:.3f}, {rec['dr_max']:.3f}] | {rec['dr_out']} | {rec['zero_gap']:+.4f} | {rec['n_gem_c']:.0f} | {rec['tv']:.3f} | {rec['n_ctrl_min']}–{rec['n_ctrl_max']} |")
d0 = df[df.scheme == "pooled"]; P(f"\nTarget cells' own gemgroups: median {df.n_gem_t.median():.0f} per target (min {df.n_gem_t.min()}, max {df.n_gem_t.max()}); "
  f"pooled-scheme per-target depth ratio IQR [{d0.depth_ratio.quantile(.25):.3f}, {d0.depth_ratio.quantile(.75):.3f}].")
# does the per-target change track any diagnostic?
dd = df[df.scheme == "gemdepth"].groupby("target").agg(n=("n_real_tc", "mean"), dr=("depth_ratio", "mean"), cv=("ctrl_umi_cv", "mean"), zg=("zero_gap", "mean"))
dp = df[df.scheme == "pooled"].groupby("target").agg(n=("n_real_tc", "mean"), dr=("depth_ratio", "mean"), zg=("zero_gap", "mean"))
chg = np.log2((dd.n + 1) / (dp.n + 1))
P(f"Spearman of per-target log2 change (gem+depth vs pooled) with the pooled depth ratio: {spearmanr(chg, dp.dr).statistic:+.3f}; "
  f"with the pooled zero-fraction gap: {spearmanr(chg, dp.zg).statistic:+.3f}; with pooled n_real: {spearmanr(chg, dp.n).statistic:+.3f}.")

# ---- 5. overlap on the five best-covered targets (+ DCP2, the largest set), seed 0
P("\n## DE-set overlap across schemes, seed 0 — five best-covered targets by cells (+ DCP2, the largest set)\n")
five = ref.sort_values("cells", ascending=False).index[:5].tolist() + ["DCP2"]
P("| target | cells | pooled | gem | gem+depth | ∩ all three | J(pooled, gem) | J(pooled, gem+depth) | new in gem+depth (not pooled) | their pooled p-rank median [max] | their pooled q median | lost from pooled | ρ(−log p) pooled vs gem+depth |")
P("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|")
def bh(p):
    m = len(p); o = np.argsort(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)
ov = {}
for t in five:
    sets = {}; qs = {}; ps = {}
    for s in SCH:
        p = T[f"{s}_s0_{t}_p"]; q = bh(p); mask = genes != t
        sets[s] = set(np.flatnonzero((q < 0.05) & mask)); qs[s] = q; ps[s] = p
    a, b, c = sets["pooled"], sets["gem"], sets["gemdepth"]
    inter = a & b & c
    J = lambda x, y: len(x & y) / len(x | y) if (x | y) else 1.0
    new = sorted(c - a); lost = a - c
    rank_a = np.empty(len(genes)); rank_a[np.argsort(ps["pooled"], kind="stable")] = np.arange(1, len(genes) + 1)
    nr = rank_a[new] if new else np.array([np.nan]); nq = qs["pooled"][new] if new else np.array([np.nan])
    rho = spearmanr(-np.log10(ps["pooled"] + 1e-300), -np.log10(ps["gemdepth"] + 1e-300)).statistic
    ov[t] = dict(pooled=len(a), gem=len(b), gemdepth=len(c), inter=len(inter), J_ab=J(a, b), J_ac=J(a, c), new=len(new), new_rank_med=float(np.nanmedian(nr)),
                 new_rank_max=float(np.nanmax(nr)), new_q_med=float(np.nanmedian(nq)), lost=len(lost), rho=float(rho),
                 new_genes=[genes[i] for i in new][:25])
    P(f"| {t} | {int(ref.loc[t].cells)} | {len(a)} | {len(b)} | {len(c)} | {len(inter)} | {J(a,b):.2f} | {J(a,c):.2f} | {len(new)} | "
      f"{np.nanmedian(nr):.0f} [{np.nanmax(nr):.0f}] | {np.nanmedian(nq):.3f} | {len(lost)} | {rho:+.3f} |")
P("\n(A pooled p-rank of k means the gene was the k-th smallest pooled p-value among the 8,248 tested, target gene included; "
  "'same genes at a looser bar' would read as ranks just beyond the pooled set size with pooled q just above 0.05.)\n")
for t in five:
    if ov[t]["new"]: P(f"- {t} new in gem+depth: {', '.join(ov[t]['new_genes'])}{' …' if ov[t]['new'] > 25 else ''}")

json.dump(dict(ka_exact=exact, summary=summ, diagnostics=diag, overlap=ov), open(W / "analysis266.json", "w"), indent=1, default=float)
(W / "tables266.md").write_text("\n".join(out) + "\n")
print("\n".join(out))
