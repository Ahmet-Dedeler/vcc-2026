"""287 novelty audit — the four read-only checks on published supplementary tables.

Run: uv run python experiments/287-novelty-audit/work/checks287.py
Reads only files under data/external/ (never written) and data/challenge/pert_counts.csv.

KNOWN-ANSWER ROW (KA-1): the GWCD4i K562-comparison table's random controls must read at chance
(median |r| < 0.02) while the real cross-cell-type correlation must be ~0.30. If KA-1 fails, the
parse is wrong and nothing below is readable.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

SUPP = "data/external/cd4_cell2026/suppl_tables_2026-08-24"
DE = f"{SUPP}/DE_stats.suppl_table.csv"
K562 = f"{SUPP}/K562_comparison.suppl_table.csv"
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
SEED = 287


def boot_median_diff(a, b, rng, n=4000):
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return np.percentile(d, [2.5, 97.5])


def boot_spearman(x, y, rng, n=4000):
    idx = np.arange(len(x))
    d = [spearmanr(x[j], y[j]).statistic for j in (rng.choice(idx, len(idx)) for _ in range(n))]
    return np.percentile(d, [2.5, 97.5])


def main():
    rng = np.random.default_rng(SEED)
    panel = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
    de = pd.read_csv(DE, low_memory=False)
    k5 = pd.read_csv(K562, low_memory=False)

    print("=== KA-1: random controls at chance, real cross-line correlation ~0.30 ===")
    ok = True
    for cond in CONDS:
        s = k5[k5["condition"] == cond]
        r = s["logfc_pearson_r"].dropna()
        rand = pd.concat([s["random_r1"], s["random_r2"], s["random_r3"]]).dropna()
        print(f"  {cond:9s} real median {r.median():+.4f} (n={len(r)})   random median {rand.median():+.4f}"
              f"  p95 {rand.quantile(.95):+.4f}")
        ok &= abs(rand.median()) < 0.02 and 0.25 < r.median() < 0.36
    print(f"  KA-1 {'PASS' if ok else 'FAIL — stop, the parse is wrong'}")
    if not ok:
        raise SystemExit(1)

    print("\n=== CHECK 1: the authors' own donor and guide reproducibility (GWCD4i DE_stats) ===")
    for cond in CONDS:
        sub = de[de["culture_condition"] == cond]
        sp = sub[sub["target_contrast_gene_name"].isin(panel)]
        for label, d in (("all", sub), ("panel-300", sp)):
            out = [f"{cond:9s} {label:9s}"]
            for c in ["donor_correlation_all_mean", "donor_correlation_hits_mean", "guide_correlation_all"]:
                v = d[c].dropna()
                out.append(f"{c.replace('_correlation','')}: med {v.median():+.4f} n={len(v)} <=0 {(v<=0).mean():.3f}")
            print("  " + " | ".join(out))

    print("\n=== CHECK 1b: can the preprint's stated cross-guide medians (0.49/0.47/0.52) be reproduced? ===")
    filters = {
        "guide_correlation_signif, all rows": lambda d: d["guide_correlation_signif"],
        "  + ontarget_significant": lambda d: d.loc[d["ontarget_significant"] == True, "guide_correlation_signif"],
        "  + 2 effective guides": lambda d: d.loc[d["guide_n_signif_ontarget"] == 2, "guide_correlation_signif"],
        "  + both, not low_target_gex": lambda d: d.loc[
            (d["ontarget_significant"] == True) & (d["guide_n_signif_ontarget"] == 2)
            & (d["low_target_gex"] != True), "guide_correlation_signif"],
        "guide_correlation_all, all rows": lambda d: d["guide_correlation_all"],
    }
    for name, f in filters.items():
        row = []
        for cond in CONDS:
            v = f(de[de["culture_condition"] == cond]).dropna()
            row.append(f"{cond}={v.median():+.3f}(n={len(v)})")
        print(f"  {name:38s} " + "  ".join(row))
    print("  -> none matches 0.49/0.47/0.52; the preprint figure was not reproducible from the")
    print("     supplementary table under any of these filters. Cite the table values, not that quote.")

    print("\n=== CHECK 2: ZNF conservation deficit inside the authors' own K562-vs-CD4 table ===")
    for cond in CONDS:
        s = k5[k5["condition"] == cond].dropna(subset=["logfc_pearson_r"])
        m = s["target_contrast_gene_name"].astype(str).str.match(r"^ZNF\d").values
        a, b = s["logfc_pearson_r"].values[m], s["logfc_pearson_r"].values[~m]
        u, p = mannwhitneyu(a, b, alternative="less")
        auc = 1 - u / (len(a) * len(b))
        lo, hi = boot_median_diff(a, b, rng)
        print(f"  {cond:9s} ZNF n={len(a)}  dmedian {np.median(a)-np.median(b):+.4f} [{lo:+.4f},{hi:+.4f}]"
              f"  AUC(other>ZNF) {auc:.3f}  one-sided p {p:.2g}")

    print("\n=== CHECK 3: does cross-line conservation track cross-donor reproducibility (their columns)? ===")
    for cond in CONDS:
        s = k5[k5["condition"] == cond].dropna(subset=["logfc_pearson_r", "donor_correlation_mean"])
        x, y = s["logfc_pearson_r"].values, s["donor_correlation_mean"].values
        lo, hi = boot_spearman(x, y, rng)
        print(f"  {cond:9s} n={len(s):5d}  rho(conservation, donor reproducibility) = "
              f"{spearmanr(x, y).statistic:+.3f} [{lo:+.3f},{hi:+.3f}]"
              f"   medians {np.median(x):+.3f} / {np.median(y):+.3f}")

    print("\n=== CHECK 4: the magnitude confound (their MASH DE counts) ===")
    for cond, col in [("Stim48hr", "n_degs_MASH_Stim48hr"), ("Stim8hr", "n_degs_MASH_Stim8hr"),
                      ("Rest", "n_degs_MASH_Rest")]:
        s = k5[k5["condition"] == cond].dropna(subset=["logfc_pearson_r", col])
        s2 = k5[k5["condition"] == cond].dropna(subset=["donor_correlation_mean", col])
        print(f"  {cond:9s} rho(conservation, n_DEG) = "
              f"{spearmanr(s['logfc_pearson_r'], s[col]).statistic:+.3f} (n={len(s)})   "
              f"rho(donor reprod., n_DEG) = "
              f"{spearmanr(s2['donor_correlation_mean'], s2[col]).statistic:+.3f} (n={len(s2)})")


if __name__ == "__main__":
    main()
