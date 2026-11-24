"""Gate 3 substitute -- the question the mandate actually asks, measured on the surface that
can carry it.

  Does a PER-CELL-TYPE promoter-accessibility feature carry line-specific information that
  basal expression does not?

Object: the LINE-SPECIFIC deviation, i.e. each matrix centred per gene across lines -- the same
shape of object 268/289 predict on Broad-16. The LINE is the experimental unit throughout.

Readouts
  Q0  known-answer rows (the matrix is the thing it claims to be)
  Q1  redundancy: how much of the accessibility deviation is basal expression
  Q2  line similarity: does accessibility re-rank donors relative to expression
  Q3  reliability: is the expression-ORTHOGONAL part reproducible between two independent
      DNase-seq experiments of the same line, or is it assay noise
  Q4  the same on the 300 VCC panel targets alone
  Q5  depth conditioning (programme rule 9)

Controls throughout: shuffled gene pairing, shuffled line labels, and the quantile-normalised
(depth-free) rendering.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parent
RNG = np.random.default_rng(315)
OUT: dict = {}


def centre(M: np.ndarray) -> np.ndarray:
    return M - M.mean(axis=0, keepdims=True)


def colcorr(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """per-column Pearson r between two (n_lines, n_genes) matrices."""
    a, b = centre(A), centre(B)
    na, nb = np.linalg.norm(a, axis=0), np.linalg.norm(b, axis=0)
    ok = (na > 1e-9) & (nb > 1e-9)
    r = np.full(A.shape[1], np.nan)
    r[ok] = (a[:, ok] * b[:, ok]).sum(axis=0) / (na[ok] * nb[ok])
    return r


def deflate(A: np.ndarray, E: np.ndarray) -> np.ndarray:
    """remove, per gene, the component of the accessibility deviation along the expression
    deviation. One coefficient per gene, fitted on the same lines (df cost is controlled for by
    the shuffled arms, which pay exactly the same cost)."""
    a, e = centre(A), centre(E)
    den = (e * e).sum(axis=0)
    beta = np.divide((a * e).sum(axis=0), den, out=np.zeros(A.shape[1]), where=den > 1e-9)
    return a - beta * e


def flatcorr(A: np.ndarray, B: np.ndarray) -> float:
    a, b = centre(A).ravel(), centre(B).ravel()
    return float(np.corrcoef(a, b)[0, 1])


def qnorm(M: np.ndarray) -> np.ndarray:
    """quantile-normalise every line to the panel's mean distribution: removes any per-experiment
    depth / peak-count scale, keeping only the within-line ORDER of promoters."""
    ref = np.sort(M, axis=1).mean(axis=0)
    return np.asarray([ref[rankdata(row, method="ordinal").astype(int) - 1] for row in M])


d = np.load(WORK / "matrix315.npz", allow_pickle=False)
A0, E, A1 = d["A0"].astype(np.float64), d["E"].astype(np.float64), d["A1"].astype(np.float64)
lines, rep1_lines, genes = [str(x) for x in d["lines"]], [str(x) for x in d["rep1_lines"]], \
    np.asarray([str(x) for x in d["genes"]])
L, G = A0.shape
print(f"[data] {L} lines x {G} genes; {len(rep1_lines)} lines with a second DNase experiment")

# ------------------------------------------------------------------- Q0  known-answer rows
within = np.asarray([spearmanr(A0[i], E[i]).statistic for i in range(L)])
short = [l.split("_")[0] for l in lines]
anchors = {a: (a in short) for a in ("K562", "HCT116", "JURKAT", "HELA", "RPMI7951")}
OUT["Q0_known_answer"] = {
    "within_line_spearman_accessibility_vs_expression": {
        "median": float(np.median(within)), "min": float(within.min()),
        "max": float(within.max()), "n_lines_positive": int((within > 0).sum())},
    "anchor_lines_present": anchors,
    "open_promoter_fraction": {"median": float(np.median((A0 > 0).mean(axis=1))),
                               "min": float((A0 > 0).mean(axis=1).min()),
                               "max": float((A0 > 0).mean(axis=1).max())},
}
print(f"[Q0] within-line Spearman(promoter accessibility, basal expression): "
      f"median {np.median(within):+.3f}  range {within.min():+.3f}..{within.max():+.3f}  "
      f"{int((within > 0).sum())}/{L} positive")
print(f"[Q0] anchors present: {anchors}")

# ------------------------------------------------------------------------ Q1  redundancy
a_res, e_res = centre(A0), centre(E)
r_gene = colcorr(A0, E)
r_gene_ok = r_gene[np.isfinite(r_gene)]
resid = deflate(A0, E)
r2_total = 1.0 - float((resid ** 2).sum() / (a_res ** 2).sum())
# shuffled-gene control: pair each gene's accessibility with a RANDOM gene's expression
perm = RNG.permutation(G)
r2_shuf_gene = 1.0 - float((deflate(A0, E[:, perm]) ** 2).sum() / (a_res ** 2).sum())
# shuffled-line control: same gene, lines permuted
lperm = RNG.permutation(L)
r2_shuf_line = 1.0 - float((deflate(A0, E[lperm]) ** 2).sum() / (a_res ** 2).sum())
OUT["Q1_redundancy"] = {
    "per_gene_r_median": float(np.median(r_gene_ok)),
    "per_gene_r_q25": float(np.percentile(r_gene_ok, 25)),
    "per_gene_r_q75": float(np.percentile(r_gene_ok, 75)),
    "frac_genes_r_positive": float((r_gene_ok > 0).mean()),
    "variance_explained_total": r2_total,
    "variance_explained_shuffled_genes": r2_shuf_gene,
    "variance_explained_shuffled_lines": r2_shuf_line,
    "n_genes": int(len(r_gene_ok)),
}
print(f"[Q1] per-gene r(accessibility deviation, expression deviation): "
      f"median {np.median(r_gene_ok):+.3f} "
      f"[{np.percentile(r_gene_ok, 25):+.3f}, {np.percentile(r_gene_ok, 75):+.3f}]; "
      f"{100*(r_gene_ok > 0).mean():.1f}% positive")
print(f"[Q1] accessibility deviation variance explained by expression: {100*r2_total:.1f}% "
      f"(shuffled genes {100*r2_shuf_gene:.1f}%, shuffled lines {100*r2_shuf_line:.1f}%)")

# ------------------------------------------------------- Q2  does accessibility re-rank lines
SA, SE = np.corrcoef(a_res), np.corrcoef(e_res)
SR = np.corrcoef(centre(resid))
off = ~np.eye(L, dtype=bool)
mantel_ae = float(spearmanr(SA[off], SE[off]).statistic)
mantel_re = float(spearmanr(SR[off], SE[off]).statistic)
nnA = [lines[np.argsort(-np.where(np.eye(L, dtype=bool)[i], -np.inf, SA[i]))[0]] for i in range(L)]
nnE = [lines[np.argsort(-np.where(np.eye(L, dtype=bool)[i], -np.inf, SE[i]))[0]] for i in range(L)]
nnR = [lines[np.argsort(-np.where(np.eye(L, dtype=bool)[i], -np.inf, SR[i]))[0]] for i in range(L)]
OUT["Q2_line_similarity"] = {
    "mantel_spearman_accessibility_vs_expression": mantel_ae,
    "mantel_spearman_orthogonal_vs_expression": mantel_re,
    "nearest_neighbour_agreement_access_vs_expr": float(np.mean(
        [a == b for a, b in zip(nnA, nnE)])),
    "nearest_neighbour_agreement_orth_vs_expr": float(np.mean(
        [a == b for a, b in zip(nnR, nnE)])),
    "nn_table": [{"line": lines[i], "nn_expression": nnE[i], "nn_accessibility": nnA[i],
                  "nn_orthogonal": nnR[i]} for i in range(L)],
}
print(f"[Q2] Mantel Spearman(line-similarity from accessibility, from expression) = {mantel_ae:+.3f}; "
      f"same nearest neighbour for {100*np.mean([a == b for a, b in zip(nnA, nnE)]):.0f}% of lines")

# ------------------------------------ Q3  is the expression-orthogonal part reproducible at all
sub = [lines.index(l) for l in rep1_lines]
A0s, A1s, Es = A0[sub], A1, E[sub]
q3 = {"n_lines": len(sub), "lines": rep1_lines}
q3["reliability_raw"] = flatcorr(A0s, A1s)
q3["reliability_orthogonal"] = flatcorr(deflate(A0s, Es), deflate(A1s, Es))
# controls that pay the SAME degrees of freedom as the deflation
pg, pl = RNG.permutation(G), RNG.permutation(len(sub))
q3["reliability_deflated_by_shuffled_genes"] = flatcorr(
    deflate(A0s, Es[:, pg]), deflate(A1s, Es[:, pg]))
q3["reliability_deflated_by_shuffled_lines"] = flatcorr(
    deflate(A0s, Es[pl]), deflate(A1s, Es[pl]))
# how much line-specific accessibility survives at all, per line
q3["per_line_reliability"] = {
    rep1_lines[i]: float(np.corrcoef(centre(A0s)[i], centre(A1s)[i])[0, 1]) for i in range(len(sub))}
OUT["Q3_reliability"] = q3
print(f"[Q3] split-experiment reliability of the line-specific accessibility deviation "
      f"({len(sub)} lines): raw {q3['reliability_raw']:+.3f}, "
      f"expression-orthogonal {q3['reliability_orthogonal']:+.3f} "
      f"(shuffled-gene deflation {q3['reliability_deflated_by_shuffled_genes']:+.3f}, "
      f"shuffled-line deflation {q3['reliability_deflated_by_shuffled_lines']:+.3f})")

# ------------------------------------------------------------------ Q4  the 300 panel targets
panel = pd.read_csv(ROOT / "data/challenge/pert_counts.csv")
pcol = [c for c in panel.columns if "target" in c.lower() or "gene" in c.lower()][0]
targets = set(panel[pcol].astype(str))
mask = np.isin(genes, list(targets))
q4 = {"n_panel_targets_on_axis": int(mask.sum()), "n_panel_targets_total": len(targets)}
if mask.sum() >= 20:
    rp = colcorr(A0[:, mask], E[:, mask])
    rp = rp[np.isfinite(rp)]
    q4["per_gene_r_median"] = float(np.median(rp))
    q4["variance_explained_total"] = 1.0 - float(
        (deflate(A0[:, mask], E[:, mask]) ** 2).sum() / (centre(A0[:, mask]) ** 2).sum())
    mp = np.isin(genes, list(targets))
    q4["reliability_raw"] = flatcorr(A0s[:, mp], A1s[:, mp])
    q4["reliability_orthogonal"] = flatcorr(
        deflate(A0s[:, mp], Es[:, mp]), deflate(A1s[:, mp], Es[:, mp]))
OUT["Q4_panel_targets"] = q4
print(f"[Q4] {q4['n_panel_targets_on_axis']}/{q4['n_panel_targets_total']} VCC panel targets on "
      f"the axis; per-gene r median {q4.get('per_gene_r_median', float('nan')):+.3f}; "
      f"variance explained {100*q4.get('variance_explained_total', float('nan')):.1f}%; "
      f"orthogonal reliability {q4.get('reliability_orthogonal', float('nan')):+.3f}")

# ----------------------------------------------------- Q5  depth conditioning (programme rule 9)
n_open = (A0 > 0).sum(axis=1)
A0q, A1q = qnorm(A0), qnorm(A1)
q5 = {
    "spearman_open_promoters_vs_residual_norm": float(spearmanr(
        n_open, np.linalg.norm(a_res, axis=1)).statistic),
    "n_open_min": int(n_open.min()), "n_open_max": int(n_open.max()),
    "qnorm_variance_explained_total": 1.0 - float(
        (deflate(A0q, E) ** 2).sum() / (centre(A0q) ** 2).sum()),
    "qnorm_reliability_raw": flatcorr(A0q[sub], A1q),
    "qnorm_reliability_orthogonal": flatcorr(deflate(A0q[sub], Es), deflate(A1q, Es)),
    "qnorm_mantel_vs_expression": float(spearmanr(
        np.corrcoef(centre(A0q))[off], SE[off]).statistic),
}
OUT["Q5_depth"] = q5
print(f"[Q5] open promoters {n_open.min()}..{n_open.max()}; Spearman with residual norm "
      f"{q5['spearman_open_promoters_vs_residual_norm']:+.3f}; after quantile normalisation: "
      f"variance explained {100*q5['qnorm_variance_explained_total']:.1f}%, "
      f"reliability raw {q5['qnorm_reliability_raw']:+.3f} / "
      f"orthogonal {q5['qnorm_reliability_orthogonal']:+.3f}")

(WORK / "analysis315.json").write_text(json.dumps(OUT, indent=2) + "\n")
pd.DataFrame(OUT["Q2_line_similarity"]["nn_table"]).to_csv(WORK / "nn315.csv", index=False)
print(f"\nwrote {WORK / 'analysis315.json'}")
