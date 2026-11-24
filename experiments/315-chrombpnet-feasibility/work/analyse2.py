"""Gate 3 substitute, primary rendering.

`analyse.py` established (Q5) that raw narrowPeak `signalValue` is NOT comparable across ENCODE
experiments: the split-experiment reliability of the raw line-specific deviation is -0.020, and
rises to +0.425 once every line is quantile-normalised to a common distribution. Peak-caller
score scale is a per-experiment nuisance exactly as programme rule 9 anticipates, so the PRIMARY
rendering here is quantile-normalised (within-line promoter ORDER only, which is also what a
per-cell-type accessibility feature would mean); raw is reported alongside.

Two QC rules, both data-quality and both declared before any statistic below was read:
  * drop a line whose peak call leaves < 25 % of promoters open (a failed or ultra-shallow
    experiment, not a closed genome);
  * drop a line whose second experiment fails the same rule, for the reliability arm only.

The experimental unit is the LINE. Every interval is a bootstrap over lines.
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
MIN_OPEN_FRAC = 0.25
NBOOT = 2000
OUT: dict = {}


def centre(M):
    return M - M.mean(axis=0, keepdims=True)


def qnorm(M):
    ref = np.sort(M, axis=1).mean(axis=0)
    return np.asarray([ref[rankdata(r, method="ordinal").astype(int) - 1] for r in M])


def colcorr(A, B):
    a, b = centre(A), centre(B)
    na, nb = np.linalg.norm(a, axis=0), np.linalg.norm(b, axis=0)
    ok = (na > 1e-9) & (nb > 1e-9)
    r = np.full(A.shape[1], np.nan)
    r[ok] = (a[:, ok] * b[:, ok]).sum(axis=0) / (na[ok] * nb[ok])
    return r


def deflate(A, E):
    a, e = centre(A), centre(E)
    den = (e * e).sum(axis=0)
    beta = np.divide((a * e).sum(axis=0), den, out=np.zeros(A.shape[1]), where=den > 1e-9)
    return a - beta * e


def var_explained(A, E):
    return 1.0 - float((deflate(A, E) ** 2).sum() / (centre(A) ** 2).sum())


def paired_reliability(A0, A1):
    """mean same-line agreement of the line-specific deviation, and its wrong-recipient control
    (rep0 of line i against rep1 of line j != i)."""
    a0, a1 = centre(A0), centre(A1)
    a0 = a0 / np.maximum(np.linalg.norm(a0, axis=1, keepdims=True), 1e-12)
    a1 = a1 / np.maximum(np.linalg.norm(a1, axis=1, keepdims=True), 1e-12)
    M = a0 @ a1.T
    n = M.shape[0]
    same = np.diag(M).copy()
    diff = M[~np.eye(n, dtype=bool)]
    return same, float(same.mean()), float(diff.mean())


def boot_lines(fn, n, nboot=NBOOT):
    vals = []
    for _ in range(nboot):
        idx = RNG.integers(0, n, n)
        try:
            vals.append(fn(idx))
        except Exception:  # noqa: BLE001
            pass
    v = np.asarray(vals)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


d = np.load(WORK / "matrix315.npz", allow_pickle=False)
A0, E, A1 = d["A0"].astype(np.float64), d["E"].astype(np.float64), d["A1"].astype(np.float64)
lines = [str(x) for x in d["lines"]]
rep1_lines = [str(x) for x in d["rep1_lines"]]
genes = np.asarray([str(x) for x in d["genes"]])

# --------------------------------------------------------------------------------- QC filter
open0 = (A0 > 0).mean(axis=1)
keep = open0 >= MIN_OPEN_FRAC
dropped = [(lines[i], float(open0[i])) for i in range(len(lines)) if not keep[i]]
A0, E = A0[keep], E[keep]
lines = [l for l, k in zip(lines, keep) if k]
open1 = (A1 > 0).mean(axis=1)
keep1 = np.asarray([(l in lines) and (open1[i] >= MIN_OPEN_FRAC)
                    for i, l in enumerate(rep1_lines)])
A1, rep1_lines = A1[keep1], [l for l, k in zip(rep1_lines, keep1) if k]
L, G = A0.shape
sub = [lines.index(l) for l in rep1_lines]
print(f"[qc] dropped {len(dropped)} lines below {MIN_OPEN_FRAC:.0%} open promoters: {dropped}")
print(f"[data] {L} lines x {G} genes; reliability arm on {len(rep1_lines)} lines")
OUT["qc"] = {"min_open_fraction": MIN_OPEN_FRAC, "dropped": dropped, "n_lines": L,
             "n_genes": G, "reliability_lines": rep1_lines}

RENDER = {"qnorm": (qnorm(A0), qnorm(A1)), "raw": (A0, A1)}

for name, (Aq, Aq1) in RENDER.items():
    r = {}
    Es = E[sub]
    Aqs = Aq[sub]

    # --- Q1 redundancy ---------------------------------------------------------------------
    rg = colcorr(Aq, E)
    rg = rg[np.isfinite(rg)]
    pg, pl = RNG.permutation(G), RNG.permutation(L)
    r["Q1"] = {
        "per_gene_r_median": float(np.median(rg)),
        "per_gene_r_q25": float(np.percentile(rg, 25)),
        "per_gene_r_q75": float(np.percentile(rg, 75)),
        "frac_positive": float((rg > 0).mean()),
        "variance_explained": var_explained(Aq, E),
        "control_shuffled_genes": var_explained(Aq, E[:, pg]),
        "control_shuffled_lines": var_explained(Aq, E[pl]),
    }
    r["Q1"]["ci_variance_explained"] = boot_lines(
        lambda i: var_explained(Aq[i], E[i]), L, 400)

    # --- Q2 line similarity ----------------------------------------------------------------
    SA, SE = np.corrcoef(centre(Aq)), np.corrcoef(centre(E))
    SR = np.corrcoef(centre(deflate(Aq, E)))
    off = ~np.eye(L, dtype=bool)
    nn = lambda S: [lines[np.argsort(-np.where(np.eye(L, dtype=bool)[i], -np.inf, S[i]))[0]]  # noqa: E731
                    for i in range(L)]
    nnA, nnE, nnR = nn(SA), nn(SE), nn(SR)
    r["Q2"] = {
        "mantel_access_vs_expr": float(spearmanr(SA[off], SE[off]).statistic),
        "mantel_orthogonal_vs_expr": float(spearmanr(SR[off], SE[off]).statistic),
        "nn_agree_access_vs_expr": float(np.mean([a == b for a, b in zip(nnA, nnE)])),
        "nn_agree_orth_vs_expr": float(np.mean([a == b for a, b in zip(nnR, nnE)])),
        "nn_same_lineage_expr": float(np.mean(
            [a.split("_", 1)[1] == b.split("_", 1)[1] for a, b in zip(lines, nnE)])),
        "nn_same_lineage_access": float(np.mean(
            [a.split("_", 1)[1] == b.split("_", 1)[1] for a, b in zip(lines, nnA)])),
        "nn_same_lineage_orth": float(np.mean(
            [a.split("_", 1)[1] == b.split("_", 1)[1] for a, b in zip(lines, nnR)])),
    }
    if name == "qnorm":
        pd.DataFrame({"line": lines, "nn_expression": nnE, "nn_accessibility": nnA,
                      "nn_orthogonal": nnR}).to_csv(WORK / "nn315.csv", index=False)

    # --- Q3 reliability + wrong-recipient control ------------------------------------------
    same_raw, m_raw, m_raw_wrong = paired_reliability(Aqs, Aq1)
    d0, d1 = deflate(Aqs, Es), deflate(Aq1, Es)
    same_orth, m_orth, m_orth_wrong = paired_reliability(d0, d1)
    pgs, pls = RNG.permutation(G), RNG.permutation(len(sub))
    _, m_sg, _ = paired_reliability(deflate(Aqs, Es[:, pgs]), deflate(Aq1, Es[:, pgs]))
    _, m_sl, _ = paired_reliability(deflate(Aqs, Es[pls]), deflate(Aq1, Es[pls]))
    n1 = len(sub)
    r["Q3"] = {
        "n_lines": n1,
        "same_line_raw": m_raw, "wrong_recipient_raw": m_raw_wrong,
        "same_line_orthogonal": m_orth, "wrong_recipient_orthogonal": m_orth_wrong,
        "control_deflate_by_shuffled_genes": m_sg,
        "control_deflate_by_shuffled_lines": m_sl,
        "ci_same_line_raw": boot_lines(lambda i: float(same_raw[i].mean()), n1),
        "ci_same_line_orthogonal": boot_lines(lambda i: float(same_orth[i].mean()), n1),
        "per_line_orthogonal": {rep1_lines[i]: float(same_orth[i]) for i in range(n1)},
        "n_lines_orthogonal_positive": int((same_orth > 0).sum()),
    }

    # --- Q4 the 300 VCC panel targets ------------------------------------------------------
    panel = pd.read_csv(ROOT / "data/challenge/pert_counts.csv")
    pcol = [c for c in panel.columns if "target" in c.lower() or "gene" in c.lower()][0]
    m = np.isin(genes, list(set(panel[pcol].astype(str))))
    so, mo, mw = paired_reliability(deflate(Aqs[:, m], Es[:, m]), deflate(Aq1[:, m], Es[:, m]))
    sr, mr, mrw = paired_reliability(Aqs[:, m], Aq1[:, m])
    r["Q4"] = {
        "n_targets_on_axis": int(m.sum()),
        "variance_explained": var_explained(Aq[:, m], E[:, m]),
        "same_line_raw": mr, "wrong_recipient_raw": mrw,
        "same_line_orthogonal": mo, "wrong_recipient_orthogonal": mw,
        "ci_same_line_orthogonal": boot_lines(lambda i: float(so[i].mean()), len(sub)),
    }
    OUT[name] = r

    print(f"\n===== rendering: {name} =====")
    print(f"[Q1] per-gene r median {r['Q1']['per_gene_r_median']:+.3f} "
          f"[{r['Q1']['per_gene_r_q25']:+.3f}, {r['Q1']['per_gene_r_q75']:+.3f}]; "
          f"variance explained {100*r['Q1']['variance_explained']:.1f}% "
          f"[{100*r['Q1']['ci_variance_explained'][0]:.1f}, "
          f"{100*r['Q1']['ci_variance_explained'][1]:.1f}] "
          f"(shuf-gene {100*r['Q1']['control_shuffled_genes']:.1f}%, "
          f"shuf-line {100*r['Q1']['control_shuffled_lines']:.1f}%)")
    print(f"[Q2] Mantel access-vs-expr {r['Q2']['mantel_access_vs_expr']:+.3f}, "
          f"orth-vs-expr {r['Q2']['mantel_orthogonal_vs_expr']:+.3f}; nn agreement "
          f"{100*r['Q2']['nn_agree_access_vs_expr']:.0f}% / "
          f"{100*r['Q2']['nn_agree_orth_vs_expr']:.0f}%; nn same lineage "
          f"expr {100*r['Q2']['nn_same_lineage_expr']:.0f}% "
          f"access {100*r['Q2']['nn_same_lineage_access']:.0f}% "
          f"orth {100*r['Q2']['nn_same_lineage_orth']:.0f}%")
    print(f"[Q3] n={r['Q3']['n_lines']} raw same-line {r['Q3']['same_line_raw']:+.3f} "
          f"{tuple(round(x,3) for x in r['Q3']['ci_same_line_raw'])} vs wrong-recipient "
          f"{r['Q3']['wrong_recipient_raw']:+.3f}; ORTHOGONAL same-line "
          f"{r['Q3']['same_line_orthogonal']:+.3f} "
          f"{tuple(round(x,3) for x in r['Q3']['ci_same_line_orthogonal'])} vs wrong-recipient "
          f"{r['Q3']['wrong_recipient_orthogonal']:+.3f} "
          f"(shuf-gene {r['Q3']['control_deflate_by_shuffled_genes']:+.3f}, "
          f"shuf-line {r['Q3']['control_deflate_by_shuffled_lines']:+.3f}); "
          f"{r['Q3']['n_lines_orthogonal_positive']}/{r['Q3']['n_lines']} lines positive")
    print(f"[Q4] {r['Q4']['n_targets_on_axis']} panel targets: variance explained "
          f"{100*r['Q4']['variance_explained']:.1f}%; raw same {r['Q4']['same_line_raw']:+.3f} "
          f"vs wrong {r['Q4']['wrong_recipient_raw']:+.3f}; orthogonal same "
          f"{r['Q4']['same_line_orthogonal']:+.3f} "
          f"{tuple(round(x,3) for x in r['Q4']['ci_same_line_orthogonal'])} "
          f"vs wrong {r['Q4']['wrong_recipient_orthogonal']:+.3f}")

(WORK / "analysis315_primary.json").write_text(json.dumps(OUT, indent=2) + "\n")
print(f"\nwrote {WORK / 'analysis315_primary.json'}")
