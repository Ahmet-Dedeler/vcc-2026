"""218 §Phase 0 — the candidate source registry and its panel census.

Every table on disk that is shaped like a donor bank (`delta.npy`, `b_ctrl.npy`,
`perts.txt`, `genes.txt`) and could in principle serve as a post-whitening
un-whitened residual.  A source can only serve if it carries enough of the 300-gene
panel that the MEDIAN row-norm sizing rule is non-vacuous (200/207's trap), so the
census computes panel coverage per PLOLO fold, not just globally.

Metadata columns are the independence axes the mandate names: platform / chemistry,
lineage, transformation state, and the lab.  They are recorded here once, by hand,
from each source's own provenance, so no downstream script has to guess.

    uv run python sources218.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

PANEL = tuple(
    line.split(",")[0]
    for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
    if line.strip()
)
assert len(PANEL) == 300

# name -> (path, platform, chemistry, lineage, transformation, lab, value semantics)
SOURCES: dict[str, dict] = {
    "k562_gwps": dict(
        path="data/donors/k562_gwps", platform="Replogle/Weissman CRISPRi",
        chemistry="10x 3' v3 (Perturb-seq)", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Weissman 2022",
        semantics="log-pseudobulk delta", role="PRODUCTION donor"),
    "orion_hct116": dict(
        path="data/donors/orion_hct116", platform="X-Atlas/Orion CRISPRi",
        chemistry="Xaira in-house barcoded (X-Atlas)", lineage="epithelial (colorectal)",
        transformation="immortalised cancer line", lab="Xaira 2026",
        semantics="log-pseudobulk delta", role="PRODUCTION donor"),
    "orion_hek293t": dict(
        path="data/donors/orion_hek293t", platform="X-Atlas/Orion CRISPRi",
        chemistry="Xaira in-house barcoded (X-Atlas)", lineage="epithelial (kidney)",
        transformation="immortalised, SV40/adenovirus-transformed", lab="Xaira 2026",
        semantics="log-pseudobulk delta", role="PRODUCTION donor"),
    "kolf2": dict(
        path="data/external/kolf2", platform="KOLF2.1J atlas CRISPRi",
        chemistry="10x (published z-score matrix)", lineage="pluripotent (iPSC)",
        transformation="NON-transformed iPSC", lab="Y-Doctor/JAX 2025",
        semantics="per-gene z-score (NOT log-fold-change)", role="candidate"),
    "feng_ipsc": dict(
        path="experiments/099-donors/work/feng_ipsc", platform="Feng iPSC CRISPRi atlas",
        chemistry="10x, 24-34 HipSci donor lines pooled",
        lineage="pluripotent (iPSC, multi-donor)",
        transformation="NON-transformed iPSC", lab="Feng 2025 (Cell Genomics)",
        semantics="pooled target-level lfc", role="candidate"),
    "gse295214_h1": dict(
        path="experiments/115-gse295214-h1/work/compact-donor",
        platform="Pan 2025 hESC CRISPRi", chemistry="10x",
        lineage="pluripotent (hESC)", transformation="NON-transformed hESC",
        lab="Pan 2025", semantics="log-pseudobulk delta", role="candidate"),
    "telohaec": dict(
        path="data/external/igvf-telohaec", platform="IGVF TeloHAEC CRISPRi",
        chemistry="10x", lineage="endothelial (aortic)",
        transformation="hTERT-immortalised primary", lab="IGVF consortium",
        semantics="log-pseudobulk delta", role="candidate"),
    "myllia_a549": dict(
        path="experiments/120-myllia-a549/work/compact-donor",
        platform="Myllia A549 screen", chemistry="10x",
        lineage="epithelial (lung adeno)", transformation="immortalised cancer line",
        lab="Myllia", semantics="log-pseudobulk delta", role="candidate"),
    "h1_2025": dict(
        path="data/donors/h1_2025", platform="VCC2025 H1 CRISPRi",
        chemistry="10x Flex", lineage="pluripotent (hESC)",
        transformation="NON-transformed hESC", lab="Arc 2025",
        semantics="log-pseudobulk delta", role="candidate (leaks h1-test)"),
    "nadig_jurkat": dict(
        path="data/donors/nadig_jurkat", platform="Nadig CRISPRi (essential library)",
        chemistry="10x", lineage="lymphoid T (Jurkat)",
        transformation="immortalised cancer line", lab="Nadig 2025",
        semantics="log-pseudobulk delta", role="training-only (0/300)"),
    "nadig_hepg2": dict(
        path="experiments/166-donor-acquisition/donors/nadig_hepg2",
        platform="Nadig CRISPRi (essential library)", chemistry="10x",
        lineage="hepatocyte (hepatoblastoma)",
        transformation="immortalised cancer line", lab="Nadig 2025",
        semantics="log-pseudobulk delta", role="training-only (0/300)"),
    "rpe1": dict(
        path="data/donors/rpe1", platform="Replogle/Weissman CRISPRi",
        chemistry="10x 3' v3", lineage="epithelial (retinal pigment)",
        transformation="hTERT-immortalised primary", lab="Weissman 2022",
        semantics="log-pseudobulk delta", role="training-only (0/300)"),
    "k562_essential": dict(
        path="data/donors/k562_essential", platform="Replogle/Weissman CRISPRi",
        chemistry="10x 3' v3", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Weissman 2022",
        semantics="log-pseudobulk delta", role="training-only (0/300)"),
    "frangieh_melanoma": dict(
        path="experiments/166-donor-acquisition/donors/frangieh_melanoma",
        platform="Perturb-CITE-seq Cas9 KO", chemistry="10x, MOI>1",
        lineage="melanocyte (patient melanoma)",
        transformation="patient-derived cancer line", lab="Frangieh 2021",
        semantics="log-pseudobulk delta", role="training-only (5/300)"),
    "perturbmatch": dict(
        path="data/external/perturbmatch", platform="PerturbMatch K562 reprocessing",
        chemistry="10x 3' v3 (same cells as k562_gwps)", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="reprocessing of Weissman 2022",
        semantics="log-pseudobulk delta", role="already blended into k562 slot"),
}

# CD4 is three banks, handled by the production `prepare_cd4_mean` path, so it is
# described here but its coverage is read from the compact banks directly.
CD4_META = dict(
    path="data/external/cd4t/cloud-compact/cd4_compact",
    platform="Primary CD4+ T CRISPRi (GSE314342)",
    chemistry="10x, 4 human donors x 3 stimulation states",
    lineage="lymphoid T (PRIMARY human cells)",
    transformation="NON-transformed primary cells", lab="Schmidt/Marson-adjacent 2026",
    semantics="log2FC x ln2", role="PRODUCTION residual (the arm this mandate follows)")

FOLD_TARGET_CACHE = HERE / "work/plolo"


def fold_targets() -> dict[str, list[str]]:
    out = {}
    for f in sorted(FOLD_TARGET_CACHE.glob("fold-*.npz")):
        held = f.stem[len("fold-"):]
        with np.load(f, allow_pickle=True) as z:
            out[held] = [str(t) for t in z["targets"]]
    return out


def main() -> None:
    panel = set(PANEL)
    folds = fold_targets()
    rows = []
    for name, meta in SOURCES.items():
        p = ROOT / meta["path"]
        perts = [x for x in (p / "perts.txt").read_text().splitlines() if x.strip()]
        genes = [x for x in (p / "genes.txt").read_text().splitlines() if x.strip()]
        delta = np.load(p / "delta.npy", mmap_mode="r")
        have = set(perts)
        cov = {h: sum(t in have for t in ts) for h, ts in folds.items()}
        frac = {h: cov[h] / len(ts) for h, ts in folds.items()}
        rows.append(dict(
            name=name, **{k: v for k, v in meta.items() if k != "path"},
            path=meta["path"], n_targets=len(perts), n_genes=len(genes),
            shape=list(map(int, delta.shape)),
            panel=len(have & panel), panel_pct=round(100 * len(have & panel) / 300, 1),
            fold_coverage=cov, fold_coverage_frac={k: round(v, 4) for k, v in frac.items()},
            median_nonvacuous={h: bool(frac[h] > 0.5) for h in folds},
        ))
    # CD4: three banks, a target counts only if all three carry it (production rule).
    cd4root = ROOT / CD4_META["path"]
    banks = {}
    for cond in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
        banks[cond] = set(
            x for x in (cd4root / cond / "perts.txt").read_text().splitlines() if x.strip())
    complete = set.intersection(*banks.values())
    cov = {h: sum(t in complete for t in ts) for h, ts in folds.items()}
    frac = {h: cov[h] / len(ts) for h, ts in folds.items()}
    rows.append(dict(
        name="cd4", **{k: v for k, v in CD4_META.items() if k != "path"},
        path=CD4_META["path"], n_targets=len(complete),
        n_genes=10113, shape=[len(complete), 10113],
        panel=len(complete & panel), panel_pct=round(100 * len(complete & panel) / 300, 1),
        fold_coverage=cov, fold_coverage_frac={k: round(v, 4) for k, v in frac.items()},
        median_nonvacuous={h: bool(frac[h] > 0.5) for h in folds},
        per_condition={k: len(v) for k, v in banks.items()}))

    (WORK / "sources218.json").write_text(json.dumps(rows, indent=2))
    hdr = f"{'source':18s} {'panel':>6s} {'targets':>8s} {'genes':>6s}  " \
          f"{'k562':>5s} {'hct':>5s} {'hek':>5s} {'kolf':>5s}  transformation"
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda r: -r["panel"]):
        c = r["fold_coverage"]
        print(f"{r['name']:18s} {r['panel']:6d} {r['n_targets']:8d} {r['n_genes']:6d}  "
              f"{c.get('k562_gwps',0):5d} {c.get('orion_hct116',0):5d} "
              f"{c.get('orion_hek293t',0):5d} {c.get('kolf2',0):5d}  {r['transformation']}")
    print("\nfold sizes:", {h: len(t) for h, t in folds.items()})


if __name__ == "__main__":
    main()
