"""278 item 1 -- assemble the per-(panel target, line, source) DEPTH table from disk.

Every number here is a cell count read from an obs table, a supplementary table shipped with the
dataset, or a QC table shipped with the dataset -- never a count matrix.  Sources with no cell
count on disk (published DE / lfc / z-score matrices) are carried with cells = NaN and a
`depth_source` of "not-on-disk" so they are visible as gaps rather than silently absent.

Writes work/depth278.csv (long: one row per (target, source)) and work/depth278_summary.json.

    PYTHONDONTWRITEBYTECODE=1 uv run python build278.py
"""
from __future__ import annotations

import json
import resource
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/278-depth-inventory"
WORK = HERE / "work"
H5 = WORK / "h5ad_counts"

PANEL = [l.strip() for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
PANEL_SET = set(PANEL)
assert len(PANEL) == 300, len(PANEL)

# state: "somatic" = a differentiated / somatic cell (line or primary); "pluripotent" = hESC / iPSC.
# `derived` marks a somatic cell type differentiated from a pluripotent source.
ROWS: list[dict] = []


def add(target: str, source: str, line: str, state: str, chem: str, unit: str,
        cells, ctrl=np.nan, guides=np.nan, batches=np.nan, depth_source: str = "", note: str = "") -> None:
    ROWS.append(dict(target=target, source=source, line=line, state=state, chemistry=chem,
                     unit=unit, cells=cells, control_cells=ctrl, guides=guides, batches=batches,
                     depth_source=depth_source, note=note))


# ---------------------------------------------------------------- 1. the three PLOLO folds (267)
FOLDS = {
    "k562_gwps": ("K562", "somatic", "10x 3' Perturb-seq (CRISPRi)", 75_328, "R267 §1 inv267.py"),
    "orion_hct116": ("HCT116", "somatic", "X-Atlas/Orion scRNA (CRISPRi)", 165_931, "R267 §1 inv267.py"),
    "orion_hek293t": ("HEK293T", "somatic", "X-Atlas/Orion scRNA (CRISPRi)", 219_071, "R267 §1 inv267.py"),
}
for key, (line, state, chem, nctrl, prov) in FOLDS.items():
    df = pd.read_csv(ROOT / f"experiments/267-replicate-validated-oracle/work/inv/{key}-targets.csv")
    for r in df.itertuples():
        if r.target not in PANEL_SET:
            continue
        add(r.target, key, line, state, chem, "cells", int(r.cells), nctrl, int(r.guides),
            int(r.batches), prov, "" if r.in_bank else "below the builder's 20-cell floor")

# ---------------------------------------------------------------- 2. CD4 T (GWCD4i, Cell 2026)
cd4 = pd.read_csv(ROOT / "data/external/cd4t/DE_stats.suppl_table.csv")
CD4_COND = {"Rest": "cd4_rest", "Stim8hr": "cd4_stim8hr", "Stim48hr": "cd4_stim48hr"}
for r in cd4.itertuples():
    t = r.target_contrast_gene_name
    if t not in PANEL_SET or r.culture_condition not in CD4_COND:
        continue
    add(t, CD4_COND[r.culture_condition], f"CD4 T ({r.culture_condition})", "somatic",
        "10x CRISPRi (primary T cells, pooled donors)", "cells", float(r.n_cells_target),
        np.nan, np.nan, np.nan, "cd4t/DE_stats.suppl_table.csv n_cells_target")

# ---------------------------------------------------------------- 3. IGVF TeloHAEC
tel = pd.read_csv(WORK / "igvf_telohaec_panel.csv")
for r in tel.itertuples():
    add(r.clean_target, "igvf_telohaec", "TeloHAEC", "somatic", "10x CRISPRi", "cells",
        float(r.target_elig), np.nan, float(r.guides), np.nan,
        "igvf-telohaec/guide_qc.csv target_eligible_cells")

# ---------------------------------------------------------------- 4. IGVF cardiomyocyte
cm = pd.read_csv(WORK / "igvf_cardiomyocyte_panel.csv")
for r in cm.itertuples():
    add(r.intended_target_name, "igvf_cardiomyocyte", "cardiomyocyte (iPSC-derived)", "derived",
        "10x CRISPRi", "cells", float(r.num_cell), np.nan, np.nan, np.nan,
        "igvf-cardiomyocyte IGVFFI0830FXFI.csv.gz num_cell",
        "on-disk effect is a SIGNIFICANT-HITS-ONLY DE table (median 6 genes/panel row)")

# ---------------------------------------------------------------- 5. Altair A549 (pseudobulk h5ad with n_cells)
alt = pd.read_csv(WORK / "altair_a549_panel.csv")
for r in alt.itertuples():
    add(r.gene, "altair_a549", "A549", "somatic", "10x CRISPRi", "cells", float(r.n_cells),
        45_320, float(r.n_guides), np.nan, "altair-a549/altair-scrna.h5ad obs n_cells")

# ---------------------------------------------------------------- 6. h5ad obs scans
def obs_rows(key: str, source: str, line: str, state: str, chem: str, note: str = "",
             ctrl_labels: tuple[str, ...] = ("control", "non-targeting", "Non-Targeting", "NonTarget",
                                             "unperturbed", "NT", "neg")) -> None:
    p = H5 / f"{key}.json"
    if not p.exists():
        return
    d = json.loads(p.read_text())
    counts = d.get("counts", {})
    nctrl = sum(int(v) for k, v in counts.items() if k in ctrl_labels)
    gp = d.get("guides_per_panel_target", {})
    for t, n in d.get("panel_counts", {}).items():
        add(t, source, line, state, chem, "cells", float(n), float(nctrl) if nctrl else np.nan,
            float(gp[t]) if t in gp else np.nan, float(d.get("n_batchs", np.nan) or np.nan),
            f"obs scan of {d['file']}", note)


obs_rows("myllia_a549_ch1", "myllia_a549_ch1", "A549", "somatic", "10x CRISPRi (channel 1)")
obs_rows("myllia_a549_ch2", "myllia_a549_ch2", "A549", "somatic", "10x CRISPRi (channel 2)")
SCP = {
    "SunshineHein2023": ("sunshine_calu3", "Calu-3", "somatic", "10x CRISPR-KO"),
    "TianKampmann2021_CRISPRi": ("tian_ipsc_neuron", "iPSC-neuron", "derived", "10x CRISPRi"),
    "TianKampmann2021_CRISPRa": ("tian_ipsc_neuron_crispra", "iPSC-neuron", "derived", "10x CRISPRa"),
    "TianKampmann2019_iPSC": ("tian2019_ipsc", "iPSC", "pluripotent", "10x CRISPRi"),
    "TianKampmann2019_day7neuron": ("tian2019_neuron", "iPSC-neuron d7", "derived", "10x CRISPRi"),
    "XuCao2023": ("xucao_hek293", "HEK293", "somatic", "10x CRISPRi"),
    "FrangiehIzar2021_RNA": ("frangieh_melanoma", "melanoma (patient)", "somatic", "10x CRISPR-KO"),
    "NormanWeissman2019_filtered": ("norman_k562", "K562", "somatic", "10x CRISPRa"),
    "GasperiniShendure2019_atscale": ("gasperini_k562_atscale", "K562", "somatic", "10x CRISPRi"),
    "GasperiniShendure2019_lowMOI": ("gasperini_k562_lowmoi", "K562", "somatic", "10x CRISPRi"),
    "WesselsSatija2023": ("wessels_thp1", "THP-1", "somatic", "10x CRISPR"),
    "PapalexiSatija2021_eccite_arrayed_RNA": ("papalexi_thp1", "THP-1", "somatic", "10x ECCITE"),
    "ReplogleWeissman2022_K562_gwps": ("k562_gwps_raw", "K562", "somatic", "10x 3' Perturb-seq"),
}
for stem, (src, line, state, chem) in SCP.items():
    obs_rows(f"scperturb__{stem}", src, line, state, chem)

# ---------------------------------------------------------------- 7. Feng HipSci iPSC (cell metadata TSVs)
for name, f, screen in [("feng_targeted", "TargetedScreen_Cell-Metadata.tsv.gz", "targeted"),
                        ("feng_gw_nonfitness", "GenomeWideScreen_NonFitnessGenes_Cell-Metadata.tsv.gz", "gw non-fitness"),
                        ("feng_gw_fitness", "GenomeWideScreen_FitnessGenes_Cell-Metadata.tsv.gz", "gw fitness")]:
    d = pd.read_csv(ROOT / "data/external/feng-ipsc" / f, sep="\t")
    d["gene"] = d.Guide_Call.astype(str).str.split("_").str[0]
    nctrl = int((d.gene == "NonTarget").sum())
    unass = int((d.Guide_Call.astype(str) == "unassigned").sum())
    # the two genome-wide files carry almost no NonTarget_ calls (12 / 36 cells) and 71 % / 58 %
    # unassigned cells: their control set is not identifiable in this metadata, so do not pretend it is
    ctrl = float(nctrl) if nctrl > 1000 else np.nan
    s = d[d.gene.isin(PANEL_SET)]
    per = s.groupby("gene").agg(cells=("Cell_ID", "size"), lines=("Cell_Line", "nunique"),
                                batches=("Batch", "nunique"), guides=("Guide_Call", "nunique"))
    for t, r in per.iterrows():
        add(t, name, "iPSC (HipSci, %d donor lines pooled)" % d.Cell_Line.nunique(), "pluripotent",
            "10x CRISPRi", "cells (pooled over donor lines)", float(r.cells), ctrl,
            float(r.guides), float(r.batches), f"feng-ipsc {f}",
            f"{screen}; {int(r.lines)} donor lines; file has {unass:,} unassigned cells "
            f"and {nctrl} NonTarget_ cells")
    per.to_csv(WORK / f"{name}_panel.csv")

# ---------------------------------------------------------------- 8. banks with NO cell count on disk
NO_CELLS = [
    ("kolf2", "KOLF2.1J iPSC", "pluripotent", "scRNA CRISPRi", "published per-target z-score matrix (int8)"),
    ("gse295214_h1", "H1 hESC", "pluripotent", "10x CRISPRi", "pseudobulk database h5ad; no per-target cell count in obs"),
    ("perturbmatch", "K562", "somatic", "10x 3' CRISPRi", "limma-voom lfc matrix; SAME cells as k562_gwps (GSE146194)"),
    ("dld1_lowmoi", "DLD-1", "somatic", "scRNA CRISPRi", "GSE337988 lfc matrices, two low-MOI replicates"),
    ("igvf_hcasmc", "HCASMC", "somatic", "scRNA CRISPRi", "sceptre DE table, no cell counts"),
] + [(f"jiang_{ln}", ln, "somatic", "CRISPR-KO bulk-like DE", "per-line log2FC tables")
     for ln in ("A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7")]
inv276 = pd.read_csv(ROOT / "experiments/276-conservation-map/work/inventory276_sources.csv")
for key, line, state, chem, note in NO_CELLS:
    src = key
    path = ROOT / str(inv276.loc[inv276.key == key, "path"].iloc[0]) if (inv276.key == key).any() else None
    if path is not None and (path / "perts.txt").exists():
        perts = [l.strip() for l in (path / "perts.txt").read_text().splitlines() if l.strip()]
        tg = [p for p in perts if p in PANEL_SET]
    else:
        tg = []
    for t in tg:
        add(t, src, line, state, chem, "no cell count on disk", np.nan, np.nan, np.nan, np.nan,
            "not-on-disk", note)

# h1_2025: the 2025 VCC harness, exactly 400 cells per target by construction
h1 = json.loads((ROOT / "data/donors/h1_2025/meta.json").read_text())
perts = [l.strip() for l in (ROOT / "data/donors/h1_2025/perts.txt").read_text().splitlines() if l.strip()]
for t in [p for p in perts if p in PANEL_SET]:
    add(t, "h1_2025", "H1 hESC", "pluripotent", "10x Flex (VCC 2025 reference regime)", "cells",
        400.0, float(h1["n_control_cells"]), np.nan, np.nan,
        "data/donors/h1_2025/meta.json cells_per_pert 400/400/400", "2025 VCC reference: exactly 400 cells")

# ---------------------------------------------------------------- write
df = pd.DataFrame(ROWS)
# rows that measure the SAME cells as another row: never counted twice, kept for traceability
DUP = {"k562_gwps_raw": "k562_gwps (this experiment's own obs scan of the same h5ad; KA-1)",
       "perturbmatch": "k562_gwps (a reanalysis of GSE146194, the same cells)"}
df["duplicate_of"] = df.source.map(DUP).fillna("")
df = df.sort_values(["target", "cells"], ascending=[True, False], na_position="last")
df.to_csv(WORK / "depth278.csv", index=False)

meas = df[df.cells.notna() & (df.duplicate_of == "")].copy()
SOM = meas[meas.state == "somatic"]
DER = meas[meas.state.isin(["somatic", "derived"])]


def best(frame: pd.DataFrame) -> pd.Series:
    return frame.groupby("target").cells.max().reindex(PANEL, fill_value=0.0).fillna(0.0)


b_som, b_der, b_any = best(SOM), best(DER), best(meas)
summary = {
    "n_rows": int(len(df)), "n_measured_rows": int(len(meas)),
    "n_sources": int(df.source.nunique()),
    "panel_targets_with_any_measured_depth": int((b_any > 0).sum()),
    "thresholds": {},
    "per_source": {},
}
for name, b in [("somatic_lines_only", b_som), ("somatic_or_derived", b_der), ("any_line", b_any)]:
    summary["thresholds"][name] = {f"ge_{k}": int((b >= k).sum()) for k in (400, 300, 200, 150, 100, 50, 20)}
    summary["thresholds"][name]["median_deepest"] = float(np.median(b))
for s, g in meas.groupby("source"):
    summary["per_source"][s] = dict(
        line=g.line.iloc[0], state=g.state.iloc[0], n_panel=int(g.target.nunique()),
        median_cells=float(g.cells.median()), max_cells=float(g.cells.max()),
        ge400=int((g.cells >= 400).sum()), ge200=int((g.cells >= 200).sum()), ge100=int((g.cells >= 100).sum()))
for s, g in df[df.cells.isna()].groupby("source"):
    summary["per_source"][s] = dict(line=g.line.iloc[0], state=g.state.iloc[0],
                                    n_panel=int(g.target.nunique()), median_cells=None,
                                    note="no per-target cell count on disk")
(WORK / "depth278_summary.json").write_text(json.dumps(summary, indent=1))

# per-target deepest table
deep = pd.DataFrame({"target": PANEL})
for name, sub in [("somatic", SOM), ("somatic_or_derived", DER), ("any", meas)]:
    idx = sub.groupby("target").cells.idxmax()
    top = sub.loc[idx].set_index("target")
    deep[f"deepest_{name}_cells"] = deep.target.map(top.cells)
    deep[f"deepest_{name}_source"] = deep.target.map(top.source)
    deep[f"deepest_{name}_line"] = deep.target.map(top.line)
deep["n_sources_measured"] = deep.target.map(meas.groupby("target").source.nunique())
deep["n_somatic_lines_ge200"] = deep.target.map(
    SOM[SOM.cells >= 200].groupby("target").line.nunique())
deep = deep.fillna({"n_sources_measured": 0, "n_somatic_lines_ge200": 0})
deep.to_csv(WORK / "depth278_per_target.csv", index=False)

print(json.dumps(summary["thresholds"], indent=1))
print("rows", len(df), "sources", df.source.nunique(), "peak RSS GiB",
      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
