"""292 item 1+3 -- the inventory over EVERY perturbation dataset on disk, USED / UNUSED with the
evidence, and the ranked reserve.

Pure bookkeeping over metadata artefacts already written by 276 / 278 / 292 -- no h5ad is opened.

Output: work/inventory292.csv, work/sha256.txt, work/reserve292.json
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/292-reserve-declaration/work"

scan = {r["file"][:-5]: r for r in json.loads((WORK / "scan292.json").read_text())}
perts = json.loads((WORK / "perts292.json").read_text())
src276 = pd.read_csv(ROOT / "experiments/276-conservation-map/work/inventory276_sources.csv")
MAP28 = set(json.loads((ROOT / "experiments/279-znf-state-dependence/work/build279.json")
                       .read_text())["sources"])
MAP_LINES = set(json.loads((ROOT / "experiments/279-znf-state-dependence/work/build279.json")
                           .read_text())["lines"])

# --------------------------------------------------------------- sha256
SHA: dict[str, str] = {}
for f in ("sha256_scperturb.txt", "sha256_other.txt"):
    for line in (WORK / f).read_text().splitlines():
        m = re.match(r"SHA2-256\((.+)\)= ([0-9a-f]{64})", line.strip())
        if m:
            p = m.group(1)
            p = p[len(str(ROOT)) + 1:] if p.startswith(str(ROOT)) else p
            SHA[p] = m.group(2)

# --------------------------------------------------------------- scPerturb rows
# every scPerturb stem that FEEDS a bank used by the programme, with the code that reads it
USED_SCP = {
    "ReplogleWeissman2022_K562_gwps": "bank k562_gwps (map) -- scripts/build_donor_table.py:40",
    "ReplogleWeissman2022_K562_essential": "bank k562_essential (map) -- scripts/build_donor_table.py:46",
    "ReplogleWeissman2022_rpe1": "bank rpe1 (map) -- scripts/build_donor_table.py:52",
    "FrangiehIzar2021_RNA": "bank frangieh_melanoma (map) -- E166/build_extra.py",
    "SunshineHein2023": "bank sunshine_calu3 (map) -- E225/build225.py:63",
    "TianKampmann2021_CRISPRi": "bank tian_ipsc_neuron (map) -- E225/build225.py:60",
    "TianKampmann2021_CRISPRa": "bank tian_ipsc_neuron_crispra (276 inventory) -- E225/build225.py:61",
    "TianKampmann2019_iPSC": "bank tian2019_ipsc (map) -- E225/build225.py:67",
    "XuCao2023": "bank xucao_hek293 (map) -- E225/build225.py:62",
    "NormanWeissman2019_filtered": "bank norman_k562 (276 inventory) -- E225/build225.py:64",
    "GasperiniShendure2019_lowMOI": "bank gasperini_k562 (map) -- E225/build225.py:66",
    "WesselsSatija2023": "bank wessels_thp1 (map) -- E225/build225.py:65",
}

# UNUSED classification: why a set is or is not reserve-eligible
TIER = {
    "AdamsonWeissman2016_GSM2406681_10X010": ("R3-core", "K562 CRISPRi, 82 genes >=50 cells, new study"),
    "DatlingerBock2017": ("R3-core", "Jurkat CROP-seq Cas9 KO, 3 guides/gene"),
    "DatlingerBock2021": ("R3-core", "Jurkat scifi Cas9 KO, 1.7k cells/gene, 384 samples"),
    "DixitRegev2016": ("R3-core", "K562 Perturb-seq, 2.5 guides/gene, deep"),
    "PapalexiSatija2021_eccite_RNA": ("R3-core", "THP-1 ECCITE Cas9 KO, 4 guides/gene"),
    "PapalexiSatija2021_eccite_arrayed_RNA": ("R3-core", "THP-1 arrayed (one gene per well)"),
    "ShifrutMarson2018": ("R3-core", "PRIMARY human CD8 T cells, Cas9 RNP, 2 donors -- the only new context"),
    "SchraivogelSteinmetz2020_TAP_SCREEN__chromosome_11_screen":
        ("excluded", "K562 CRISPRi but TARGETED readout: 3,185 measured genes cover only 30 of the "
                     "map's 2,739 response-axis genes -- no comparable response vector"),
    "SchraivogelSteinmetz2020_TAP_SCREEN__chromosome_8_screen":
        ("excluded", "K562 CRISPRi but TARGETED readout: 4,191 measured genes cover only 28 of the "
                     "map's 2,739 response-axis genes -- no comparable response vector"),
    "TianKampmann2019_day7neuron":
        ("R3-restricted", "iPSC-neuron d7 CRISPRi -- SAME STUDY as the map bank tian2019_ipsc"),
    "AdamsonWeissman2016_GSM2406675_10X001": ("excluded", "7 genes, pilot lane of the same study"),
    "AdamsonWeissman2016_GSM2406677_10X005": ("excluded", "median 3 cells per gene (3-guide combos)"),
    "FrangiehIzar2021_protein": ("excluded", "SAME CELLS as the used RNA bank; 24 protein features"),
    "PapalexiSatija2021_eccite_protein": ("excluded", "4 protein features, not a transcriptome"),
    "PapalexiSatija2021_eccite_arrayed_protein": ("excluded", "4 protein features"),
    "GasperiniShendure2019_atscale": ("excluded", "high-MOI enhancer screen; label = guide combination"),
    "GasperiniShendure2019_highMOI": ("excluded", "high-MOI enhancer screen; label = guide combination"),
    "XieHon2017": ("excluded", "enhancer coordinates, no gene targets"),
    "JoungZhang2023_atlas": ("excluded-recoverable",
                             "1.15 M hESC cells, 3,369 TF ORFs -- ids NOT resolvable from on-disk "
                             "metadata; modality is ORF OVEREXPRESSION, not knockdown"),
    "JoungZhang2023_combinatorial": ("excluded", "combinatorial ORF pairs only"),
    "McFarlandTsherniak2020": ("excluded", "1 gene (GPX4) + drugs; 209 lines but chemical"),
    "CuiHacohen2023": ("excluded", "MOUSE, cytokine treatment (symbols match by accident)"),
    "GehringPachter2019": ("excluded", "MOUSE, BMP4 treatment"),
    "LaraAstiasoHuntly2023_exvivo": ("excluded", "MOUSE"),
    "LaraAstiasoHuntly2023_invivo": ("excluded", "MOUSE"),
    "LaraAstiasoHuntly2023_leukemia": ("excluded", "MOUSE"),
    "LiangWang2023": ("excluded", "MOUSE organoid"),
    "SantinhaPlatt2023": ("excluded", "MOUSE in vivo"),
    "SchiebingerLander2019_GSE106340": ("excluded", "MOUSE, chemical reprogramming"),
    "SchiebingerLander2019_GSE115943": ("excluded", "MOUSE, chemical reprogramming"),
    "WeinrebKlein2020": ("excluded", "MOUSE, cytokine"),
    "AissaBenevolenskaya2021": ("excluded", "drug"),
    "ChangYe2021": ("excluded", "drug"),
    "LotfollahiTheis2023": ("excluded", "drug combinations"),
    "SrivatsanTrapnell2020_sciplex2": ("excluded", "drug"),
    "SrivatsanTrapnell2020_sciplex3": ("excluded", "drug"),
    "SrivatsanTrapnell2020_sciplex4": ("excluded", "drug"),
    "ZhaoSims2021": ("excluded", "drug, primary glioblastoma"),
}

rows = []
for stem in sorted(scan):
    s, p = scan[stem], perts[stem]
    d = s.get("desc", {})

    def top1(k):
        v = d.get(k)
        return v["top"][0][0] if v and v["top"] else ""

    used = stem in USED_SCP
    tier, note = TIER.get(stem, ("", ""))
    rows.append(dict(
        collection="scperturb", key=stem, path=s["path"], size_bytes=s["size_bytes"],
        sha256=SHA.get(s["path"], ""),
        organism=top1("organism"), line=top1("cell_line") or top1("celltype"),
        modality=top1("perturbation_type"), n_cells=s.get("n_cells"),
        n_var=s.get("n_var"), n_labels=s.get("n_labels"),
        n_gene_perts=p["n_gene_perts"], n_ge50=p["n_ge50"],
        n_in_held_map=p["n_in_held"], n_in_panel=p["n_in_panel"],
        median_cells_per_pert=p["median_cells"], q10_cells_per_pert=p["q10_cells"],
        control_cells=p["control_cells"], median_guides_per_gene=p["median_guides_per_gene"],
        batch_col=s.get("batch_col"), n_batches=s.get("n_batchs"),
        guide_col=s.get("guide_col"), n_guides=s.get("n_guides"),
        status="USED" if used else "UNUSED",
        evidence=USED_SCP.get(stem, "no code reads it; only obs/var censuses (E166/shared.py, "
                                    "E218/census166_copy.py, E225/census225.py, E270 obs scan, "
                                    "E278/scan278_h5ad.py, this experiment)"),
        reserve_tier=tier, reserve_note=note, grammar=p["grammar"]))

# --------------------------------------------------------------- non-scPerturb rows
for r in src276.itertuples():
    key = r.key
    in_map = key in MAP28
    rows.append(dict(
        collection="donor-bank", key=key, path=r.path, size_bytes="", sha256="",
        organism="mouse" if key == "hspc_crispri500" else "human",
        line=r.line, modality=r.perturbation, n_cells="", n_var=r.n_genes,
        n_labels=r.n_perts, n_gene_perts=r.n_perts, n_ge50="",
        n_in_held_map="", n_in_panel=r.n_panel,
        median_cells_per_pert="", q10_cells_per_pert="", control_cells="",
        median_guides_per_gene="", batch_col="", n_batches="", guide_col="", n_guides="",
        status="UNUSED" if key == "hspc_crispri500" else "USED",
        evidence=("276 inventory row, present=False (mouse; no bank built)"
                  if key == "hspc_crispri500" else
                  ("one of the 28 conservation-map banks (E279/build279.json)" if in_map
                   else "E276/work/inventory276_sources.csv (map inventory / production donor)")),
        reserve_tier="excluded" if key == "hspc_crispri500" else "",
        reserve_note="MOUSE HSPC DE table" if key == "hspc_crispri500" else "", grammar=""))

OTHER = [
    ("sealed", "hcasmc_igvf", "data/external/hcasmc-igvf/IGVFFI0206LUDV.tsv.gz", "human",
     "HCASMC (coronary artery smooth muscle)", "CRISPRi", 27,
     "UNUSED", "SEALED by 126 / docs/DEF-RUNBOOK.md:294-320 -- no decompression, no header read, "
               "no value-dependent statistic without a new hash-bound authorisation",
     "sealed-elsewhere", "already a reserve under an older, competition-track rule; NOT part of R3"),
    ("sealed", "hcasmc_igvf_aux", "data/external/hcasmc-igvf/IGVFFI9546UECF.tsv.gz", "human",
     "HCASMC", "CRISPRi", "", "UNUSED", "sibling of the sealed file (76,842 B)",
     "sealed-elsewhere", "not opened"),
    ("no-data", "gse249595", "data/external/gse249595", "human", "-", "CRISPR (guide library)", 0,
     "UNUSED", "guide FASTA + one supplementary xlsx + a GEO file listing; no expression matrix on disk",
     "excluded", "no data to score"),
    ("reference", "cmap", "data/external/cmap", "human", "38 lines", "shRNA / CRISPR-KO (L1000)",
     148, "USED", "270 computed SS_cmap features; conservation276.csv carries C_cmap_sh / C_cmap_xpr",
     "", "outcomes inspected"),
    ("reference", "broad16line", "data/external/broad16line", "human", "16 cancer lines",
     "Cas9 KO", 0, "USED",
     "268 / 283 / 289 / 291 built and scored responses; the four reserved lines are SPENT (C2)",
     "", "outcomes inspected"),
    ("reference", "tahoe", "data/external/tahoe", "human", "50 lines", "drug", 0, "USED",
     "270 §4 audit; chemical only -- cannot carry a genetic-perturbation claim", "", ""),
    ("reference", "cd4_cell2026", "data/external/cd4_cell2026", "human", "CD4 T", "CRISPRi", 300,
     "USED", "287 and 291 (claim C4) read K562_comparison.suppl_table.csv", "", ""),
    ("reference", "vcc2025", "data/external/vcc2025", "human", "H1", "CRISPRi", 20, "USED",
     "bank h1_2025 (map)", "", ""),
    ("reference", "orion", "data/external/orion", "human", "HCT116 / HEK293T", "CRISPRi", 294,
     "USED", "banks orion_hct116 / orion_hek293t (map); HEK293T half-B reserve C1 SPENT", "", ""),
    ("annotation", "depmap26q1 / ccle_broad_2025 / regulon / omnipath-dorothea / slim / sigrecon",
     "data/external/{depmap26q1,ccle_broad_2025,regulon,omnipath-dorothea,slim,sigrecon}", "human",
     "-", "annotation / embeddings / ranked lists", "", "USED",
     "feature sources for 281 / 284 / 290; not perturbation measurements", "", ""),
    ("model", "state-pretrained / state-vcc / xcell / perturbation-decomposition",
     "data/external/{state-pretrained,state-vcc,xcell,perturbation-decomposition}", "-", "-",
     "model checkpoints / code", "", "USED", "baselines and method references, not data", "", ""),
]
for coll, key, path, org, line, mod, npanel, status, ev, tier, note in OTHER:
    rows.append(dict(collection=coll, key=key, path=path,
                     size_bytes=(ROOT / path).stat().st_size if (ROOT / path).is_file() else "",
                     sha256=SHA.get(path, ""), organism=org, line=line, modality=mod,
                     n_cells="", n_var="", n_labels="", n_gene_perts="", n_ge50="",
                     n_in_held_map="", n_in_panel=npanel, median_cells_per_pert="",
                     q10_cells_per_pert="", control_cells="", median_guides_per_gene="",
                     batch_col="", n_batches="", guide_col="", n_guides="", status=status,
                     evidence=ev, reserve_tier=tier, reserve_note=note, grammar=""))

cols = list(rows[0].keys())
with (WORK / "inventory292.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow(r)

with (WORK / "sha256.txt").open("w") as f:
    for p, h in sorted(SHA.items()):
        f.write(f"{h}  {p}\n")

n_used = sum(1 for r in rows if r["status"] == "USED")
print(f"inventory292.csv: {len(rows)} rows, USED {n_used}, UNUSED {len(rows) - n_used}")
print(f"sha256.txt: {len(SHA)} digests")
core = [k for k, (t, _) in TIER.items() if t == "R3-core"]
restr = [k for k, (t, _) in TIER.items() if t == "R3-restricted"]
print("R3 core:", len(core), "restricted:", len(restr))
res = {"core": core, "restricted": restr,
       "files": {k: {"path": scan[k]["path"], "sha256": SHA[scan[k]["path"]],
                     "size_bytes": scan[k]["size_bytes"]} for k in core + restr}}
(WORK / "reserve292.json").write_text(json.dumps(res, indent=1))
