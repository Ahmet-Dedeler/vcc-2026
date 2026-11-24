"""811 step 1 -- inventory: which public banks carry which fold panels' targets (A/B/C and the lab's WORKING panels L1-L6).

Reads only target lists (perts.txt / meta json / the CD4 h5ad's target categories / 813's TeloHAEC coverage record); never
an expression matrix. Never reads the sealed panels (L7, L8). Writes work/inventory811.json and prints a table.

    uv run python experiments/811-lolo-model/inv811.py
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work" / "inventory811.json"
PANELS = {"ABC": ROOT / "data/challenge/pert_counts.csv",
          **{f"L{k}": ROOT / f"experiments/801-lab-folds/work/L{k}/bundle/pert_counts.csv" for k in range(1, 7)}}
T223 = ROOT / "experiments/223-silent-zero-audit/work/tables"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
# name: (target-list source, line, lab / assay note, independent-line group)
BANKS = {
    "k562_gwps": (ROOT / "data/donors/k562_gwps/perts.txt", "K562", "Weissman 2022 GWPS, CRISPRi 10x 3'", "K562"),
    "viperturb_k562": (ROOT / "data/donors/viperturb_k562/perts.txt", "K562", "VIPerturb-seq 2026, CRISPRi Flex, ~48 cells", "K562"),
    "perturbmatch_k562": (ROOT / "data/external/perturbmatch/perts.txt", "K562", "GWPS re-analysis (same cells as k562_gwps)", "K562"),
    "orion_hct116": (ROOT / "data/donors/orion_hct116/perts.txt", "HCT116", "Xaira X-Atlas/Orion, CRISPRi Fix-Cryo", "HCT116"),
    "orion_hek293t": (ROOT / "data/donors/orion_hek293t/perts.txt", "HEK293T", "Xaira X-Atlas/Orion, CRISPRi Fix-Cryo", "HEK293T"),
    "cd4_gwcd4i": ("CD4H5AD", "CD4 T (primary)", "Marson GWCD4i, CRISPRi Flex, 3 conditions", "CD4"),
    "h1_pan2025": (ROOT / "experiments/115-gse295214-h1/work/compact-donor/perts.txt", "H1 hESC", "Pan 2025 GSE295214 CRISPRi", "H1"),
    "h1_vcc2025": (ROOT / "data/donors/h1_2025_full/perts.txt", "H1 hESC", "Arc VCC 2025, CRISPRi Flex (300 targets)", "H1"),
    "kolf2_local": (ROOT / "data/external/kolf2/perts.txt", "KOLF2.1J iPSC", "Mali atlas z-scores, local A/B/C rows only", "KOLF2"),
    "kolf2_official": ("KOLF2META", "KOLF2.1J iPSC", "Mali atlas z-scores, full official list (not on disk as rows)", "KOLF2"),
    "telohaec_pilot": (ROOT / "data/external/igvf-telohaec/perts.txt", "TeloHAEC", "IGVF pilot IGVFFI8645YMGX, CRISPRi", "TeloHAEC"),
    "telohaec_gw": ("TELO813", "TeloHAEC", "IGVF genome-wide 2026-10-12 (813; not on disk)", "TeloHAEC"),
    "dld1_lowmoi": (T223 / "dld1_lowmoi/perts.txt", "DLD-1", "GSE337988 author LFC (low MOI)", "DLD1"),
    "igvf_cardiomyocyte": (T225 / "igvf_cardiomyocyte/perts.txt", "iPSC-cardiomyocyte", "IGVF flat DE", "CM"),
    "altair_a549": (T223 / "altair_a549/perts.txt", "A549", "Altair optical-pooled scRNA", "A549"),
    "cmap_xpr": (T223 / "cmap_xpr/perts.txt", "many (L1000)", "LINCS L1000 CRISPR KO consensus (inferred genes)", "CMAP"),
    "cmap_sh_cgs": (T223 / "cmap_sh_cgs/perts.txt", "many (L1000)", "LINCS L1000 shRNA CGS (inferred genes)", "CMAP"),
    "k562_essential": (ROOT / "data/donors/k562_essential/perts.txt", "K562", "Weissman essential", "K562"),
    "rpe1": (ROOT / "data/donors/rpe1/perts.txt", "RPE1", "Weissman essential", "RPE1"),
    "nadig_jurkat": (ROOT / "data/donors/nadig_jurkat/perts.txt", "Jurkat", "Nadig 2025 essential", "Jurkat"),
}
for j in sorted(T225.glob("jiang_*")):
    if j.name.startswith("jiang_") and j.name not in ("jiang_6line", "jiang_5line_noK562"):
        BANKS[j.name] = (j / "perts.txt", j.name.split("_", 1)[1], "Jiang 2025 pathway screens", j.name.split("_", 1)[1])


def panel(p: Path) -> list[str]:
    with open(p) as f:
        return [r["target_gene"] for r in csv.DictReader(f)]


def targets_of(src) -> set[str]:
    if src == "CD4H5AD":
        with h5py.File(ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad", "r") as h:
            node = h["obs"]["target_contrast_gene_name"]
            vals = node["categories"][:] if isinstance(node, h5py.Group) else node[:]
            return {v.decode() if isinstance(v, bytes) else str(v) for v in vals}
    if src == "KOLF2META":
        return set(json.loads((ROOT / "data/external/kolf2/official/meta.json").read_text())["perts"])
    if src == "TELO813":
        return set()
    return {ln.strip() for ln in Path(src).read_text().splitlines() if ln.strip()}


def main():
    pans = {k: panel(p) for k, p in PANELS.items()}
    telo = json.loads((ROOT / "experiments/813-public-data/work/meta/telohaec_rep1_coverage.json").read_text())
    rows = {}
    for name, (src, line, note, grp) in BANKS.items():
        T = targets_of(src)
        cov = {k: len(set(v) & T) for k, v in pans.items()}
        if src == "TELO813":
            cov = {"ABC": telo["panel_in_library"], **{f"L{k}": telo[f"L{k}_in_library"] for k in range(1, 7)}}
            n = telo["n_target_genes"]
        else:
            n = len(T)
        rows[name] = dict(line=line, group=grp, note=note, n_targets=n, coverage=cov,
                          source=str(src if isinstance(src, str) else Path(src).relative_to(ROOT)),
                          sha256=(hashlib.sha256(Path(src).read_bytes()).hexdigest() if isinstance(src, Path) else None))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dict(panels={k: hashlib.sha256("\n".join(sorted(v)).encode()).hexdigest() for k, v in pans.items()},
                                   banks=rows), indent=1))
    hdr = f"{'bank':20s} {'line':18s} {'n_targets':>9s} " + " ".join(f"{k:>4s}" for k in pans)
    print(hdr)
    for name, r in rows.items():
        print(f"{name:20s} {r['line'][:18]:18s} {r['n_targets']:9d} " + " ".join(f"{r['coverage'][k]:4d}" for k in pans))


if __name__ == "__main__":
    main()
