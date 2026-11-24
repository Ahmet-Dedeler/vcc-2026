"""276 item 1 -- INVENTORY.

For every compact-donor bank on disk (the `delta.npy` / `b_ctrl.npy` / `genes.txt` / `perts.txt`
object `scripts/build_donor_table.py` writes, which is the object production uses everywhere),
plus the two LINCS/CMap per-(target, line) signature banks 270 extracted, report:

  * the (dataset, line) identity and the perturbation type / unit
  * how many of the 300 panel targets it carries
  * how many of the challenge's 18,533 gene symbols it measures
  * the leak surface against the challenge lines (A=JURKAT, B=HELA, C=CAL33; AGENTS.md 2.1)
    and against our own PLOLO folds (K562, HCT116, HEK293T, KOLF2.1J)

Writes work/inventory276.json and work/inventory276.csv.  Read-only; nothing outside work/.

    PYTHONDONTWRITEBYTECODE=1 uv run python inv276.py
"""
from __future__ import annotations

import json
import resource
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/276-conservation-map"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

PANEL = [l.strip() for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
OFFICIAL = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
PANEL_SET, OFFICIAL_SET = set(PANEL), set(OFFICIAL)

# (key, path, dataset, cell line, assay, perturbation, unit of the effect row)
BANKS: list[tuple[str, str, str, str, str, str, str]] = [
    # --- production donors (the three PLOLO folds' own lines) ---
    ("k562_gwps", "data/donors/k562_gwps", "Replogle 2022 GWPS", "K562", "10x 3' Perturb-seq", "CRISPRi", "cells (267 gem groups)"),
    ("orion_hct116", "data/donors/orion_hct116", "X-Atlas/Orion", "HCT116", "Orion scRNA", "CRISPRi", "cells (batches)"),
    ("orion_hek293t", "data/donors/orion_hek293t", "X-Atlas/Orion", "HEK293T", "Orion scRNA", "CRISPRi", "cells (batches)"),
    # --- other single-cell genetic screens with a compact-donor table on disk ---
    ("k562_essential", "data/donors/k562_essential", "Replogle 2022 essential", "K562", "10x 3' Perturb-seq", "CRISPRi", "cells"),
    ("rpe1", "data/donors/rpe1", "Replogle 2022 RPE1", "RPE1", "10x 3' Perturb-seq", "CRISPRi", "cells"),
    ("nadig_jurkat", "data/donors/nadig_jurkat", "Nadig 2025", "Jurkat", "10x Perturb-seq", "CRISPRi", "cells (gem groups x guides)"),
    ("nadig_hepg2", "experiments/166-donor-acquisition/donors/nadig_hepg2", "Nadig 2025", "HepG2", "10x Perturb-seq", "CRISPRi", "cells (gem groups x guides)"),
    ("kolf2", "data/external/kolf2", "KOLF2.1J screen", "KOLF2.1J iPSC", "scRNA", "CRISPRi", "z-scores"),
    ("feng_ipsc", "experiments/099-donors/work/feng_ipsc", "Feng 2025 HipSci", "iPSC (34 donor lines pooled)", "10x", "CRISPRi", "cells (median 22/target)"),
    ("gse295214_h1", "experiments/115-gse295214-h1/work/compact-donor", "GSE295214", "H1 hESC", "10x", "CRISPRi", "cells"),
    ("h1_2025", "data/donors/h1_2025", "VCC 2025 H1", "H1 hESC", "10x Flex", "CRISPRi", "cells"),
    ("perturbmatch", "data/external/perturbmatch", "PerturbMatch K562 reanalysis", "K562", "10x 3'", "CRISPRi", "cells"),
    ("frangieh_melanoma", "experiments/166-donor-acquisition/donors/frangieh_melanoma", "Frangieh 2021", "melanoma (patient)", "10x", "CRISPR-KO", "cells"),
    ("igvf_telohaec", "data/external/igvf-telohaec", "IGVF", "TeloHAEC", "bulk/sc DE", "CRISPRi", "DE table"),
    ("igvf_cardiomyocyte", "experiments/225-source-sweep/work/tables/igvf_cardiomyocyte", "IGVF", "cardiomyocyte", "DE", "CRISPRi", "DE table"),
    ("myllia_a549", "experiments/120-myllia-a549/work/compact-donor", "Myllia", "A549", "scRNA", "CRISPRi", "cells (2 channels)"),
    ("altair_a549", "experiments/223-silent-zero-audit/work/tables/altair_a549", "Altair", "A549", "scRNA", "CRISPRi", "cells"),
    ("dld1_lowmoi", "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi", "DLD-1 low-MOI", "DLD-1", "scRNA", "CRISPRi", "lfc (2 replicates)"),
    ("dld1_low1", "experiments/223-silent-zero-audit/work/tables/dld1_low1", "DLD-1 low-MOI rep 1", "DLD-1", "scRNA", "CRISPRi", "lfc"),
    ("dld1_low2", "experiments/223-silent-zero-audit/work/tables/dld1_low2", "DLD-1 low-MOI rep 2", "DLD-1", "scRNA", "CRISPRi", "lfc"),
    ("sunshine_calu3", "experiments/225-source-sweep/work/tables/sunshine_calu3", "Sunshine 2023", "Calu-3", "10x", "CRISPR", "cells"),
    ("tian_ipsc_neuron", "experiments/225-source-sweep/work/tables/tian_ipsc_neuron", "Tian 2021", "iPSC-neuron", "10x", "CRISPRi", "cells"),
    ("tian_ipsc_neuron_crispra", "experiments/225-source-sweep/work/tables/tian_ipsc_neuron_crispra", "Tian 2021", "iPSC-neuron", "10x", "CRISPRa", "cells"),
    ("tian2019_ipsc", "experiments/225-source-sweep/work/tables/tian2019_ipsc", "Tian 2019", "iPSC", "10x", "CRISPRi", "cells"),
    ("xucao_hek293", "experiments/225-source-sweep/work/tables/xucao_hek293", "Xu 2023", "HEK293", "10x", "CRISPRi", "cells"),
    ("norman_k562", "experiments/225-source-sweep/work/tables/norman_k562", "Norman 2019", "K562", "10x", "CRISPRa", "cells"),
    ("gasperini_k562", "experiments/225-source-sweep/work/tables/gasperini_k562", "Gasperini 2019", "K562", "10x", "CRISPRi", "cells"),
    ("wessels_thp1", "experiments/225-source-sweep/work/tables/wessels_thp1", "Wessels", "THP-1", "10x", "CRISPR", "cells"),
    # --- Jiang 2025, per line (6 lines, unperturbed cytokine state) ---
    ("jiang_A549", "experiments/225-source-sweep/work/tables/jiang_A549", "Jiang 2025", "A549", "DE table", "CRISPR-KO", "log2FC"),
    ("jiang_BXPC3", "experiments/225-source-sweep/work/tables/jiang_BXPC3", "Jiang 2025", "BxPC-3", "DE table", "CRISPR-KO", "log2FC"),
    ("jiang_HAP1", "experiments/225-source-sweep/work/tables/jiang_HAP1", "Jiang 2025", "HAP1", "DE table", "CRISPR-KO", "log2FC"),
    ("jiang_HT29", "experiments/225-source-sweep/work/tables/jiang_HT29", "Jiang 2025", "HT-29", "DE table", "CRISPR-KO", "log2FC"),
    ("jiang_K562", "experiments/225-source-sweep/work/tables/jiang_K562", "Jiang 2025", "K562", "DE table", "CRISPR-KO", "log2FC"),
    ("jiang_MCF7", "experiments/225-source-sweep/work/tables/jiang_MCF7", "Jiang 2025", "MCF7", "DE table", "CRISPR-KO", "log2FC"),
    # --- CD4 T cells, three activation states (the incumbent's residual source) ---
    ("cd4_rest", "data/external/cd4t/cloud-compact/cd4_compact/cd4_rest", "GWCD4i", "CD4 T (rest)", "10x", "CRISPRi", "cells (pooled donors)"),
    ("cd4_stim8hr", "data/external/cd4t/cloud-compact/cd4_compact/cd4_stim8hr", "GWCD4i", "CD4 T (stim 8 h)", "10x", "CRISPRi", "cells (pooled donors)"),
    ("cd4_stim48hr", "data/external/cd4t/cloud-compact/cd4_compact/cd4_stim48hr", "GWCD4i", "CD4 T (stim 48 h)", "10x", "CRISPRi", "cells (pooled donors)"),
    # --- aggregates over several lines: inventoried, EXCLUDED from the map (not one line) ---
    ("jiang_6line", "experiments/225-source-sweep/work/tables/jiang_6line", "Jiang 2025 (6-line mean)", "AGGREGATE", "DE table", "CRISPR-KO", "log2FC"),
    ("jiang_5line_noK562", "experiments/225-source-sweep/work/tables/jiang_5line_noK562", "Jiang 2025 (5-line mean)", "AGGREGATE", "DE table", "CRISPR-KO", "log2FC"),
    ("cmap_sh_cgs", "experiments/223-silent-zero-audit/work/tables/cmap_sh_cgs", "LINCS/CMap (line-collapsed)", "AGGREGATE", "L1000", "shRNA", "z"),
    ("cmap_xpr", "experiments/223-silent-zero-audit/work/tables/cmap_xpr", "LINCS/CMap (line-collapsed)", "AGGREGATE", "L1000", "CRISPR-KO", "z"),
    ("cmap_sh_cgs_lm", "experiments/223-silent-zero-audit/work/tables/cmap_sh_cgs_lm", "LINCS/CMap landmark (line-collapsed)", "AGGREGATE", "L1000", "shRNA", "z"),
    ("hspc_crispri500", "data/external/hspc_crispri500", "mouse HSPC CRISPRi", "MOUSE", "DE", "CRISPRi", "DE table"),
]

# the lines that matter for leakage
CHALLENGE_LINES = {"JURKAT", "HELA", "CAL33"}      # A / B / C, exp 041 CCLE match (AGENTS.md 2.1)
FOLD_LINES = {"K562", "HCT116", "HEK293T", "KOLF2.1J"}  # our own PLOLO folds -- also SOURCES here


def scan_bank(path: Path) -> dict | None:
    if not (path / "delta.npy").exists():
        return None
    perts = [p.strip() for p in (path / "perts.txt").read_text().splitlines() if p.strip()]
    genes = [g.strip() for g in (path / "genes.txt").read_text().splitlines() if g.strip()]
    d = np.load(path / "delta.npy", mmap_mode="r")
    meta = {}
    if (path / "meta.json").exists():
        try:
            meta = json.loads((path / "meta.json").read_text())
        except Exception:
            meta = {}
    hit = sorted(PANEL_SET.intersection(perts))
    return dict(
        n_perts=len(perts), n_genes=len(genes), shape=list(d.shape),
        n_panel=len(hit), panel=hit,
        n_genes_official=len(OFFICIAL_SET.intersection(genes)),
        cells_per_pert=meta.get("cells_per_pert"), source=meta.get("source"),
    )


def main() -> None:
    rows = []
    for key, rel, dataset, line, assay, pert, unit in BANKS:
        p = ROOT / rel
        s = scan_bank(p)
        if s is None:
            rows.append(dict(key=key, path=rel, dataset=dataset, line=line, assay=assay,
                             perturbation=pert, unit=unit, present=False))
            print(f"MISSING {key:26s} {rel}")
            continue
        up = line.upper().replace("-", "").replace(" ", "")
        leak_challenge = any(c in up for c in CHALLENGE_LINES)
        leak_fold = any(c.replace(".", "") in up.replace(".", "") for c in FOLD_LINES)
        rows.append(dict(key=key, path=rel, dataset=dataset, line=line, assay=assay,
                         perturbation=pert, unit=unit, present=True,
                         leak_challenge_line=bool(leak_challenge), is_fold_line=bool(leak_fold),
                         **{k: v for k, v in s.items() if k != "panel"}))
        rows[-1]["panel"] = s["panel"]
        print(f"{key:26s} perts {s['n_perts']:>6d}  genes {s['n_genes']:>6d} "
              f"(official {s['n_genes_official']:>5d})  PANEL {s['n_panel']:>3d}   {line}")

    # ---- CMap, per (target, line) ---------------------------------------------------------
    cmap = {}
    for tag, f in (("cmap_sh", "cmap270_trt_sh_cgs.npz"), ("cmap_xpr", "cmap270_trt_xpr.npz")):
        z = np.load(ROOT / "experiments/270-data-reconnaissance/work" / f, allow_pickle=True)
        tgt = np.asarray(z["target"], dtype=object)
        lin = np.asarray(z["line"], dtype=object)
        lines = sorted(set(lin))
        per_line = {L: sorted(set(tgt[lin == L])) for L in lines}
        cmap[tag] = dict(n_sigs=int(len(tgt)), n_targets=len(set(tgt)), n_lines=len(lines),
                         lines=lines, n_genes=int(z["Z"].shape[1]),
                         targets_per_line={L: len(v) for L, v in per_line.items()})
        print(f"{tag:26s} sigs {len(tgt):>6d}  targets {len(set(tgt)):>3d}  lines {len(lines):>3d}")

    # ---- targets x number of independent lines / datasets ---------------------------------
    #  a source counts if it is a single line (not AGGREGATE / MOUSE) and carries the target
    per_target_lines: dict[str, set[str]] = {t: set() for t in PANEL}
    per_target_sets: dict[str, set[str]] = {t: set() for t in PANEL}
    for r in rows:
        if not r.get("present") or r["line"] in ("AGGREGATE", "MOUSE"):
            continue
        for t in r["panel"]:
            per_target_lines[t].add(r["line"])
            per_target_sets[t].add(r["key"])
    cmap_lines: dict[str, set[str]] = {t: set() for t in PANEL}
    for tag, f in (("cmap_sh", "cmap270_trt_sh_cgs.npz"), ("cmap_xpr", "cmap270_trt_xpr.npz")):
        z = np.load(ROOT / "experiments/270-data-reconnaissance/work" / f, allow_pickle=True)
        tgt = np.asarray(z["target"], dtype=object); lin = np.asarray(z["line"], dtype=object)
        for t, L in zip(tgt, lin):
            if t in cmap_lines:
                cmap_lines[t].add(L)

    tbl = []
    for t in PANEL:
        sc_lines = per_target_lines[t]
        tbl.append(dict(target=t, n_sc_sources=len(per_target_sets[t]), n_sc_lines=len(sc_lines),
                        n_cmap_lines=len(cmap_lines[t]),
                        n_all_lines=len(sc_lines | cmap_lines[t]),
                        sc_sources=";".join(sorted(per_target_sets[t])),
                        cmap_lines=";".join(sorted(cmap_lines[t]))))

    import csv
    with (WORK / "inventory276_targets.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tbl[0]))
        w.writeheader(); w.writerows(tbl)
    with (WORK / "inventory276_sources.csv").open("w", newline="") as fh:
        cols = ["key", "path", "dataset", "line", "assay", "perturbation", "unit", "present",
                "n_perts", "n_genes", "n_genes_official", "n_panel", "leak_challenge_line", "is_fold_line"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    (WORK / "inventory276.json").write_text(json.dumps(
        dict(panel=len(PANEL), official_genes=len(OFFICIAL), sources=rows, cmap=cmap),
        indent=1, default=str))

    n_sc = np.array([r["n_sc_lines"] for r in tbl])
    n_all = np.array([r["n_all_lines"] for r in tbl])
    print("\n--- targets by number of independent LINES ---")
    for lo in (0, 1, 2, 3, 5, 8, 12):
        print(f"  >= {lo:2d} single-cell/DE lines: {int((n_sc >= lo).sum()):3d}      "
              f">= {lo:2d} lines incl. CMap: {int((n_all >= lo).sum()):3d}")
    print(f"  median lines per target: sc {np.median(n_sc):.1f}  all {np.median(n_all):.1f}")
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.2f} GiB")


if __name__ == "__main__":
    main()
