"""276 — shared loading for the conservation map.  Read-only outside work/.

Builds, for every included source, the panel-target rows on the FROZEN common gene axis, in the
residual space of `work/DEFINITION-276.md` §2 (each source centred on its own across-target mean).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/276-conservation-map"
WORK = HERE / "work"

PANEL = [l.strip() for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
OFFICIAL = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
PANEL_IDX = {t: i for i, t in enumerate(PANEL)}

# key -> (relative path, cell line, dataset label).  Single-line sources with >= 5 panel targets,
# CRISPRi/CRISPR-KO only (DEFINITION-276 §1).
SOURCES: dict[str, tuple[str, str, str]] = {
    "k562_gwps":          ("data/donors/k562_gwps", "K562", "Replogle22-GWPS"),
    "perturbmatch":       ("data/external/perturbmatch", "K562", "PerturbMatch"),
    "jiang_K562":         ("experiments/225-source-sweep/work/tables/jiang_K562", "K562", "Jiang25"),
    "orion_hct116":       ("data/donors/orion_hct116", "HCT116", "X-Atlas/Orion"),
    "orion_hek293t":      ("data/donors/orion_hek293t", "HEK293T", "X-Atlas/Orion"),
    "xucao_hek293":       ("experiments/225-source-sweep/work/tables/xucao_hek293", "HEK293", "Xu23"),
    "kolf2":              ("data/external/kolf2", "KOLF2.1J", "KOLF2.1J screen"),
    "feng_ipsc":          ("experiments/099-donors/work/feng_ipsc", "iPSC-HipSci", "Feng25"),
    "tian_ipsc_neuron":   ("experiments/225-source-sweep/work/tables/tian_ipsc_neuron", "iPSC-neuron", "Tian21"),
    "gse295214_h1":       ("experiments/115-gse295214-h1/work/compact-donor", "H1", "GSE295214"),
    "h1_2025":            ("data/donors/h1_2025", "H1", "VCC2025-H1"),
    "frangieh_melanoma":  ("experiments/166-donor-acquisition/donors/frangieh_melanoma", "melanoma", "Frangieh21"),
    "igvf_telohaec":      ("data/external/igvf-telohaec", "TeloHAEC", "IGVF"),
    "igvf_cardiomyocyte": ("experiments/225-source-sweep/work/tables/igvf_cardiomyocyte", "cardiomyocyte", "IGVF"),
    "myllia_a549":        ("experiments/120-myllia-a549/work/compact-donor", "A549", "Myllia"),
    "altair_a549":        ("experiments/223-silent-zero-audit/work/tables/altair_a549", "A549", "Altair"),
    "jiang_A549":         ("experiments/225-source-sweep/work/tables/jiang_A549", "A549", "Jiang25"),
    "dld1_lowmoi":        ("experiments/223-silent-zero-audit/work/tables/dld1_lowmoi", "DLD-1", "DLD-1 low-MOI"),
    "sunshine_calu3":     ("experiments/225-source-sweep/work/tables/sunshine_calu3", "Calu-3", "Sunshine23"),
    "jiang_BXPC3":        ("experiments/225-source-sweep/work/tables/jiang_BXPC3", "BxPC-3", "Jiang25"),
    "jiang_HAP1":         ("experiments/225-source-sweep/work/tables/jiang_HAP1", "HAP1", "Jiang25"),
    "jiang_HT29":         ("experiments/225-source-sweep/work/tables/jiang_HT29", "HT-29", "Jiang25"),
    "jiang_MCF7":         ("experiments/225-source-sweep/work/tables/jiang_MCF7", "MCF7", "Jiang25"),
    "cd4_rest":           ("data/external/cd4t/cloud-compact/cd4_compact/cd4_rest", "CD4T", "GWCD4i rest"),
    "cd4_stim8hr":        ("data/external/cd4t/cloud-compact/cd4_compact/cd4_stim8hr", "CD4T", "GWCD4i stim8h"),
    "cd4_stim48hr":       ("data/external/cd4t/cloud-compact/cd4_compact/cd4_stim48hr", "CD4T", "GWCD4i stim48h"),
}

# the source(s) that must be dropped when a given line is HELD (leave-the-held-line-out)
HELD_LINE_SOURCES: dict[str, tuple[str, ...]] = {
    "K562": ("k562_gwps", "perturbmatch", "jiang_K562"),
    "HCT116": ("orion_hct116",),
    "HEK293T": ("orion_hek293t", "xucao_hek293"),   # HEK293 is the parental line of HEK293T
    "KOLF2.1J": ("kolf2", "feng_ipsc"),             # Feng's HipSci pool contains kolf_2
    "H1": ("gse295214_h1", "h1_2025"),
}
CMAP_HELD_LINE: dict[str, tuple[str, ...]] = {
    "K562": ("K562",), "HCT116": ("HCT116",), "HEK293T": ("HEK293T", "HEK293"),
    "KOLF2.1J": (), "H1": (),
}


def read_bank(rel: str) -> tuple[list[str], list[str], Path]:
    p = ROOT / rel
    perts = [x.strip() for x in (p / "perts.txt").read_text().splitlines() if x.strip()]
    genes = [x.strip() for x in (p / "genes.txt").read_text().splitlines() if x.strip()]
    return perts, genes, p


def common_axis_paths(rels: list[str]) -> list[str]:
    axis = set(OFFICIAL)
    for rel in rels:
        _, genes, _ = read_bank(rel)
        axis &= set(genes)
    # deterministic order: the challenge's own gene order
    return [g for g in OFFICIAL if g in axis]


def common_axis(keys: list[str]) -> list[str]:
    return common_axis_paths([SOURCES[k][0] for k in keys])


def load_source_path(rel: str, axis: list[str],
                     targets: list[str] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (rows, avail): rows (n_targets, |axis|) float32 residual, avail (n_targets,) bool.

    The residual is the source's row minus the mean over exactly the rows it contributes to
    THIS map (DEFINITION-276 s2)."""
    tg = PANEL if targets is None else targets
    tidx = {t: i for i, t in enumerate(tg)}
    perts, genes, p = read_bank(rel)
    gidx = {g: i for i, g in enumerate(genes)}
    cols = np.asarray([gidx[g] for g in axis], dtype=np.int64)
    pidx = {t: i for i, t in enumerate(perts)}
    have = [(t, pidx[t]) for t in tg if t in pidx]
    out = np.zeros((len(tg), len(axis)), dtype=np.float32)
    avail = np.zeros(len(tg), dtype=bool)
    if not have:
        return out, avail
    D = np.load(p / "delta.npy", mmap_mode="r")
    order = np.argsort([r for _, r in have])
    rows_sorted = np.asarray([have[i][1] for i in order], dtype=np.int64)
    blk = np.asarray(D[rows_sorted][:, cols], dtype=np.float64)
    blk -= blk.mean(axis=0, keepdims=True)          # centre on this source's across-target mean
    for i, oi in enumerate(order):
        t = have[oi][0]
        out[tidx[t]] = blk[i].astype(np.float32)
        avail[tidx[t]] = True
    return out, avail


def load_source(key: str, axis: list[str],
                targets: list[str] | None = None) -> tuple[np.ndarray, np.ndarray]:
    return load_source_path(SOURCES[key][0], axis, targets)


def pair_list(keys: list[str], cross_line_only: bool = True) -> list[tuple[int, int]]:
    lines = [SOURCES[k][1] for k in keys]
    out = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            same = lines[i] == lines[j]
            if cross_line_only and same:
                continue
            if (not cross_line_only) and (not same):
                continue
            out.append((i, j))
    return out


def cosine_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    d = na * nb
    out = np.full(A.shape[0], np.nan)
    ok = d > 0
    out[ok] = (A[ok] * B[ok]).sum(axis=1) / d[ok]
    return out
