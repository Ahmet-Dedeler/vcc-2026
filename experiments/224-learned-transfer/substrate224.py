"""224 Phase 0 — the leave-one-LINE-out substrate, and its census.

Enumerates every donor-shaped table on disk (`delta.npy` + `perts.txt` + `genes.txt`),
records provenance and value semantics, and computes the pairwise shared-target matrix
that decides which surfaces a leave-one-LINE-out study can actually be run on.

Nothing here fits anything.  It writes `work/substrate224.json` and a cached common-axis
tensor `work/surface-<name>.npz` for each surface the ladder will use.

CONVENTIONS THAT ARE LOAD-BEARING (all inherited from 159 §5, which this reproduces):
  * every line is MEAN-REMOVED over its own targets first.  The across-target mean IS the
    metric's zero (`b` in AGENTS.md §4), so no predictor may be credited for it.
  * panel target coordinates are dropped from the gene axis (`pds_cosine` drops all 300).
  * rows are unit-normalised before averaging (AGENTS.md §4.2: equal donor weights win).
  * k562_gwps and k562_essential are the SAME cell line, two libraries.  They may never
    appear on opposite sides of a leave-one-LINE-out split.

    uv run python substrate224.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

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
PANEL_SET = set(PANEL)

# name -> (path, cell-line identity group, lineage, semantics, note)
# `group` is the LINE identity: two entries sharing a group are the same cell line and must
# never straddle a LOLO split.
LINES: dict[str, dict] = {
    "k562_gwps": dict(
        path="data/donors/k562_gwps", group="K562", lineage="myeloid (CML)",
        semantics="log-pseudobulk delta", lab="Weissman 2022", role="production donor"),
    "k562_essential": dict(
        path="data/donors/k562_essential", group="K562", lineage="myeloid (CML)",
        semantics="log-pseudobulk delta", lab="Weissman 2022", role="training-only (0/300)"),
    "orion_hct116": dict(
        path="data/donors/orion_hct116", group="HCT116", lineage="epithelial (colorectal)",
        semantics="log-pseudobulk delta", lab="Xaira 2026", role="production donor"),
    "orion_hek293t": dict(
        path="data/donors/orion_hek293t", group="HEK293T", lineage="epithelial (kidney)",
        semantics="log-pseudobulk delta", lab="Xaira 2026", role="production donor"),
    "rpe1": dict(
        path="data/donors/rpe1", group="RPE1", lineage="epithelial (retinal pigment)",
        semantics="log-pseudobulk delta", lab="Weissman 2022", role="training-only (0/300)"),
    "nadig_jurkat": dict(
        path="data/donors/nadig_jurkat", group="JURKAT", lineage="lymphoid T",
        semantics="log-pseudobulk delta", lab="Nadig 2025", role="training-only (0/300)"),
    "nadig_hepg2": dict(
        path="experiments/166-donor-acquisition/donors/nadig_hepg2", group="HEPG2",
        lineage="hepatocyte (hepatoblastoma)", semantics="log-pseudobulk delta",
        lab="Nadig 2025", role="training-only (0/300)"),
    "h1_2025": dict(
        path="data/donors/h1_2025", group="H1", lineage="pluripotent (hESC)",
        semantics="log-pseudobulk delta", lab="Arc 2025",
        role="candidate — LEAKS h1-test (85/85) and is inside 222's union"),
    "frangieh_melanoma": dict(
        path="experiments/166-donor-acquisition/donors/frangieh_melanoma", group="MELANOMA",
        lineage="melanocyte (patient melanoma)", semantics="log-pseudobulk delta",
        lab="Frangieh 2021", role="training-only (5/300)"),
    "dld1_lowmoi": dict(
        path="experiments/223-silent-zero-audit/work/tables/dld1_lowmoi", group="DLD1",
        lineage="epithelial (colorectal)", semantics="lfc matrix (223 build)",
        lab="GSE337988", role="223's live source (67/300)"),
    "telohaec": dict(
        path="data/external/igvf-telohaec", group="TELOHAEC",
        lineage="endothelial (aortic)", semantics="log-pseudobulk delta",
        lab="IGVF", role="candidate — 105 and 223 both NO-GO"),
    "kolf2": dict(
        path="data/external/kolf2", group="KOLF2", lineage="pluripotent (iPSC)",
        semantics="per-gene z-score (NOT lfc)", lab="Y-Doctor/JAX 2025",
        role="candidate — 159 §4 anti-aligned, carries its own blob"),
    "feng_ipsc": dict(
        path="experiments/099-donors/work/feng_ipsc", group="FENG_IPSC",
        lineage="pluripotent (iPSC, multi-donor pool)", semantics="pooled target-level lfc",
        lab="Feng 2025", role="candidate — 223 fails the row-matched control at every dose"),
    "myllia_a549": dict(
        path="experiments/120-myllia-a549/work/compact-donor", group="A549",
        lineage="epithelial (lung adeno)", semantics="log-pseudobulk delta",
        lab="Myllia", role="candidate"),
    "gse295214_h1": dict(
        path="experiments/115-gse295214-h1/work/compact-donor", group="H1",
        lineage="pluripotent (hESC)", semantics="log-pseudobulk delta",
        lab="Pan 2025", role="candidate — same LINE as h1_2025"),
}

# CD4 is three stimulation-state banks of PRIMARY human T cells from 4 donors.  Production
# treats them as one source; 159 §4 treats rest/stim48hr as two rows of one line.
CD4_ROOT = "data/external/cd4t/cloud-compact/cd4_compact"
CD4_CONDS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")


def load_meta(rel: str) -> tuple[list[str], list[str], tuple[int, int]]:
    p = ROOT / rel
    perts = [x.strip() for x in (p / "perts.txt").read_text().splitlines() if x.strip()]
    genes = [x.strip() for x in (p / "genes.txt").read_text().splitlines() if x.strip()]
    shape = tuple(map(int, np.load(p / "delta.npy", mmap_mode="r").shape))
    return perts, genes, shape


def unit(A: np.ndarray) -> np.ndarray:
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def load_effects(rel: str, targets: list[str], genes: list[str]) -> np.ndarray:
    """Rows for `targets`, columns for `genes`, from a memmapped delta.npy."""
    p = ROOT / rel
    pl = [x.strip() for x in (p / "perts.txt").read_text().splitlines() if x.strip()]
    gl = [x.strip() for x in (p / "genes.txt").read_text().splitlines() if x.strip()]
    ip = {x: i for i, x in enumerate(pl)}
    ig = {x: i for i, x in enumerate(gl)}
    d = np.load(p / "delta.npy", mmap_mode="r")
    ri = np.asarray([ip[t] for t in targets])
    ci = np.asarray([ig[g] for g in genes])
    # row-chunked: fancy-indexing a 17,767 x 38,606 memmap in one go materialises the whole
    # row block twice (float32 then float64) and took the crosspanel run to 9.6 GiB
    out = np.empty((len(ri), len(ci)), dtype=np.float64)
    for s in range(0, len(ri), 256):
        out[s:s + 256] = np.asarray(d[ri[s:s + 256]][:, ci], dtype=np.float64)
    return np.nan_to_num(out, copy=False)


def load_basal(rel: str, genes: list[str]) -> np.ndarray | None:
    p = ROOT / rel / "b_ctrl.npy"
    if not p.exists():
        return None
    gl = [x.strip() for x in (ROOT / rel / "genes.txt").read_text().splitlines() if x.strip()]
    ig = {x: i for i, x in enumerate(gl)}
    b = np.nan_to_num(np.load(p).astype(np.float64)).ravel()
    if b.shape[0] != len(gl):
        return None
    return b[np.asarray([ig[g] for g in genes])]


def main() -> None:
    inv = {}
    for name, meta in LINES.items():
        try:
            perts, genes, shape = load_meta(meta["path"])
        except Exception as exc:                                     # noqa: BLE001
            inv[name] = dict(**meta, error=str(exc))
            continue
        inv[name] = dict(
            **meta, n_targets=len(perts), n_genes=len(genes), shape=list(shape),
            panel=len(set(perts) & PANEL_SET),
            has_basal=(ROOT / meta["path"] / "b_ctrl.npy").exists())

    # CD4: production rule is a target counts only if all three banks carry it.
    banks = {c: set(x.strip() for x in
                    (ROOT / CD4_ROOT / c / "perts.txt").read_text().splitlines() if x.strip())
             for c in CD4_CONDS}
    cd4_complete = set.intersection(*banks.values())
    cd4_genes = [x.strip() for x in
                 (ROOT / CD4_ROOT / "cd4_rest/genes.txt").read_text().splitlines() if x.strip()]
    inv["cd4_rest"] = dict(
        path=f"{CD4_ROOT}/cd4_rest", group="CD4T", lineage="lymphoid T (PRIMARY)",
        semantics="log2FC x ln2", lab="GSE314342", role="production residual source",
        n_targets=len(banks["cd4_rest"]), n_genes=len(cd4_genes),
        shape=[len(banks["cd4_rest"]), len(cd4_genes)],
        panel=len(banks["cd4_rest"] & PANEL_SET), has_basal=True,
        cd4_complete_all_three=len(cd4_complete),
        cd4_complete_panel=len(cd4_complete & PANEL_SET))
    inv["cd4_stim48hr"] = dict(
        path=f"{CD4_ROOT}/cd4_stim48hr", group="CD4T", lineage="lymphoid T (PRIMARY)",
        semantics="log2FC x ln2", lab="GSE314342", role="production residual source",
        n_targets=len(banks["cd4_stim48hr"]), n_genes=len(cd4_genes),
        shape=[len(banks["cd4_stim48hr"]), len(cd4_genes)],
        panel=len(banks["cd4_stim48hr"] & PANEL_SET), has_basal=True)

    names = [n for n in inv if "error" not in inv[n]]
    pertsets, genesets = {}, {}
    for n in names:
        p = ROOT / inv[n]["path"]
        pertsets[n] = set(x.strip() for x in (p / "perts.txt").read_text().splitlines()
                          if x.strip())
        genesets[n] = set(x.strip() for x in (p / "genes.txt").read_text().splitlines()
                          if x.strip())

    print(f"{'line':18s} {'group':10s} {'targets':>8s} {'genes':>7s} {'panel':>6s} "
          f"{'basal':>6s}  role")
    print("-" * 100)
    for n in sorted(names, key=lambda x: -inv[x]["panel"]):
        r = inv[n]
        print(f"{n:18s} {r['group']:10s} {r['n_targets']:8d} {r['n_genes']:7d} "
              f"{r['panel']:6d} {str(r['has_basal']):>6s}  {r['role']}")

    print("\nPAIRWISE SHARED TARGETS (upper) / SHARED GENES (lower)")
    hdr = f"{'':18s}" + "".join(f"{n[:9]:>10s}" for n in names)
    print(hdr)
    pair = {}
    for a in names:
        row = f"{a:18s}"
        for b in names:
            if a == b:
                row += f"{'-':>10s}"
                continue
            st = len(pertsets[a] & pertsets[b])
            sg = len(genesets[a] & genesets[b])
            pair[f"{a}|{b}"] = dict(targets=st, genes=sg)
            row += f"{st:>10d}" if names.index(a) < names.index(b) else f"{sg:>10d}"
        print(row)

    out = dict(inventory=inv, pairwise=pair, panel=len(PANEL))
    (WORK / "substrate224.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {WORK / 'substrate224.json'}")


if __name__ == "__main__":
    main()
