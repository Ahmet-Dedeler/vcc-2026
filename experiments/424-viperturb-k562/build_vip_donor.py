"""424 — build data/donors/viperturb_k562 with the PRODUCTION donor-table code, unchanged.

`scripts/build_donor_table.build()` is the function that made data/donors/k562_gwps. It reads an h5ad (backed) for the
per-cell labels and gene axis, and `group_sums()` streams the count matrix. The VIPerturb bins are Seurat RDS, so the
two readers are replaced by stubs that serve the SAME inputs from convert_bin.py's exact per-label sums:
  * `ad.read_h5ad`  -> a stub with .obs['perturbation'] = the per-cell label (meta.data$gene; the 7,949 NO-TARGET cells
                       appear once, they are the identical cells in all three bins) and var_names = the 19,068 features;
  * `group_sums`    -> for each builder group, the label's summed raw counts (the builder's own group order).
Everything after that — min 20 cells, b = log1p(5e4 * P / sum P), delta = b_p - b_ctrl on the donor's own gene axis,
float32, sorted labels, meta.json — is the production code, byte for byte.

KNOWN ANSWERS (abort on failure):
  KA-stub   the stubbed group_sums reproduces the builder's contract: per group, the sum over exactly the cells the
            builder assigned to it (checked via cell counts per group against convert_bin's counts).
  KA-bins   NO-TARGET: identical barcodes and identical summed counts in binA/B/C; every target label in exactly one
            bin; per-label cell counts == the Zenodo manifest's cell_count for every label.

    uv run python build_vip_donor.py
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
OUT = ROOT / "data/donors/viperturb_k562"
BINS = ("binA", "binB", "binC")
CONTROL = "NO-TARGET"
sys.path.insert(0, str(ROOT / "scripts"))
import build_donor_table as bdt  # noqa: E402


def load_bin(b):
    d = WORK / "conv" / b
    return dict(sums=np.load(d / "sums.npy", mmap_mode="r"),
                groups=(d / "groups.txt").read_text().split("\n"),
                ncells=np.asarray([int(x) for x in (d / "ncells.txt").read_text().split("\n")]),
                features=(d / "features.txt").read_text().split("\n"),
                cells=pd.read_csv(d / "cells.tsv.gz", sep="\t", usecols=["barcode", "gene"]),
                info=json.loads((d / "info.json").read_text()))


class StubAnnData:
    def __init__(self, labels, genes):
        self.obs = pd.DataFrame({"perturbation": labels})
        self.var_names = pd.Index(genes)
        self.n_obs, self.n_vars = len(labels), len(genes)


def main():
    bins = {b: load_bin(b) for b in BINS}
    feats = bins["binA"]["features"]
    assert all(bins[b]["features"] == feats for b in BINS), "feature axes differ between bins"
    # ---- KA-bins
    ntc = {b: bins[b]["groups"].index(CONTROL) for b in BINS}
    bc = {b: bins[b]["cells"].loc[bins[b]["cells"].gene == CONTROL, "barcode"].to_numpy() for b in BINS}
    for b in BINS[1:]:
        assert np.array_equal(bc[b], bc["binA"]), f"NO-TARGET barcodes differ in {b}"
        assert np.array_equal(bins[b]["sums"][ntc[b]], bins["binA"]["sums"][ntc["binA"]]), f"NO-TARGET sums differ in {b}"
    label_sum, label_n, label_bin = {}, {}, {}
    for b in BINS:
        for k, g in enumerate(bins[b]["groups"]):
            if g == CONTROL and b != "binA":
                continue
            assert g not in label_sum, f"label {g} in two bins"
            label_sum[g] = (b, k)
            label_n[g] = int(bins[b]["ncells"][k])
            label_bin[g] = b
    man = pd.read_csv(WORK / "dl/genome_wide_manifest.txt", sep=r"\s+")
    mcount = man.groupby("gene").cell_count.first().to_dict()
    bad = {g: (label_n.get(g), mcount.get(g)) for g in set(mcount) | set(label_n) if label_n.get(g) != mcount.get(g)}
    print(f"KA-bins: NO-TARGET identical in 3 bins ({len(bc['binA'])} cells); {len(label_n)} labels; "
          f"labels whose cell count != manifest: {len(bad)} {list(bad.items())[:5]}")
    assert not bad, "cell counts do not reproduce the manifest"

    # per-cell labels in bin order, NO-TARGET once
    labels = np.concatenate([bins[b]["cells"].gene.to_numpy().astype(str)[
        (bins[b]["cells"].gene.to_numpy() != CONTROL) | (b == "binA")] for b in BINS])
    assert len(labels) == sum(label_n.values())

    # The builder holds the full (groups x genes) float64 sums plus ~3 temporaries of that size: ~9.3 GB for
    # 17.5k x 19k. Rows are independent (b_ctrl comes from the control alone), so build() runs on four contiguous
    # ranges of the sorted labels, each with the control, and the partial tables are concatenated in order.
    uniq = np.unique(labels)
    targets = [g for g in uniq if g != CONTROL]
    ranges = np.array_split(np.asarray(targets), 4)
    checked = {}
    parts = []
    for ci, rg in enumerate(ranges):
        keep = set(rg.tolist()) | {CONTROL}
        lab_c = labels[np.isin(labels, list(keep))]
        stub = StubAnnData(lab_c, feats)

        def group_sums(path, codes, n_groups, n_genes, block=16_384, lab_c=lab_c):
            """Stub for bdt.group_sums: the summed raw counts of each builder group, from convert_bin's exact sums."""
            assert n_genes == len(feats)
            out = np.zeros((n_groups, n_genes), dtype=np.float64)
            first = np.full(n_groups, -1)
            seen = codes >= 0
            idx = np.flatnonzero(seen)
            first[codes[idx[::-1]]] = idx[::-1]
            cnt = np.bincount(codes[seen], minlength=n_groups)
            for k in range(n_groups):
                g = lab_c[first[k]]
                b, j = label_sum[g]
                out[k] = bins[b]["sums"][j]
                assert cnt[k] == label_n[g], f"KA-stub: group {g} has {cnt[k]} cells, convert_bin counted {label_n[g]}"
                checked[g] = int(cnt[k])
            return out

        bdt.group_sums = group_sums
        bdt.ad.read_h5ad = lambda path, backed=None, stub=stub: stub
        d = WORK / "build_chunks" / f"c{ci}"
        m = bdt.build(Path(f"chunk{ci}"), CONTROL, d, 20, label_col="perturbation", gene_col=None)
        parts.append((d, m))
    print(f"KA-stub: {len(checked)} builder groups served, cell counts match convert_bin for all")
    # ---- assemble (the builder's own file formats)
    b0 = np.load(parts[0][0] / "b_ctrl.npy")
    g0 = (parts[0][0] / "genes.txt").read_text()
    for d, _ in parts[1:]:
        assert np.array_equal(np.load(d / "b_ctrl.npy"), b0) and (d / "genes.txt").read_text() == g0, "chunk mismatch"
    perts = sum(((d / "perts.txt").read_text().split("\n") for d, _ in parts), [])
    assert perts == sorted(perts) and len(set(perts)) == len(perts)
    delta = np.concatenate([np.load(d / "delta.npy") for d, _ in parts])
    OUT.mkdir(parents=True, exist_ok=True)
    np.save(OUT / "delta.npy", delta)
    np.save(OUT / "b_ctrl.npy", b0)
    (OUT / "perts.txt").write_text("\n".join(perts))
    (OUT / "genes.txt").write_text(g0)
    norms = np.linalg.norm(delta, axis=1)
    cells = [label_n[p] for p in perts]
    meta = dict(
        source="zenodo:18460279 genome_wide_bin{A,B,C}.RDS (VIPerturb-seq K562; see provenance.json)",
        control_label=CONTROL, target_sum=bdt.TS,
        n_perturbations=len(perts), n_genes=int(delta.shape[1]),
        n_control_cells=label_n[CONTROL],
        min_cells=20, skipped_low_cells=int(sum(m["skipped_low_cells"] for _, m in parts)),
        cells_per_pert=dict(min=int(min(cells)), median=int(np.median(cells)), max=int(max(cells))),
        effect_norm=dict(min=float(norms.min()), median=float(np.median(norms)), max=float(norms.max())),
    )
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"wrote {OUT}: delta {delta.shape}; meta {json.dumps(meta)}")
    del delta

    # ---- provenance sidecar
    sha = json.loads((WORK / "dl/sha256.json").read_text()) if (WORK / "dl/sha256.json").exists() else {}
    cells = pd.read_csv(WORK / "conv/cells_all_summary.csv") if (WORK / "conv/cells_all_summary.csv").exists() else None
    perts = (OUT / "perts.txt").read_text().split("\n")
    panel = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
    k562 = set((ROOT / "data/donors/k562_gwps/perts.txt").read_text().split("\n"))
    prov = dict(
        dataset="VIPerturb-seq genome-wide K562 CRISPRi screen (Bradu ... Satija, bioRxiv 2026.02.12.705613)",
        zenodo_record="https://zenodo.org/records/18460279", license="CC BY 4.0",
        files={f"genome_wide_{b}.RDS": dict(url=f"https://zenodo.org/api/records/18460279/files/genome_wide_{b}.RDS/content",
                                            bytes=bins[b]["info"]["input_bytes"], md5_zenodo=sha.get(b, {}).get("md5"),
                                            sha256=sha.get(b, {}).get("sha256"), n_cells=bins[b]["info"]["n_cells"],
                                            nnz=bins[b]["info"]["nnz"]) for b in BINS},
        manifest=dict(url="https://zenodo.org/api/records/18460279/files/genome_wide_manifest.txt/content",
                      sha256=hashlib.sha256((WORK / "dl/genome_wide_manifest.txt").read_bytes()).hexdigest()),
        reader="experiments/424-viperturb-k562/rds_stream.py + convert_bin.py (streaming XDR parser; group sums == R "
               "readRDS on binC for NO-TARGET + 4 targets, |d| = 0 on all 19,068 features; per-cell colSums == "
               "meta.data$nCount_RNA and nnz == nFeature_RNA for every cell in every bin)",
        builder="scripts/build_donor_table.py build() (production code, unmodified; its two readers stubbed by "
                "experiments/424-viperturb-k562/build_vip_donor.py)",
        builder_sha256=hashlib.sha256((ROOT / "scripts/build_donor_table.py").read_bytes()).hexdigest(),
        counts_layer="RNA assay, layer 'counts' (raw UMI; all values integer, audited)",
        control=dict(label=CONTROL, n_cells=label_n[CONTROL],
                     note="the same 7,949 cells (identical barcodes and counts) are shipped in every bin; used once"),
        qc=dict(
            cells="none added: the bins' cells are the authors' post-QC, single-guide-assigned cells (meta.data$gene)",
            min_cells=20, knockdown_filter="none (mirrors K562-GWPS, whose table keeps every perturbation with >= 20 "
            "cells regardless of knockdown); per-target own-gene knockdown is measured and reported in "
            "experiments/424-viperturb-k562/RESULTS.md and knockdown.csv",
            guides="pooled per target gene (1-3 sgRNAs per gene, GuEST-List); no per-guide filter, as K562-GWPS pools "
                   "its dual-guide constructs per gene"),
        counts=dict(n_labels=len(label_n) - 1, n_perturbations=len(perts),
                    panel_ge20=int(len(set(panel) & set(perts))),
                    panel_ge20_not_in_k562_gwps=sorted(set(panel) & set(perts) - k562)),
        cells_per_target={g: n for g, n in sorted(label_n.items()) if g != CONTROL},
        target_bin={g: b for g, b in sorted(label_bin.items()) if g != CONTROL},
        meta=meta)
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=1))
    print(json.dumps({k: v for k, v in prov["counts"].items()}, indent=1))


if __name__ == "__main__":
    main()
