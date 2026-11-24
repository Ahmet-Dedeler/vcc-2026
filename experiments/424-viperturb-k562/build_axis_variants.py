"""424 — VIPerturb donor tables on RESTRICTED gene axes, built by the same production code (build_vip_donor's stubs
around scripts/build_donor_table.build), low-memory (8 label chunks, file-offset row reads, streamed .npy writes).

Why: K562-GWPS's table is built on the gene axis of its scPerturb file, i.e. Replogle's 8,248 K562-EXPRESSED genes
(normalisation total over those genes only). The primary viperturb_k562 table keeps the whole 19,068-probe Flex axis.
"Same gene axis rules" can be read either way, so the axis is priced, not assumed:
  full     the primary table rebuilt here with 8 chunks (KNOWN ANSWER: byte-identical delta/b_ctrl/perts/genes to
           data/donors/viperturb_k562, i.e. the chunking cannot change a byte)
  gwpsaxis VIPerturb features ∩ K562-GWPS's gene axis (7,539 genes), totals over those genes (the K562 donor's own axis)
  floor    VIPerturb features whose NO-TARGET fraction >= K562-GWPS's smallest control fraction (5.97e-6; 9,209 genes),
           i.e. K562-GWPS's expression floor applied to VIPerturb's own controls
Variants are written under ./work/donors/<name> (not data/donors). The table on disk in data/donors stays the full axis.

    uv run python build_axis_variants.py full gwpsaxis floor
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import hashlib  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import build_vip_donor as B  # noqa: E402  (main() is not run on import)

ROOT, WORK, CONTROL, BINS, bdt = B.ROOT, B.WORK, B.CONTROL, B.BINS, B.bdt


class RowReader:
    """Row reads from a float64 C-order .npy by file offset (no mmap, so nothing lingers in RSS)."""

    def __init__(self, path):
        self.f = open(path, "rb")
        v = np.lib.format.read_magic(self.f)
        rh = np.lib.format.read_array_header_1_0 if v == (1, 0) else np.lib.format.read_array_header_2_0
        self.shape, fo, dt = rh(self.f)
        assert not fo and dt == np.float64
        self.off = self.f.tell()

    def row(self, j):
        self.f.seek(self.off + 8 * self.shape[1] * j)
        return np.fromfile(self.f, dtype=np.float64, count=self.shape[1])


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def build(name, n_chunks=8):
    feats = (WORK / "conv/binA/features.txt").read_text().split("\n")
    rd = {b: RowReader(WORK / f"conv/{b}/sums.npy") for b in BINS}
    groups = {b: (WORK / f"conv/{b}/groups.txt").read_text().split("\n") for b in BINS}
    ncell = {b: [int(x) for x in (WORK / f"conv/{b}/ncells.txt").read_text().split("\n")] for b in BINS}
    label_sum, label_n = {}, {}
    for b in BINS:
        for k, g in enumerate(groups[b]):
            if g == CONTROL and b != "binA":
                continue
            label_sum[g], label_n[g] = (b, k), ncell[b][k]
    if name == "full":
        cols = np.arange(len(feats))
    elif name == "gwpsaxis":
        gg = set((ROOT / "data/donors/k562_gwps/genes.txt").read_text().split("\n"))
        cols = np.asarray([i for i, g in enumerate(feats) if g in gg])
    elif name == "floor":
        bg = np.load(ROOT / "data/donors/k562_gwps/b_ctrl.npy").astype(np.float64)
        floor = float((np.expm1(bg) / 5e4).min())
        c = rd["binA"].row(label_sum[CONTROL][1])
        cols = np.flatnonzero(c / c.sum() >= floor)
    else:
        raise ValueError(name)
    gsub = [feats[i] for i in cols]
    import pandas as pd
    labels = np.concatenate([pd.read_csv(WORK / f"conv/{b}/cells.tsv.gz", sep="\t", usecols=["gene"]).gene.to_numpy()
                             .astype(str) for b in BINS])
    # NO-TARGET appears in all three bins: keep binA's copy only (bins are concatenated A, B, C)
    nA = sum(ncell["binA"])
    labels = np.concatenate([labels[:nA], labels[nA:][labels[nA:] != CONTROL]])
    assert len(labels) == sum(label_n.values())
    targets = [g for g in np.unique(labels) if g != CONTROL]
    out = WORK / "donors" / (f"viperturb_k562_{name}")
    tmp = WORK / f"build_chunks_{name}"
    out.mkdir(parents=True, exist_ok=True)
    perts, b0, g0, skipped, norms_all = [], None, None, 0, []
    rows_written = 0
    fh = open(out / "delta.npy.part", "wb")
    for ci, rg in enumerate(np.array_split(np.asarray(targets), n_chunks)):
        keep = set(rg.tolist()) | {CONTROL}
        lab_c = labels[np.isin(labels, list(keep))]
        stub = B.StubAnnData(lab_c, gsub)

        def group_sums(path, codes, n_groups, n_genes, block=16_384, lab_c=lab_c):
            assert n_genes == len(gsub)
            o = np.zeros((n_groups, n_genes), dtype=np.float64)
            first = np.full(n_groups, -1)
            seen = codes >= 0
            idx = np.flatnonzero(seen)
            first[codes[idx[::-1]]] = idx[::-1]
            cnt = np.bincount(codes[seen], minlength=n_groups)
            for k in range(n_groups):
                g = lab_c[first[k]]
                b, j = label_sum[g]
                o[k] = rd[b].row(j)[cols]
                assert cnt[k] == label_n[g]
            return o

        bdt.group_sums = group_sums
        bdt.ad.read_h5ad = lambda path, backed=None, stub=stub: stub
        m = bdt.build(Path(f"chunk{ci}"), CONTROL, tmp / f"c{ci}", 20, label_col="perturbation", gene_col=None)
        skipped += m["skipped_low_cells"]
        d = np.load(tmp / f"c{ci}" / "delta.npy")
        bc = np.load(tmp / f"c{ci}" / "b_ctrl.npy")
        gt = (tmp / f"c{ci}" / "genes.txt").read_text()
        if b0 is None:
            b0, g0 = bc, gt
        assert np.array_equal(bc, b0) and gt == g0
        perts += (tmp / f"c{ci}" / "perts.txt").read_text().split("\n")
        assert d.dtype == np.float32 and d.flags.c_contiguous
        fh.write(d.tobytes())
        rows_written += d.shape[0]
        norms_all.append(np.linalg.norm(d, axis=1))
        del d
        shutil.rmtree(tmp / f"c{ci}")
    fh.close()
    assert perts == sorted(perts)
    # prepend the exact header np.save would write for this array
    hdr = {"descr": np.lib.format.dtype_to_descr(np.dtype("<f4")), "fortran_order": False,
           "shape": (rows_written, len(gsub))}
    with open(out / "delta.npy", "wb") as f:
        np.lib.format.write_array_header_1_0(f, hdr)
        with open(out / "delta.npy.part", "rb") as src:
            shutil.copyfileobj(src, f, length=1 << 24)
    (out / "delta.npy.part").unlink()
    np.save(out / "b_ctrl.npy", b0)
    (out / "perts.txt").write_text("\n".join(perts))
    (out / "genes.txt").write_text(g0)
    norms = np.concatenate(norms_all)
    cells = [label_n[p] for p in perts]
    meta = dict(source=f"zenodo:18460279 VIPerturb-seq K562, axis variant '{name}' (experiments/424)", control_label=CONTROL,
                target_sum=bdt.TS, n_perturbations=len(perts), n_genes=len(gsub), n_control_cells=label_n[CONTROL],
                min_cells=20, skipped_low_cells=int(skipped),
                cells_per_pert=dict(min=int(min(cells)), median=int(np.median(cells)), max=int(max(cells))),
                effect_norm=dict(min=float(norms.min()), median=float(np.median(norms)), max=float(norms.max())))
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"[{name}] wrote {out}: {rows_written} x {len(gsub)}; effect norm median {meta['effect_norm']['median']:.3f}", flush=True)
    return out, meta


def main(names):
    for name in names:
        out, meta = build(name)
        if name == "full":
            prim = ROOT / "data/donors/viperturb_k562"
            same = {f: sha(out / f) == sha(prim / f) for f in ("delta.npy", "b_ctrl.npy", "perts.txt", "genes.txt")}
            pm = json.loads((prim / "meta.json").read_text())
            same["meta"] = all(pm[k] == meta[k] for k in meta if k != "source")
            print(f"KA-CHUNK full axis rebuilt with 8 chunks vs data/donors/viperturb_k562 (4 chunks): {same}", flush=True)
            assert all(same.values()), "chunking changed the table"
            shutil.rmtree(out)                     # identical to the primary; do not keep a duplicate
    import resource
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30:.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1:])
