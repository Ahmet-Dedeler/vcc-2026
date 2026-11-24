"""216 Phase 1b — read the OTHER published objects' structure without downloading them.

`GWCD4i.pseudobulk_merged.h5ad` is 41.5 GiB and `GWCD4i.DE_stats.by_guide.h5mu` is 27.4 GiB.
The Phase-1 questions about them -- how many genes, how many panel targets, how many guides
per target -- are answered by `.obs` and `.var`, which are a few MB inside those files.  So
this opens each one over HTTP range requests through a file-like object handed to `h5py`,
and downloads only the bytes HDF5 actually asks for.  Bytes fetched are counted and
reported, so the download budget is auditable.

    uv run python remote216.py
"""
from __future__ import annotations

import io
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BUCKET = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com"

PANEL = tuple(line.split(",")[0]
              for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
              if line.strip())
OFFICIAL = tuple(line.strip()
                 for line in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
                 if line.strip())


class HTTPRangeFile(io.RawIOBase):
    """Minimal seekable read-only file over S3 byte ranges, with a block cache."""

    BLOCK = 1 << 20

    def __init__(self, key: str, size: int):
        self.url = f"{BUCKET}/{key}"
        self.size = size
        self.pos = 0
        self.cache: dict[int, bytes] = {}
        self.bytes_fetched = 0
        self.requests = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, offset, whence=io.SEEK_SET):
        self.pos = (offset if whence == io.SEEK_SET
                    else self.pos + offset if whence == io.SEEK_CUR
                    else self.size + offset)
        return self.pos

    def tell(self):
        return self.pos

    def _block(self, idx: int) -> bytes:
        if idx not in self.cache:
            start = idx * self.BLOCK
            end = min(start + self.BLOCK, self.size) - 1
            out = subprocess.run(["curl", "-sS", "-r", f"{start}-{end}", self.url],
                                 capture_output=True, check=True)
            self.cache[idx] = out.stdout
            self.bytes_fetched += len(out.stdout)
            self.requests += 1
            if len(self.cache) > 4096:            # bounded memory: 4 GiB worst case
                self.cache.pop(next(iter(self.cache)))
        return self.cache[idx]

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        buf = bytearray()
        p = self.pos
        while len(buf) < n:
            idx, off = divmod(p, self.BLOCK)
            chunk = self._block(idx)[off:]
            take = min(len(chunk), n - len(buf))
            buf += chunk[:take]
            p += take
        b[:n] = buf[:n]
        self.pos += n
        return n


def decode(node):
    """Decode an AnnData-encoded 1-D obs/var column (plain array or categorical group)."""
    if isinstance(node, h5py.Group):
        cats = node["categories"][:]
        codes = node["codes"][:]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else c for c in cats])
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
        return out
    v = node[:]
    if v.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])
    return v


def inspect(key: str, size: int, obs_cols, var_path="var", obs_path="obs", root="") -> dict:
    f = HTTPRangeFile(key, size)
    rep: dict[str, object] = {"key": key, "size": size}
    with h5py.File(f, "r") as h:
        base = h[root] if root else h
        rep["top_level"] = sorted(base.keys())
        v = base[var_path]
        rep["var_columns"] = sorted(v.keys())
        gname = "gene_name" if "gene_name" in v else "_index"
        genes = decode(v[gname])
        rep["n_var"] = int(len(genes))
        rep["genes_on_official_axis"] = int(len(set(genes) & set(OFFICIAL)))
        o = base[obs_path]
        rep["obs_columns"] = sorted(o.keys())
        cols = {}
        for c in obs_cols:
            if c in o:
                cols[c] = decode(o[c])
        rep["n_obs"] = int(len(next(iter(cols.values())))) if cols else None
        if "perturbed_gene_name" in cols:
            pg = cols["perturbed_gene_name"]
            rep["n_unique_perturbed_genes"] = int(len(set(pg.tolist())))
            rep["panel_targets_present"] = int(len(set(pg.tolist()) & set(PANEL)))
        if "guide_id" in cols and "perturbed_gene_name" in cols:
            import collections
            per = collections.defaultdict(set)
            for g, gid in zip(cols["perturbed_gene_name"], cols["guide_id"]):
                per[g].add(gid)
            panel_per = {g: len(s) for g, s in per.items() if g in set(PANEL)}
            rep["panel_guides_per_target_median"] = (
                float(np.median(list(panel_per.values()))) if panel_per else None)
        for c in ("culture_condition", "donor_id", "guide_type", "keep_for_DE",
                  "keep_test_genes", "10xrun_id"):
            if c in cols:
                u, n = np.unique(cols[c].astype(str), return_counts=True)
                rep[f"obs_{c}"] = {str(a): int(b) for a, b in zip(u, n)}
        # panel coverage under the file's own DE-eligibility flags
        if "perturbed_gene_name" in cols and "keep_for_DE" in cols:
            keep = cols["keep_for_DE"].astype(str)
            mask = np.isin(keep, ("True", "true", "1"))
            rep["panel_targets_keep_for_DE"] = int(
                len(set(cols["perturbed_gene_name"][mask].tolist()) & set(PANEL)))
            if "culture_condition" in cols:
                cc = cols["culture_condition"].astype(str)
                per_cond = {}
                for cond in np.unique(cc):
                    sel = mask & (cc == cond)
                    per_cond[str(cond)] = int(
                        len(set(cols["perturbed_gene_name"][sel].tolist()) & set(PANEL)))
                rep["panel_targets_keep_for_DE_per_condition"] = per_cond
        rep["_genes"] = genes.tolist()
    rep["bytes_fetched"] = f.bytes_fetched
    rep["http_requests"] = f.requests
    print(f"[{key}] fetched {f.bytes_fetched/1e6:.1f} MB in {f.requests} requests", flush=True)
    return rep


def main() -> None:
    inv = {r["key"]: r for r in json.loads((WORK / "s3_inventory.json").read_text())}
    out: dict[str, object] = {}

    k = "marson2025_data/GWCD4i.pseudobulk_merged.h5ad"
    out["pseudobulk"] = inspect(
        k, inv[k]["size"],
        obs_cols=("perturbed_gene_name", "culture_condition", "donor_id", "guide_id",
                  "guide_type", "keep_for_DE", "keep_test_genes", "n_cells", "10xrun_id"))

    k2 = "marson2025_data/GWCD4i.DE_stats.by_guide.h5mu"
    try:
        out["by_guide_guide_1"] = inspect(
            k2, inv[k2]["size"],
            obs_cols=("target_contrast_gene_name", "culture_condition", "n_guides"),
            root="mod/guide_1")
    except Exception as exc:                                   # noqa: BLE001
        out["by_guide_guide_1"] = {"error": repr(exc)}

    genes_pb = set(out["pseudobulk"].get("_genes", []))
    with h5py.File(ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad", "r") as h:
        gn = h["var/gene_name"]
        de_genes = set(decode(gn).tolist())
    out["gene_axis_comparison"] = dict(
        pseudobulk_n_var=len(genes_pb),
        de_table_n_var=len(de_genes),
        pseudobulk_on_official=len(genes_pb & set(OFFICIAL)),
        de_table_on_official=len(de_genes & set(OFFICIAL)),
        extra_official_genes_in_pseudobulk=len((genes_pb - de_genes) & set(OFFICIAL)),
        de_genes_not_in_pseudobulk=len(de_genes - genes_pb))
    extra = sorted((genes_pb - de_genes) & set(OFFICIAL))
    (WORK / "extra_official_genes_pseudobulk.json").write_text(json.dumps(extra))

    for v in out.values():
        if isinstance(v, dict):
            v.pop("_genes", None)
    out["total_bytes_fetched"] = sum(
        v.get("bytes_fetched", 0) for v in out.values() if isinstance(v, dict))
    (WORK / "remote216.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str)[:6000])


if __name__ == "__main__":
    main()
