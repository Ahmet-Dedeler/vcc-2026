"""597 addendum (2026-10-02) — save the real-Flex NTC pseudobulks themselves (597's score run kept only the scores).

Same RDS, parser and selection as `score597.py score <rds> <spill> hash.ID Gene_assignment '^NO\\.TARGET$'`. Per line
(K562 / HAP1 / HEK293): the RNA assay's raw `counts` summed over that line's non-targeting cells. Data: Satija lab
VIPerturb-seq pilot, Zenodo 18460279, CC BY 4.0 (10x Flex v1, custom probes, 1 % formaldehyde 10 min fixation; HEK293,
not HEK293T).

Known answer (else nothing is written): the RDS sha256, the cell counts and the per-line totals equal 597's run
(work/score597.json, work/score597.log: 498 / 478 / 486 cells; 2,623,903 / 3,103,451 / 5,555,320 counts).

Writes work/pseudobulk597.npz:
  lines (3,), n_cells (3,), genes (F,) the object's RNA feature names, counts (3, F) integer-valued float64 sums,
  log1p_cp50k (3, F) = log1p(counts / row total * 5e4) over the object's own features,
  official_genes (18,533,) = data/challenge/gene_names.csv order, counts_official (3, 18,533) with NaN where the probe set
  lacks the gene (exactly what 597 passed to 590's matcher), meta (a JSON string)
and work/pseudobulk597.json (the npz sha256 + the same meta). Light: < 1 GB RSS, seconds.

  uv run python experiments/597-flex-hek293-floor/pseudobulk597.py <rds> <spill dir>
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from score597 import LINES, OFFICIAL, column, dimnames, log, named, parse  # noqa: E402

WANT = {"rds_sha256": "3218665cfdc93b3aa98a2ef78c67ebf98320a82ccfb10ecf0d01242c8eee5217",
        "cells": {"K562": 498, "HAP1": 478, "HEK293": 486},
        "totals": {"K562": 2623903, "HAP1": 3103451, "HEK293": 5555320}}
LINE_COL, GUIDE_COL, NT_REGEX = "hash.ID", "Gene_assignment", r"^NO\.TARGET$"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main(rds: Path, spill: Path) -> int:
    rds_sha = sha256(rds)
    obj = parse(rds, spill)
    rna = named(obj.attrs["assays"], "RNA")
    cnt = named(rna.attrs["layers"], "counts")
    assert list(cnt.attrs["class"]) == ["dgCMatrix"], cnt.attrs["class"]
    features = dimnames(rna.attrs["features"])[0]
    meta = obj.attrs["meta.data"]
    cells = list(meta.attrs["row.names"])
    assert [int(v) for v in cnt.attrs["Dim"]] == [len(features), len(cells)]
    p = np.asarray(cnt.attrs["p"], dtype=np.int64)
    ispill, xspill = cnt.attrs["i"], cnt.attrs["x"]
    line = column(meta, LINE_COL).astype(str)
    guide = column(meta, GUIDE_COL).astype(str)
    nt = np.asarray([bool(re.search(NT_REGEX, g)) for g in guide])
    lab = np.full(len(cells), "", dtype=object)
    for key in LINES:                       # score597's mapping, unchanged
        lab[[key.lower() in str(x).lower().replace("-", "") for x in line]] = key
    groups = [k for k in LINES if ((lab == k) & nt).sum() > 0]
    code = np.full(len(cells), -1, np.int64)
    for gi, k in enumerate(groups):
        code[np.flatnonzero((lab == k) & nt)] = gi
    acc = np.zeros((len(groups), len(features)), dtype=np.float64)
    integer = True
    with open(ispill.path, "rb") as fi, open(xspill.path, "rb") as fx:
        for c0 in range(0, len(cells), 4000):
            c1 = min(c0 + 4000, len(cells))
            a, z = int(p[c0]), int(p[c1])
            fi.seek(4 * a); fx.seek(4 * a)
            ii = np.fromfile(fi, dtype=np.int32, count=z - a)
            xx = np.fromfile(fx, dtype=np.float32, count=z - a).astype(np.float64)
            integer &= bool(np.all(xx == np.round(xx)))
            local = np.repeat(np.arange(c1 - c0), np.diff(p[c0:c1 + 1]))
            g = code[c0:c1][local]
            m = g >= 0
            np.add.at(acc, (g[m], ii[m]), xx[m])
    got = {"rds_sha256": rds_sha, "cells": {k: int((code == i).sum()) for i, k in enumerate(groups)},
           "totals": {k: int(round(acc[i].sum())) for i, k in enumerate(groups)}}
    ka = got == WANT and integer
    log(f"known answer {'PASS' if ka else 'FAIL'}: {got} integer={integer}")
    if not ka:
        return 1
    pos = {g: j for j, g in enumerate(features)}
    off = np.full((len(groups), len(OFFICIAL)), np.nan)
    idx = np.array([pos.get(g, -1) for g in OFFICIAL])
    have = idx >= 0
    off[:, have] = acc[:, idx[have]]
    cp = np.log1p(acc / acc.sum(axis=1, keepdims=True) * 5e4)
    meta_out = {"source": "Satija lab VIPerturb-seq pilot, multimodal_cell_line_mixing_pilot.rds, Zenodo 18460279, CC BY 4.0",
                "chemistry": "10x Flex v1, custom probes, 1% formaldehyde 10 min; KRAB-dCas9-MeCP2 lines pooled, split by hashing",
                "selection": {"line_col": LINE_COL, "guide_col": GUIDE_COL, "nt_regex": NT_REGEX,
                              "assay": "RNA", "layer": "counts"},
                "known_answer": got, "n_features": len(features), "n_official_present": int(have.sum()),
                "script_sha256": sha256(Path(__file__)), "caveats": "HEK293 (not HEK293T); ~480-500 NT cells per line"}
    out = HERE / "work/pseudobulk597.npz"
    np.savez_compressed(out, lines=np.asarray(groups), n_cells=np.asarray([got["cells"][k] for k in groups]),
                        genes=np.asarray(features), counts=acc, log1p_cp50k=cp, official_genes=np.asarray(OFFICIAL),
                        counts_official=off, meta=np.asarray(json.dumps(meta_out)))
    (HERE / "work/pseudobulk597.json").write_text(json.dumps({"npz": str(out.relative_to(HERE.parents[1])),
                                                             "npz_sha256": sha256(out), **meta_out}, indent=1) + "\n")
    log(f"wrote {out.name} ({out.stat().st_size} bytes): {len(features)} native genes, {int(have.sum())} on the official axis")
    shutil.rmtree(spill, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
