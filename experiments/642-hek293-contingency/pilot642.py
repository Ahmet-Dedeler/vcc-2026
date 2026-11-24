"""642 step 1: real-Flex non-targeting pseudobulks (K562 / HAP1 / HEK293) from Satija's VIPerturb pilot, onto the
official axis, saved small (work/flexpilot642.npz). Same cells and arithmetic as 597's `score` (hash.ID, Gene_assignment
== NO.TARGET, raw RNA counts); only the matcher step is left out. The RDS is read in place (597's download, Ahmet's OK
2026-10-01) and never copied into the repo.

    uv run python experiments/642-hek293-contingency/pilot642.py <rds> <spill dir>
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/597-flex-hek293-floor"))
import score597 as S  # noqa: E402

LINES = ("K562", "HAP1", "HEK293")


def main(rds: Path, spill: Path) -> None:
    obj = S.parse(rds, spill)
    rna = S.named(obj.attrs["assays"], "RNA")
    cnt = S.named(rna.attrs["layers"], "counts")
    features = S.dimnames(rna.attrs["features"])[0]
    meta = obj.attrs["meta.data"]
    cells = list(meta.attrs["row.names"])
    p = np.asarray(cnt.attrs["p"], dtype=np.int64)
    ispill, xspill = cnt.attrs["i"], cnt.attrs["x"]
    line = S.column(meta, "hash.ID").astype(str)
    guide = S.column(meta, "Gene_assignment").astype(str)
    nt = np.asarray([bool(re.search(r"^NO\.TARGET$", g)) for g in guide])
    lab = np.full(len(cells), "", dtype=object)
    for key in LINES:
        lab[[key.lower() in str(x).lower().replace("-", "") for x in line]] = key
    code = np.full(len(cells), -1, np.int64)
    for gi, k in enumerate(LINES):
        code[(lab == k) & nt] = gi
    acc = np.zeros((len(LINES), len(features)), dtype=np.float64)
    with open(ispill.path, "rb") as fi, open(xspill.path, "rb") as fx:
        for c0 in range(0, len(cells), 4000):
            c1 = min(c0 + 4000, len(cells))
            a, z = int(p[c0]), int(p[c1])
            fi.seek(4 * a); fx.seek(4 * a)
            ii = np.fromfile(fi, dtype=np.int32, count=z - a)
            xx = np.fromfile(fx, dtype=np.float32, count=z - a).astype(np.float64)
            assert np.all(xx == np.round(xx))
            local = np.repeat(np.arange(c1 - c0), np.diff(p[c0:c1 + 1]))
            g = code[c0:c1][local]
            m = g >= 0
            np.add.at(acc, (g[m], ii[m]), xx[m])
    pos = {g: j for j, g in enumerate(features)}
    off = np.full((len(LINES), len(S.OFFICIAL)), np.nan)
    for j, gname in enumerate(S.OFFICIAL):
        q = pos.get(gname)
        if q is not None:
            off[:, j] = acc[:, q]
    n = {k: int(((lab == k) & nt).sum()) for k in LINES}
    tot = {k: float(acc[i].sum()) for i, k in enumerate(LINES)}
    # 597's recorded cells and totals are the known answer
    assert n == {"K562": 498, "HAP1": 478, "HEK293": 486}, n
    assert {k: round(v) for k, v in tot.items()} == {"K562": 2623903, "HAP1": 3103451, "HEK293": 5555320}, tot
    h = hashlib.md5()
    with open(rds, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    meta_out = dict(source="Zenodo 18460279 multimodal_cell_line_mixing_pilot.rds (Satija lab, VIPerturb-seq pilot), CC BY 4.0",
                    chemistry="10x Flex v1 (custom probes); KRAB-dCas9-MeCP2 lines, hashed", rds_md5=h.hexdigest(),
                    n_cells=n, totals=tot, n_official_present=int(np.isfinite(off[0]).sum()))
    assert meta_out["rds_md5"] == "3f910a17baa164df2bd66bfac3197360"
    np.savez(HERE / "work/flexpilot642.npz", counts=off, lines=np.asarray(LINES), meta=json.dumps(meta_out))
    print(json.dumps(meta_out, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
