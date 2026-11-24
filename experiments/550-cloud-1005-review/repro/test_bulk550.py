#!/usr/bin/env python3
"""Synthetic h5ads (and a tiny .vcc) for bulk_identity550.py. No challenge data."""
from __future__ import annotations

import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import subprocess

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
CHECK = HERE / "bulk_identity550.py"


def write_h5ad(path: Path, bump: bool) -> None:
    genes = [f"g{i}" for i in range(6)]
    rows = []
    obs = []
    for ctx in ("A", "B"):
        for tgt in ("T1", "T2"):
            for cell in range(4):
                v = np.zeros(6, dtype=np.float32)
                v[0] = 2
                v[1] = cell + 1
                if bump and ctx == "B" and tgt == "T2" and cell == 0:
                    v[3] = 7
                rows.append(v)
                obs.append({"context": ctx, "target_gene": tgt})
    X = sparse.csr_matrix(np.vstack(rows))
    a = ad.AnnData(X=X, obs=pd.DataFrame(obs), var=pd.DataFrame(index=genes))
    a.write_h5ad(path)


def write_vcc(path: Path, h5ad: Path) -> None:
    blob = subprocess.check_output(["zstd", "-q", "-c"], input=h5ad.read_bytes())
    with tarfile.open(path, "w") as archive:
        info = tarfile.TarInfo(name="pred.h5ad.zst")
        info.size = len(blob)
        import io
        archive.addfile(info, io.BytesIO(blob))


def run(left: Path, right: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(CHECK), str(left), str(right)], capture_output=True, text=True)


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="bulk550-"))
    same_a = tmp / "a.h5ad"
    same_b = tmp / "b.h5ad"
    bumped = tmp / "bumped.h5ad"
    write_h5ad(same_a, bump=False)
    write_h5ad(same_b, bump=False)
    write_h5ad(bumped, bump=True)
    ok = run(same_a, same_b)
    print("identical rc", ok.returncode)
    print(ok.stdout.strip())
    print(ok.stderr.strip())
    bad = run(same_a, bumped)
    print("bumped rc", bad.returncode)
    print(bad.stdout.strip())
    vcc = tmp / "a.vcc"
    write_vcc(vcc, same_a)
    via = run(vcc, same_b)
    print("vcc rc", via.returncode)
    print(via.stdout.strip())
    print(via.stderr.strip())
    if ok.returncode != 0 or bad.returncode != 1 or via.returncode != 0:
        return 1
    if "DIFF B T2" not in bad.stdout:
        print("expected the bumped group to be named")
        return 1
    print("BULK CHECK PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
