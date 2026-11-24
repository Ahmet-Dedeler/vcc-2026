"""645: cache HPA's cell-line RNA atlas (rna_celline.tsv, v23 2023-05-24; downloaded 2026-10-02 with Ahmet's standing OK)
as a (cell line x gene) log2(TPM + 1) matrix -- the same units as def_identity's CCLE cache (DepMap log2(TPM + 1)).

    uv run python experiments/645-hpa-references/load_hpa645.py   -> data/external/hpa/hpa_cache645.npz
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TSV = ROOT / "data/external/hpa/rna_celline.tsv"
OUT = ROOT / "data/external/hpa/hpa_cache645.npz"


def main():
    df = pd.read_csv(TSV, sep="\t", usecols=["Gene name", "Cell line", "TPM"],
                     dtype={"Gene name": "category", "Cell line": "category", "TPM": np.float32})
    lines = df["Cell line"].cat.categories.astype(str).to_numpy()
    genes = df["Gene name"].cat.categories.astype(str).to_numpy()
    X = np.zeros((len(lines), len(genes)), dtype=np.float32)
    # duplicated gene names (several Ensembl ids) are summed in TPM space, as a symbol-level profile
    np.add.at(X, (df["Cell line"].cat.codes.to_numpy(), df["Gene name"].cat.codes.to_numpy()), df["TPM"].to_numpy())
    X = np.log2(1.0 + X).astype(np.float32)
    np.savez(OUT, X=X, lines=np.asarray(lines, dtype="U"), genes=np.asarray(genes, dtype="U"))
    meta = {"source": "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip", "zip_sha256":
            "e5b11779308c2720cfe89d9cb19932ad4c0c3fea0fa37e8c61a14d501995852a", "zip_bytes": 205861141,
            "tsv_bytes": TSV.stat().st_size, "n_lines": len(lines), "n_genes": len(genes), "units": "log2(1 + TPM), symbol-summed"}
    (OUT.with_suffix(".json")).write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))
    pat = ("HEK", "293", "RPE", "TERT", "HaCaT", "BJ", "MCF10", "HME", "RPTEC", "HBEC", "TIME", "HUVEC", "BEAS", "fHDF",
           "K-562", "HeLa", "JURKAT", "Jurkat", "HCT", "A549", "MCF7", "Hep G2", "CAL", "THP", "HAP1", "U-2 OS")
    print([l for l in lines if any(p.lower() in l.lower() for p in pat)])


if __name__ == "__main__":
    main()
