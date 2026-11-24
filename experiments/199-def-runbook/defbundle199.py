"""Build a synthetic D/E/F bundle with a NEW 300-target panel, for the Oct-22 rehearsal.

180's rehearsal reused the real `pert_counts.csv`, so `panel_support` never had to move and
the donor coverage never changed — the two things a new panel actually breaks. This one draws
a fresh 300, so the census, the `panel_support` re-pin and the preflight are exercised for
real.

The panel is drawn from the donors' own coverage, matching AGENTS.md §2.1's selection
criterion ("adequate representation in existing public datasets"), and DISJOINT from the
current 300 so nothing can pass by carrying an A/B/C answer forward.

Also writes a deliberately CORRUPTED copy of the bundle (two gene symbols transposed in
`gene_names.csv`) as the preflight's negative control: a checker that passes on that file is
not checking the gene axis.

    uv run python experiments/199-def-runbook/defbundle199.py
"""

from __future__ import annotations

import hashlib
import json
import shutil
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data/challenge"
OUT = Path(__file__).resolve().parent / "work/def-bundle"
BAD = Path(__file__).resolve().parent / "work/def-bundle-badaxis"
RENAME = {"A": "D", "B": "E", "C": "F"}
SEED = 199


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def draw_panel(n: int = 300) -> list[str]:
    genes = [str(g) for g in pd.read_csv(SRC / "gene_names.csv")["gene_name"]]
    current = {str(t) for t in pd.read_csv(SRC / "pert_counts.csv")["target_gene"]}
    covered: set[str] = set()
    for donor in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        covered |= {ln.strip() for ln in
                    (ROOT / f"data/donors/{donor}/perts.txt").read_text().splitlines() if ln.strip()}
    eligible = sorted(set(genes) & covered - current)
    rng = np.random.default_rng(SEED)
    return sorted(rng.choice(np.asarray(eligible), size=n, replace=False).tolist())


def relabel(src: Path, dst: Path, new: str) -> dict:
    """Copy a control file and rewrite obs['context'] in place — never loads X."""
    shutil.copyfile(src, dst)
    with h5py.File(dst, "r+") as h:
        node = h["obs"]["context"]
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in node["categories"][:]]
        if len(cats) != 1:
            raise SystemExit(f"{src.name}: expected one context category, got {cats}")
        attrs = dict(node["categories"].attrs)
        del node["categories"]
        node.create_dataset("categories", data=np.asarray([new], dtype=object),
                            dtype=h5py.string_dtype())
        for k, v in attrs.items():
            node["categories"].attrs[k] = v
        rows = int(h["X"].attrs["shape"][0])
        nnz = int(h["X"]["data"].shape[0])
    return {"from": cats[0], "to": new, "rows": rows, "nnz": nnz, "sha256": sha256(dst)}


def main() -> None:
    t0 = time.time()
    for d in (OUT, BAD):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    panel = draw_panel()
    current = {str(t) for t in pd.read_csv(SRC / "pert_counts.csv")["target_gene"]}
    assert len(panel) == 300 and not (set(panel) & current)
    pd.DataFrame({"target_gene": panel}).to_csv(OUT / "pert_counts.csv", index=False)
    shutil.copyfile(SRC / "gene_names.csv", OUT / "gene_names.csv")

    manifest = json.loads((SRC / "manifest.json").read_text())
    manifest["partition"] = "test"
    manifest["panel_id"] = "vcc2026-synthetic-def-199"
    manifest["contexts"] = [RENAME[c] for c in manifest["contexts"]]
    manifest["per_context"] = {RENAME[k]: v for k, v in manifest["per_context"].items()}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    files = {}
    for old, new in RENAME.items():
        rec = relabel(SRC / f"context_{old}.h5ad", OUT / f"context_{new}.h5ad", new)
        files[f"context_{new}.h5ad"] = rec
        print(f"context_{old}.h5ad -> context_{new}.h5ad  obs.context ['{new}']  "
              f"{rec['rows']} rows  nnz {rec['nnz']}  sha {rec['sha256'][:16]}...", flush=True)

    # negative control: same bundle, two gene symbols transposed on the axis
    for f in OUT.iterdir():
        shutil.copyfile(f, BAD / f.name)
    g = pd.read_csv(BAD / "gene_names.csv")
    v = list(g["gene_name"])
    v[100], v[101] = v[101], v[100]
    pd.DataFrame({"gene_name": v}).to_csv(BAD / "gene_names.csv", index=False)

    out = {
        "seed": SEED,
        "panel_id": manifest["panel_id"],
        "contexts": manifest["contexts"],
        "panel_first5": panel[:5], "panel_last5": panel[-5:],
        "panel_overlap_with_current_300": 0,
        "files": {**{k: v["sha256"] for k, v in files.items()},
                  "manifest.json": sha256(OUT / "manifest.json"),
                  "gene_names.csv": sha256(OUT / "gene_names.csv"),
                  "pert_counts.csv": sha256(OUT / "pert_counts.csv")},
        "bad_axis_bundle": str(BAD.relative_to(ROOT)),
        "bad_axis_transposition": [v[101], v[100]],
        "seconds": round(time.time() - t0, 1),
    }
    (Path(__file__).resolve().parent / "work/def-bundle199.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("panel_id", "contexts", "panel_first5", "seconds")}, indent=2))


if __name__ == "__main__":
    main()
