"""180 — build a synthetic D/E/F bundle: A/B/C renamed, nothing else changed.

The final panel lands on Oct 22 and the builder has never seen a bundle that is not A/B/C.
The only way to know the context adapter works before then is to make one. This copies the
real controls, renames the files AND the `obs['context']` label to D / E / F, and rewrites
`manifest.json` -- so the bundle is byte-different from the pinned one in exactly the way a
real new round would be, while every number the builder computes from it stays comparable to
the A/B/C build.

The `obs['context']` categories dataset is rewritten in place through h5py rather than by a
read-modify-write through anndata, so the whole job stays under a few hundred MB of RSS
instead of loading three 110M-nonzero matrices.

    uv run python defbundle180.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "data/challenge"
OUT = HERE / "work/def-bundle"
RENAME = {"A": "D", "B": "E", "C": "F"}


def sha256(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block):
            digest.update(chunk)
    return digest.hexdigest()


def relabel_context(path: Path, old: str, new: str) -> None:
    with h5py.File(path, "r+") as handle:
        group = handle["obs"]["context"]
        categories = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["categories"][:]
        ]
        if categories != [old]:
            raise ValueError(f"{path}: expected categories ['{old}'], got {categories}")
        attrs = dict(group["categories"].attrs)
        del group["categories"]
        dataset = group.create_dataset(
            "categories",
            data=np.asarray([new], dtype=object),
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
        for key, value in attrs.items():
            dataset.attrs[key] = value


def main() -> int:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    for name in ("gene_names.csv", "pert_counts.csv"):
        shutil.copyfile(SRC / name, OUT / name)

    manifest = json.loads((SRC / "manifest.json").read_text())
    manifest["panel_id"] = "vcc2026-synthetic-def-1"
    manifest["partition"] = "synthetic-def"
    manifest["contexts"] = [RENAME[c] for c in manifest["contexts"]]
    manifest["per_context"] = {
        RENAME[key]: value for key, value in manifest["per_context"].items()
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    for old, new in RENAME.items():
        target = OUT / f"context_{new}.h5ad"
        shutil.copyfile(SRC / f"context_{old}.h5ad", target)
        relabel_context(target, old, new)
        with h5py.File(target, "r") as handle:
            categories = [
                v.decode() if isinstance(v, bytes) else str(v)
                for v in handle["obs"]["context"]["categories"][:]
            ]
            rows = int(handle["X"].attrs["shape"][0])
            nnz = int(handle["X"]["data"].shape[0])
        source_sha = sha256(SRC / f"context_{old}.h5ad")
        target_sha = sha256(target)
        print(
            f"[def] context_{old}.h5ad -> context_{new}.h5ad  obs.context {categories}  "
            f"{rows} rows  nnz {nnz}\n"
            f"[def]   source sha {source_sha}\n[def]   target sha {target_sha}",
            flush=True,
        )
        if source_sha == target_sha:
            raise AssertionError("the relabelled file is byte-identical to its source")

    digests = {
        f"bundle/{p.name}": sha256(p) for p in sorted(OUT.iterdir()) if p.is_file()
    }
    (HERE / "work/def-bundle-hashes.json").write_text(json.dumps(digests, indent=1) + "\n")
    print(json.dumps(digests, indent=1))
    print(f"[def] bundle at {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
