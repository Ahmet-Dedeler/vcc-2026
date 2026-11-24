"""180 gate 3 — the D/E/F rehearsal delivered a23's counts, relabelled and nothing else.

The synthetic bundle IS A/B/C with the files, `obs['context']` and the manifest renamed, so a
correct context adapter must emit the same 2.38e9 counts in the same order and differ from
`submissions/a23-massfix.h5ad` in exactly one place: the three context labels.

That is a far stronger statement than "the build did not crash". If the adapter had, say,
paired context D with context_E's controls, or reordered the emission loop, or re-seeded a
stream, `X/data` would move and this would say so.

Streamed in chunks, so peak RSS stays a few hundred MB on two 19 GB files.

    uv run python defidentity180.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
A23 = REPO / "submissions/a23-massfix.h5ad"
DEF = WORK / "def-rehearsal.h5ad"
CHUNK = 1 << 26          # 64M elements per read


def dataset_digest(path: Path, key: str) -> tuple[str, int]:
    digest = hashlib.sha256()
    with h5py.File(path, "r") as handle:
        dataset = handle[key]
        total = int(dataset.shape[0])
        for start in range(0, total, CHUNK):
            block = dataset[start:min(start + CHUNK, total)]
            digest.update(np.ascontiguousarray(block).tobytes())
    return digest.hexdigest(), total


def categorical(path: Path, column: str):
    with h5py.File(path, "r") as handle:
        group = handle["obs"][column]
        cats = [v.decode() if isinstance(v, bytes) else str(v)
                for v in group["categories"][:]]
        codes = group["codes"][:]
    return cats, codes


def main() -> int:
    if not DEF.is_file():
        raise SystemExit(f"{DEF} is missing; run the D/E/F rehearsal build first")
    out = {}
    for key in ("X/data", "X/indices", "X/indptr"):
        a, na = dataset_digest(A23, key)
        b, nb = dataset_digest(DEF, key)
        out[key] = {"a23": a, "def": b, "equal": a == b, "n_a23": na, "n_def": nb}
        print(f"[id] {key:<12} n {na:>13,} vs {nb:>13,}   "
              f"{'IDENTICAL' if a == b else 'DIFFERENT'}   {a[:16]} / {b[:16]}", flush=True)

    for column in ("context", "target_gene"):
        ca, ka = categorical(A23, column)
        cb, kb = categorical(DEF, column)
        same_codes = np.array_equal(ka, kb)
        out[f"obs/{column}"] = {
            "a23_categories": ca if len(ca) < 20 else f"{len(ca)} categories",
            "def_categories": cb if len(cb) < 20 else f"{len(cb)} categories",
            "codes_identical": bool(same_codes),
            "categories_identical": ca == cb,
        }
        print(f"[id] obs/{column:<11} categories "
              f"{ca if len(ca) < 20 else str(len(ca)) + ' cats'} -> "
              f"{cb if len(cb) < 20 else str(len(cb)) + ' cats'}; codes "
              f"{'IDENTICAL' if same_codes else 'DIFFERENT'}", flush=True)

    with h5py.File(A23, "r") as ha, h5py.File(DEF, "r") as hb:
        var_a = [v.decode() if isinstance(v, bytes) else str(v) for v in ha["var"]["_index"][:]]
        var_b = [v.decode() if isinstance(v, bytes) else str(v) for v in hb["var"]["_index"][:]]
    out["var_index_identical"] = var_a == var_b
    print(f"[id] var index identical: {var_a == var_b} ({len(var_a)} genes)")

    counts_identical = all(out[k]["equal"] for k in ("X/data", "X/indices", "X/indptr"))
    labels_moved = not out["obs/context"]["categories_identical"]
    codes_same = out["obs/context"]["codes_identical"]
    verdict = (
        "PASS" if counts_identical and labels_moved and codes_same
        and out["obs/target_gene"]["codes_identical"] and out["var_index_identical"]
        else "FAIL"
    )
    print(f"\n[gate3] {verdict}: counts identical={counts_identical}, "
          f"context labels moved={labels_moved}, context codes identical={codes_same}")
    out["verdict"] = verdict
    (WORK / "defidentity180.json").write_text(json.dumps(out, indent=1))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
