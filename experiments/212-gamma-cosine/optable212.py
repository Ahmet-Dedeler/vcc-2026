"""212 — rule 11's operator table: every whitening operator on disk, with path, sha256,
max eigenvalue, condition number and shape.

Rule 11 exists because "eight distinct whitening operators exist across `experiments/`,
several with near-identical names, and on 2026-08-31 a whole evening's (gamma, lam) screen
was run against the wrong one". 170's gate 1 failed once for exactly this reason. This
script resolves every candidate to a file and a digest so nothing in RESULTS.md is a
description.

    uv run python optable212.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
PROD_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> None:
    seen, rows = set(), []
    cands = sorted(set(
        list((ROOT / "data/donors").glob("whitener*.npz"))
        + list(ROOT.glob("experiments/*/out*/whitener*.npz"))
        + list(ROOT.glob("experiments/*/cloud-output/whiteners/*.npz"))
        + list(ROOT.glob("experiments/*/work/operators/*.npz"))
        + list(ROOT.glob("experiments/*/*/whitener*.npz"))))
    for p in cands:
        try:
            with np.load(p) as z:
                if "evals" not in z:
                    continue
                ev = np.asarray(z["evals"], dtype=np.float64)
                shape = tuple(int(x) for x in np.asarray(z["evecs"]).shape) \
                    if "evecs" in z else None
        except Exception as exc:                        # noqa: BLE001
            rows.append(dict(path=str(p.relative_to(ROOT)), error=str(exc)))
            continue
        sha = sha256_file(p)
        if sha in seen:
            dup = True
        else:
            seen.add(sha)
            dup = False
        meta = p.with_suffix(".json")
        basis = None
        if meta.exists():
            try:
                basis = json.loads(meta.read_text()).get("basis")
            except Exception:                            # noqa: BLE001
                basis = None
        rows.append(dict(
            path=str(p.relative_to(ROOT)), sha256=sha, sha16=sha[:16],
            n=int(ev.size), evecs_shape=shape, basis=basis,
            max_eig=float(ev.max()), min_eig=float(ev.min()),
            condition=float(ev.max() / ev.min()) if ev.min() > 0 else float("inf"),
            is_production=(sha == PROD_SHA), duplicate_content=dup,
            bytes=int(p.stat().st_size)))
    rows.sort(key=lambda r: (not r.get("is_production", False), r["path"]))
    with open(WORK / "optable212.json", "w") as fh:
        json.dump(rows, fh, indent=1)
    print(f"{'path':<74} {'sha16':<17} {'max eig':>10} {'cond':>12} {'n':>7}")
    for r in rows:
        if "error" in r:
            print(f"{r['path']:<74} ERROR {r['error'][:40]}")
            continue
        mark = " <= PRODUCTION" if r["is_production"] else ""
        print(f"{r['path']:<74} {r['sha16']:<17} {r['max_eig']:10.4f} "
              f"{r['condition']:12.4g} {r['n']:7d}{mark}")
    n_prod = sum(1 for r in rows if r.get("is_production"))
    print(f"\n{len(rows)} operator files, {len(seen)} distinct digests; "
          f"{n_prod} of them ARE the production operator (same bytes, different paths).")


if __name__ == "__main__":
    main()
