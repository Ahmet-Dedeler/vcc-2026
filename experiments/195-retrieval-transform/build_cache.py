"""195 — build a memory-lean cache of the 019 LODO harness objects.

Experiment 019's `interact.py` is the harness that produced AGENTS.md 4.2's headline
(incumbent 0.7378 -> whitened 0.7873 -> equal-weight 0.7978).  It loads three donor banks
whole and keeps two float32 copies of each deviation matrix; for the two Orion lines that is
~16 GB, over this run's 12 GB ceiling.  Nothing in the harness needs the full 38,606-gene
axis: every quantity it computes lives on `feat`, the 7,760 genes shared by all three lines
minus the 259 VCC panel targets they all measure.

So we stream each donor's `delta.npy` from a memmap in row blocks, subtract the perturbation
mean, and keep only the `feat` columns.  `gene_norm`'s per-gene divisor is defined on the
donor's FULL axis (its 25th-percentile floor is taken over all genes), so that is accumulated
over the full axis first and only then restricted -- restricting before flooring would change
the floor and silently move every number downstream.

Written per donor:
    s_raw   (n_perts, 7760) float32   deviations delta - mean(delta), feat columns
    scale   (7760,)         float64   _norm_scale over the FULL axis, feat columns
    perts   (n_perts,)      str
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import _norm_scale  # noqa: E402

OUT = Path(__file__).resolve().parent / "work" / "cache"
LODO = ["k562_gwps", "orion_hct116", "orion_hek293t"]
BLOCK = 1024


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_arr(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    banks = {n: ROOT / "data" / "donors" / n for n in LODO}
    genes = {n: (p / "genes.txt").read_text().split("\n") for n, p in banks.items()}
    perts = {n: (p / "perts.txt").read_text().split("\n") for n, p in banks.items()}

    shared = sorted(set.intersection(*[set(g) for g in genes.values()]))
    vcc = list(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"])
    idx = {n: {p: i for i, p in enumerate(perts[n])} for n in LODO}
    common = [g for g in vcc if all(g in idx[n] for n in LODO)]
    feat = [g for g in shared if g not in set(common)]
    print(f"shared genes {len(shared)}  common panel perts {len(common)}  feat {len(feat)}")

    manifest = {"feat_n": len(feat), "common_n": len(common), "donors": {}}
    for n in LODO:
        src = banks[n] / "delta.npy"
        gpos = {g: i for i, g in enumerate(genes[n])}
        col = np.fromiter((gpos[g] for g in feat), np.int64, len(feat))
        A = np.load(src, mmap_mode="r")
        nrow, ncol = A.shape
        assert ncol == len(genes[n])

        # pass 1: perturbation mean over the full axis
        acc = np.zeros(ncol, np.float64)
        for i in range(0, nrow, BLOCK):
            acc += np.asarray(A[i:i + BLOCK], np.float64).sum(0)
        m = (acc / nrow).astype(np.float32)

        # pass 2: full-axis mean |s| (for _norm_scale's floor) and the feat-column deviations
        absacc = np.zeros(ncol, np.float64)
        s_raw = np.empty((nrow, len(feat)), np.float32)
        for i in range(0, nrow, BLOCK):
            blk = np.asarray(A[i:i + BLOCK], np.float32) - m
            absacc += np.abs(blk).sum(0, dtype=np.float64)
            s_raw[i:i + BLOCK] = blk[:, col]
        del A

        # _norm_scale(s) with s the full-axis deviation matrix, then restricted to feat.
        full_scale = _norm_scale_from_mean(absacc / nrow)
        scale = full_scale[col]

        np.save(OUT / f"{n}.s_raw.npy", s_raw)
        np.save(OUT / f"{n}.scale.npy", scale)
        manifest["donors"][n] = {
            "delta_path": str(src.relative_to(ROOT)),
            "delta_sha256": sha256_file(src),
            "n_perts": int(nrow), "n_genes_full": int(ncol),
            "s_raw_sha256": sha256_arr(s_raw),
            "scale_sha256": sha256_arr(scale),
        }
        print(f"{n:16} {s_raw.shape} cached  sha {manifest['donors'][n]['s_raw_sha256'][:16]}")
        del s_raw

    (OUT / "feat.txt").write_text("\n".join(feat))
    (OUT / "common.txt").write_text("\n".join(common))
    manifest["feat_sha256"] = hashlib.sha256("\n".join(feat).encode()).hexdigest()
    manifest["common_sha256"] = hashlib.sha256("\n".join(common).encode()).hexdigest()
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print("wrote", OUT / "manifest.json")


def _norm_scale_from_mean(sc: np.ndarray, floor_pct: float = 25.0) -> np.ndarray:
    """`_norm_scale` re-expressed on a precomputed mean-|s| vector (identical output)."""
    live = sc[sc > 0]
    floor = np.percentile(live, floor_pct) if live.size else 1e-9
    return np.maximum(sc, floor)


if __name__ == "__main__":
    main()
