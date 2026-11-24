#!/usr/bin/env python3
"""537 — does the H1 harness's missing zero repair matter? (RESULTS §3.1)

446bj2z builds with `VCC_DE_ZERO_REPAIR=-5` (428: repair gated genes whose block total is zero, delivered lfc ~ -35); the
real-emitter H1 harness (435 → 528 → 537) runs without it. The repair acts only on a gated (target, gene) pair with a zero
block total that is not the target's own gene, so it is a no-op wherever no such pair exists. This counts them in every
H1 stats file priced here (this directory's, 528's, 435's; pre-zs and zs). < 0.5 GB, seconds.
usage: uv run python experiments/537-arm-pipeline/zerocheck537.py
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def main():
    ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    keep = ref["keep"].astype(bool)
    genes, targets = ref["genes"].astype(str), ref["targets"].astype(str)
    own = np.equal.outer(targets, genes)
    pats = ("experiments/537-arm-pipeline/work/h1_*/stats-*.npz", "experiments/528-callprop-block/work/h1_*/stats-*.npz",
            "experiments/435-zs-on-coarse/work/h1/stats-base_s*_bps_prodem_mde0.1cap0.1min0.3coarse*.npz")
    files = sorted(f for p in pats for f in glob.glob(str(ROOT / p)))
    n_files, non_own, own_pairs = 0, 0, set()
    for f in files:
        s = np.asarray(np.load(f, allow_pickle=True)["sums"])
        z = (s == 0) & keep[None, :]
        n_files += 1
        non_own += int((z & ~own).sum())
        own_pairs |= {(str(targets[i]), str(genes[j])) for i, j in zip(*np.nonzero(z & own))}
    print(f"{n_files} stats files; zero-total gated pairs that are NOT the target's own gene: {non_own}; own-gene zero pairs: "
          f"{sorted(own_pairs)} -> the zero repair is {'a no-op on every file' if non_own == 0 else 'NOT a no-op'}")


if __name__ == "__main__":
    main()
