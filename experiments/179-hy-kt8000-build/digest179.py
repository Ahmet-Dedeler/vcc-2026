"""179 gate 4 — the DELIVERED pseudobulk identity, on the two real artifacts.

WHAT THIS IS FOR
----------------
`pds` (BRIEF section 1) and the profile half of `mse` are computed from the per-(context,
perturbation) pseudobulk, i.e. from each gene's PANEL TOTAL over the 400 emitted cells.
`emit_counts` splices the parametric and resampled halves and only THEN repairs each gene's
column total to `floor(n * L * profile)`; `L = int(median(lib_sizes))` and the picks come from
the `pick_rng` stream, which the emission law never touches. So the claim is that the two
builds deliver bit-identical panel totals and `Delta pds = Delta mse_profile = 0` is an
IDENTITY rather than an estimate.

This checks the claim on the artifacts themselves rather than on the algebra. Experiment 156's
first a23 verification was an algebraic tautology that read 1.00000000 for a build with the fix
and for one without it; 171 section 6 is the non-vacuous form and this is that form applied to
the two packages: one sha256 of the delivered `int64` panel totals per block.

WHAT IT READS IF THE THING IT TESTS IS ABSENT
---------------------------------------------
Without the post-splice repair the Poisson/rounding residual lands in the totals and every
block differs (177 measured exactly this: the arms separate the moment the repair is removed).
The script also prints the digest of the raw NONZERO STRUCTURE per block, which MUST differ on
900/900 -- if that column matched too, the two files would simply be the same file and the
totals identity would be vacuous.

    uv run python digest179.py <built.h5ad> <a23.h5ad>
"""
from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path

import h5py
import numpy as np

N_GENES = 18533


def read_obs_blocks(path: Path):
    """(context, target_gene) per row block, straight off the categorical encoding."""
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        out = {}
        for name in ("context", "target_gene"):
            group = obs[name]
            cats = np.asarray(group["categories"][:]).astype(str)
            codes = np.asarray(group["codes"][:])
            out[name] = cats[codes]
        indptr = np.asarray(handle["X"]["indptr"][:])
        shape = tuple(handle["X"].attrs["shape"])
    labels = list(zip(out["context"], out["target_gene"]))
    blocks, start = [], 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            blocks.append((labels[start], start, i))
            start = i
    return blocks, indptr, shape


def block_digests(path: Path, blocks, indptr, tag: str):
    totals_sha, struct_sha, totals_sum = {}, {}, {}
    t0 = time.time()
    with h5py.File(path, "r") as handle:
        data = handle["X"]["data"]
        indices = handle["X"]["indices"]
        for k, (key, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            cols = indices[lo:hi].astype(np.int64)
            vals = data[lo:hi]
            totals = np.zeros(N_GENES, dtype=np.int64)
            np.add.at(totals, cols, vals.astype(np.int64))
            totals_sha[key] = hashlib.sha256(totals.tobytes()).hexdigest()[:16]
            totals_sum[key] = int(totals.sum())
            # the nonzero PATTERN, which the emission law does move
            struct_sha[key] = hashlib.sha256(
                np.bincount(cols, minlength=N_GENES).astype(np.int64).tobytes()
            ).hexdigest()[:16]
            if (k + 1) % 100 == 0:
                print(f"  [{tag}] {k + 1}/{len(blocks)} [{time.time() - t0:.0f}s]", flush=True)
    return totals_sha, struct_sha, totals_sum


def main() -> None:
    built, incumbent = Path(sys.argv[1]), Path(sys.argv[2])
    b_blocks, b_indptr, b_shape = read_obs_blocks(built)
    a_blocks, a_indptr, a_shape = read_obs_blocks(incumbent)
    print(f"built      {built}  shape {b_shape}  {len(b_blocks)} blocks")
    print(f"incumbent  {incumbent}  shape {a_shape}  {len(a_blocks)} blocks")
    b_keys = [k for k, _, _ in b_blocks]
    a_keys = [k for k, _, _ in a_blocks]
    if b_keys != a_keys:
        raise SystemExit(
            "block (context, target) ORDER differs between the two files -- the digest "
            f"comparison would be meaningless. first mismatch at "
            f"{next(i for i, (x, y) in enumerate(zip(b_keys, a_keys)) if x != y)}"
        )
    sizes = {hi - lo for _, lo, hi in b_blocks} | {hi - lo for _, lo, hi in a_blocks}
    print(f"block sizes {sizes} (must be {{400}})")

    b_tot, b_str, b_sum = block_digests(built, b_blocks, b_indptr, "built")
    a_tot, a_str, a_sum = block_digests(incumbent, a_blocks, a_indptr, "a23")

    same_tot = sum(int(b_tot[k] == a_tot[k]) for k in b_keys)
    same_str = sum(int(b_str[k] == a_str[k]) for k in b_keys)
    same_sum = sum(int(b_sum[k] == a_sum[k]) for k in b_keys)
    print()
    print(f"per-gene PANEL TOTAL digest identical on   {same_tot}/{len(b_keys)}")
    print(f"panel total SUM identical on               {same_sum}/{len(b_keys)}")
    print(f"nonzero-STRUCTURE digest identical on      {same_str}/{len(b_keys)}"
          "   (must be far below 900, or the two files are the same file)")
    bad = [k for k in b_keys if b_tot[k] != a_tot[k]]
    if bad:
        print(f"FIRST MISMATCHES: {bad[:5]}")
    print("\nsample rows:")
    for k in b_keys[:3] + b_keys[-2:]:
        print(f"  {k}  totals {b_tot[k]} vs {a_tot[k]}   struct {b_str[k]} vs {a_str[k]}")
    print(f"\nGATE 4 {'PASS' if same_tot == len(b_keys) else 'FAIL'}: "
          f"{same_tot}/{len(b_keys)} delivered pseudobulk digests identical")


if __name__ == "__main__":
    main()
