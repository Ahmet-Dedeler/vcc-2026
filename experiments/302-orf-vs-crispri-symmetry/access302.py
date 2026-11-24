"""302 — THE SINGLE ACCESS to reserve R4.  One streaming pass over `JoungZhang2023_atlas.h5ad`'s `X`.

Collects, for every scored ORF code and for every (scored symbol × cell shard), the summed raw
counts on the frozen 15,710-gene axis and the summed library size over ALL 37,528 columns, plus
the two known-answer rows that need the cell-level values (KA-E library size, KA-H integrality).
Everything downstream (`analyse302.py`) reads only the sums this writes — `X` is opened once.

    PYTHONDONTWRITEBYTECODE=1 uv run python access302.py            # the reserve
    PYTHONDONTWRITEBYTECODE=1 uv run python access302.py --atlas <mock.h5ad> --tag mock

Writes work/sums302[_tag].npz and work/access302[_tag].json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib302 as L  # noqa: E402

BLOCK = 20_000


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default=str(L.ATLAS))
    ap.add_argument("--tag", default="")
    ap.add_argument("--block", type=int, default=BLOCK)
    args = ap.parse_args()
    atlas = Path(args.atlas)
    tag = f"_{args.tag}" if args.tag else ""
    real = atlas.resolve() == L.ATLAS.resolve()

    t0 = time.time()
    if real:
        d = sha256(atlas)
        print(f"KA-0  atlas sha256 {d[:16]}… -> {'MATCH' if d == L.ATLAS_SHA else 'MISMATCH'}",
              flush=True)
        assert d == L.ATLAS_SHA, "atlas digest does not match the declared reserve"
        d2 = sha256(L.ORFMAP)
        assert d2 == L.ORFMAP_SHA, "301's ORF map digest changed"
        print(f"KA-D  ORF map sha256 {d2[:16]}… MATCH  ({time.time() - t0:.0f}s for both digests)",
              flush=True)

    D = L.build_design(atlas)
    nA = len(D["axis"])
    nS, nC = D["n_sym_groups"], D["n_code_groups"]
    print(f"design: {D['n_cells']} cells, axis {nA} genes, scored symbols {len(D['scored'])} "
          f"(primary {len(D['primary'])}, secondary {len(D['secondary'])}), "
          f"{len(D['used_codes'])} ORF codes, {nS} symbol groups × {L.N_SHARDS} shards, "
          f"{nC} code groups", flush=True)

    axis_pos = np.full(D["n_var"], -1, np.int32)
    axis_pos[D["axis_var"]] = np.arange(nA, dtype=np.int32)

    g_sym, g_code, shard = D["g_sym"], D["g_code"], D["shard"]
    g_ss = np.where(g_sym >= 0, g_sym * L.N_SHARDS + shard, -1).astype(np.int32)

    acc_sym = np.zeros(nS * L.N_SHARDS * nA, np.float64)
    acc_code = np.zeros(nC * nA, np.float64)
    lib_sym = np.zeros(nS * L.N_SHARDS, np.float64)
    lib_code = np.zeros(nC, np.float64)
    cells_sym = np.zeros(nS * L.N_SHARDS, np.int64)
    cells_code = np.zeros(nC, np.int64)

    ka_e_max = 0.0
    ka_h_max = 0.0
    n_nnz = 0
    with h5py.File(atlas, "r") as f:
        Xd, Xi, Xp = f["X/data"], f["X/indices"], f["X/indptr"]
        indptr = Xp[:].astype(np.int64)
        ncounts = D["ncounts"]
        nblk = (D["n_cells"] + args.block - 1) // args.block
        for bi in range(nblk):
            lo = bi * args.block
            hi = min(lo + args.block, D["n_cells"])
            s0, s1 = int(indptr[lo]), int(indptr[hi])
            if s1 == s0:
                continue
            data = Xd[s0:s1]
            idx = Xi[s0:s1]
            n_nnz += s1 - s0
            cnt = np.diff(indptr[lo:hi + 1])
            rowl = np.repeat(np.arange(hi - lo, dtype=np.int32), cnt)

            # KA-H: integrality
            fr = float(np.abs(data - np.rint(data)).max()) if data.size else 0.0
            ka_h_max = max(ka_h_max, fr)
            # KA-E: computed library size == obs['ncounts']
            lib = np.bincount(rowl, weights=data.astype(np.float64), minlength=hi - lo)
            ka_e_max = max(ka_e_max, float(np.abs(lib - ncounts[lo:hi]).max()))

            gss_b = g_ss[lo:hi]
            gc_b = g_code[lo:hi]
            use = gss_b >= 0
            if use.any():
                np.add.at(lib_sym, gss_b[use], lib[use])
                np.add.at(cells_sym, gss_b[use], 1)
                np.add.at(lib_code, gc_b[use], lib[use])
                np.add.at(cells_code, gc_b[use], 1)

                col = axis_pos[idx]
                keep = (col >= 0) & use[rowl]
                if keep.any():
                    v = data[keep].astype(np.float64)
                    c = col[keep].astype(np.int64)
                    r = rowl[keep]
                    acc_sym += np.bincount(gss_b[r].astype(np.int64) * nA + c, weights=v,
                                           minlength=acc_sym.size)
                    acc_code += np.bincount(gc_b[r].astype(np.int64) * nA + c, weights=v,
                                            minlength=acc_code.size)
                del col, keep
            del data, idx, rowl, cnt, lib
            if bi % 10 == 0 or bi == nblk - 1:
                rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
                print(f"  block {bi + 1}/{nblk} rows {lo}-{hi} nnz {n_nnz/1e9:.3f}e9  "
                      f"{time.time() - t0:.0f}s  peak RSS {rss:.2f} GiB", flush=True)

    acc_sym = acc_sym.reshape(nS, L.N_SHARDS, nA)
    acc_code = acc_code.reshape(nC, nA)
    print(f"KA-E  max |computed library size − obs['ncounts']| over all cells: {ka_e_max:.3g} "
          f"-> {'PASS' if ka_e_max == 0 else 'FAIL'}", flush=True)
    print(f"KA-H  max fractional part of X: {ka_h_max:.3g} "
          f"-> {'PASS' if ka_h_max == 0 else 'FAIL'}", flush=True)
    print(f"cells accumulated: symbols {int(cells_sym.sum())}, codes {int(cells_code.sum())}",
          flush=True)

    np.savez_compressed(
        L.WORK / f"sums302{tag}.npz",
        axis=np.asarray(D["axis"]), scored=np.asarray(D["scored"]),
        primary=np.asarray(D["primary"]), secondary=np.asarray(D["secondary"]),
        used_codes=np.asarray(D["used_codes"]), acc_sym=acc_sym, acc_code=acc_code,
        lib_sym=lib_sym.reshape(nS, L.N_SHARDS), lib_code=lib_code,
        cells_sym=cells_sym.reshape(nS, L.N_SHARDS), cells_code=cells_code,
        draws=np.asarray([np.concatenate(d) for d in D["draws"]]),
        g_gfp=D["g_gfp"], g_mch=D["g_mch"], c_gfp=D["c_gfp"], c_mch=D["c_mch"],
        axis_h1=D["axis_h1"])
    out = dict(atlas=str(atlas), real=real, n_cells=D["n_cells"], n_axis=nA, n_nnz=int(n_nnz),
               n_scored=len(D["scored"]), n_primary=len(D["primary"]),
               n_secondary=len(D["secondary"]), n_codes=len(D["used_codes"]),
               ka_e_max=ka_e_max, ka_e_pass=bool(ka_e_max == 0),
               ka_h_max=ka_h_max, ka_h_pass=bool(ka_h_max == 0),
               cells_sym=int(cells_sym.sum()), cells_code=int(cells_code.sum()),
               seconds=time.time() - t0,
               rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (L.WORK / f"access302{tag}.json").write_text(json.dumps(out, indent=1))
    print(f"\ndone in {out['seconds']:.0f}s, peak RSS {out['rss_gib']:.2f} GiB "
          f"-> work/sums302{tag}.npz")


if __name__ == "__main__":
    sys.exit(main())
