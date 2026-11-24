"""305 step 2 -- the reusable sufficient statistics: exact per-construct, per-half count sums.

One split per invocation, streamed SERIALLY in CSR row blocks:

    PYTHONDONTWRITEBYTECODE=1 uv run python sums305.py <train|validation|test>

For every construct of the split (targeting AND non-targeting), the split's cells are partitioned into
two deterministic, batch-balanced halves (PREREGISTRATION sec 5) and the exact integer count sum over
ALL 18,080 gene columns is accumulated for each half.  The full library total per group is therefore
kept, so the frozen response formula's denominator is exact no matter which columns are read later.

Writes `work/sums305_<split>.npz` (int64 sums, names, halves, cell counts, UMI totals) and
`work/sums305_<split>.json` (metadata, per-half lane composition, hashes, KA-5's direct-slice check).
Nothing is written outside this directory; no old script is invoked.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402


def direct_slice_sum(path: Path, rows: np.ndarray, n_genes: int) -> np.ndarray:
    """KA-5: sum the given rows by reading each row's CSR slice directly -- no chunk loop, no one-hot."""
    out = np.zeros(n_genes, dtype=np.int64)
    with h5py.File(path, "r") as h:
        indptr = h["X/indptr"][:]
        for r in np.sort(rows):
            s, e = int(indptr[r]), int(indptr[r + 1])
            np.add.at(out, h["X/indices"][s:e], h["X/data"][s:e].astype(np.int64))
    return out


def main() -> None:
    split = sys.argv[1]
    path = X.SPLITS[split]
    t0 = time.time()
    P = X.L.Logger(X.HERE / f"sums305_{split}.log")
    P(f"305 sums305.py {split} -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    gid, tgt, bat, genes, shape = X.read_obs(path)
    n_cells, n_genes = shape
    cons = np.asarray([X.construct_of(g) for g in gid])
    allc = sorted(set(cons.tolist()))
    cpos = {c: i for i, c in enumerate(allc)}
    is_ctrl = tgt == X.CTRL
    ntc_cons = sorted(set(cons[is_ctrl].tolist()))
    ntc_par = {c: i % 2 for i, c in enumerate(ntc_cons)}          # 293 run293.py:92

    # ---- deterministic, batch-balanced halves (PREREGISTRATION sec 5); no RNG
    half = np.full(n_cells, -1, dtype=np.int8)
    for c in allc:
        ix = np.flatnonzero(cons == c)
        half[ix] = X.halves(ix, bat[ix])
    assert (half >= 0).all()
    groups = 2 * np.asarray([cpos[c] for c in cons], dtype=np.int64) + half.astype(np.int64)
    n_groups = 2 * len(allc)
    P(f"  {n_cells:,} cells x {n_genes} genes; {len(allc)} constructs "
      f"({len(allc) - len(ntc_cons)} targeting + {len(ntc_cons)} NTC); {n_groups} half-groups; "
      f"accumulator {n_groups * n_genes * 8 / 2**20:.0f} MiB")

    acc, cnt, g_chk = X.stream_group_sums(path, groups, n_groups, log=P)
    assert g_chk == n_genes
    S = np.rint(acc).astype(np.int64)
    del acc
    umi = S.sum(1)
    P(f"  streamed in {time.time() - t0:.0f}s, peak {X.peak_gib():.2f} GiB; "
      f"total UMIs {umi.sum():,}; median UMI/cell {umi.sum() / n_cells:,.0f}")

    names = np.asarray([c for c in allc for _ in (0, 1)])
    hids = np.asarray([h for _ in allc for h in (0, 1)], dtype=np.int8)

    # ---- KA-5: chunked sums vs a direct per-row slice, on the smallest construct under 1,500 cells
    sizes = {c: int((cons == c).sum()) for c in allc}
    ka5c = min((c for c in allc if sizes[c] < 1500), key=lambda c: sizes[c])
    ka5rows = np.flatnonzero(cons == ka5c)
    direct = direct_slice_sum(path, ka5rows, n_genes)
    chunked = S[2 * cpos[ka5c]] + S[2 * cpos[ka5c] + 1]
    ka5_ok = bool(np.array_equal(direct, chunked))
    P(f"  KA-5 direct-slice check on {ka5c} ({ka5rows.size} cells, {int(direct.sum()):,} UMIs): "
      f"{'PASS' if ka5_ok else 'FAIL'} (max |diff| {int(np.abs(direct - chunked).max())})")

    lane_comp = {}
    for c in allc:
        for h in (0, 1):
            m = (cons == c) & (half == h)
            lane_comp[f"{c}|{h}"] = {b: int(((bat == b) & m).sum()) for b in sorted(set(bat[m].tolist()))}

    meta = dict(
        split=split, file=str(path.relative_to(X.ROOT)), sha256=X.sha256_file(path),
        n_cells=int(n_cells), n_genes=int(n_genes), constructs=allc,
        ntc_constructs=ntc_cons, ntc_parity=ntc_par,
        target_of_construct={c: str(tgt[np.flatnonzero(cons == c)[0]]) for c in allc},
        cells_per_half={f"{names[i]}|{hids[i]}": int(cnt[i]) for i in range(n_groups)},
        umi_per_half={f"{names[i]}|{hids[i]}": int(umi[i]) for i in range(n_groups)},
        lane_composition=lane_comp,
        ka5=dict(construct=ka5c, n_cells=int(ka5rows.size), ok=ka5_ok,
                 max_abs_diff=int(np.abs(direct - chunked).max()), umis=int(direct.sum())),
        seconds=round(time.time() - t0, 1), peak_rss_gib=round(X.peak_gib(), 3))
    np.savez_compressed(X.HERE / f"sums305_{split}.npz", sums=S, names=names, halves=hids,
                        cells=cnt.astype(np.int64), umis=umi.astype(np.int64), genes=genes)
    X.dump(meta, X.HERE / f"sums305_{split}.json")
    P(f"  wrote sums305_{split}.npz "
      f"({(X.HERE / f'sums305_{split}.npz').stat().st_size / 2**20:.1f} MiB), sha256 "
      f"{X.sha256_file(X.HERE / f'sums305_{split}.npz')[:16]}")
    P(f"done {meta['seconds']}s, peak {meta['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
