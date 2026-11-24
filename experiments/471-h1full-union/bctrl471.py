"""471 — `b_ctrl.npy` for the denoised H1 bank (`data/donors/h1_2025_full`), so it can join the builder's H1 union.

The union's compact-donor format needs {delta, b_ctrl, genes, perts}; 458d wrote only delta/genes/perts. b_ctrl is the
bank's own control profile in the same frame as delta: log1p(5e4 * P_ctrl / sum P_ctrl) over ALL control cells, each
counted once across the three VCC 2025 files (build_donor_table.py's convention; float32).

It re-runs `denoise458d.py`'s summation loop verbatim (same files, chunking, duplicate-control rule and seeded halves, so
the float64 accumulation order is identical), then:
  KA-1: every perturbation's recomputed delta equals the committed `delta.npy` (26662c52...) bit for bit;
  KA-2: the cell counts equal the official pert_counts CSVs (458d's check).
Only then it writes `b_ctrl.npy` next to the bank and prints its sha256 for the builder's UNION_SOURCE_REGISTRY.
Nothing else in the bank changes.   usage: uv run python experiments/471-h1full-union/bctrl471.py [--write]
"""
from __future__ import annotations

import csv
import hashlib
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
FILES = [ROOT / "data/external/vcc2025/train/adata_Training.h5ad", ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
         ROOT / "data/external/vcc2025/test/adata_Test.h5ad"]
COUNTS = [ROOT / "data/external/vcc2025/train/pert_counts_Training.csv",
          ROOT / "data/external/vcc2025/validation/pert_counts_Validation.csv",
          ROOT / "data/external/vcc2025/test/pert_counts_Test.csv"]
BANK = ROOT / "data/donors/h1_2025_full"
DELTA_SHA = "26662c527164ddd434ae19a13076e915d1f82de9683c75e27e0af1f5b28591c0"
TS, CHUNK = 5e4, 4000
t0 = time.time()
log = lambda *x: print("[471]", *x, f"({time.time() - t0:.0f}s)", flush=True)
rng = np.random.default_rng(4580)          # 458d's seed: the halves decide the float64 accumulation order


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 24), b""):
            h.update(c)
    return h.hexdigest()


def cat(n):
    return n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)


def main():
    assert sha(BANK / "delta.npy") == DELTA_SHA, "the bank's delta is not 458d's 26662c52"
    genes = None
    sums: dict[str, np.ndarray] = {}
    ncell: dict[str, int] = {}
    seen_ctrl: set[str] = set()
    for path in FILES:                     # --- denoise458d.py's loop, verbatim apart from the log prefix ---
        with h5py.File(path, "r") as f:
            g = f["var/_index"]
            gv = (g[:] if not isinstance(g, h5py.Group) else g["values"][:]).astype(str)
            if genes is None:
                genes = gv
            assert np.array_equal(genes, gv), f"{path.name}: gene axis differs"
            tv = cat(f["obs/target_gene"])
            ix = f["obs/_index"]
            ids = (ix[:] if not isinstance(ix, h5py.Group) else ix["values"][:]).astype(str)
            is_ctrl = tv == "non-targeting"
            dup = is_ctrl & np.isin(ids, list(seen_ctrl)) if seen_ctrl else np.zeros(tv.size, bool)
            seen_ctrl |= set(ids[is_ctrl])
            half = rng.integers(0, 2, tv.size)
            keys = np.where(is_ctrl, "__ctrl__", tv)
            ukeys, kidx = np.unique(keys, return_inverse=True)
            for k in ukeys:
                sums.setdefault(k, np.zeros((3, len(genes))))
            kk, kc = np.unique(keys[~dup], return_counts=True)
            for k, c in zip(kk, kc):
                ncell[k] = ncell.get(k, 0) + int(c)
            n = tv.size
            ip = f["X/indptr"]
            log(f"{path.name}: {n} cells, {int(is_ctrl.sum())} controls ({int(dup.sum())} already counted)")
            for s in range(0, n, CHUNK):
                e = min(n, s + CHUNK)
                a, b = int(ip[s]), int(ip[e])
                M = sp.csr_matrix((f["X/data"][a:b], f["X/indices"][a:b], ip[s:e + 1] - a), shape=(e - s, len(genes)))
                keep = ~dup[s:e]
                for hv in (0, 1):
                    sel = keep & (half[s:e] == hv)
                    if not sel.any():
                        continue
                    G = sp.csr_matrix((np.ones(int(sel.sum())), (kidx[s:e][sel], np.flatnonzero(sel))),
                                      shape=(len(ukeys), e - s))
                    part = np.asarray((G @ M).todense())
                    for j in np.flatnonzero(np.asarray(G.sum(1)).ravel() > 0):
                        sums[ukeys[j]][0] += part[j]
                        sums[ukeys[j]][1 + hv] += part[j]
    want = {}
    for c in COUNTS:
        for r_ in csv.DictReader(open(c)):
            want[r_["target_gene"]] = int(r_["n_cells"])
    bad = {k: (ncell.get(k), v) for k, v in want.items() if ncell.get(k) != v}
    assert not bad and set(want) == set(sums) - {"__ctrl__"}, f"KA-2 cell counts: {list(bad.items())[:5]}"
    log(f"KA-2: all {len(want)} perturbations' cell counts equal the official pert_counts CSVs")
    lp = lambda c: np.log1p(TS * c / c.sum())
    ctrl = sums.pop("__ctrl__")
    perts = sorted(sums)
    cF = lp(ctrl[0])
    delta = np.stack([lp(sums[p][0]) - cF for p in perts]).astype(np.float32)
    old = np.load(BANK / "delta.npy")
    assert [l.strip() for l in open(BANK / "perts.txt")] == perts and [l.strip() for l in open(BANK / "genes.txt")] == list(genes)
    same = np.array_equal(old, delta)
    log(f"KA-1: recomputed delta == delta.npy bit for bit: {same} (max |d| {float(np.abs(old - delta).max()):.3g})")
    assert same, "KA-1 failed"
    b_ctrl = cF.astype(np.float32)
    log(f"b_ctrl: {b_ctrl.shape}, controls {len(seen_ctrl)} unique cells, expressed (CPM > 10) {int((np.expm1(b_ctrl) / TS * 1e6 > 10).sum())}")
    if "--write" in sys.argv:
        np.save(BANK / "b_ctrl.npy", b_ctrl)
        log(f"wrote {BANK / 'b_ctrl.npy'} sha256 {sha(BANK / 'b_ctrl.npy')}")
        for name in ("delta.npy", "genes.txt", "perts.txt"):
            log(f"registry sha256 {name}: {sha(BANK / name)}")


if __name__ == "__main__":
    main()
