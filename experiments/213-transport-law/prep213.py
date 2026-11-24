"""213 Phase 1 — donor-side sufficient statistics under an arbitrary transport law.

`E062.prepare_main_source` does exactly three things to a donor bank, in this order:

  1. `mean_all`   = across-ALL-perturbations column mean of `delta[:, all_source_columns]`
  2. `raw_median` = median over rows of `|| delta[row, all_source_columns] - mean_all ||_2`
  3. `effects[t]` = (delta[row_t, unique_source_columns] - unique_mean) * (EFFECT_SCALE/raw_median)

Under a transport law the object that is centred, scaled and averaged is not `delta` but the
law-space effect `g = law.to_law(delta, p_ctrl_donor)`.  Steps 1-3 are otherwise byte-identical,
including the 256-row chunking and the float64 accumulation order, so that `logc20` reproduces
`prepare_main_source` EXACTLY rather than approximately.

WHAT IS CACHED: only the two full-matrix streams (`mean_on_axis`, `raw_median_norm`) --
148 kB per (donor, law).  The 300 panel rows are re-gathered on demand because they are
5 M elements and cost nothing, and caching them in float32 would cost the identity row.

PERTURBMATCH IS NOT TRANSFORMED.  Its `delta.npy` is a limma-voom log-fold-change
(`data/external/perturbmatch/meta.json`), not `log1p(TS p_pert) - log1p(TS p_ctrl)`, so there
is no `p^pert` to re-express and no law-space form of it that is not an invention.  It enters
every arm exactly as it enters production, which is also what keeps the identity row exact.
`plolo213.py --no-pm` prices what that impurity is worth.

    uv run python prep213.py --laws all
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
PREP = WORK / "prep"
PREP.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"

from plolo import E062, E064, OFFICIAL, PANEL, PANEL_SET, PM_PATH, DONOR_PATHS  # noqa: E402
import law213 as L  # noqa: E402

CHUNK = 256
EFFECT_SCALE = plolo.EFFECT_SCALE
SOURCES = ("k562_gwps", "orion_hct116", "orion_hek293t")


def bank_for(name):
    if name == "perturbmatch":
        return E064.load_perturbmatch(PM_PATH)
    return E062.Bank.load(name, DONOR_PATHS[name])


class DonorBlock:
    """A donor's mapped columns held once in RAM, so a law sweep costs no extra I/O."""

    def __init__(self, name: str, load_block: bool = True):
        self.name = name
        self.bank = bank_for(name)
        self.amap = E062.map_bank_to_axis(self.bank, OFFICIAL)
        cols = self.amap.all_source_columns
        d = np.load(self.bank.delta_path, mmap_mode="r")
        self.n = d.shape[0]
        self.block = None
        if load_block:
            self.block = np.empty((self.n, cols.size), dtype=np.float32)
            for lo in range(0, self.n, 2048):
                hi = min(lo + 2048, self.n)
                self.block[lo:hi] = d[lo:hi, cols]
        del d
        # the donor's own basal abundance on those same columns
        b = np.asarray(np.load(self.bank.b_ctrl_path), dtype=np.float64)
        self.p_all = np.expm1(b[cols]) / L.TS
        # position of each unique source column inside `all_source_columns`
        pos = {int(s): i for i, s in enumerate(cols.tolist())}
        self.uniq_in_all = np.asarray(
            [pos[int(s)] for s in self.amap.unique_source_columns], dtype=np.int64)
        self.p_uniq = self.p_all[self.uniq_in_all]

    def ensure_block(self) -> None:
        if self.block is not None:
            return
        cols = self.amap.all_source_columns
        d = np.load(self.bank.delta_path, mmap_mode="r")
        self.block = np.empty((self.n, cols.size), dtype=np.float32)
        for lo in range(0, self.n, 2048):
            hi = min(lo + 2048, self.n)
            self.block[lo:hi] = d[lo:hi, cols]
        del d

    def transform(self, rows: slice | np.ndarray, law, p) -> np.ndarray:
        # `p` stays 1-D: every law broadcasts (n_genes,) against (n_rows, n_genes), and the
        # rank laws key their cached frame on the 1-D control vector.
        return law.to_law(np.asarray(self.block[rows], dtype=np.float64), p)

    def stats(self, law) -> tuple[np.ndarray, float]:
        """`mean_all` over all rows and the median centred row norm, in law space.

        At the incumbent's own parameter the law is the identity on `delta`, so the two
        streams are delegated to `E062`'s own functions on the memmap.  That is not a
        special case in the science -- it is what "the incumbent is an interior point of the
        family" means operationally -- and it makes the identity row BIT-exact instead of
        1.3e-15 off.  The general path is measured against it in `gate213.py` (falsifier 1):
        the two agree to 1.3e-15 on `pre`, an allocation-alignment artefact of
        `np.linalg.norm` on differently-aligned buffers, 30,000x below the 4.13e-08 float32
        round-trip 204 itself asserts.
        """
        if not law.is_incumbent and self.block is None:
            self.ensure_block()
        if law.is_incumbent:
            mean_all = E062.stream_column_mean(
                self.bank, self.amap.all_source_columns, chunk_rows=CHUNK, max_rss_gib=40.0)
            raw_median, _, _ = E062.centered_statistics(
                self.bank, self.amap.all_source_columns, mean_all, need_gene_scale=False,
                chunk_rows=CHUNK, max_rss_gib=40.0)
            return mean_all, float(raw_median)
        total = np.zeros(self.block.shape[1], dtype=np.float64)
        for lo in range(0, self.n, CHUNK):
            hi = min(lo + CHUNK, self.n)
            total += self.transform(slice(lo, hi), law, self.p_all).sum(axis=0, dtype=np.float64)
        mean_all = total / self.n
        norms = np.empty(self.n, dtype=np.float64)
        for lo in range(0, self.n, CHUNK):
            hi = min(lo + CHUNK, self.n)
            norms[lo:hi] = np.linalg.norm(
                self.transform(slice(lo, hi), law, self.p_all) - mean_all, axis=1)
        return mean_all, float(np.median(norms)) or 1.0


def cache_path(donor: str, law_key: str) -> Path:
    return PREP / f"{donor}--{law_key}.npz"


def build(donor: str, law, block: DonorBlock | None = None) -> dict:
    path = cache_path(donor, law.key)
    if path.exists():
        z = np.load(path)
        return dict(mean_all=z["mean_all"], raw_median=float(z["raw_median"]))
    own = block is None
    if own:
        block = DonorBlock(donor)
    t0 = time.time()
    mean_all, raw_median = block.stats(law)
    np.savez(path, mean_all=mean_all, raw_median=np.float64(raw_median))
    print(f"  [{donor}/{law.key}] raw_median {raw_median:.6f} ({time.time()-t0:.0f}s)",
          flush=True)
    if own:
        del block
    return dict(mean_all=mean_all, raw_median=raw_median)


def prepared(donor: str, law, targets, panel_positions, block: DonorBlock | None = None):
    """The law-space analogue of `E062.PreparedSource`, for `targets` in the given order."""
    st = build(donor, law, block)
    blk = block if block is not None else DonorBlock(donor)
    mean_all, raw_median = st["mean_all"], st["raw_median"]
    unique_mean = mean_all[blk.uniq_in_all]
    mean_on_axis = np.zeros(len(OFFICIAL), dtype=np.float64)
    mean_on_axis[blk.amap.unique_target_positions] = unique_mean

    rows, avail = [], np.zeros(len(targets), dtype=bool)
    idx = []
    for i, t in enumerate(targets):
        r = blk.bank.pert_pos.get(t)
        if r is not None:
            idx.append((i, r))
            avail[i] = True
    effects = np.zeros((len(targets), len(OFFICIAL)), dtype=np.float64)
    if idx:
        take = np.asarray([r for _, r in idx], dtype=np.int64)
        if blk.block is not None:
            h = np.asarray(blk.block[take][:, blk.uniq_in_all], dtype=np.float64)
        else:                                     # stream the 300 rows straight off disk
            d = np.load(blk.bank.delta_path, mmap_mode="r")
            usc = blk.amap.unique_source_columns
            h = np.empty((take.size, usc.size), dtype=np.float64)
            for k, r in enumerate(take.tolist()):
                h[k] = np.asarray(d[r, usc], dtype=np.float64)
            del d
        g = law.to_law(h, blk.p_uniq)
        vals = (g - unique_mean) * (EFFECT_SCALE / raw_median)
        for k, (i, _) in enumerate(idx):
            effects[i, blk.amap.unique_target_positions] = vals[k]
    effects[:, panel_positions] = 0.0
    mean_on_axis[panel_positions] = 0.0
    del rows
    return E062.PreparedSource(
        name=blk.bank.name, effects=effects, target_available=avail,
        mean_on_axis=mean_on_axis, raw_median_norm=raw_median,
        mapped_gene_count=int(blk.amap.unique_target_positions.size),
        duplicate_symbols=blk.amap.duplicate_symbols)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["all"])
    ap.add_argument("--donors", nargs="*", default=list(SOURCES))
    args = ap.parse_args()
    laws = L.all_laws()
    keys = list(laws) if args.laws == ["all"] else args.laws
    done = {}
    for donor in args.donors:
        todo = [k for k in keys if not cache_path(donor, k).exists()]
        if not todo:
            print(f"[{donor}] all {len(keys)} laws cached", flush=True)
            continue
        print(f"[{donor}] loading block for {len(todo)} laws", flush=True)
        blk = DonorBlock(donor)
        print(f"[{donor}] block {blk.block.shape} "
              f"{blk.block.nbytes/2**30:.2f} GiB  uniq {blk.uniq_in_all.size}", flush=True)
        for k in todo:
            build(donor, laws[k], blk)
        done[donor] = len(todo)
        del blk
        import gc
        gc.collect()
    (WORK / "prep213.json").write_text(json.dumps(
        {"laws": {k: laws[k].asdict() for k in keys}, "built": done}, indent=1))


if __name__ == "__main__":
    main()
