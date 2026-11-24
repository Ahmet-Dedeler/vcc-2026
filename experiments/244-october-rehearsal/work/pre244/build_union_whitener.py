"""Rebuild the production union-axis whitening operator (`data/donors/whitener_union.npz`).

Provenance
----------
The shipped operator was **not** produced by a checked-in script. It was fitted on
2026-08-22 19:09 (+03) by an inline heredoc, recovered verbatim from the Claude Code
session transcript `d207ce74-...jsonl` (line 1973, ts 2026-08-22T15:59:42Z); its stdout
is still on disk as `data/donors/fit_union.log`:

    union axis: 17942 genes (intersection was 7,312)
    fitting (33787, 17942)
    saved

That heredoc was::

    spec  = ChallengeSpec.load("data/challenge")
    model = MultiDonorModel({n: Donor.load(p) for n, p in paths.items()}, spec.genes,
                            gene_norm=False, eb=False)
    union = set()
    for n in model.names:
        union |= set(model.target_genes[model.tcols[n]])
    union -= set(spec.perts)
    ti    = {g: i for i, g in enumerate(model.target_genes)}
    wcols = np.array(sorted(ti[g] for g in union), dtype=np.int64)
    rng   = np.random.default_rng(0)
    rows  = []
    for n in model.names:
        S = np.zeros((model.s[n].shape[0], len(model.target_genes)), np.float32)
        S[:, model.tcols[n]] = model.s[n][:, model.dcols[n]]
        M = S[:, wcols]
        if M.shape[0] > 12000:
            M = M[np.sort(rng.choice(M.shape[0], 12000, False))]
        rows.append(model._renorm[n] * M.astype(np.float64)); del S, M
    B = np.vstack(rows); del rows
    W = Whitener.fit(B); W.save("data/donors/whitener_union.npz")

This module reproduces that computation **bit-for-bit** while keeping peak RSS around
10 GiB instead of the ~30 GiB the original needed. Every arithmetic step is preserved
exactly; only the *storage* is streamed:

* donor `delta` matrices are `mmap`ed, never all resident;
* `s = delta - m` is never materialised -- the basis is assembled a row-block at a time,
  and `(delta - m)[:, dcols]` is evaluated as `delta[:, dcols] - m[dcols]`, which is the
  same float32 subtraction on the same pairs of numbers;
* row norms are reduced per row, so blocking rows cannot change any sum;
* the float64 basis lives in a scratch `.npy` memmap, and `B.T @ B` is handed that memmap.
  NumPy dispatches `X.T @ X` to a symmetric-rank-k update whose result is bit-identical
  whether the operand is file-backed or resident (verified).

Pruned variant (`--basis-prune essential`)
------------------------------------------
Same rule as `experiments/031-somatic/run.py`: drop every basis row whose perturbation
name appears in the common-essential list (`data/donors/k562_essential/perts.txt`).
Pruning is applied *after* the 12,000-row-per-donor cap, so the pruned basis is a strict
subset of the production basis and the two operators differ in nothing but those rows.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import os
import json
import resource
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vcc2026.whiten import Whitener  # noqa: E402

DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")


# --------------------------------------------------------------------------------------
# memory discipline
# --------------------------------------------------------------------------------------
def rss_gib() -> float:
    """*Current* RSS in GiB, via `ps`.

    Deliberately not `ru_maxrss`: paging a mmap'ed donor matrix through the process
    inflates the high-water mark with file pages the kernel is free to evict, so the
    peak figure would trip the guard on memory we are not really holding.
    """
    import subprocess
    try:
        kb = subprocess.run(["ps", "-o", "rss=", "-p", str(os.getpid())],
                            capture_output=True, text=True, check=True).stdout.strip()
        return int(kb) / 2**20
    except Exception:
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / 2**30 if sys.platform == "darwin" else peak / 2**20


def peak_gib() -> float:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 2**30 if sys.platform == "darwin" else peak / 2**20


def guard(limit: float, where: str) -> None:
    used = rss_gib()
    if used > limit:
        raise MemoryError(f"RSS {used:.2f} GiB exceeded the {limit:.1f} GiB budget "
                          f"at {where}")


def drop_pages(mm: np.memmap | None) -> None:
    """Return a memmap's file-backed pages to the OS so they stop counting against RSS."""
    if mm is None:
        return
    try:
        mm._mmap.close()  # type: ignore[union-attr]
    except Exception:
        pass


# --------------------------------------------------------------------------------------
# donor bookkeeping -- the parts of MultiDonorModel the operator actually uses
# --------------------------------------------------------------------------------------
class DonorView:
    """`MultiDonorModel`'s per-donor state, computed without materialising `s`.

    `gene_norm=False` and `eb=False` in the production call, so `s` is exactly
    `delta - delta.mean(axis=0)` and nothing else touches it.
    """

    def __init__(self, name: str, path: Path, target_genes: np.ndarray,
                 tindex: dict[str, int], chunk_rows: int, max_rss: float):
        self.name = name
        self.path = path
        self.delta_path = path / "delta.npy"
        self.genes = np.asarray((path / "genes.txt").read_text().split("\n"))
        self.perts = np.asarray((path / "perts.txt").read_text().split("\n"))

        delta = np.load(self.delta_path, mmap_mode="r")
        self.n_rows, self.n_genes = delta.shape

        # `Donor.decompose()`: m = delta.mean(axis=0), a float32 strided reduction over
        # rows. Run it straight on the memmap -- the iteration order, and therefore every
        # rounding decision, is identical to the resident case.
        self.m = np.asarray(delta.mean(axis=0))
        guard(max_rss, f"{name} mean")

        # Column map onto the official axis, in donor-gene order (duplicates: last wins,
        # exactly as MultiDonorModel's scatter does).
        tcols, dcols = [], []
        for j, g in enumerate(self.genes):
            t = tindex.get(g)
            if t is not None:
                tcols.append(t)
                dcols.append(j)
        self.tcols = np.asarray(tcols, dtype=np.int64)
        self.dcols = np.asarray(dcols, dtype=np.int64)

        # `_scale[n] = median(||s[:, dcols]||_2)`. `np.linalg.norm(..., axis=1)` reduces
        # inside each row, so blocking over rows leaves every sum untouched.
        md = self.m[self.dcols]
        norms = np.empty(self.n_rows, dtype=np.float32)
        for lo in range(0, self.n_rows, chunk_rows):
            hi = min(self.n_rows, lo + chunk_rows)
            blk = np.asarray(delta[lo:hi][:, self.dcols], dtype=np.float32) - md
            norms[lo:hi] = np.linalg.norm(blk, axis=1)
            del blk
        drop_pages(delta)
        del delta
        gc.collect()
        self.scale = float(np.median(norms))
        self.target_symbols = target_genes[self.tcols]
        guard(max_rss, f"{name} scale")

    def block(self, rows: np.ndarray, n_target: int, wcols: np.ndarray,
              renorm: float, out: np.ndarray, cursor: int,
              chunk_rows: int, max_rss: float) -> None:
        """Write `renorm * S[rows][:, wcols]` (float64) into `out[cursor:...]`.

        `S` is the donor's effect matrix scattered onto the 18,533-gene official axis and
        then gathered onto the operator axis -- the same two indexing steps the original
        performed, done a row-block at a time.
        """
        delta = np.load(self.delta_path, mmap_mode="r")
        md = self.m[self.dcols]
        try:
            for lo in range(0, len(rows), chunk_rows):
                hi = min(len(rows), lo + chunk_rows)
                idx = rows[lo:hi]
                s_blk = np.asarray(delta[idx][:, self.dcols], dtype=np.float32) - md
                S = np.zeros((hi - lo, n_target), dtype=np.float32)
                S[:, self.tcols] = s_blk
                del s_blk
                out[cursor + lo: cursor + hi] = renorm * S[:, wcols].astype(np.float64)
                del S
                guard(max_rss, f"{self.name} basis")
        finally:
            drop_pages(delta)
            del delta
            gc.collect()


# --------------------------------------------------------------------------------------
def sha256_file(path: Path, bufsize: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def axis_sha256(symbols) -> str:
    """sha256 over the operator's ORDERED column/gene axis (experiment 235).

    THE FAILURE THIS EXISTS TO MAKE LOUD (`DEF-RUNBOOK` section 0.1, priced by 212 at
    **-0.0429 overall on `pds` alone**).  The axis below is
    `sorted_by_official_index(union(donor coverage) - panel)`, so it is a FUNCTION OF THE
    PANEL.  Point an operator fitted for one panel at another and every check that exists
    passes -- `n_cols` is **17,942 on all four panels on disk** (235's `axis235.json`), the
    eigenvector shape matches, the file sha256 matches -- while **98.8%** of the coordinates
    pair with the wrong gene.

    The digest is over the ordered SYMBOLS, not over the integer column indices the builder
    already records in `meta['operator']['columns_sha256']`.  Measured, that matters: 199's
    `def-bundle-badaxis` (two symbols transposed in `gene_names.csv`) has a column-index
    digest IDENTICAL to `def-bundle`'s and a symbol digest that differs.  An index digest
    cannot see a re-ordered gene axis; this can.

    Canonical form: newline-joined UTF-8, no trailing newline.  `make_submission_pm_cd4.py`
    carries the same one-line implementation as `operator_axis_sha256` -- a second copy on
    purpose, because the builder verifies the sha256 of every helper it imports and importing
    this module would couple two independently-pinned files.  235's known-answer row 3
    measures the two implementations against each other on a real refit, and its test suite
    checks them on six inputs including a permutation and a non-ASCII symbol.
    """
    return hashlib.sha256("\n".join(map(str, symbols)).encode("utf-8")).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", default="data/challenge")
    ap.add_argument("--donor-root", default="data/donors")
    ap.add_argument("--donors", default=",".join(DONOR_ORDER),
                    help="donor directory names, in the order the basis stacks them")
    ap.add_argument("--basis-per-donor", type=int, default=12000,
                    help="row cap per donor. 12000 is what production used: k562_gwps "
                         "has 9,787 perturbations and is taken whole, the two Orion lines "
                         "are subsampled to 12,000 each -> 33,787 basis rows.")
    ap.add_argument("--seed", type=int, default=0,
                    help="seeds the single Generator that draws every donor's rows, in "
                         "donor order. Donors under the cap consume no random numbers.")
    ap.add_argument("--basis-prune", default="none", choices=["none", "essential"],
                    help="'essential' drops basis rows whose perturbation is in the "
                         "common-essential list (experiment 031 semantics), applied after "
                         "the row cap so the basis stays a subset of the production one.")
    ap.add_argument("--essential-list", default="data/donors/k562_essential/perts.txt")
    ap.add_argument("--out", default=None,
                    help="default: whitener_union.npz, or whitener_union_pruned.npz "
                         "when --basis-prune is not 'none'")
    ap.add_argument("--scratch", default=None, help="directory for the float64 basis memmap")
    ap.add_argument("--chunk-rows", type=int, default=2048)
    ap.add_argument("--max-rss-gib", type=float, default=20.0)
    ap.add_argument("--expect-sha256", default=None,
                    help="abort non-zero unless the written file hashes to this")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the axis, scales and row selection, then stop")
    args = ap.parse_args()

    root = Path(args.donor_root)
    names = [n.strip() for n in args.donors.split(",") if n.strip()]
    out = Path(args.out) if args.out else root / (
        "whitener_union.npz" if args.basis_prune == "none" else "whitener_union_pruned.npz")

    # ---- the official axis ------------------------------------------------------------
    bundle = Path(args.bundle)
    target_genes = pd.read_csv(bundle / "gene_names.csv")["gene_name"].to_numpy()
    panel = pd.read_csv(bundle / "pert_counts.csv")["target_gene"].to_numpy()
    tindex = {g: i for i, g in enumerate(target_genes)}
    print(f"official axis {len(target_genes)} genes | panel {len(panel)} targets", flush=True)

    donors = {}
    for n in names:
        d = DonorView(n, root / n, target_genes, tindex, args.chunk_rows, args.max_rss_gib)
        donors[n] = d
        print(f"  {n:14} {d.n_rows:6d} perts x {d.n_genes:6d} genes | "
              f"{len(d.tcols):6d} on the official axis | scale {d.scale:.6f}", flush=True)

    # ---- donor magnitude equalisation (MultiDonorModel._renorm) ------------------------
    ref = float(np.median([donors[n].scale for n in names]))
    renorm = {n: ref / (donors[n].scale or 1.0) for n in names}
    print(f"reference scale {ref:.6f} | renorm " +
          " ".join(f"{n}={renorm[n]:.6f}" for n in names), flush=True)

    # ---- operator axis: union of donor coverage, minus the panel targets ---------------
    union: set[str] = set()
    for n in names:
        union |= set(donors[n].target_symbols)
    union -= set(panel)
    wcols = np.array(sorted(tindex[g] for g in union), dtype=np.int64)
    axis_symbols = [str(gene) for gene in target_genes[wcols]]
    axis_digest = axis_sha256(axis_symbols)
    panel_digest = axis_sha256(map(str, panel))
    print(f"union axis: {len(wcols)} genes (intersection was 7,312)", flush=True)
    # `n_cols` is 17,942 on every panel measured, so it identifies nothing. This does.
    print(f"axis sha256 {axis_digest}\npanel sha256 {panel_digest}", flush=True)

    # ---- row selection: one Generator, consumed in donor order -------------------------
    rng = np.random.default_rng(args.seed)
    selection: dict[str, np.ndarray] = {}
    for n in names:
        n_rows = donors[n].n_rows
        if n_rows > args.basis_per_donor:
            selection[n] = np.sort(rng.choice(n_rows, args.basis_per_donor, False))
        else:
            selection[n] = np.arange(n_rows)

    # ---- optional prune (experiment 031 semantics) ------------------------------------
    essential: set[str] = set()
    dropped: dict[str, int] = {n: 0 for n in names}
    dropped_names: dict[str, list[str]] = {n: [] for n in names}
    if args.basis_prune != "none":
        essential = {ln.strip() for ln in Path(args.essential_list).read_text().splitlines()
                     if ln.strip()}
        print(f"prune list: {len(essential)} common-essential perturbations", flush=True)
        for n in names:
            # exp031 recovers row->name through pert_index, so a duplicated symbol keeps
            # only its last row; mirror that exactly.
            pert_index = {p: i for i, p in enumerate(donors[n].perts)}
            row_name: list[str | None] = [None] * donors[n].n_rows
            for pert, idx in pert_index.items():
                if 0 <= idx < len(row_name):
                    row_name[idx] = pert
            rows = selection[n]
            is_ess = np.array([row_name[i] in essential for i in rows], dtype=bool)
            dropped[n] = int(is_ess.sum())
            dropped_names[n] = [row_name[i] for i in rows[is_ess]]  # type: ignore[misc]
            selection[n] = rows[~is_ess]

    counts = {n: len(selection[n]) for n in names}
    total = sum(counts.values())
    print("basis rows " + " + ".join(f"{n}={counts[n]}" for n in names) + f" = {total}",
          flush=True)
    if args.basis_prune != "none":
        print("dropped     " + " + ".join(f"{n}={dropped[n]}" for n in names) +
              f" = {sum(dropped.values())}", flush=True)

    if args.dry_run:
        print("dry run: stopping before the eigendecomposition")
        return

    # ---- assemble the float64 basis into a scratch memmap ------------------------------
    scratch = Path(args.scratch) if args.scratch else out.parent / ".whitener_scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    basis_path = scratch / f"basis_{out.stem}.npy"
    B = np.lib.format.open_memmap(basis_path, mode="w+", dtype=np.float64,
                                  shape=(total, len(wcols)))
    cursor = 0
    for n in names:
        donors[n].block(selection[n], len(target_genes), wcols, renorm[n],
                        B, cursor, args.chunk_rows, args.max_rss_gib)
        cursor += counts[n]
        print(f"  basis: {n} written ({counts[n]} rows), RSS {rss_gib():.2f} GiB",
              flush=True)
    B.flush()

    # ---- Whitener.fit, inlined so the basis can be released before the eigensolve -------
    print(f"fitting {B.shape}", flush=True)
    C = (B.T @ B) / B.shape[0]
    guard(args.max_rss_gib, "covariance")
    drop_pages(B)
    del B
    gc.collect()
    basis_path.unlink(missing_ok=True)
    print(f"  covariance {C.shape} {C.dtype}, RSS {rss_gib():.2f} GiB", flush=True)

    evals, evecs = np.linalg.eigh(C)
    del C
    gc.collect()
    evals = np.maximum(evals, evals.max() * 1e-12)
    guard(args.max_rss_gib, "eigendecomposition")
    print(f"  eigh done, RSS {rss_gib():.2f} GiB", flush=True)

    Whitener(evals, evecs).save(out)
    del evals, evecs
    gc.collect()
    try:
        scratch.rmdir()
    except OSError:
        pass

    digest = sha256_file(out)
    meta = {
        "n_cols": int(len(wcols)),
        # The panel this operator was FITTED for, as a hash of the axis it produced. The
        # builder recomputes the same digest from the panel it is building against and
        # refuses to run if the two differ (`load_production_operator`, experiment 235).
        "axis_sha256": axis_digest,
        "panel_sha256": panel_digest,
        "basis": [int(total), int(len(wcols))],
        "donors": names,
    }
    if args.basis_prune != "none":
        meta.update({
            "basis_prune": args.basis_prune,
            "essential_list": args.essential_list,
            "essential_list_size": len(essential),
            "basis_per_donor_cap": args.basis_per_donor,
            "seed": args.seed,
            "rows_kept_per_donor": counts,
            "rows_dropped_per_donor": dropped,
            "rows_dropped_total": int(sum(dropped.values())),
            "dropped_in_panel": sorted(
                set().union(*(set(v) for v in dropped_names.values())) & set(panel.tolist())
            ),
            "sha256": digest,
        })
    out.with_suffix(".json").write_text(json.dumps(meta))
    print(f"saved {out}\n  sha256 {digest}\n  RSS {rss_gib():.2f} GiB "
          f"(process high-water {peak_gib():.2f} GiB)", flush=True)

    if args.expect_sha256:
        ok = digest == args.expect_sha256
        print(("MATCH" if ok else "MISMATCH") + f": expected {args.expect_sha256}")
        raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
