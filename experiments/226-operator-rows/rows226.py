"""226 Phase 1 — the whitening operator's COSINE-vs-FITTING-ROWS curve.

The shipped operator (`data/donors/whitener_union_pruned.npz`, sha `e2f65ecf`) is fitted on
**29,394** rows x **17,942** columns -- n/p = 1.64.  193 asked whether that is enough and
answered on `pds` using the three operators that happened to be on disk; those three differ in
*composition* as well as size (the prune, the 12,000-row cap), so its size axis is confounded,
and it never reported the **cosine** channel at all.  212 then proved cosine and `pds` move in
OPPOSITE directions along the operator's other axis (gamma), so a flat `pds` size axis is not
evidence of a flat cosine size axis.

This script measures the clean thing: hold composition EXACTLY fixed (stratify the subsample
by donor, so every point has production's donor mix to within one row), vary only `n`, and
read both channels.

Design
------
* Pool = production's own 29,394 rows, reproduced by `scripts/build_union_whitener.py`'s
  selection logic (seed 0, 12,000/donor cap, then 091's common-essential prune).
* The float64 basis is assembled ONCE into a scratch memmap; each point re-accumulates
  `C = B[idx].T @ B[idx] / |idx|` in row blocks, so no subset is ever materialised.
* `f = 1.0` consumes the whole pool at any seed, so it must reproduce the shipped operator
  BIT-FOR-BIT.  That is gate KA-3 and it is checked against the shipped npz, not against a
  number someone wrote down.
* Scoring is 193's `phase2_basis.score()` verbatim in structure: a15 config, `kolf2` fold,
  three challenge frames.  `kolf2` is the ONLY leakage-free fold for a production-composition
  operator (162's rule; all three production donors are in the basis).

Endpoints
---------
`pds`      -- macro challenge-frame PDS, the board's 5/6 retrieval channel.
`cos`      -- macro per-target transfer cosine (212's functional).
`cos_agg`  -- energy-weighted aggregate cosine `sqrt(1 - sse)` (197's functional).
AGENTS.md's 2026-09-04 warning: these two are NOT interchangeable and neither converts to the
board by a scalar.  Both are reported; nothing here is multiplied by a transfer factor.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import plolo                                                       # noqa: E402
from plolo import build_fold, deliver, endpoints, solve_fraction   # noqa: E402
from run_arms import CONFIGS                                       # noqa: E402
from build_union_whitener import DonorView, drop_pages             # noqa: E402

SPAN_PDS = 0.452          # 162's board span for pds (193 uses the same)
PROD_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
PROD_OP = ROOT / "data/donors/whitener_union_pruned.npz"
DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")

# Phase 2 sources that contribute ZERO rows to the shipped operator (222 falsifier 1 proved
# this for H1: the operator's sha was unchanged when H1 joined the stack, because H1 enters
# as a RESIDUAL source).  `kolf2` is deliberately absent and must stay absent -- it is the
# evaluation fold, and putting it in the basis is the leak this experiment exists to avoid.
EXTRA_ROOTS = {
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    "rpe1": ROOT / "data/donors/rpe1",
    "h1_2025": ROOT / "data/donors/h1_2025",
    "igvf-telohaec": ROOT / "data/external/igvf-telohaec",
    "perturbmatch": ROOT / "data/external/perturbmatch",
}
FORBIDDEN_IN_BASIS = {"kolf2"}


def rss_gib() -> float:
    kb = subprocess.run(["ps", "-o", "rss=", "-p", str(os.getpid())],
                        capture_output=True, text=True, check=True).stdout.strip()
    return int(kb) / 2**20


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def sha256_file(path: Path, buf: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(buf):
            h.update(chunk)
    return h.hexdigest()


def cd4_cached(targets):
    """193's cache, re-pointed at this experiment's work dir."""
    key = hashlib.sha256("\n".join(targets).encode()).hexdigest()[:16]
    path = HERE / "work" / f"cd4-{key}.npz"
    if path.exists():
        z = np.load(path)
        return plolo.E062.CD4Prepared(
            mode=str(z["mode"]), mean_effect=z["mean_effect"], gene_votes=z["gene_votes"],
            target_condition_count=z["tcc"], condition_details=())
    c = plolo.cd4_for(targets)
    np.savez(path, mode=c.mode, mean_effect=np.asarray(c.mean_effect),
             gene_votes=np.asarray(c.gene_votes),
             tcc=np.asarray(c.target_condition_count))
    return c


def frames_for(fold):
    out = []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        out.append((f"chal{ctx}", b, (cpm > plolo.MIN_CPM) & (b > 0)))
    return out


# ---------------------------------------------------------------- the row pool
class Pool:
    """Production's basis rows, assembled once into a float64 scratch memmap."""

    def __init__(self, bundle: Path, donor_root: Path, prune: str, cap: int, seed: int,
                 scratch: Path, chunk_rows: int, max_rss: float,
                 extra_donors: tuple[str, ...] = (), extra_roots: dict | None = None,
                 pin_axis: bool = True, pin_renorm: bool = True):
        target_genes = pd.read_csv(bundle / "gene_names.csv")["gene_name"].to_numpy()
        panel = pd.read_csv(bundle / "pert_counts.csv")["target_gene"].to_numpy()
        tindex = {g: i for i, g in enumerate(target_genes)}
        self.target_genes = target_genes
        self.panel = panel
        names = list(DONOR_ORDER) + list(extra_donors)
        self.names = names

        roots = dict(extra_roots or {})
        donors = {}
        for n in names:
            donors[n] = DonorView(n, roots.get(n, donor_root / n), target_genes, tindex,
                                  chunk_rows, max_rss)
            print(f"  {n:16} {donors[n].n_rows:6d} x {donors[n].n_genes:6d} | "
                  f"{len(donors[n].tcols):6d} on axis | scale {donors[n].scale:.6f}",
                  flush=True)
        self.donors = donors

        # `MultiDonorModel._renorm`: the reference is the MEDIAN donor scale, so adding a
        # donor MOVES it and rescales every production row too.  That would confound "more
        # rows" with "all the old rows changed", so `pin_renorm` holds the reference at the
        # three production donors' median -- the extra rows enter on production's own scale
        # and nothing already in the basis moves by a single bit.
        ref_names = list(DONOR_ORDER) if pin_renorm else names
        ref = float(np.median([donors[n].scale for n in ref_names]))
        self.renorm = {n: ref / (donors[n].scale or 1.0) for n in names}
        self.ref_scale = ref
        print(f"  reference scale {ref:.6f} (over {ref_names}) | renorm " +
              " ".join(f"{n}={self.renorm[n]:.6f}" for n in names), flush=True)

        # Same argument for the axis: `union(donor coverage) - panel` GROWS when a donor is
        # added, and a wider operator is a different operator, not a better-sampled one.
        # `pin_axis` freezes the axis at production's 17,942, so an extra donor contributes
        # ROWS only; genes it does not cover enter that row as exact zeros.
        axis_names = list(DONOR_ORDER) if pin_axis else names
        union: set[str] = set()
        for n in axis_names:
            union |= set(donors[n].target_symbols)
        union -= set(panel)
        self.wcols = np.array(sorted(tindex[g] for g in union), dtype=np.int64)
        self.pin_axis, self.pin_renorm = pin_axis, pin_renorm
        # Only the bit-gate run takes production's single-syrk path (which needs the basis
        # mmap'ed and costs 4.2 GiB of resident file pages).  Every other run streams, and
        # the f = 1.0 point is re-measured through the streaming path so the blocked-vs-syrk
        # arithmetic difference is itself measured rather than assumed negligible.
        self.exact_full = False
        print(f"  operator axis: {len(self.wcols)} genes (over {axis_names})", flush=True)
        for n in names:
            if n not in axis_names:
                cov = len(set(donors[n].target_symbols) & set(target_genes[self.wcols]))
                print(f"    {n}: covers {cov}/{len(self.wcols)} of the pinned axis "
                      f"({100*cov/len(self.wcols):.1f}%)", flush=True)

        rng = np.random.default_rng(seed)
        sel: dict[str, np.ndarray] = {}
        for n in names:
            nr = donors[n].n_rows
            sel[n] = (np.sort(rng.choice(nr, cap, False)) if nr > cap else np.arange(nr))
        if prune == "essential":
            ess = {ln.strip() for ln in
                   (donor_root / "k562_essential/perts.txt").read_text().splitlines()
                   if ln.strip()}
            for n in names:
                pidx = {p: i for i, p in enumerate(donors[n].perts)}
                row_name: list[str | None] = [None] * donors[n].n_rows
                for pert, i in pidx.items():
                    if 0 <= i < len(row_name):
                        row_name[i] = pert
                rows = sel[n]
                is_ess = np.array([row_name[i] in ess for i in rows], dtype=bool)
                sel[n] = rows[~is_ess]
        self.selection = sel
        self.counts = {n: len(sel[n]) for n in names}
        self.total = sum(self.counts.values())
        print("  pool rows " + " + ".join(f"{n}={self.counts[n]}" for n in names) +
              f" = {self.total}", flush=True)

        # donor label per pool row, for stratified subsampling
        self.donor_of_row = np.concatenate(
            [np.full(self.counts[n], i, dtype=np.int8) for i, n in enumerate(names)])

        scratch.mkdir(parents=True, exist_ok=True)
        self.basis_path = scratch / "pool_basis.npy"
        B = np.lib.format.open_memmap(self.basis_path, mode="w+", dtype=np.float64,
                                      shape=(self.total, len(self.wcols)))
        cur = 0
        for n in names:
            donors[n].block(sel[n], len(target_genes), self.wcols, self.renorm[n],
                            B, cur, chunk_rows, max_rss)
            cur += self.counts[n]
            print(f"  basis: {n} written ({self.counts[n]}), RSS {rss_gib():.2f} GiB",
                  flush=True)
        B.flush()
        drop_pages(B)
        del B
        gc.collect()

    def strat_subsample(self, frac: float, seed: int) -> np.ndarray:
        """Rows for a subsample that keeps each donor's SHARE of the pool fixed."""
        rng = np.random.default_rng(seed)
        out, off = [], 0
        for n in self.names:
            k = self.counts[n]
            take = int(round(frac * k))
            take = max(1, min(k, take))
            out.append(off + np.sort(rng.choice(k, take, False)))
            off += k
        return np.concatenate(out)

    def covariance(self, idx: np.ndarray, chunk: int, permute_genes: int | None = None,
                   max_rss: float = 11.0) -> np.ndarray:
        """`B[idx].T @ B[idx] / |idx|`, accumulated in row blocks.

        `permute_genes` is falsifier 2: draw the same number of rows but permute each row's
        gene labels independently, destroying the covariance structure while preserving the
        per-row marginal distribution and the sample size exactly.
        """
        p = len(self.wcols)
        if (self.exact_full and permute_genes is None and len(idx) == self.total
                and np.array_equal(idx, np.arange(self.total))):
            # The whole pool: take production's own arithmetic (one syrk over the memmap),
            # so the f = 1.0 point is bit-comparable to the shipped operator rather than
            # merely close.  Blocked accumulation would reorder the sums.
            B = np.load(self.basis_path, mmap_mode="r")
            try:
                C = (B.T @ B) / B.shape[0]
            finally:
                drop_pages(B)
                del B
                gc.collect()
            return C
        # Deliberately NOT mmap here.  Fancy-indexing a memmap pulls every touched page into
        # the process's RSS and they stay there until the kernel evicts them, so a blocked
        # pass over a 4.2 GiB basis adds 4.2 GiB of file pages on top of the 2.6 GiB `C` --
        # which is exactly how the first run of this script hit 12.3 GiB.  Plain `read` goes
        # through the page cache instead, where it does not count against the process.
        C = np.zeros((p, p), dtype=np.float64)
        prng = np.random.default_rng(permute_genes) if permute_genes is not None else None
        rowbytes = p * 8
        with open(self.basis_path, "rb") as fh:
            np.lib.format.read_magic(fh)
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(fh)
            assert shape == (self.total, p) and dtype == np.dtype("float64") and not fortran
            data0 = fh.tell()
            buf = np.empty((chunk, p), dtype=np.float64)
            for lo in range(0, len(idx), chunk):
                rows = idx[lo:lo + chunk]
                for j, r in enumerate(rows):
                    fh.seek(data0 + int(r) * rowbytes)
                    fh.readinto(memoryview(buf[j]).cast("B"))
                blk = buf[:len(rows)]
                if prng is not None:
                    for r in range(blk.shape[0]):
                        prng.shuffle(blk[r])
                C += blk.T @ blk
                if rss_gib() > max_rss:
                    raise MemoryError(f"RSS {rss_gib():.2f} GiB at covariance")
            del buf
        gc.collect()
        C /= len(idx)
        return C


def eigensystem(C: np.ndarray, floor_rel: float = 1e-12):
    """`Whitener.fit`'s tail, inlined so `C` can be released as early as possible."""
    ev, V = np.linalg.eigh(C)
    ev = np.maximum(ev, ev.max() * floor_rel)
    return ev, V


# ------------------------------------------------------- per-target detail (F3, F4)
def _unit(A):
    n = np.linalg.norm(A, axis=1, keepdims=True)
    return A / np.where(n > 0, n, 1.0)


def detail(P, T, derange: np.ndarray):
    """Per-target cosine and per-target pds, plus the WRONG-RECIPIENT cosine.

    `derange` is a fixed derangement of the target index (no fixed points), so
    `cos_wrong[i] = cos(prediction_i, truth_{derange[i]})` is the same functional evaluated
    on a deliberately mis-paired recipient.  Falsifier 3 (184's null): if `cos_wrong` tracks
    `cos_own` as `n` grows, the operator is moving the whole space, not the pairing.
    """
    A, B = _unit(P), _unit(T)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    pds_t = (beaten + 0.5 * ties) / (n - 1)
    wrong = S[np.arange(n), derange].copy()
    return own, pds_t.astype(np.float64), wrong


# ------------------------------------------------------------------- the scorer
def score(fold, cfg, evals, evecs, cd4, frames, truth, wcols, derange=None):
    """193 `phase2_basis.score`, with the operator passed in memory instead of by path."""
    assert evecs.shape[0] == len(wcols) == len(evals), (evecs.shape, len(wcols), len(evals))
    pre = np.asarray(fold["pre"], dtype=np.float64) * cfg["alpha"]
    sub = pre[:, wcols]
    before = float(np.median(np.linalg.norm(sub, axis=1))) or 1.0
    w = (evals + cfg["lam"]) ** (-cfg["gamma"] / 2.0)
    t = ((sub @ evecs) * w) @ evecs.T
    after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
    whitened = pre.copy()
    whitened[:, wcols] = t * (before / after)
    del t, sub
    if "cd4_f" in cfg:
        frac, _ = solve_fraction(whitened, cd4, cfg["cd4_f"])
    else:
        frac = cfg["cd4_frac"]
    final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, wcols)
    final[:, fold["panel_positions"]] = 0.0
    del whitened
    per = {}
    for fname, b, gate in frames:
        mat = deliver(final, b, gate, cfg["mass"])
        sub = mat[:, fold["scored"]]
        e = endpoints(sub, truth)
        row = dict(pds=e["pds"], cos=e["cos"], sgn_w=e["sgn_w"], sse=e["sse"],
                   cos_agg=float(np.sqrt(max(0.0, 1.0 - e["sse"]))))
        if derange is not None:
            own, pds_t, wrong = detail(sub, truth, derange)
            row["cos_per_target"] = own.tolist()
            row["pds_per_target"] = pds_t.tolist()
            row["cos_wrong"] = float(wrong.mean())
            row["cos_wrong_per_target"] = wrong.tolist()
        per[fname] = row
        del mat, sub
    ei = endpoints(final[:, fold["scored"]], truth)
    del final
    gc.collect()
    return dict(
        chal_pds=float(np.mean([v["pds"] for v in per.values()])),
        chal_cos=float(np.mean([v["cos"] for v in per.values()])),
        chal_cos_agg=float(np.mean([v["cos_agg"] for v in per.values()])),
        chal_cos_wrong=(float(np.mean([v["cos_wrong"] for v in per.values()]))
                        if derange is not None else None),
        intended_pds=ei["pds"], intended_cos=ei["cos"],
        intended_cos_agg=float(np.sqrt(max(0.0, 1.0 - ei["sse"]))),
        per_frame=per, cd4_fraction=float(frac), max_ev=float(evals.max()),
        n_floored=int((evals <= evals.max() * 1e-12 * (1 + 1e-9)).sum()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", default="kolf2")
    ap.add_argument("--config", default="a15")
    ap.add_argument("--fracs", default="0.25,0.375,0.5,0.625,0.75,0.875,1.0")
    ap.add_argument("--seeds", default="1")
    ap.add_argument("--extra-seed-fracs", default="",
                    help="fractions to repeat at every seed in --seeds (noise band)")
    ap.add_argument("--prune", default="essential", choices=["none", "essential"])
    ap.add_argument("--cap", type=int, default=12000)
    ap.add_argument("--pool-seed", type=int, default=0)
    ap.add_argument("--extra-donors", default="",
                    help="Phase 2: sources that contribute NO row to the shipped operator")
    ap.add_argument("--no-pin-axis", action="store_true",
                    help="let an extra donor widen the operator axis (breaks the fold guard)")
    ap.add_argument("--no-pin-renorm", action="store_true",
                    help="let an extra donor move the median-scale reference")
    ap.add_argument("--permute-null", action="store_true",
                    help="falsifier 2: gene-label-permuted refit at every fraction")
    ap.add_argument("--chunk-rows", type=int, default=1024)
    ap.add_argument("--max-rss-gib", type=float, default=11.0)
    ap.add_argument("--bit-gate", action="store_true",
                    help="assert the frac=1.0 eigensystem is bit-identical to the shipped npz")
    ap.add_argument("--out", default="work/rows226.json")
    args = ap.parse_args()

    t0 = time.time()
    scratch = HERE / "work" / "scratch226"
    extra = tuple(x for x in args.extra_donors.split(",") if x)
    # falsifier 1, enforced in code: the held fold can never enter the fitting basis.
    bad = (set(extra) | set(DONOR_ORDER)) & (FORBIDDEN_IN_BASIS | {args.fold})
    if bad:
        raise SystemExit(f"LEAK ABORT: {sorted(bad)} would be both fitted on and scored on")
    unknown = set(extra) - set(EXTRA_ROOTS)
    if unknown:
        raise SystemExit(f"unknown extra donors {sorted(unknown)}")
    print("=== pool ===", flush=True)
    pool = Pool(ROOT / "data/challenge", ROOT / "data/donors", args.prune, args.cap,
                args.pool_seed, scratch, args.chunk_rows, args.max_rss_gib, extra,
                extra_roots={n: p for n, p in EXTRA_ROOTS.items() if n in extra},
                pin_axis=not args.no_pin_axis, pin_renorm=not args.no_pin_renorm)
    pool.exact_full = bool(args.bit_gate)

    print("\n=== fold ===", flush=True)
    fold = build_fold(args.fold)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    wcols = np.asarray(fold["wcols"])
    print(f"  fold {args.fold}: {len(fold['targets'])} targets, "
          f"{len(fold['scored'])} scored genes, wcols {len(wcols)}", flush=True)
    if len(wcols) != len(pool.wcols) or not np.array_equal(wcols, pool.wcols):
        raise SystemExit(
            f"PANEL-IDENTITY ABORT: fold axis ({len(wcols)}) != operator axis "
            f"({len(pool.wcols)}) -- see falsifier 5")
    cd4 = cd4_cached(fold["targets"])
    frames = frames_for(fold)
    cfg = CONFIGS[args.config]

    # falsifier 3's fixed derangement of the target index -- the SAME permutation at every
    # row count, so `cos_wrong(n)` is a curve on identical pairings and only `n` differs.
    nt = len(fold["targets"])
    drng = np.random.default_rng(20260904)
    derange = np.arange(nt)
    while True:
        derange = drng.permutation(nt)
        if not np.any(derange == np.arange(nt)):
            break
    assert not np.any(derange == np.arange(nt))

    fracs = [float(x) for x in args.fracs.split(",") if x]
    seeds = [int(x) for x in args.seeds.split(",") if x]
    extra_fracs = [float(x) for x in args.extra_seed_fracs.split(",") if x]

    points = []
    for f in fracs:
        points.append((f, seeds[0], False))
    for f in extra_fracs:
        for s in seeds[1:]:
            points.append((f, s, False))
    if args.permute_null:
        for f in fracs:
            points.append((f, seeds[0], True))

    out = {
        "fold": args.fold, "config": args.config, "cfg": {k: v for k, v in cfg.items()},
        "pool": {"total": pool.total, "counts": pool.counts, "p": int(len(pool.wcols)),
                 "prune": args.prune, "cap": args.cap, "pool_seed": args.pool_seed,
                 "extra_donors": list(extra), "renorm": pool.renorm,
                 "ref_scale": pool.ref_scale, "pin_axis": pool.pin_axis,
                 "pin_renorm": pool.pin_renorm},
        "n_targets": len(fold["targets"]), "scored_genes": int(len(fold["scored"])),
        "targets": list(fold["targets"]), "derangement": derange.tolist(),
        "fold_audit": fold["audit"], "rows": [],
    }
    outp = HERE / args.out

    for frac, seed, perm in points:
        idx = pool.strat_subsample(frac, seed)
        n = len(idx)
        tag = f"f{frac:g}_s{seed}" + ("_perm" if perm else "")
        ts = time.time()
        C = pool.covariance(idx, args.chunk_rows,
                            permute_genes=(seed * 1000 + int(frac * 1000)) if perm else None,
                            max_rss=args.max_rss_gib)
        tc = time.time() - ts
        te = time.time()
        evals, evecs = eigensystem(C)
        del C
        gc.collect()
        teig = time.time() - te

        gate = None
        if args.bit_gate and frac == 1.0 and not perm:
            with np.load(PROD_OP) as z:
                pe = np.asarray(z["evals"], dtype=np.float64)
                pv = np.asarray(z["evecs"])       # float32 on disk
            same_ev = bool(np.array_equal(pe, evals))
            same_vec = bool(np.array_equal(pv, evecs.astype(np.float32)))
            gate = dict(evals_bit_identical=same_ev, evecs_bit_identical=same_vec,
                        max_abs_ev=float(np.max(np.abs(pe - evals))),
                        max_abs_vec=float(np.max(np.abs(pv - evecs.astype(np.float32)))))
            del pe, pv
            gc.collect()
            print(f"  KA-3 bit gate: evals {same_ev}  evecs {same_vec}", flush=True)

        r = score(fold, cfg, evals, evecs, cd4, frames, truth, wcols, derange=derange)
        del evals, evecs
        gc.collect()
        r.update(frac=frac, seed=seed, permuted=perm, n_rows=int(n), tag=tag,
                 t_cov_s=round(tc, 1), t_eigh_s=round(teig, 1),
                 rss_peak_gib=round(peak_gib(), 2))
        if gate:
            r["bit_gate"] = gate
        out["rows"].append(r)
        print(f"[{tag}] n={n:6d}  n/p={n/len(wcols):.3f}  pds {r['chal_pds']:.6f}  "
              f"cos {r['chal_cos']:.7f}  cos_agg {r['chal_cos_agg']:.7f}  "
              f"floored {r['n_floored']}  ({tc:.0f}s cov + {teig:.0f}s eigh, "
              f"peak {peak_gib():.2f} GiB)", flush=True)
        outp.parent.mkdir(parents=True, exist_ok=True)
        outp.write_text(json.dumps(out, indent=1))

    try:
        pool.basis_path.unlink(missing_ok=True)
        scratch.rmdir()
    except OSError:
        pass
    print(f"\nwrote {outp}  ({time.time()-t0:.0f}s total, peak {peak_gib():.2f} GiB)")


if __name__ == "__main__":
    main()
