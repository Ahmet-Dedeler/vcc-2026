"""213 — refit the whitening operator IN A LAW'S OWN SPACE, and score with it.

WHY.  The shipped operator's covariance `C = B^T B / n` was estimated from donor effects
expressed as `log1p(TS p_pert) - log1p(TS p_ctrl)`.  A different transport law re-expresses
those same effects, so it changes what the operator is whitening.  Holding the operator fixed
while changing the law is therefore a confounded comparison in one direction, and refitting it
is a confounded comparison in the other (the refit operator is estimated from a different, and
possibly better- or worse-conditioned, basis).  Both are reported.

THE CONSTRUCTION mirrors `scripts/build_union_whitener.py`, which itself reproduces the
2026-08-22 heredoc that made `data/donors/whitener_union_pruned.npz`:

    per donor:  s = g - mean_rows(g)                     (g = the law-space effect)
                map onto the official axis, take `wcols`
                cap at 12,000 rows, seed-0 `rng.choice`, sorted
                drop rows whose perturbation is common-essential (`--basis-prune essential`)
                multiply by `renorm[n] = median_n(scale) / scale[n]`
    stack, C = B^T B / n_rows, eigh, floor eigenvalues at max*1e-12

It is a MIRROR, not a byte-copy of the shipped artifact: `scale[n]` here is the median centred
row norm on `all_source_columns` (what `prep213` already caches per law) rather than
`MultiDonorModel._scale` on `dcols`, and the row cap draws from this script's own generator.
So the identity check is BEHAVIOURAL, not a sha: refitting at the INCUMBENT law and scoring
the fold with it must land on the shipped operator's own PLOLO numbers.  If it does not, the
refit column is not comparable and says so.

    uv run python refit213.py --laws logc20 logc40 --fold kolf2
"""
from __future__ import annotations

import argparse
import gc
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
OPS = WORK / "ops"
OPS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
import score204  # noqa: E402
score204.WORK = plolo.WORK.parent

from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for  # noqa: E402
import law213 as L  # noqa: E402
import prep213 as P  # noqa: E402
import plolo213 as S  # noqa: E402

CAP = 12000
SEED = 0
ESSENTIAL = ROOT / "data/donors/k562_essential/perts.txt"


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30


def build_operator(held: str, law, wcols: np.ndarray, prune: bool = True) -> Path:
    out = OPS / f"op_{held}_{law.key}{'_pruned' if prune else ''}.npz"
    if out.exists():
        return out
    names = FOLDS[held]["sources"]
    ess = {ln.strip() for ln in ESSENTIAL.read_text().splitlines() if ln.strip()} if prune \
        else set()
    blocks, scales, sel = {}, {}, {}
    rng = np.random.default_rng(SEED)
    for n in names:
        blk = P.DonorBlock(n, load_block=False)
        st = P.build(n, law, blk)
        scales[n] = st["raw_median"]
        blocks[n] = blk
        rows = (np.sort(rng.choice(blk.n, CAP, False)) if blk.n > CAP else np.arange(blk.n))
        if prune:
            keep = np.asarray([blk.bank.perts[i] not in ess for i in rows], dtype=bool)
            rows = rows[keep]
        sel[n] = rows
    ref = float(np.median([scales[n] for n in names]))
    renorm = {n: ref / (scales[n] or 1.0) for n in names}
    total = int(sum(len(sel[n]) for n in names))
    print(f"[refit {held}/{law.key}] rows " +
          " + ".join(f"{n}={len(sel[n])}" for n in names) + f" = {total}; "
          f"renorm " + " ".join(f"{n}={renorm[n]:.5f}" for n in names), flush=True)

    scratch = WORK / f".basis_{held}_{law.key}.npy"
    B = np.lib.format.open_memmap(scratch, mode="w+", dtype=np.float64,
                                  shape=(total, wcols.size))
    # official-axis position of each of `wcols`, per donor
    cur = 0
    for n in names:
        blk = blocks[n]
        blk.ensure_block()
        st = P.build(n, law, blk)
        mean_all = st["mean_all"]
        uniq_mean = mean_all[blk.uniq_in_all]
        tpos = blk.amap.unique_target_positions
        # which of `wcols` this donor covers, and where they sit in the donor's unique list
        wpos = {int(c): i for i, c in enumerate(wcols.tolist())}
        dst, src = [], []
        for j, t in enumerate(tpos.tolist()):
            k = wpos.get(int(t))
            if k is not None:
                dst.append(k)
                src.append(j)
        dst = np.asarray(dst, dtype=np.int64)
        src = np.asarray(src, dtype=np.int64)
        rows = sel[n]
        for lo in range(0, rows.size, 1024):
            hi = min(lo + 1024, rows.size)
            g = law.to_law(np.asarray(blk.block[rows[lo:hi]], dtype=np.float64), blk.p_all)
            g = g[:, blk.uniq_in_all] - uniq_mean
            chunk = np.zeros((hi - lo, wcols.size), dtype=np.float64)
            chunk[:, dst] = g[:, src]
            B[cur + lo:cur + hi] = chunk * renorm[n]
        cur += rows.size
        blk.block = None
        gc.collect()
        print(f"    basis {n} written ({rows.size} rows), RSS {rss():.2f} GiB", flush=True)
    B.flush()
    del blocks
    gc.collect()
    t0 = time.time()
    C = (B.T @ B) / B.shape[0]
    del B
    gc.collect()
    scratch.unlink(missing_ok=True)
    print(f"    covariance {C.shape} in {time.time()-t0:.0f}s, RSS {rss():.2f} GiB",
          flush=True)
    t0 = time.time()
    evals, evecs = np.linalg.eigh(C)
    del C
    gc.collect()
    evals = np.maximum(evals, evals.max() * 1e-12)
    print(f"    eigh in {time.time()-t0:.0f}s, eig [{evals.min():.3e}, {evals.max():.4f}], "
          f"RSS {rss():.2f} GiB", flush=True)
    np.savez(out, evals=evals.astype(np.float64), evecs=evecs.astype(np.float32))
    del evals, evecs
    gc.collect()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["logc20", "logc40"])
    ap.add_argument("--fold", default="kolf2")
    ap.add_argument("--no-prune", action="store_true")
    args = ap.parse_args()

    laws = L.all_laws()
    held = args.fold
    fold = build_fold(held)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    cd4 = cd4_for(fold["targets"])
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    shipped_path = Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                        else PROD_OPS / "whitener_union_pruned.npz")

    out_path = WORK / "refit213.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["law"], r["operator"]) for r in rows}

    for key in args.laws:
        law = laws[key]
        srcs = S.law_sources(held, fold, law)
        pre, _ = score204.weighted_baseline(srcs, fold["targets"], "equal",
                                            {n: 1.0 for n, _ in srcs})
        del srcs
        gc.collect()
        for which in ("shipped", "refit"):
            if (held, key, which) in done:
                continue
            path = (shipped_path if which == "shipped"
                    else build_operator(held, law, wcols, prune=not args.no_prune))
            op = Operator(path)
            t0 = time.time()
            e = S.score_law(fold, pre, op, cd4, fold["scored"], truth, law)
            rows.append(dict(fold=held, law=key, operator=which,
                             op_path=str(Path(path).relative_to(ROOT)),
                             max_eig=float(op.evals.max()), min_eig=float(op.evals.min()),
                             n=int(op.evals.size), **e,
                             seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=1))
            print(f"[{held}/{key}/{which}] a15cos {e['a15_cos']:.7f} "
                  f"a15pds {e['a15_pds']:.6f} a23cos {e['a23_cos']:.7f} "
                  f"a23pds {e['a23_pds']:.6f}  max_eig {op.evals.max():.4f}", flush=True)
            del op
            gc.collect()
        del pre
        gc.collect()
    print(f"\npeak RSS {rss():.2f} GiB")


if __name__ == "__main__":
    main()
