"""195 Stage F — same test, with training queries that match the panel query exactly.

Stages B/C/E train on single-donor queries while the panel query is an equal-weight average of
two donors.  Averaging is not a neutral operation here -- experiment 019's whole result is that
whitening changes the sign of the donor-weighting effect -- so a correction fitted for
single-donor queries could in principle be the wrong correction for averaged ones, and that
would be a mundane explanation for the transfer failure rather than a line-boundary one.

So build the matched problem.  With four training lines per fold, a query can be the average of
two of them and the gallery a third, which is structurally identical to the panel
(two donors in, one unseen line as the gallery).  The inner leave-one-training-LINE-out then
holds the gallery line out of the fit entirely.  If the correction still fails to cross the
line boundary here, the mismatch explanation is closed.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fitlib import Adapter, BatchStats, batch_forward, make_batches  # noqa: E402
from foldio import Fold  # noqa: E402
from harness import LODO, Bank, load_banks  # noqa: E402
from stage_e import fit  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREP = HERE / "work" / "prep_multi"


def direction_perts(banks, tgt: str, pred: str, vcc: set) -> list[str]:
    """Reproduce prep_multi's `shared` list for one direction (target-bank order)."""
    return [p for p in banks[tgt].perts if p in banks[pred].idx and p not in vcc]


def build(f: Fold, banks, vcc, tgt: str, preds: tuple[str, str], gamma: float, k: int,
          batch: int, rows_split, seed: int) -> list[BatchStats]:
    """Query = mean of the two predictor lines' effects; gallery = the target line's truth."""
    lists = {p: direction_perts(banks, tgt, p, vcc) for p in preds}
    pos = {p: {q: i for i, q in enumerate(lists[p])} for p in preds}
    keep = [q for q in lists[preds[0]] if q in pos[preds[1]]]
    if len(keep) < batch:
        return []
    idx = {p: np.fromiter((pos[p][q] for q in keep), np.int64, len(keep)) for p in preds}
    Y = [np.load(PREP / f"{f.held}.{tgt}_from_{p}.Yq.npy", mmap_mode="r") for p in preds]
    Tsrc = np.load(PREP / f"{f.held}.{tgt}_from_{preds[0]}.Yt.npy", mmap_mode="r")
    w = (f.evals ** (-gamma / 2.0)).astype(np.float32)
    rng = np.random.default_rng(seed)
    sel = rows_split(len(keep))
    out = []
    for b in make_batches(len(sel), batch, rng):
        rr = np.sort(sel[b])
        Q = 0.5 * (np.asarray(Y[0][idx[preds[0]][rr]], np.float64)
                   + np.asarray(Y[1][idx[preds[1]][rr]], np.float64))
        Z = Q * w
        T = np.asarray(Tsrc[idx[preds[0]][rr]], np.float64)
        bs = BatchStats.__new__(BatchStats)
        bs.Zp, bs.Tb = Z[:, -k:].copy(), T[:, -k:].copy()
        bs.R = Z @ T.T - bs.Zp @ bs.Tb.T
        bs.c = np.einsum("ij,ij->i", Z, Z) - np.einsum("ij,ij->i", bs.Zp, bs.Zp)
        bs.tn = np.linalg.norm(T, axis=1)
        bs.B = len(rr)
        out.append(bs)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=256)
    ap.add_argument("--rank", type=int, default=8)
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--out", default=str(HERE / "work" / "stage_f.json"))
    args = ap.parse_args()
    k, r, g = args.k, args.rank, 0.65
    vcc = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"])
    res = {"cfg": vars(args), "folds": {}}

    for held in LODO:
        f = Fold(held)
        meta = json.loads((PREP / f"{held}.meta.json").read_text())
        lines = meta["lines"]
        banks = load_banks()
        for n in lines:
            if n not in banks:
                banks[n] = Bank(n)

        def split(n, frac=0.2, s=195):
            perm = np.random.default_rng(s).permutation(n)
            nd = int(round(frac * n))
            return np.sort(perm[nd:]), np.sort(perm[:nd])

        rows = {}
        for X in lines:
            seen = [n for n in lines if n != X]
            fit_b, dev_b = [], []
            for A in seen:
                for P in itertools.combinations([n for n in seen if n != A], 2):
                    fit_b += build(f, banks, vcc, A, P, g, k, args.batch,
                                   lambda n: split(n)[0], 31)
                    dev_b += build(f, banks, vcc, A, P, g, k, args.batch,
                                   lambda n: split(n)[1], 32)
            uns_b = []
            for P in itertools.combinations(seen, 2):
                uns_b += build(f, banks, vcc, X, P, g, k, args.batch, lambda n: np.arange(n), 33)
            if not fit_b or not uns_b:
                continue
            M = fit(fit_b, k, r, args.epochs, 0.01, 1e-3, 0.05, 195)
            row = {"n_fit_q": int(sum(b.B for b in fit_b)),
                   "n_unseen_q": int(sum(b.B for b in uns_b)),
                   "d_train": f.batch_pds(fit_b, M, k) - f.batch_pds(fit_b, None, k),
                   "d_dev": (f.batch_pds(dev_b, M, k) - f.batch_pds(dev_b, None, k))
                   if dev_b else None,
                   "d_unseen_line": f.batch_pds(uns_b, M, k) - f.batch_pds(uns_b, None, k),
                   "d_panel": f.panel_scores(g, M, k)[0] - f.panel_scores(g)[0]}
            rows[X] = row
            print(f"  {held:14} hold gallery line {X:14} fitq={row['n_fit_q']:6d} "
                  f"unsq={row['n_unseen_q']:5d}  d_train {row['d_train']:+.4f}  "
                  f"d_dev {row['d_dev']:+.4f}  d_unseenLINE {row['d_unseen_line']:+.4f}  "
                  f"d_panel {row['d_panel']:+.4f}", flush=True)
        res["folds"][held] = rows
        Path(args.out).write_text(json.dumps(res, indent=2))

    cells = [(h, X) for h in LODO for X in res["folds"][h]]
    macro = {kk: float(np.mean([res["folds"][h][X][kk] for h, X in cells]))
             for kk in ["d_train", "d_dev", "d_unseen_line", "d_panel"]}
    res["macro"] = macro
    res["n_cells"] = len(cells)
    print(f"\n== two-predictor queries, leave-one-gallery-LINE-out, macro over {len(cells)} cells ==")
    for kk, v in macro.items():
        print(f"   {kk:15} {v:+.5f}")
    Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
