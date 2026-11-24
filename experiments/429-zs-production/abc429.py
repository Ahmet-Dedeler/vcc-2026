"""429: the production zs call on the real A/B/C blocks of the 428a package (426a + zero repair; the carriers zs would
see in production), with each context's own scorer-order control block built exactly as the patched builder builds it.
423 section 8 item 5: the A/B/C mass-holder census must look like H1's before the arm is priced.

Per block: K from the ballast libraries (rows > 300k counts; asserted to be the last K rows), then
src_rank_mag.apply(block, 400 - K, scorer_control, "dn_k250_zs"). Its in-pass asserts (per-gene totals, libraries,
U per moved gene in the evaluator's arithmetic, ballast rows, no sign flip) are the production asserts. Recorded: genes
moved, counts moved, the delivered lfc change on moved genes (sum of |Δl|, median |l| before/after), the feasible
fraction, and the pure mass-holder census (gated down genes whose nonzero carriers all sit above the control max).
The same summary is computed on H1's dumped 429 blocks for comparison (h1 rows).
No truth is read: this prices capacity, not sign accuracy.
usage: uv run python abc429.py <pred.h5ad> [--every N]"""
import argparse
import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(ROOT / "src"), str(HERE), str(ROOT / "experiments/323-current-family-validation")):
    sys.path.insert(0, _p)
import src_rank_mag as RM  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def summarise(X, n_car, cc, tpos):
    """apply() plus the capacity summary for one block."""
    keep = np.asarray(cc.keep, bool)
    l0 = RM.block_lfc(X, cc)
    X2, d = RM.apply(X, n_car, cc, "dn_k250_zs")
    l1 = RM.block_lfc(X2, cc)
    gi = np.flatnonzero(keep)
    notself = gi != tpos
    mv = (l1 != l0) & notself
    # pure mass-holders: gated, down, every nonzero carrier above the control max
    s = X.sum(1)[:n_car] / 1e6
    Xk = X[:n_car][:, gi]
    v = np.log1p(Xk / s[:, None])
    cmax = cc.cs[:, -1]
    nz = Xk > 0
    holder = nz.any(0) & np.all(~nz | (v > cmax[None, :]), axis=0) & (l0 < 0) & notself
    return X2, dict(n_car=int(n_car), K=int(X.shape[0] - n_car), n_moved=int(mv.sum()),
                    counts_moved=float(d["zs_counts_moved"]), dl_abs_sum=float(np.abs(l1[mv] - l0[mv]).sum()),
                    l_before_med=float(np.median(np.abs(l0[mv]))) if mv.any() else float("nan"),
                    l_after_med=float(np.median(np.abs(l1[mv]))) if mv.any() else float("nan"),
                    frac_massw=float(d["zs_frac_massw"]), holders=int(holder.sum()),
                    holders_moved=int((holder & mv).sum()), sub_exhausted=float(d["zs_sub_budget_exhausted_cells"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pred")
    ap.add_argument("--every", type=int, default=1, help="every N-th target per context")
    ap.add_argument("--bundle", default=str(ROOT / "data/challenge"))
    ap.add_argument("--gate", action="store_true", help="also build the emitter's rc_control and compare the gates (+1.6 GB)")
    a = ap.parse_args()
    out = {"pred": a.pred, "rows": [], "gate": {}}
    t0 = time.monotonic()
    with h5py.File(a.pred, "r") as h:
        genes = str_col(h, "var/_index")
        ctxs = str_col(h, "obs/context")
        tg = str_col(h, "obs/target_gene")
        ip = np.asarray(h["X"]["indptr"][:])
        gpos = {g: i for i, g in enumerate(genes)}
        for ctx in sorted(set(ctxs)):
            controls = ad.read_h5ad(Path(a.bundle) / f"context_{ctx}.h5ad")
            cm = controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
            assert list(map(str, controls.var_names)) == list(genes), "gene order"
            cc = RM.scorer_control_block(cm, min_cpm=5.0)
            rc = rc_emit.control_block_from_counts(cm, min_cpm=5.0) if a.gate else None
            gd = int((cc.keep != rc.keep).sum()) if a.gate else None
            out["gate"][ctx] = {"scorer_gate": int(cc.keep.sum()), "differs_from_emitter_gate": gd}
            print(f"[{ctx}] scorer-order control block {int(cc.keep.sum())} genes; differs from the emitter gate on {gd} "
                  f"({time.monotonic() - t0:.0f}s)", flush=True)
            del controls, cm, rc
            rows = np.flatnonzero(ctxs == ctx)
            starts = rows[::400]
            gi = np.flatnonzero(cc.keep)
            kpos = {int(g): j for j, g in enumerate(gi)}
            for k, st in enumerate(starts):
                if k % a.every:
                    continue
                st = int(st)
                assert tg[st + 399] == tg[st] and ctxs[st + 399] == ctx
                lo, hi = int(ip[st]), int(ip[st + 400])
                X = sp.csr_matrix((np.asarray(h["X"]["data"][lo:hi]), np.asarray(h["X"]["indices"][lo:hi]),
                                   ip[st:st + 401] - lo), shape=(400, genes.size)).toarray().astype(np.int64)
                lib = X.sum(1)
                ball = lib > 300_000
                K = int(ball.sum())
                assert K >= 1 and ball[400 - K:].all(), f"{ctx}/{tg[st]}: ballast rows are not the last K"
                tpos = kpos.get(gpos.get(str(tg[st]), -1), -1)
                _, rec = summarise(X, 400 - K, cc, tpos)
                rec.update(context=str(ctx), target=str(tg[st]))
                out["rows"].append(rec)
                if len(out["rows"]) % 50 == 0:
                    print(f"{len(out['rows'])} blocks ({time.monotonic() - t0:.0f}s) last {ctx}/{tg[st]} K={K} "
                          f"moved {rec['n_moved']}", flush=True)
            del cc
    # H1 comparison on the dumped seed-0 blocks (MDE + no repair; H1 has no zeros)
    sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))
    import evaluate323 as E  # noqa: PLC0415
    r, cc, _ = E.load_reference()
    gi = np.flatnonzero(cc.keep)
    kpos = {int(g): j for j, g in enumerate(gi)}
    hg = {str(g): i for i, g in enumerate(r["genes"])}
    for p in sorted((HERE / "work/blocks").glob("base_s0_bps_*.npz")):
        z = np.load(p)
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
        _, rec = summarise(X, int(z["n_car"]), cc, kpos.get(hg.get(str(z["target"]), -1), -1))
        rec.update(context="h1_mde" if "mde" in p.stem else "h1_base", target=str(z["target"]))
        out["rows"].append(rec)
    (HERE / "work/abc/abc429.json").write_text(json.dumps(out, indent=1))
    keys = ("n_moved", "counts_moved", "dl_abs_sum", "l_before_med", "l_after_med", "frac_massw", "holders", "holders_moved",
            "K", "sub_exhausted")
    print("context blocks " + " ".join(keys))
    for ctx in sorted({r_["context"] for r_ in out["rows"]}):
        R = [r_ for r_ in out["rows"] if r_["context"] == ctx]
        print(ctx, len(R), " ".join(f"{np.nanmedian([r_[k] for r_ in R]):.4g}" for k in keys), "(medians)")
        print(ctx, len(R), " ".join(f"{np.nanmean([r_[k] for r_ in R]):.4g}" for k in keys), "(means)")
    print(f"done {time.monotonic() - t0:.0f}s")


if __name__ == "__main__":
    main()
