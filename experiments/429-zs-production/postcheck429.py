"""429 post-build check (423 section 8 item 4) for a `VCC_DE_MAG` arm against its pass-off parent built with the same
seed (e.g. 429a = 428a + VCC_DE_MAG=dn_k250_zs vs 428a). Both restored prediction h5ads (scripts/restore_prediction_h5ad.py).
Per (context, target) block, streamed:
  * identical obs order, gene axis, and number of rows;
  * identical per-gene block totals (so pds and mse's profile half are bit-identical) and identical libraries;
  * identical ballast rows (rows > 300k counts, the last K) and identical K;
  * every differing gene is in the scorer gate, is down (arm lfc < parent lfc < 0), and keeps its exact rank sum U
    against the context's controls in the evaluator's arithmetic (count / (lib / 1e6));
  * outside the gate, only sub-gate compensation moves.
Prints a census (blocks with moves, genes moved per block, Σ|Δl|) and exits 1 on any violation.
usage: uv run python postcheck429.py <arm.h5ad> <parent.h5ad> [--bundle data/challenge] [--every N]
       uv run python postcheck429.py --selftest <parent.h5ad>   (applies the pass in memory to every 50th block)"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import src_rank_mag as RM  # noqa: E402


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def block(h, ip, st, g):
    lo, hi = int(ip[st]), int(ip[st + 400])
    return sp.csr_matrix((np.asarray(h["X"]["data"][lo:hi]), np.asarray(h["X"]["indices"][lo:hi]), ip[st:st + 401] - lo),
                         shape=(400, g)).toarray().astype(np.int64)


def compare(A, B, cc):
    """Violations (list of str) and (genes moved, Σ|Δl|) for arm block A vs parent block B."""
    bad = []
    if not np.array_equal(A.sum(0), B.sum(0)):
        bad.append("per-gene totals")
    if not np.array_equal(A.sum(1), B.sum(1)):
        bad.append("libraries")
    ka, kb = int((A.sum(1) > 300_000).sum()), int((B.sum(1) > 300_000).sum())
    if ka != kb or not np.array_equal(A[400 - ka:], B[400 - kb:]):
        bad.append("ballast rows")
    keep = np.asarray(cc.keep, bool)
    gi = np.flatnonzero(keep)
    diff = np.flatnonzero((A != B).any(0))
    dk = np.isin(diff, gi)
    la, lb = RM.block_lfc(A, cc), RM.block_lfc(B, cc)
    kpos = {int(g): j for j, g in enumerate(gi)}
    s_a, s_b = A.sum(1) / 1e6, B.sum(1) / 1e6
    for g in diff[dk]:
        j = kpos[int(g)]
        if not (la[j] < lb[j] < 0):
            bad.append(f"gene {g}: not a down magnitude move ({lb[j]:.3f} -> {la[j]:.3f})")
        if int(RM.ranks_of(A[:, g], s_a, cc.cs[j]).sum()) != int(RM.ranks_of(B[:, g], s_b, cc.cs[j]).sum()):
            bad.append(f"gene {g}: U moved")
    mv = la != lb
    return bad, int(dk.sum()), float(np.abs(la[mv] - lb[mv]).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", nargs="?")
    ap.add_argument("parent")
    ap.add_argument("--bundle", default=str(ROOT / "data/challenge"))
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        a.every = max(a.every, 50)
    fails, moved, dls, nblk = 0, [], [], 0
    hp = h5py.File(a.parent, "r")
    ha = hp if a.selftest else h5py.File(a.arm, "r")
    genes = str_col(hp, "var/_index")
    ctxs, tg = str_col(hp, "obs/context"), str_col(hp, "obs/target_gene")
    if not a.selftest:
        if not (np.array_equal(genes, str_col(ha, "var/_index")) and np.array_equal(ctxs, str_col(ha, "obs/context"))
                and np.array_equal(tg, str_col(ha, "obs/target_gene"))):
            print("[FAIL] obs/var axes differ")
            sys.exit(1)
    ipp, ipa = np.asarray(hp["X"]["indptr"][:]), np.asarray(ha["X"]["indptr"][:])
    for ctx in sorted(set(ctxs)):
        controls = ad.read_h5ad(Path(a.bundle) / f"context_{ctx}.h5ad")
        cm = controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
        cc = RM.scorer_control_block(cm, min_cpm=5.0)
        del controls, cm
        starts = np.flatnonzero(ctxs == ctx)[::400]
        for k, st in enumerate(starts):
            if k % a.every:
                continue
            B = block(hp, ipp, int(st), genes.size)
            if a.selftest:
                K = int((B.sum(1) > 300_000).sum())
                A, _ = RM.apply(B, 400 - K, cc, "dn_k250_zs")
            else:
                A = block(ha, ipa, int(st), genes.size)
            bad, n, dl = compare(A, B, cc)
            nblk += 1
            moved.append(n)
            dls.append(dl)
            if bad:
                fails += 1
                print(f"[FAIL] {ctx}/{tg[st]}: {'; '.join(bad[:5])}")
        del cc
        print(f"[{ctx}] checked; blocks so far {nblk}, failures {fails}", flush=True)
    print(f"{nblk} blocks, {fails} with violations; blocks with moves {int(np.sum(np.array(moved) > 0))}; "
          f"genes moved per block median {np.median(moved):.0f} (max {max(moved)}); Σ|Δl| per block mean {np.mean(dls):.2f}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
