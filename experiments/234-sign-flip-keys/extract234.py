"""234 — extract one PLOLO fold's truth / pool / delivered-sign / KEY tensors.

One fold per process (232 measured 5.06 GiB for a single full arm rebuild; the 17,942^2
whitener is the cost).  Everything downstream runs on the small `.npz` this writes, so the
screen itself is a few hundred MB.

WHAT IS REPRODUCED EXACTLY.  The stack, the render and the pool are 232's, unchanged:
  * `final` = shipped `h1union` stack (whiten @ alpha/gamma/lam, CD4 @ 0.34, H1 union @ 0.28,
    panel zeroed) -- 232 `common232.FoldBundle.final(stack="base222")`.
  * `H_A`   = `deliver(final, b_ctrl(A), gate(A), None)` -- 188's frame, a15 render.  This is
    the object 232 measures `q_w` on, so the flip screen must measure it on the same one.
  * pools   = 188's `wald89` / `lfc89` / `h89`, top-89 per row.
The KA row printed at the end is 232 §2b's delivered `q_w` per fold.

    uv run python extract234.py k562_gwps
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"                 # 234's OWN cache; nothing outside is written
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import ALPHA, GAMMA, LAM, source_residual  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, topn  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/220-ordering-keys-retry"))
import common220 as C220  # noqa: E402

from common234 import KA232_QW_DELIVERED, qw, prep_weights, row_rank_frac  # noqa: E402

CD4_F, UNION_G = 0.34, 0.28
UNION = ("gse295214_h1", "h1_2025")
SRC232 = ROOT / "experiments/232-sign-vote/work"        # READ ONLY


def ctx_gate(ctx: str):
    b = plolo.challenge_bctrl(ctx)
    e = np.expm1(b)
    return b, (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (b > 0)


def build(held: str) -> dict:
    fold = build_fold(held, limit_gib=7.5)
    sc = np.asarray(fold["scored"])
    wcols = np.asarray(fold["wcols"])
    truth_h = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_h[:, sc]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    targets = fold["targets"]

    op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                  if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))
    wh = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    del op
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    frac_cd4, _ = solve_fraction(wh, cd4, CD4_F)
    frac_u, _ = solve_fraction_carrying(wh, h1, UNION_G)
    res = merged(cd4, frac_cd4, h1, frac_u)
    final, _ = E062.add_post_whiten_residual(wh, res, 1.0, wcols)
    final[:, fold["panel_positions"]] = 0.0
    del wh, cd4, h1, res

    bA, gA = ctx_gate("A")
    H_A = deliver(final, bA, gA, None)

    b = fold["b_ctrl"]
    cpm_full = 1e6 * np.expm1(b) / float(np.expm1(b).sum())
    LR = to_lfc(truth_h, b, gA)[:, sc]
    pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_full[sc], 0)), POOL_N),
             "lfc89": topn(np.abs(LR), POOL_N),
             "h89": topn(np.abs(truth), POOL_N)}
    return dict(fold=fold, sc=sc, wcols=wcols, truth=truth, pre=pre, final=final,
                H_A=H_A, pools=pools, b_ctrl=b, cpm=cpm_full, targets=targets)


def donor_keys(held: str, sc: np.ndarray, T: int) -> dict:
    """Per-(target, gene) donor statistics, read from 232's frozen `sources-*.npz`."""
    z = np.load(SRC232 / f"sources-{held}.npz", allow_pickle=True)
    E = np.asarray(z["effects"], dtype=np.float64)[:, :, sc]
    avail = np.asarray(z["available"], dtype=bool)
    supp = np.asarray(z["support"], dtype=bool)[:, sc]
    names = [str(x) for x in z["names"]]
    cov = avail[:, :, None] & supp[:, None, :]
    S = E.shape[0]

    s = (np.sign(E) * cov).astype(np.int8)
    pos = (s > 0).sum(0).astype(np.int16)
    neg = (s < 0).sum(0).astype(np.int16)
    ncov = cov.sum(0).astype(np.int8)

    zw = np.zeros_like(E)
    for k in range(S):
        rows = avail[k]
        sd = E[k][rows].std(0)
        zw[k] = np.abs(E[k]) / np.where(sd > 0, sd, np.inf)[None, :]
    zw *= cov
    zsum = zw.sum(0)
    zscore = (s * zw).sum(0)

    absum = (np.abs(E) * cov).sum(0)
    ssum = (E * cov).sum(0)
    n = np.maximum(ncov, 1).astype(np.float64)
    mean = ssum / n
    var = ((E - mean[None]) ** 2 * cov).sum(0) / np.maximum(n - 1, 1)
    sem = np.sqrt(np.maximum(var, 0) / n)

    zpos = np.where(cov, zw, np.inf)
    out = dict(
        ncov=ncov.astype(np.float32),
        margin_maj=(np.abs(pos.astype(np.float64) - neg)
                    / np.maximum(pos + neg, 1)).astype(np.float32),
        margin_z=(np.abs(zscore) / np.maximum(zsum, 1e-300)).astype(np.float32),
        coh=(np.abs(ssum) / np.maximum(absum, 1e-300)).astype(np.float32),
        tdonor=np.where(ncov >= 2, np.abs(mean) / np.maximum(sem, 1e-300),
                        0.0).astype(np.float32),
        minz=np.where(ncov >= 1, np.min(zpos, axis=0), 0.0).astype(np.float32),
        sumz=zsum.astype(np.float32),
        disagree=((pos > 0) & (neg > 0)),
        s_mean=np.sign(ssum).astype(np.int8))
    return out, names


def main() -> None:
    t0 = time.time()
    held = sys.argv[1]
    d = build(held)
    sc, truth, final, H_A = d["sc"], d["truth"], d["final"], d["H_A"]
    T = truth.shape[0]

    Hs = H_A[:, sc]
    Fs = final[:, sc]
    Ps = d["pre"][:, sc]
    sgn_del = np.sign(Hs).astype(np.int8)
    sgn_fin = np.sign(Fs).astype(np.int8)
    correct = ((sgn_del != 0) & (sgn_del == np.sign(truth)))

    theta = C220.theta_on([plolo.OFFICIAL[j] for j in sc])
    sd_g = Fs.std(0)
    keys = dict(
        absH=np.abs(Hs).astype(np.float32),
        absfin=np.abs(Fs).astype(np.float32),
        abspre=np.abs(Ps).astype(np.float32),
        key189=(np.abs(Hs) * np.exp(np.clip(0.5 * theta, -50, 50))[None, :]).astype(np.float32),
        theta=np.broadcast_to(theta[None, :], Hs.shape).astype(np.float32).copy(),
        expr=np.broadcast_to(d["b_ctrl"][sc][None, :], Hs.shape).astype(np.float32).copy(),
        cpm=np.broadcast_to(np.log1p(d["cpm"][sc])[None, :],
                            Hs.shape).astype(np.float32).copy(),
        zpred=(np.abs(Fs) / np.where(sd_g > 0, sd_g, np.inf)[None, :]).astype(np.float32),
        absH_rank=row_rank_frac(np.abs(Hs)),
        whitenkept=(np.sign(Ps) == np.sign(Fs)).astype(np.float32),
        rownorm=np.broadcast_to(
            np.linalg.norm(H_A[:, d["wcols"]], axis=1)[:, None],
            Hs.shape).astype(np.float32).copy())
    dk, names = donor_keys(held, sc, T)
    disagree = dk.pop("disagree")
    s_mean = dk.pop("s_mean")
    keys.update({k: v for k, v in dk.items()})

    # KA row: 232 §2b's delivered |t|-weighted sign accuracy on `wald89`
    w, den = prep_weights(truth, d["pools"]["wald89"])
    q = qw(sgn_del, truth, w, den)
    exp = KA232_QW_DELIVERED.get(held)
    ok = exp is None or abs(q - exp) < 1e-12
    print(f"[{held}] delivered q_w(wald89) = {q:.15f}" +
          (f"   232 BASE: {exp:.15f}   D {abs(q-exp):.1e}   {'PASS' if ok else 'FAIL'}"
           if exp else ""), flush=True)
    assert ok, f"KA row failed for {held}: q_w {q} vs 232's {exp}"

    out = WORK / f"ex234-{held}.npz"
    np.savez_compressed(
        out, truth=truth.astype(np.float32), sgn_del=sgn_del, sgn_fin=sgn_fin,
        correct=correct, disagree=disagree, s_mean=s_mean, sc=sc,
        targets=np.asarray(d["targets"]),
        **{f"pool_{k}": v for k, v in d["pools"].items()},
        **{f"key_{k}": v for k, v in keys.items()})
    meta = dict(fold=held, n_targets=int(T), n_scored=int(sc.size), donors=names,
                q_w_delivered=q, q_w_232=exp, keys=sorted(keys),
                seconds=time.time() - t0)
    (WORK / f"ex234-{held}.json").write_text(json.dumps(meta, indent=1))
    print(f"  -> {out.name} ({out.stat().st_size/1e6:.0f} MB, {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
