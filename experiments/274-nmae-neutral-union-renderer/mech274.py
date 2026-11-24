"""274 -- the mechanism, per fold: where on the row does the union's transfer-cosine gain sit, and what does the
renderer remove?  Per target on frame A (a23 render, scored columns): the inner product with the truth and the
prediction energy split into the truth's wald89 pool / the incumbent's own top-89 / top-400 by |P| / the rest, for
the incumbent, the six-source union at g = 0.35 (R0) and its two key renderings (wald89 ratio = oracle, all rank = R2b).
One fold per process (~5.5 GiB).  Writes work/mech274-<fold>.json.

    uv run python experiments/274-nmae-neutral-union-renderer/mech274.py --fold k562_gwps
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import render274 as R  # noqa: E402  (sets the builder env, imports plolo with 274's caches)
import numpy as np  # noqa: E402

from plolo import FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator, build_fold, cd4_for, deliver  # noqa: E402
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, topn  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    MOD, _ = R.load_builder()
    fold = build_fold(held)
    targets = fold["targets"]
    truth_h = np.asarray(fold["truth"], dtype=np.float64)
    pre = np.asarray(fold["pre"], dtype=np.float64)
    cd4 = cd4_for(targets)
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    panel_positions = np.asarray(fold["panel_positions"], dtype=np.int64)
    unions = {sel: MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0, sources=R.sources_of(MOD, R.SEL[sel]))[0]
              for sel in ("pair", "stack6")}
    op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production" else PROD_OPS / "whitener_union_pruned.npz")
    op = Operator(op_path)
    wh = op.apply(pre * R.ALPHA, wcols, R.GAMMA, R.LAM)
    del op
    gc.collect()
    ctx = R.challenge_contexts()
    bA, gA = ctx["A"]
    sc = np.asarray(fold["scored"])
    Ts = truth_h[:, sc]
    cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
    LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
    wald = topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N)
    ok_cols = np.ones(R.N_GENES, dtype=bool)
    ok_cols[panel_positions] = False

    final_inc, info_inc = R.compose(MOD, wh, cd4, wcols, panel_positions, 5, 0.34, unions["pair"], 0.28)
    X_ref = final_inc / info_inc["preserve_scale"]
    absr = np.where(ok_cols[None, :], np.abs(X_ref), -1.0)
    own_sets = {K: [np.argpartition(-absr[i], K - 1)[:K] for i in range(X_ref.shape[0])] for K in (89, 400, 2000)}
    wald_cols = [sc[np.flatnonzero(wald[i])] for i in range(X_ref.shape[0])]
    fb, ib = R.compose(MOD, wh, cd4, wcols, panel_positions, 5, 0.34, unions["stack6"], 0.35)
    X6, changed = R.recover(fb, ib["preserve_scale"], X_ref)
    arms = {"cd4svd5": final_inc.copy(), "stack6 g=0.35": fb.copy()}
    for rset, rmode, K in (("wald89", "ratio", None), ("all", "rank", None), ("own", "ratio", 89), ("all", "norm", None)):
        Xo, _, _ = R.render(X6, X_ref, rset, rmode, K, ok_cols, own_sets, wald_cols)
        arms[R.rname("stack6 g=0.35", rset, rmode, K)], _ = R.preserve(Xo, wh, wcols)
    Xe, _ = R.render_exact_pool(X6, final_inc, wald_cols, changed, wh, wcols)
    arms[R.rname("stack6 g=0.35", "wald89", "exact")], _ = R.preserve(Xe, wh, wcols)
    # the incumbent's own top-K by |P| on the SCORED columns (the sets the renderer can see), mapped to scored positions
    sc_pos = {g: j for j, g in enumerate(sc)}
    own_sc = {}
    for K in (89, 400):
        own_sc[K] = [np.array([sc_pos[g] for g in own_sets[K][i] if g in sc_pos]) for i in range(X_ref.shape[0])]
    out = dict(fold=held, n_targets=len(targets), n_carried=int(changed.sum()), targets=list(targets), carried=changed.tolist(), arms={})
    Tn = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-300)
    t_energy_pool = np.array([float((Tn[i, wald[i]] ** 2).sum()) for i in range(Ts.shape[0])])
    out["truth_energy_share_wald89"] = t_energy_pool.tolist()
    for name, final in arms.items():
        final[:, panel_positions] = 0.0
        P = deliver(final, bA, gA, "minnorm")[:, sc]
        Pn = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
        rec = {}
        ip_all = np.einsum("ij,ij->i", Pn, Tn)                        # = own cosine
        rec["own_cos"] = ip_all.tolist()
        for tag, sets in (("wald89", [np.flatnonzero(wald[i]) for i in range(Ts.shape[0])]), ("own89", own_sc[89]), ("own400", own_sc[400])):
            ip_pool = np.array([float((Pn[i, s] * Tn[i, s]).sum()) for i, s in enumerate(sets)])
            en_pool = np.array([float((Pn[i, s] ** 2).sum()) for i, s in enumerate(sets)])
            l1_pool = np.array([float(np.abs(P[i, s]).sum()) for i, s in enumerate(sets)])
            sgn = np.array([float(np.where(np.sign(P[i, s]) == np.sign(Ts[i, s]), np.abs(Ts[i, s]), 0.0).sum() / max(np.abs(Ts[i, s]).sum(), 1e-300)) for i, s in enumerate(sets)])
            rec[f"ip_{tag}"] = ip_pool.tolist()
            rec[f"energy_{tag}"] = en_pool.tolist()
            rec[f"l1_{tag}"] = l1_pool.tolist()
            rec[f"signacc_{tag}"] = sgn.tolist()
        rec["l1_all"] = np.abs(P).sum(1).tolist()
        out["arms"][name] = rec
        print(f"[{held}] {name:36s} own cos {ip_all.mean():.5f}  ip wald89 {np.mean(rec['ip_wald89']):.5f}  energy wald89 {np.mean(rec['energy_wald89']):.4f}  ({time.time()-t0:.0f}s)", flush=True)
        del P, Pn
        gc.collect()
    (R.WORK / f"mech274-{held}.json").write_text(json.dumps(out))
    print(f"[{held}] wrote mech274-{held}.json ({time.time()-t0:.0f}s, peak {R.peak_gib():.2f} GiB)", flush=True)


if __name__ == "__main__":
    main()
