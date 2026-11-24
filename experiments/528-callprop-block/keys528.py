#!/usr/bin/env python3
"""528 — the call-propensity block keys for the H1 real-emitter screen, and their capture / jac emulation on the exact
production call set (head = top 400 by |key|, CPM block 2,200 among the rest; lib521's reconstruction of 435).

Features (from the recipient's OWN controls only; feat528_H1.npz): lc = log10 of the emitter's ref_lin (mean per-cell CPM),
x = excess dropout (517). Training truth never includes the held-out target:
  * `tv`  : fitted on h1-trainval's 126 perturbations (truth528: the scorer's calls vs h1-test's controls), applied to all
            85 h1-test targets (one key; the deployment form: coefficients from a public Flex screen, features from the
            recipient's controls);
  * `2f`  : 517's 2-fold by target on h1-test (halves from default_rng(517)); a target's key never saw its own truth.
Models: `p` = 517's logistic on [lc, x, lc x, lc^2]; `l` = the linear logit on [lc, x], used as its linear key lc + b x.
Controls: `rand` = the same swap size, random genes from the CPM band just below the kept block (matched expression).
Writes work/keys528.npz (per-target keys on the gated axis, reference target order) and work/keys528.json.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/keys528.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(HERE))
import lib521 as L  # noqa: E402
from callprop528 import block_by, feats, logit_fit  # noqa: E402

N_CPM = 2200


def gates_ref(S):
    return (S["rq"] < 0.05) & np.isfinite(S["rl"]) & (S["rl"] != 0) & ~S["own"]


def gates_tv(genes):
    z = np.load(HERE / "work/truth528_trainval.npz", allow_pickle=False)
    assert np.array_equal(z["genes"].astype(str), genes)
    perts = z["perts"].astype(str)
    q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
    own = np.equal.outer(perts, genes)
    return (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own, perts


def fit_keys(F2, F4, y, n):
    mp = logit_fit(F4, y, n)
    ml = logit_fit(F2, y, n)
    b = float(ml.coef_[0][1] / ml.coef_[0][0])
    return mp.predict_proba(F4)[:, 1], F2[:, 0] + b * F2[:, 1], dict(p_coef=mp.coef_[0].tolist() + [float(mp.intercept_[0])],
                                                                     l_coef=ml.coef_[0].tolist() + [float(ml.intercept_[0])],
                                                                     l_beta=b)


def main():
    S = L.load_state(0)
    T, G = S["key"].shape
    H = dict(np.load(HERE / "work/feat528_H1.npz", allow_pickle=True))
    z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    keep = np.asarray(z["keep"], bool)
    ref_lin = S["ref_lin"]
    x = H["x"][keep]
    lc = np.log10(np.maximum(ref_lin, 1e-3))
    F2, F4 = np.column_stack([lc, x]), feats(lc, x)
    gate = gates_ref(S)
    gtv, perts = gates_tv(S["genes"])
    info = {}
    keys = {}
    # tv: one key for every target
    ptv, ltv, info["tv"] = fit_keys(F2, F4, gtv.mean(0), gtv.shape[0])
    keys["ptv"] = np.tile(ptv, (T, 1))
    keys["ltv"] = np.tile(ltv, (T, 1))
    # 2f: 517's halves
    rng = np.random.default_rng(517)
    perm = rng.permutation(T)
    halves = (perm[: T // 2], perm[T // 2:])
    keys["p2f"] = np.zeros((T, G))
    keys["l2f"] = np.zeros((T, G))
    for h, (fit, ev) in enumerate((halves, halves[::-1])):
        p_, l_, info[f"2f_half{h}"] = fit_keys(F2, F4, gate[fit].mean(0), len(fit))
        keys["p2f"][ev] = p_
        keys["l2f"][ev] = l_
    # all 211 H1 perturbations (test + trainval): the A/B/C / October coefficients (NOT used on H1: it saw the truth)
    yall = np.concatenate([gate, gtv]).mean(0)
    _, _, info["all211"] = fit_keys(F2, F4, yall, T + gtv.shape[0])
    # lcap: the linear key lc + beta x with beta chosen to maximise the pooled capture of h1-trainval's gate genes by a
    # 2,200-gene block (panel and trainval genes excluded) -- selected on trainval only, so leak-free for h1-test
    excl = np.isin(S["genes"], list(perts) + list(S["targets"]))
    grid = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 50.0]
    capv = []
    for b_ in grid:
        k_ = lc + b_ * x
        bb = np.zeros(G, bool)
        bb[np.argsort(-np.where(excl, -np.inf, k_), kind="stable")[:N_CPM]] = True
        capv.append(float(gtv[:, bb].sum() / gtv.sum()))
    beta_c = grid[int(np.argmax(capv))]
    info["lcap"] = dict(beta=beta_c, grid=grid, trainval_capture=capv)
    keys["lcap"] = np.tile(lc + beta_c * x, (T, 1))
    # shuf: the falsifier -- the same beta on x shuffled among genes of similar expression (20 lc bins), so the swap has
    # lcap's size and expression profile but no gene-specific dropout information
    bins = np.digitize(lc, np.quantile(lc, np.linspace(0, 1, 21)[1:-1]))
    xs = x.copy()
    rs = np.random.default_rng(52802)
    for b_ in np.unique(bins):
        idx = np.flatnonzero(bins == b_)
        xs[idx] = x[rs.permutation(idx)]
    keys["shuf"] = np.tile(lc + beta_c * xs, (T, 1))
    info["shuf"] = dict(beta=beta_c, bins=20, seed=52802)
    # a matched random control: swap as many genes as ptv does (per-target mean), from the CPM band just below the kept part
    head = S["head"]
    blk0 = S["block"]
    b_ptv = np.array([block_by(keys["ptv"][i], head[i]) for i in range(T)])
    n_sw = int(round(np.mean([(b_ptv[i] & ~blk0[i]).sum() for i in range(T)])))
    keep_rand = N_CPM - n_sw
    u = np.random.default_rng(52801).random(G)
    keys["rand"] = np.zeros((T, G))
    for i in range(T):
        order = np.argsort(-np.where(head[i], -np.inf, ref_lin), kind="stable")
        band = order[keep_rand: keep_rand + 3 * n_sw]            # the next 3 n_sw genes by CPM after the kept part
        kr = np.full(G, -1.0)
        kr[band] = u[band]
        keys["rand"][i] = kr
    info["rand"] = dict(n_swap=n_sw, keep=keep_rand, band=f"CPM ranks [{keep_rand}, {keep_rand + 3 * n_sw}) among non-head")
    # capture / jac emulation for each key (full swap) and the partial swaps of ptv
    pool = L.pools(S, 1.0)
    rows = {}

    def summarize(name, blocks):
        cap = sum(int((blocks[i] & gate[i]).sum()) for i in range(T)) / int(gate.sum())
        per = float(np.mean([(blocks[i] & gate[i]).sum() / gate[i].sum() for i in range(T) if gate[i].sum() > 0]))
        jr = L.jac_rows(head | blocks, pool)
        sw = float(np.mean([(blocks[i] & ~blk0[i]).sum() for i in range(T)]))
        swin = np.concatenate([np.flatnonzero(blocks[i] & ~blk0[i]) for i in range(T)]) if sw else np.array([], int)
        rows[name] = dict(capture=cap, capture_per_target=per, jac_emul=float(jr.mean()), swapped_per_target=sw,
                          swapped_in_median_ref_lin=float(np.median(ref_lin[swin])) if swin.size else None,
                          swapped_in_median_x=float(np.median(x[swin])) if swin.size else None, _jr=jr)
    summarize("cpm", blk0)
    for k in ("ptv", "ltv", "p2f", "l2f", "lcap", "shuf"):
        summarize(k, np.array([block_by(keys[k][i], head[i]) for i in range(T)]))
    for kn in ("ltv", "lcap"):
        for f in (0.10, 0.25, 0.50):
            nk = int(round(N_CPM * (1 - f)))
            bl = []
            for i in range(T):
                kept = block_by(ref_lin, head[i], nk)
                bl.append(kept | block_by(keys[kn][i], head[i] | kept, N_CPM - nk))
            summarize(f"{kn}-k{nk}", np.array(bl))
    bl = []
    for i in range(T):
        kept = block_by(ref_lin, head[i], keep_rand)
        bl.append(kept | block_by(keys["rand"][i], head[i] | kept, N_CPM - keep_rand))
    summarize(f"rand-k{keep_rand}", np.array(bl))
    rng = np.random.default_rng(528)
    base = rows["cpm"]["_jr"]
    for k, r in rows.items():
        d = r.pop("_jr") - base
        bs = np.array([d[rng.integers(0, T, T)].mean() for _ in range(4000)])
        r["d_jac_emul"] = float(d.mean())
        r["d_jac_emul_ci"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
        print(f"{k:12s} capture {r['capture']:.4f} (per-target {r['capture_per_target']:.4f})  jac_emul {r['jac_emul']:.5f} "
              f"d {d.mean():+.5f} [{r['d_jac_emul_ci'][0]:+.5f}, {r['d_jac_emul_ci'][1]:+.5f}]  swapped {r['swapped_per_target']:6.1f}"
              f"  in: ref_lin {r['swapped_in_median_ref_lin']}  x {r['swapped_in_median_x']}", flush=True)
    for k, v in info.items():
        print(k, json.dumps(v))
    arrs = {k: np.ascontiguousarray(v) for k, v in keys.items()}
    np.savez(HERE / "work/keys528.npz", targets=S["targets"], **arrs)
    shas = {k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in arrs.items()}
    (HERE / "work/keys528.json").write_text(json.dumps(dict(info=info, rows=rows, key_sha256=shas), indent=1))


if __name__ == "__main__":
    main()
