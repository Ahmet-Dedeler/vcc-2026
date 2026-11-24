#!/usr/bin/env python3
"""410 BALLAST prototype on H1: the bulk (pds/mse) at the per-target SCALED design, the DE channel at 395's UNSCALED design.

Per target: 400 - K carrier cells are emitted by the shipped k3c5bq80 block from the UNSCALED effect (same spread / key /
call-set design as base395); then K ballast cells carry Q = alloc(S, profile_bulk) - T_carriers, where S = S_c + K * LB,
so the block's per-gene totals are the exact-total allocation of the SCALED profile (clipped where Q < 0). The scorer's
pseudobulk sums raw counts (ballast cells dominate it); per-cell CPM (DE: arithmetic means, Wilcoxon ranks) gives each
ballast cell 1/400 of the weight. K is chosen per target from the design ratio min_g profile_bulk/profile_de, capped.
usage: uv run python run410b.py ARM   (ARM: bal_lo03_hi3 | bal_lo03_hi3_k8 | bal_ident)
  bal_ident: ballast with w = 1 (bulk == DE design): isolates the ballast construction's own cost on every member.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import run400 as H  # noqa: E402
import run410h  # noqa: E402  (sets H.WORK to work/h1 and routes block_for_arm to k3c5bq80)

R, E = H.R, H.E
LB = 999_000          # ballast library (the per-cell cap is 1e6)
RATIO = False
DROP_ONLY = False
JKBOOST = False
SINK = False


def alloc(total: int, profile: np.ndarray) -> np.ndarray:
    raw = profile * total
    base = np.floor(raw).astype(np.int64)
    rem = int(total - base.sum())
    if rem > 0:
        idx = np.argsort(-(raw - base))[:rem]
        base[idx] += 1
    return base


TS = 50_000.0


def jk_sum(block):
    """Delete-one jackknife of the log-pseudobulk (cell-eval2's C: (n-1)/n * sum_g sum_i (v_ig - vbar_g)^2)."""
    X = np.asarray(block, dtype=np.float64)
    P = X.sum(0); S = P.sum(); L = X.sum(1)
    V = np.log1p(TS * (P[None, :] - X) / (S - L)[:, None])
    n = X.shape[0]
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())


def _balance_pair_noise(eps):
    """Make integer pair noise sum to exactly zero by REDUCING magnitudes on the heavier side, largest first (only
    shrinks |eps|, so every bound 0 <= a +- eps holds; always terminates because the heavier side's mass >= |sum|)."""
    d = int(eps.sum())
    if d == 0:
        return eps
    idx = np.flatnonzero(eps > 0) if d > 0 else np.flatnonzero(eps < 0)
    order = idx[np.argsort(-np.abs(eps[idx]), kind="stable")]
    mags = np.abs(eps[order])
    cum = np.cumsum(mags)
    stop = int(np.searchsorted(cum, abs(d)))            # entries order[:stop] are zeroed, order[stop] takes the rest
    eps[order[:stop]] = 0
    left = abs(d) - (int(cum[stop - 1]) if stop > 0 else 0)
    if left:
        eps[order[stop]] -= left if d > 0 else -left
    if int(eps.sum()) != 0:
        raise AssertionError("pair noise balance failed")
    return eps


def jk_boost(block, n_car, K, target, rng):
    """Pair-wise anti-correlated composition noise on the ballast rows (sums, hence the pseudobulk, are exact; the
    pair's contributions to the per-cell arithmetic mean cancel). rho chosen by bisection so jk >= target."""
    if K < 2:
        return block, 0.0, jk_sum(block)
    base = block[n_car:].copy()
    G = base.shape[1]
    z = np.where(rng.random(G) < 0.5, -1, 1)
    def build(rho):
        out = base.copy()
        for j in range(0, K - 1, 2):
            a, b = out[j], out[j + 1]
            m = np.minimum(a, b)
            eps = np.floor(rho * m).astype(np.int64) * z
            eps = _balance_pair_noise(eps)
            out[j] = a + eps; out[j + 1] = b - eps
        return np.vstack([block[:n_car], out])
    lo, hi = 0.0, 0.95
    cur = jk_sum(block)
    if cur >= target:
        return block, 0.0, cur
    best = build(hi); v = jk_sum(best)
    if v < target:
        return best, hi, v
    for _ in range(12):
        mid = 0.5 * (lo + hi); cand = build(mid); v = jk_sum(cand)
        if v >= target:
            hi, best = mid, cand
        else:
            lo = mid
    return best, hi, jk_sum(best)


def cap_mass(Q, K, cap=995_000):
    """Clipping returns mass to the ballast; trim any excess over K*cap proportionally from the largest entries."""
    excess = int(Q.sum()) - K * cap
    if excess <= 0:
        return Q, 0
    take = np.floor(Q * (excess / Q.sum())).astype(np.int64)
    short = excess - int(take.sum())
    if short > 0:
        take[np.argsort(-(Q - take), kind="stable")[:short]] += 1
    return Q - take, excess


def ballast_block(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, w_final=1.0, t_idx=None, b_ctrl=None):
    # K from the design ratio: need S_c / S <= min_g profile_bulk / profile_de over genes the carriers will carry
    live = profile_de > 1e-9
    ratio = float(np.min(profile_bulk[live] / profile_de[live]))
    s_c_est = float(np.sum(libs[: 400 - kmax]))
    need = s_c_est * (1.0 / max(ratio, 1e-6) - 1.0)
    K = int(min(kmax, max(1, np.ceil(need / LB))))
    carriers, dg = H.block_for_arm("k3c5bq80", profile_de, libs[: 400 - K], emit_rng, spread, key, i, emit_c)
    carriers = np.asarray(carriers)
    if DROP_ONLY:
        return carriers.astype(np.int64), dict(dg, ballast_K=K, ballast_ratio=ratio, ballast_clip_genes=0, ballast_clip_mass=0, ballast_trim=0)
    T_c = carriers.sum(0).astype(np.int64)
    S = int(T_c.sum()) + K * LB
    sink_used = 0
    if SINK and t_idx is not None:
        # SINK mode: non-target genes read back EXACTLY w * (carriers' delivered log effect); the target gene (excluded from
        # every metric) absorbs the composition's missing mass. Falls back to ratio mode where the sink would go negative.
        bc = np.log1p(TS * T_c / T_c.sum())
        x = np.clip(np.expm1(b_ctrl + w_final * (bc - b_ctrl)), 0.0, None) / TS
        x[t_idx] = 0.0
        xt = 1.0 - x.sum()
        if xt >= T_c[t_idx] / S:
            x[t_idx] = xt
            T_b = alloc(S, x / x.sum())
            sink_used = 1
    if not sink_used and RATIO:
        # ratio mode: the per-target scale is applied as a RATIO on the carriers' DELIVERED composition, so w == 1
        # reproduces the incumbent's delivered pseudobulk (emitter rounding included) instead of the design profile
        r = np.where(profile_de > 0, profile_bulk / np.maximum(profile_de, 1e-300), 1.0)
        tgt = (T_c / T_c.sum()) * r
        T_b = alloc(S, tgt / tgt.sum())
    elif not sink_used:
        T_b = alloc(S, profile_bulk)
    Q = T_b - T_c
    clipped = Q < 0
    clip_mass = int(-Q[clipped].sum())
    Q = np.where(clipped, 0, Q)
    Q, trim = cap_mass(Q, K)
    # split Q evenly over K cells (remainder to the first cells), each <= 1e6
    per = Q // K
    ball = np.tile(per, (K, 1))
    rem = Q - per * K
    off = np.arange(Q.size) % K          # rotate the remainders across cells so no cell collects them all
    for k in range(K):
        ball[k] += (rem > ((k - off) % K)).astype(np.int64)
    assert np.array_equal(ball.sum(0), Q) and ball.max() >= 0 and ball.sum(1).max() <= 1_000_000, ball.sum(1)
    block = np.vstack([carriers.astype(np.int64), ball]).astype(np.int64)
    jk_rho, jk_after, jk_target = 0.0, float("nan"), float("nan")
    if JKBOOST:
        jk_target = jk_sum(carriers)          # the carriers' own jackknife: what the block had before ballast
        block, jk_rho, jk_after = jk_boost(block, carriers.shape[0], K, jk_target, np.random.default_rng(1000 + i))
        assert np.array_equal(block.sum(0), np.vstack([carriers, ball]).sum(0)), "jk boost moved the totals"
        assert block.min() >= 0 and block.sum(1).max() <= 1_000_000
    dg = dict(dg)
    dg.update(jk_rho=jk_rho, jk_after=jk_after, jk_target=jk_target, sink_used=sink_used)
    dg.update(ballast_K=K, ballast_ratio=ratio, ballast_clip_genes=int(clipped.sum()), ballast_clip_mass=clip_mass, ballast_trim=trim,
              ballast_libs=[int(x) for x in ball.sum(1)], carrier_total=int(T_c.sum()))
    return block, dg


def main():
    arm = sys.argv[1]
    global RATIO, DROP_ONLY, JKBOOST, SINK
    RATIO = "_ratio" in arm
    SINK = "_sink" in arm
    JKBOOST = "_jk" in arm
    DROP_ONLY = arm.startswith("drop")
    kmax = 8 if "_k8" in arm else 16
    WORK = run410h.WORK
    if (WORK / f"metrics-{arm}.json").exists():
        raise FileExistsError(arm)
    b = H.builder389()
    eff, genes, targets = H.native_cap15(b)
    eff_de_native = eff
    w_final_all = np.ones(len(targets))
    if arm.startswith("de_x"):
        # DE-side global amplitude: carriers at k x the shipped effect, bulk pinned at the incumbent (ballast, ratio mode)
        kx = float(arm.split("_")[1][1:].replace("p", "."))
        eff_de_native = eff * kx
        eff_bulk = eff.copy()
    elif arm.startswith("de_unwhite"):
        # DE-side: carriers carry the UNWHITENED donor mean (323 stage.npz, same arm), row-norm matched to the shipped rows;
        # the bulk stays the incumbent (ballast, ratio mode) -- only the DE members can move
        with np.load(ROOT / "experiments/323-current-family-validation/work/stage.npz", allow_pickle=False) as zs:
            assert (zs["targets"].astype(str) == np.asarray(targets).astype(str)).all() and (zs["genes"] == genes).all()
            pre = np.asarray(zs["cap15"], dtype=np.float64)
        scale = np.linalg.norm(eff, axis=1) / np.maximum(np.linalg.norm(pre, axis=1), 1e-12)
        eff_de_native = pre * scale[:, None]
        eff_bulk = eff.copy()
    elif arm.startswith("de_med"):
        # DE-side per-target scale: carriers at w_de (median target 1, clip 0.5..1.8), bulk held at the incumbent design
        w, key, carried, nd = run410h.h1_weights(targets, 0.5, 1.8)
        eff_de_native = eff * w[:, None]
        eff_bulk = eff.copy()
    elif arm.startswith("bal_ident") or arm == "drop1":
        eff_bulk = eff.copy()
    else:
        lo, hi = (0.3, 1e9) if "_inf" in arm else (0.3, 3.0)
        w, key, carried, nd = run410h.h1_weights(targets, lo, hi)
        panel = set(map(str, targets))
        nonp = np.array([str(g) not in panel for g in genes])
        e2 = np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp])
        c = float(np.sqrt(e2.sum() / (w * w * e2).sum()))
        eff_bulk = eff * (w * c)[:, None]
        w_final_all = w * c
    r, cc, emit_c = E.load_reference()
    H.PRIOR_FULL, _ = H.prior_signs_full(r["genes"])
    H.TARGET_NAMES = [str(t) for t in r["targets"]]
    eff_de = R.map_axis(eff_de_native, genes, r["genes"])
    eff_bk = R.map_axis(eff_bulk, genes, r["genes"])
    pick_rng, emit_rng = b.make_rngs(0)
    sums, jks, ps, qs, lfcs, diags = [], [], [], [], [], []
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        prof_de = b.profile_from_effect(eff_de[i], r["b_ctrl"], min_cpm=10.)
        prof_bk = b.profile_from_effect(eff_bk[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(
            eff_de[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(eff_de[i], r["b_ctrl"], min_cpm=10.)[0]
        gpos = {str(g): k for k, g in enumerate(r["genes"])}
        block, dg = ballast_block(prof_de, prof_bk, r["libraries"][picks], emit_rng, spread,
                                  np.where(gate, eff_de[i], 0.), i, emit_c, kmax,
                                  w_final=float(w_final_all[i]), t_idx=gpos.get(str(target)), b_ctrl=r["b_ctrl"])
        if H.K.LAST_PRIOR_REC:
            dg = {**dg, **H.K.LAST_PRIOR_REC.pop()}
        dg = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in dg.items()}
        dg["target"] = str(target)
        diags.append(dg)
        s, j, p, q, l = E.block_stats(block, cc)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        if (i + 1) % 20 == 0:
            print(f"{arm} {i + 1}/85 ({time.monotonic() - t0:.1f}s) K={dg['ballast_K']} clip={dg['ballast_clip_genes']}", flush=True)
    raw, rows, legs = E.score(r, cc, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{arm}.npz", targets=r["targets"], sums=np.asarray(sums), jk=np.asarray(jks),
             p=np.asarray(ps), q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{arm}.csv")
    R.write_json(WORK / f"diag-{arm}.json", diags)
    out = {"arm": arm, "raw": raw, "mse_legs": legs, "seconds": time.monotonic() - t0,
           "K_hist": {int(k): int(v) for k, v in zip(*np.unique([d["ballast_K"] for d in diags], return_counts=True))},
           "clip_genes_median": float(np.median([d["ballast_clip_genes"] for d in diags]))}
    R.write_json(WORK / f"metrics-{arm}.json", out)
    print(arm, json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
