#!/usr/bin/env python3
"""423 (copied from 416's run416.py): base carriers + rank-preserving magnitude post-pass variants (rp423.py),
all variants scored from the SAME emitted base blocks within one process per seed.
416 — DE-carrier menu on H1 under ballast, seed-replicated. Copy of 410's run410b.py (ballast, ratio mode, jackknife
boost) with the CARRIER design as the one variable and the emission seed as a replicate axis; 410's files are read-only.

Per target: 400 - K carrier cells emitted by the K-class block (k3c5bq80 unless the design says otherwise) from the
design's carrier effect; K ballast cells carry alloc(S, bulk composition) - carrier totals (ratio mode: the bulk profile
is applied as a ratio on the carriers' DELIVERED composition); pair-wise anti-correlated composition noise restores the
carriers' own jackknife. The bulk effect is the shipped cap15 matrix for every design (--bulk base) or 410b's key-scaled
bulk (--bulk ps, lo 0.3 / hi 3, energy-preserving: the shipped 410b construction).

usage: uv run python run416.py DESIGN SEED [--bulk base|ps] [--kfrom RUNTAG]
  DESIGN: base | g000 | g030 | g090 | gain | signrep | prior_off | q70 | q90 | conf4
          follow-ups: pre (410's de_unwhite: gamma 0 without CD4) | w065nocd4 | gm030 (prep416b.py)
  --kfrom: force each target's ballast count K to the one recorded in work/h1/diag-RUNTAG.json (K-matched control)
known answers: `x125 0` == 410's de_x1p25_ratio_jk (carriers x1.25, cap15 bulk, K 1-3, jk boost live);
               `base 0 --bulk ps` == 410's bal_lo03_hi3_ratio_jk (the shipped 410b construction).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "h1"
WORK.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(ROOT / "experiments/416-de-carrier-menu"))   # k_emit416 (read-only)
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import run400 as H  # noqa: E402  (read-only)
import k_emit416 as K  # noqa: E402  (400's emitter + optional prior gate)

H.WORK = WORK
R, E = H.R, H.E
LB = 999_000
TS = 50_000.0
OP_H1 = ROOT / "experiments/323-current-family-validation/work/whitener_h1_pruned.npz"
STAGE = ROOT / "experiments/323-current-family-validation/work/stage.npz"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CARRIERS = HERE / "work/carriers416.npz"
CLEAN4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")

PRIOR_MODE = "q"      # "q" (k3c5bq80 family) | "off" (k3plain: K-class, unsteered)
PROTECT_Q = 0.8
GATE_FULL = None      # full-axis bool: prior may flip only here (conf4)


# ---------------------------------------------------------------- 410's ballast machinery, verbatim (run410b.py)
def alloc(total: int, profile: np.ndarray) -> np.ndarray:
    raw = profile * total
    base = np.floor(raw).astype(np.int64)
    rem = int(total - base.sum())
    if rem > 0:
        idx = np.argsort(-(raw - base))[:rem]
        base[idx] += 1
    return base


def jk_sum(block):
    X = np.asarray(block, dtype=np.float64)
    P = X.sum(0); S = P.sum(); L = X.sum(1)
    V = np.log1p(TS * (P[None, :] - X) / (S - L)[:, None])
    n = X.shape[0]
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())


def _balance_pair_noise(eps):
    d = int(eps.sum())
    if d == 0:
        return eps
    idx = np.flatnonzero(eps > 0) if d > 0 else np.flatnonzero(eps < 0)
    order = idx[np.argsort(-np.abs(eps[idx]), kind="stable")]
    mags = np.abs(eps[order])
    cum = np.cumsum(mags)
    stop = int(np.searchsorted(cum, abs(d)))
    eps[order[:stop]] = 0
    left = abs(d) - (int(cum[stop - 1]) if stop > 0 else 0)
    if left:
        eps[order[stop]] -= left if d > 0 else -left
    if int(eps.sum()) != 0:
        raise AssertionError("pair noise balance failed")
    return eps


def jk_boost(block, n_car, K_, target, rng):
    if K_ < 2:
        return block, 0.0, jk_sum(block)
    base = block[n_car:].copy()
    G = base.shape[1]
    z = np.where(rng.random(G) < 0.5, -1, 1)

    def build(rho):
        out = base.copy()
        for j in range(0, K_ - 1, 2):
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


def cap_mass(Q, K_, cap=995_000):
    excess = int(Q.sum()) - K_ * cap
    if excess <= 0:
        return Q, 0
    take = np.floor(Q * (excess / Q.sum())).astype(np.int64)
    short = excess - int(take.sum())
    if short > 0:
        take[np.argsort(-(Q - take), kind="stable")[:short]] += 1
    return Q - take, excess


# ---------------------------------------------------------------- the one variable: the carrier block
def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    cs = H.rc_emit.CallSet(*H.CALLSET)
    kw = dict(spread=spread, overdispersion=4., control=emit_c, key=key, jitter_seed=1000 + i, call_set_first=True,
              call_set=cs, emit_counts=H.EMIT_COUNTS)
    if PRIOR_MODE == "off":
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng, classes=list(H.CLASSES3), lfc_target=None, **kw)
    gated = np.asarray(emit_c.keep, dtype=bool)
    spec = {"prior": np.ascontiguousarray(H.PRIOR_FULL[gated]), "shrink_floor": 0.002, "protect_q": PROTECT_Q}
    if GATE_FULL is not None:
        spec["gate"] = np.ascontiguousarray(GATE_FULL[gated])
    return K.rc_hybrid_block_k(profile, libs, .1, emit_rng, classes=list(H.CLASSES3), lfc_target=spec, **kw)


def ballast_block(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, k_force=None):
    live = profile_de > 1e-9
    ratio = float(np.min(profile_bulk[live] / profile_de[live]))
    s_c_est = float(np.sum(libs[: 400 - kmax]))
    need = s_c_est * (1.0 / max(ratio, 1e-6) - 1.0)
    K_ = int(min(kmax, max(1, np.ceil(need / LB))))
    if k_force is not None:
        K_ = int(k_force)
    carriers, dg = carrier_block(profile_de, libs[: 400 - K_], emit_rng, spread, key, i, emit_c)
    carriers = np.asarray(carriers)
    T_c = carriers.sum(0).astype(np.int64)
    S = int(T_c.sum()) + K_ * LB
    # ratio mode (410b): the bulk/DE design ratio applied on the carriers' DELIVERED composition
    r = np.where(profile_de > 0, profile_bulk / np.maximum(profile_de, 1e-300), 1.0)
    tgt = (T_c / T_c.sum()) * r
    T_b = alloc(S, tgt / tgt.sum())
    Q = T_b - T_c
    clipped = Q < 0
    clip_mass = int(-Q[clipped].sum())
    Q = np.where(clipped, 0, Q)
    Q, trim = cap_mass(Q, K_)
    per = Q // K_
    ball = np.tile(per, (K_, 1))
    rem = Q - per * K_
    off = np.arange(Q.size) % K_
    for k in range(K_):
        ball[k] += (rem > ((k - off) % K_)).astype(np.int64)
    assert np.array_equal(ball.sum(0), Q) and ball.max() >= 0 and ball.sum(1).max() <= 1_000_000, ball.sum(1)
    block = np.vstack([carriers.astype(np.int64), ball]).astype(np.int64)
    jk_target = jk_sum(carriers)
    block, jk_rho, jk_after = jk_boost(block, carriers.shape[0], K_, jk_target, np.random.default_rng(1000 + i))
    assert np.array_equal(block.sum(0), np.vstack([carriers, ball]).sum(0)), "jk boost moved the totals"
    assert block.min() >= 0 and block.sum(1).max() <= 1_000_000
    dg = dict(dg)
    dg.update(jk_rho=jk_rho, jk_after=jk_after, jk_target=jk_target, sink_used=0)
    dg.update(ballast_K=K_, ballast_ratio=ratio, ballast_clip_genes=int(clipped.sum()), ballast_clip_mass=clip_mass,
              ballast_trim=trim, ballast_libs=[int(x) for x in ball.sum(1)], carrier_total=int(T_c.sum()))
    return block, dg


class H1Op:
    """run410h's minimal operator on the 323 H1 axis (read-only)."""
    def __init__(self):
        with np.load(OP_H1) as z:
            self.evals = np.asarray(z["evals"], float)
            self.evecs = np.asarray(z["evecs"], float)
        with np.load(STAGE, allow_pickle=False) as z:
            self.cols = np.asarray(z["columns"], dtype=np.int64)
            self.genes = [str(g) for g in z["genes"]]
        assert self.evecs.shape[0] == len(self.cols) == 18157

    def _check_axis(self, wcols):
        assert np.array_equal(np.asarray(wcols), self.cols)


def ps_bulk(eff, genes, targets):
    """410b's shipped bulk on H1: the key410 weights, clip(0.3, 3), energy-preserving over non-panel genes."""
    import key410  # noqa: PLC0415  (410's file, read-only)
    op = H1Op()
    key, carried, nd = key410.shared_energy_key(list(map(str, targets)), DONORS, op, op.cols, cols=None)
    w, med = key410.weights(key, carried, 0.3, 3.0)
    del op
    panel = set(map(str, targets))
    nonp = np.array([str(g) not in panel for g in genes])
    e2 = np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp])
    c = float(np.sqrt(e2.sum() / (w * w * e2).sum()))
    return eff * (w * c)[:, None]


def clean4_gate(genes_full):
    """conf4: True where the four CLEAN4 banks' P_down all sit on the same side of 0.5 (covered genes only)."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float)
    P = PB[[names.index(b) for b in CLEAN4]]
    unanimous = np.abs(np.sign(P - 0.5).sum(0)) == len(CLEAN4)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.zeros(len(genes_full), dtype=bool)
    out[col >= 0] = unanimous[col[col >= 0]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seed", type=int)
    ap.add_argument("--bulk", choices=["base", "ps"], default="base")
    ap.add_argument("--variants", default="", help="comma list of post-pass variants (rp423.parse_variant)")
    ap.add_argument("--dump", type=int, default=0, help="explore: save the first N base blocks and stop (no scoring)")
    a = ap.parse_args()
    import rp423  # noqa: PLC0415
    btag = f"_s{a.seed}" + ("_bps" if a.bulk == "ps" else "")
    variants = [v for v in a.variants.split(",") if v]
    specs = {v: rp423.parse_variant(v) for v in variants}
    b = H.builder389()
    eff, genes, targets = H.native_cap15(b)
    eff_de_native = eff
    eff_bulk = eff.copy() if a.bulk == "base" else ps_bulk(eff, genes, targets)
    r, cc, emit_c = E.load_reference()
    H.PRIOR_FULL, cov = H.prior_signs_full(r["genes"])
    eff_de = R.map_axis(eff_de_native, genes, r["genes"])
    eff_bk = R.map_axis(eff_bulk, genes, r["genes"])
    pick_rng, emit_rng = b.make_rngs(a.seed)
    names = ["base"] + variants
    acc = {v: dict(sums=[], jks=[], ps=[], qs=[], lfcs=[], diags=[]) for v in names}
    genes_rec = {v: [] for v in variants}
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
        block, dg = ballast_block(prof_de, prof_bk, r["libraries"][picks], emit_rng, spread,
                                  np.where(gate, eff_de[i], 0.), i, emit_c, 16)
        if K.LAST_PRIOR_REC:
            dg = {**dg, **K.LAST_PRIOR_REC.pop()}
        dg = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in dg.items()}
        dg["target"] = str(target)
        n_car = 400 - int(dg["ballast_K"])
        if a.dump:
            import scipy.sparse as sp  # noqa: PLC0415
            X = sp.csr_matrix(block.astype(np.int32))
            np.savez_compressed(HERE / f"work/explore/block{btag}_{i:02d}.npz", data=X.data, indices=X.indices,
                                indptr=X.indptr, shape=np.array(X.shape), n_car=n_car, target=str(target))
            print(f"dumped {i} {target} n_car={n_car}", flush=True)
            if i + 1 >= a.dump:
                return
            continue
        s, j, p, q, l = E.block_stats(block, cc)
        for nm, val in (("sums", s), ("jks", j), ("ps", p), ("qs", q), ("lfcs", l)):
            acc["base"][nm].append(val)
        acc["base"]["diags"].append(dg)
        R0T0 = rp423.r1_tie(block, cc) if variants else None
        for v in variants:
            blk2, pdg = rp423.postpass(block, n_car, cc, l, specs[v])
            s2, j2, p2, q2, l2 = E.block_stats(blk2, cc)
            tg, fr, lt_, dA_ = pdg.pop("_fr_by_gene")
            genes_rec[v].append(np.column_stack([np.full(tg.size, i), tg, fr, lt_, dA_]))
            pdg.update(rp423.verify(block, blk2, n_car, cc, s, p, q, l, s2, p2, q2, l2, R0T0=R0T0,
                                    touched=pdg.pop("_touched")))
            pdg["target"] = str(target)
            for nm, val in (("sums", s2), ("jks", j2), ("ps", p2), ("qs", q2), ("lfcs", l2)):
                acc[v][nm].append(val)
            acc[v]["diags"].append(pdg)
            del blk2
        if (i + 1) % 10 == 0:
            print(f"s{a.seed} {i + 1}/85 ({time.monotonic() - t0:.1f}s) K={dg['ballast_K']}", flush=True)
    for v in names:
        tag = ("base" if v == "base" else v) + btag
        A = acc[v]
        raw, rows, legs = E.score(r, cc, A["sums"], A["jks"], A["ps"], A["qs"], A["lfcs"])
        np.savez(WORK / f"stats-{tag}.npz", targets=r["targets"], sums=np.asarray(A["sums"]), jk=np.asarray(A["jks"]),
                 p=np.asarray(A["ps"]), q=np.asarray(A["qs"]), lfc=np.asarray(A["lfcs"]))
        pl.DataFrame(rows).write_csv(WORK / f"per-target-{tag}.csv")
        R.write_json(WORK / f"diag-{tag}.json", A["diags"])
        out = {"tag": tag, "variant": v, "spec": (None if v == "base" else specs[v]), "seed": a.seed, "bulk": a.bulk,
               "raw": raw, "mse_legs": legs, "seconds": time.monotonic() - t0}
        if v != "base":
            out["summary"] = rp423.summarise(A["diags"])
            np.save(WORK / f"genes-{tag}.npy", np.vstack(genes_rec[v]))   # cols: target, gene, frac, lfc, dA
        R.write_json(WORK / f"metrics-{tag}.json", out)
        print(tag, json.dumps({"raw": raw, "summary": out.get("summary")}), flush=True)


if __name__ == "__main__":
    main()
