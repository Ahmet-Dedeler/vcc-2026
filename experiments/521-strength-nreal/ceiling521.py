"""521 — the ceiling: what per-target call-set sizing is worth against TODAY's shipped call set on H1 (435 production
stats, seed 0), on all four DE members, before any predictor enters.

Delivered lfc pinned (variant b: the emitter keeps steering the dropped block genes; nmae exactly pinned) and variant a
(production emitter with a smaller CallSet: dropped block genes lose the clean4 block prior and take the design's own
unsteered lfc, read from 416's prior_off run at the same seed, plus the shipped MDE add in that sign).
Order: model order for baseline and arm alike (calls in realized order, then |key|).

Contents when a target is shrunk to K genes (all subsets of today's calls, so every call is deliverable):
  keyord  top K of today's calls in their realized p order (= |key| order: the head first)
  cpmord  top K of today's calls by control CPM (detectability, 192's membership key)
  headcpm min(K, 400) head genes by |key| then block genes by CPM
Policies: ORACLE K_p = m * n_real_p(f); negligible-only oracle (n_real < 20 -> K); call nothing (K = 0).
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/ceiling521.py
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402

FS = (1.0, 0.5, 0.35)


def content_rank(S, i, kind):
    """rank (0 = kept first) of today's calls under a content rule; non-calls get +inf."""
    sig = (S["oq"][i] < L.ALPHA) & ~S["own"][i]
    idx = np.flatnonzero(sig)
    if kind == "keyord":
        o = idx[np.lexsort((np.nan_to_num(S["op"][i][idx], nan=1.0), np.nan_to_num(S["oq"][i][idx], nan=1.0)))]
    elif kind == "cpmord":
        o = idx[np.argsort(-S["ref_lin"][idx], kind="stable")]
    elif kind == "headcpm":
        h = S["head"][i][idx]
        hk = idx[h][np.argsort(-S["absk"][i][idx[h]], kind="stable")]
        rest = idx[~h][np.argsort(-S["ref_lin"][idx[~h]], kind="stable")]
        o = np.concatenate([hk, rest])
    else:
        raise ValueError(kind)
    r = np.full(sig.size, np.inf)
    r[o] = np.arange(o.size)
    return r


def shrink_calls(S, K, kind, ranks=None):
    """delivered call set: today's calls with each target cut to its top K_p (K_p = inf keeps all). Unsuppressible
    calls (today's calls outside the prescribed set) are kept whatever K is: the emitter cannot remove them."""
    sig = (S["oq"] < L.ALPHA) & ~S["own"]
    uns = sig & ~S["C"]
    out = sig.copy()
    for i in range(sig.shape[0]):
        k = K[i]
        if not np.isfinite(k):
            continue
        r = ranks[i] if ranks is not None else content_rank(S, i, kind)
        out[i] = (r < k) | uns[i]
    return out


def lfc_variant_a(S, call):
    """variant a: block genes that leave the call set lose the block prior; they take 416 prior_off's delivered lfc (the
    unsteered split) plus the shipped MDE add (0.1/sqrt(cnt), cap 0.1, cnt >= 0.3) in that sign, times the target's
    recorded mass-balance factor (mde_fu / mde_fd)."""
    off = np.load(str(L.OFF416).format(seed=S["seed"]))
    lo = np.asarray(off["lfc"], float)
    diag = json.load(open(str(L.DIAG435).format(seed=S["seed"])))
    ref = np.load(L.REF)
    libs = np.asarray(ref["libraries"], float)
    cnt = S["ref_lin"] * float(np.median(libs)) / 1e6
    madd = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)
    out = S["ol"].copy()
    for i in range(out.shape[0]):
        drop = S["block"][i] & ~S["head"][i] & ~call[i]
        fu, fd = float(diag[i].get("mde_fu", 1.0)), float(diag[i].get("mde_fd", 1.0))
        s = np.sign(lo[i])
        out[i, drop] = lo[i, drop] + np.where(s > 0, fu, -fd)[drop] * madd[drop]
    return out


def main():
    S = L.load_state(0)
    sig = (S["oq"] < L.ALPHA) & ~S["own"]
    base_rows = {}
    res = {}
    ranks = {kind: [content_rank(S, i, kind) for i in range(85)] for kind in ("keyord", "cpmord", "headcpm")}
    # how different is the unsteered lfc on block genes from today's (sanity for variant a)
    off = np.asarray(np.load(str(L.OFF416).format(seed=0))["lfc"], float)
    blk = S["block"] & ~S["head"] & sig
    agree = (np.sign(off[blk]) == np.sign(S["ol"][blk])).mean()
    nb = ~S["C"] & ~S["own"]
    print(f"variant-a input: block genes' sign agreement unsteered vs today {agree:.3f}; non-block genes' |lfc| ratio "
          f"today/unsteered median {np.median(np.abs(S['ol'][nb]) / np.maximum(np.abs(off[nb]), 1e-9)):.3f}")
    adj = L.adjud(S)
    ok_t = (np.sign(S["ol"]) == np.sign(S["rl"]))
    ok_o = (np.sign(off) == np.sign(S["rl"]))
    print(f"  block genes, sign acc vs truth (all adjudicable): today {ok_t[blk & adj].mean():.4f}, unsteered {ok_o[blk & adj].mean():.4f}")
    for f in FS:
        pool = L.pools(S, f)
        nreal = pool.sum(1)
        m0 = L.members(S, sig, S["ol"], f, rows=True)
        base_rows[f] = m0
        print(f"\nf {f}: n_real median {np.median(nreal):.0f}, mean {nreal.mean():.0f}, zero {int((nreal == 0).sum())}, <20 "
              f"{int((nreal < 20).sum())}; BASE fid {m0['fid']:.5f} jac {m0['jac']:.5f} reach {m0['reach']:.5f} nmae {m0['nmae']:.5f}")
        for kind in ("keyord", "cpmord", "headcpm"):
            for lab, K in [(f"ORACLE m{m}", np.maximum(np.round(m * nreal), 0).astype(float)) for m in (1.0, 1.5, 2.0, 3.0, 5.0)] + \
                    [(f"NEG<20 ->K{k}", np.where(nreal < 20, float(k), np.inf)) for k in (0, 10, 25, 50, 100, 200, 400)] + \
                    [(f"ZERO-only ->0", np.where(nreal == 0, 0.0, np.inf))]:
                call = shrink_calls(S, K, kind, ranks[kind])
                mb = L.members(S, call, S["ol"], f)
                la = lfc_variant_a(S, call)
                ma = L.members(S, call, la, f)
                db = L.overall_delta(mb, m0)
                da = L.overall_delta(ma, m0)
                key = f"f{f}|{kind}|{lab}"
                res[key] = dict(b=mb, a=ma, d_b=db, d_a=da)
                print(f"  {kind:7s} {lab:14s} (b) fid {mb['fid'] - m0['fid']:+.5f} jac {mb['jac'] - m0['jac']:+.5f} "
                      f"reach {mb['reach'] - m0['reach']:+.5f} nmae {mb['nmae'] - m0['nmae']:+.5f} -> overall {db:+.5f} | "
                      f"(a) reach {ma['reach'] - m0['reach']:+.5f} nmae {ma['nmae'] - m0['nmae']:+.5f} -> {da:+.5f}")
    (L.WORK / "ceiling521.json").write_text(json.dumps(dict(res=res, base={str(k): {m: v for m, v in b.items() if m != 'rows'}
                                                                                for k, b in base_rows.items()}), indent=1, default=float))


if __name__ == "__main__":
    main()
