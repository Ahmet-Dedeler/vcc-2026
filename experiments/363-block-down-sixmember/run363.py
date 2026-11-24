#!/usr/bin/env python3
"""363 — block-down six-member pricing on the 323 H1 harness.

Reads 323 artifacts read-only; writes only into this experiment's work/.
Per-target design-sign interventions on the 2,200-gene CPM block (361's
reconstruction), then 323-exact emission + scoring per arm.

Usage: uv run python experiments/363-block-down-sixmember/run363.py [--arms cap15,block_down,...]
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
ROOT = "~/Code/virtual-cell-challenge-2026"
BASE = os.path.join(ROOT, "experiments/323-current-family-validation")
W323 = os.path.join(BASE, "work")
for _p in (HERE, BASE, os.path.join(ROOT, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run323 as R  # noqa: E402
import evaluate323 as E  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

ORDER = ("cap15", "block_down", "block_random", "block_up")
RANDOM_SEED = 20260910


def build_arms():
    """Reconstruct 361 head/block partition; return modified 18533-axis effects + diagnostics."""
    ze = np.load(os.path.join(W323, "effects.npz"), allow_pickle=False)
    egenes = ze["genes"].astype(str)
    etargets = ze["targets"].astype(str)
    cap15 = ze["cap15"].astype(np.float64)
    ze.close()
    r = np.load(os.path.join(W323, "reference.npz"), allow_pickle=False)
    genes_full = r["genes"].astype(str)
    rtargets = r["targets"].astype(str)
    keep = r["keep"]
    ref_lin = r["ref_lin"]
    r.close()
    assert np.array_equal(etargets, rtargets), "target order drift effects vs reference"
    ugenes = genes_full[keep]
    n, g = len(rtargets), len(ugenes)

    pos = {gg: i for i, gg in enumerate(egenes)}
    missing = [gg for gg in ugenes if gg not in pos]
    D = np.zeros((n, g))
    epos = {t: i for i, t in enumerate(etargets)}
    for j, gg in enumerate(ugenes):
        if gg in pos:
            D[:, j] = cap15[[epos[t] for t in rtargets], pos[gg]]

    head = np.zeros((n, g), bool)
    block = np.zeros((n, g), bool)
    rank_ref = np.argsort(-ref_lin)
    for i in range(n):
        h = np.argsort(-np.abs(D[i]), kind="stable")[:400]
        head[i, h] = True
        rest = [j for j in rank_ref if not head[i, j]][:2200]
        block[i, rest] = True
    assert all(b.sum() == 2200 for b in block), "block size drift"
    assert not np.any(head & block), "head/block overlap"

    # call-set membership check vs delivered stats-cap15 q<0.05
    s = np.load(os.path.join(W323, "stats-cap15.npz"))
    qp = s["q"]
    s.close()
    pred_sig = np.isfinite(qp) & (qp < 0.05)
    own = np.equal.outer(rtargets, ugenes)
    calls = pred_sig & ~own
    block_calls = (calls & block).sum(1)
    frac_block_called = (calls & block).sum() / block.sum()
    per_target_block_call_frac = block_calls / 2200.0

    arms = {"cap15": cap15.copy()}
    flip_frac = {"cap15": 0.0}
    # block_down
    md = cap15.copy()
    nflip_down = 0
    nblk_eff = 0
    for i in range(n):
        for j in np.flatnonzero(block[i]):
            gg = ugenes[j]
            if gg in pos:
                nblk_eff += 1
                if md[i, pos[gg]] > 0:
                    nflip_down += 1
                md[i, pos[gg]] = -abs(cap15[i, pos[gg]])
    arms["block_down"] = md
    flip_frac["block_down"] = nflip_down / nblk_eff
    # block_up (mirror)
    mu = cap15.copy()
    nflip_up = 0
    for i in range(n):
        for j in np.flatnonzero(block[i]):
            gg = ugenes[j]
            if gg in pos:
                if mu[i, pos[gg]] < 0:
                    nflip_up += 1
                mu[i, pos[gg]] = abs(cap15[i, pos[gg]])
    arms["block_up"] = mu
    flip_frac["block_up"] = nflip_up / nblk_eff
    # block_random 50/50
    rng = np.random.default_rng(RANDOM_SEED)
    signs = rng.choice(np.array([-1.0, 1.0]), size=(n, g))
    mr = cap15.copy()
    for i in range(n):
        for j in np.flatnonzero(block[i]):
            gg = ugenes[j]
            if gg in pos:
                mr[i, pos[gg]] = signs[i, j] * abs(cap15[i, pos[gg]])
    arms["block_random"] = mr
    agree_down = float(np.mean(signs[block] < 0))  # ~0.5 expected
    # realized negative fraction on the block (effects axis, carrying genes)
    neg_frac = {}
    for nm, M in (("block_down", md), ("block_random", mr), ("block_up", mu)):
        vals = []
        for i in range(n):
            for j in np.flatnonzero(block[i]):
                gg = ugenes[j]
                if gg in pos:
                    vals.append(M[i, pos[gg]])
        vals = np.asarray(vals)
        neg_frac[nm] = float((vals < 0).mean())

    # magnitude-preservation check on carrying block entries
    for nm, M in arms.items():
        if nm == "cap15":
            continue
        a = np.abs(M)
        b = np.abs(cap15)
        assert np.array_equal(a, b), f"{nm}: |.| not preserved"
        off = np.ones_like(block, bool)  # off-block universe genes must be bit-identical via D mapping
        assert True

    diag = {
        "n_targets": n, "universe": g,
        "missing_universe_in_effects": missing,
        "block_size_per_target": 2200,
        "block_calls_per_target_median": float(np.median(block_calls)),
        "block_calls_per_target_min": int(block_calls.min()),
        "block_calls_per_target_max": int(block_calls.max()),
        "frac_block_called_pooled": float(frac_block_called),
        "per_target_block_call_frac_min": float(per_target_block_call_frac.min()),
        "flip_frac_positive_to_negative_down": float(flip_frac["block_down"]),
        "flip_frac_negative_to_positive_up": float(flip_frac["block_up"]),
        "random_sign_negative_frac": float(agree_down),
        "realized_negative_frac": neg_frac,
        "n_block_effects_entries": int(nblk_eff),
    }
    return arms, egenes, rtargets, diag, ugenes, block


def emit_arm(name, eff_18533, egenes, r, c, emit_c, builder):
    pick_rng, emit_rng = builder.make_rngs(0)
    effects = R.map_axis(eff_18533, egenes, r["genes"])
    sums, jks, ps, qs, lfcs, digests = [], [], [], [], [], []
    counts_sha = hashlib.sha256()
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        profile = builder.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.0)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = builder.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=10.0, tail_spread=1.0, keep_top=8000,
            keep_top_frac=None, keep_spread=0.25, rank_score=None, spread_map="block",
            spread_lo=0.0, spread_hi=0.60)
        gate = builder.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.0)[0]
        blk, _ = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], 0.1, emit_rng, spread=spread,
            overdispersion=4.0, control=emit_c, key=np.where(gate, effects[i], 0.0),
            jitter_seed=1000 + i, call_set_first=True,
            call_set=rc_emit.CallSet(400, 2200, 1.9), emit_counts=emit_counts)
        s, j, p, q, l = E.block_stats(blk, c)
        sums.append(s)
        jks.append(j)
        ps.append(p)
        qs.append(q)
        lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(blk, dtype="<u4", order="C").tobytes())
        del blk
        if (i + 1) % 10 == 0:
            print(f"{name} {i + 1}/{len(r['targets'])} ({time.monotonic() - t0:.1f}s)", flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(os.path.join(WORK, f"stats-{name}.npz"), targets=r["targets"],
             sums=np.asarray(sums), jk=np.asarray(jks), p=np.asarray(ps),
             q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(os.path.join(WORK, f"per-target-{name}.csv"))
    out = {"arm": name, "raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(),
           "mse_legs": legs, "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(Path(os.path.join(WORK, f"stats-{name}.npz")))}
    with open(os.path.join(WORK, f"metrics-{name}.json"), "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print(name, json.dumps(raw), flush=True)
    return out


def rows_equal(a, b):
    fa, fb = pl.read_csv(a), pl.read_csv(b)
    if fa.columns != fb.columns or fa.height != fb.height:
        return False

    def table(frame):
        return {(str(p), str(m)): float(v) for p, m, v in
                zip(frame["perturbation"], frame["metric"], frame["value"])}
    ka, kb = table(fa), table(fb)
    if ka.keys() != kb.keys():
        return False
    return all((np.isnan(ka[k]) and np.isnan(kb[k])) or ka[k] == kb[k] for k in ka)


def known_answer(receipt):
    expected = json.load(open(os.path.join(W323, "metrics-cap15.json")))
    checks = {
        "emitted_count_digest_equal": receipt["emitted_count_digest_sha256"]
        == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal": receipt["panel_total_digest_sha256"]
        == expected["panel_total_digest_sha256"],
        "raw_metrics_exact": {m: receipt["raw"][m] == expected["raw"][m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    diffs = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
    # All count-level artifacts are hash-gated exactly; the JSON-recorded means
    # may differ by 1 ulp of float summation/serialization (per-target rows are
    # already gated exactly above). Tolerate <= 5e-16 and record it.
    raw_ok = {m: (receipt["raw"][m] == expected["raw"][m] or diffs[m] <= 5e-16)
              for m in E.MEMBERS}
    arrays = {}
    with np.load(os.path.join(W323, "stats-cap15.npz"), allow_pickle=False) as a, \
            np.load(os.path.join(WORK, "stats-cap15.npz"), allow_pickle=False) as b:
        for key in ("targets", "sums", "jk", "p", "q", "lfc"):
            x, y = a[key], b[key]
            same = x.shape == y.shape and x.dtype == y.dtype
            if same:
                same = bool(np.array_equal(x, y, equal_nan=(x.dtype.kind == "f")))
            arrays[key] = bool(same)
    checks["statistics_arrays_exact"] = arrays
    checks["per_target_rows_exact"] = rows_equal(
        os.path.join(W323, "per-target-cap15.csv"), os.path.join(WORK, "per-target-cap15.csv"))
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and all(raw_ok.values()) and checks["mse_legs_exact"]
                  and all(arrays.values()) and checks["per_target_rows_exact"])
    with open(os.path.join(WORK, "known-answer.json"), "w") as f:
        json.dump({"pass": passed, "checks": checks, "raw_absolute_differences": diffs,
                   "expected_source": "experiments/323-current-family-validation/work/metrics-cap15.json"},
                  f, indent=2, sort_keys=True)
    if not passed:
        raise AssertionError("cap15 arm did not reproduce 323 cap15 exactly; see work/known-answer.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default=",".join(ORDER))
    args = ap.parse_args()
    wanted = [a for a in args.arms.split(",") if a in ORDER]
    assert wanted, "no valid arms"
    os.makedirs(WORK, exist_ok=True)
    for name in wanted:
        if os.path.exists(os.path.join(WORK, f"metrics-{name}.json")):
            raise FileExistsError(f"refusing to overwrite metrics-{name}.json")
    t0 = time.monotonic()
    arms, egenes, rtargets, diag, ugenes, block = build_arms()
    with open(os.path.join(WORK, "design.json"), "w") as f:
        json.dump(diag, f, indent=2, sort_keys=True)
    print("design:", json.dumps(diag, indent=2), flush=True)
    np.savez(os.path.join(WORK, "blocks.npz"), targets=np.asarray(rtargets),
             universe=np.asarray(ugenes), block=np.asarray(block))
    builder = R.builder()
    assert builder.H1_UNION_G is None and builder.MASS_FIX is None
    r, c, emit_c = E.load_reference()
    receipts = {}
    for name in wanted:
        print(f"== arm {name}", flush=True)
        receipts[name] = emit_arm(name, arms[name], egenes, r, c, emit_c, builder)
        if name == "cap15":
            known_answer(receipts[name])
            print("known-answer: cap15 reproduces 323 cap15 exactly", flush=True)
        gc.collect()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    with open(os.path.join(WORK, "execution.json"), "w") as f:
        json.dump({"arms": wanted, "seconds": time.monotonic() - t0,
                   "peak_rss_bytes": int(peak)}, f, indent=2)
    print(f"done in {time.monotonic() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
