#!/usr/bin/env python3
"""528 — unit tests of the emitter patch (k_emit528.py = patch_emitter528.py applied to production rc_emit_k.py 1658ae2c),
on a small synthetic context (no real data, < 1 GB). Run before any H1 job:

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/test_emitter528.py

  T1  block_key unset: the patched emitter's block, diagnostics and prior record equal production's, bit for bit
      (prior branch with q80 protection, coarse MDE, 496's complement hook, three library classes).
  T2  block_key = ref_lin itself, block_keep 0: the same ranking, so the same block, bit for bit.
  T3  a random block_key with block_keep = n_cpm: the shipped CPM block, so the same block, bit for bit.
  T4  a random block_key, block_keep 0 and a partial keep: the prescribed set is head | block_from_key; per-gene panel
      totals equal production's (pass 1 untouched); the steered prior block has the same size as today's.
  T5  block_from_key edge cases: stable ties toward the lower index, NaN refused, block_keep out of range refused,
      block_key without a call set refused.
  T6  prescribed_mask with and without a key.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from vcc2026 import rc_emit  # noqa: E402
from vcc2026 import rc_emit_k as P  # noqa: E402  (production, read-only)
from vcc2026.donor_transfer import emit_counts as EMIT  # noqa: E402
import k_emit528 as K  # noqa: E402

CLASSES3 = ((1 / 3, 16078.75), (1 / 3, 20190.5), (1 / 3, 25330.5))
CS = rc_emit.CallSet(30, 150, 1.9)


def context(seed=528):
    rng = np.random.default_rng(seed)
    m = 500
    # counts per cell at a 20k library: 900 DE-universe genes, 20 very high genes that carry most of the library, and
    # 2,500 sub-gate genes (< 5 CPM) that the emitter's row repair uses as slack
    mu = np.concatenate([10 ** rng.uniform(-0.7, 1.4, 900), rng.uniform(300, 700, 20), 10 ** rng.uniform(-3.0, -1.4, 2500)])
    mu *= 20000.0 / mu.sum()
    G = mu.size
    libs = np.clip(rng.normal(20000, 4000, m), 6000, None)
    shape = rng.uniform(0.3, 5.0, G)                           # bursty and Poisson-like genes
    lam = rng.gamma(shape[None, :], (mu / shape)[None, :] * (libs[:, None] / 20000.0))
    ctrl = rng.poisson(lam).astype(np.int64)
    control = rc_emit.control_block_from_counts(ctrl)
    eff = np.where(rng.random(G) < 0.3, rng.normal(0, 0.6, G), rng.normal(0, 0.05, G))
    comp = mu * np.exp(eff)
    profile = comp / comp.sum()
    key = np.where(control.keep, eff, 0.0)
    n_keep = int(control.keep.sum())
    prior = np.where(rng.random(n_keep) < 0.7, -1.0, 1.0)
    cnt = control.ref_lin * 20000 / 1e6
    mde = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)
    hook = np.where(rng.random(n_keep) < 0.3, np.where(rng.random(n_keep) < 0.6, -1.0, 1.0), 0.0)
    spec = {"prior": prior, "shrink_floor": 0.002, "protect_q": 0.8, "mde_add": mde, "mde_scope": "all",
            "mde_coarse": True, "mde_add_sign": hook}
    return dict(control=control, profile=profile, key=key, spec=spec, n_keep=n_keep, rng=rng)


def run(mod, C, **extra):
    libs = np.clip(np.random.default_rng(7).normal(20000, 3000, 120), 8000, None)
    kw = dict(spread=0.25, overdispersion=4.0, control=C["control"], key=C["key"], jitter_seed=1000, call_set_first=True,
              call_set=CS, classes=list(CLASSES3), lfc_target=dict(C["spec"]), emit_counts=EMIT)
    del mod.LAST_PRIOR_REC[:]
    block, diag = mod.rc_hybrid_block_k(C["profile"], libs, 0.1, np.random.default_rng(11), **kw, **extra)
    rec = mod.LAST_PRIOR_REC.pop() if mod.LAST_PRIOR_REC else None
    return block, diag, rec


def same(a, b):
    (ba, da, ra), (bb, db, rb) = a, b
    return (np.array_equal(ba, bb) and set(da) == set(db)
            and all(np.array_equal(np.asarray(da[k]), np.asarray(db[k]), equal_nan=True) for k in da) and ra == rb)


def main():
    C = context()
    ref_lin = np.asarray(C["control"].ref_lin, dtype=np.float64)
    n_keep = C["n_keep"]
    print(f"synthetic context: {n_keep} genes in the 5-CPM universe; CallSet {CS}", flush=True)
    base = run(P, C)
    ok = {}
    ok["T1 unset == production"] = same(run(K, C), base)
    ok["T2 key = ref_lin, keep 0 == production"] = same(run(K, C, block_key=ref_lin, block_keep=0), base)
    rk = np.random.default_rng(5).random(n_keep)
    ok["T3 random key, keep n_cpm == production"] = same(run(K, C, block_key=rk, block_keep=CS.n_cpm), base)
    # T4: a real swap
    kabs = np.abs(np.asarray(C["key"])[C["control"].keep])
    head = np.zeros(n_keep, bool)
    head[np.argsort(-kabs, kind="stable")[:CS.n_head]] = True
    cpm_blk = np.zeros(n_keep, bool)
    cpm_blk[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[:CS.n_cpm]] = True
    for keep_n in (0, 100):
        b4, d4, r4 = run(K, C, block_key=rk, block_keep=keep_n)
        want = head | K.block_from_key(head, ref_lin, CS.n_cpm, rk, keep_n)
        ok[f"T4 keep {keep_n}: totals == production"] = bool(np.array_equal(b4.sum(0), base[0].sum(0)))
        ok[f"T4 keep {keep_n}: prescribed size"] = d4["call_set_prescribed"] == CS.size
        # the called set = the realised q < .05 genes minus the unsuppressible plus the missed; check against `want`
        z, p, q, lfc = rc_emit.de_from_block(b4, C["control"])
        dcall = q < 0.05
        ok[f"T4 keep {keep_n}: missed/unsuppressible consistent with head|block_from_key"] = (
            int((want & ~dcall).sum()) == d4["call_set_missed"] and int((dcall & ~want).sum()) == d4["call_set_unsuppressible"])
        top_keep = np.zeros(n_keep, bool)
        top_keep[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[:keep_n]] = True
        ok[f"T4 keep {keep_n}: keeps the top-{keep_n} CPM genes and swaps some"] = (
            not bool((top_keep & ~want).any()) and int((want & ~(head | cpm_blk)).sum()) > 0)
        ok[f"T4 keep {keep_n}: steered block size <= today's"] = r4["prior_n_block"] <= CS.n_cpm
        print(f"  keep {keep_n}: swapped {int((want & ~(head | cpm_blk)).sum())}, missed {d4['call_set_missed']}, "
              f"unsuppressible {d4['call_set_unsuppressible']}, prior block {r4['prior_n_block']} (base {base[2]['prior_n_block']})",
              flush=True)
    # T5: edge cases
    hd = np.array([True, False, False, False, False, False])
    rl = np.array([9.0, 5.0, 5.0, 4.0, 3.0, 1.0])
    bk = np.array([0.0, 1.0, 2.0, 2.0, 0.5, 2.0])
    b = K.block_from_key(hd, rl, 2, bk, 0)
    ok["T5 ties toward the lower index"] = b.tolist() == [False, False, True, True, False, False]
    b = K.block_from_key(hd, rl, 3, bk, 1)
    ok["T5 keep 1 then key"] = b.tolist() == [False, True, True, True, False, False]
    b = K.block_from_key(hd, rl, 3, bk, 3)
    cpm3 = np.zeros(6, bool)
    cpm3[np.argsort(-np.where(hd, -np.inf, rl), kind="stable")[:3]] = True
    ok["T5 keep n_cpm == CPM block"] = np.array_equal(b, cpm3)
    for bad, args in (("NaN key", (hd, rl, 2, np.array([0, 1, np.nan, 2, 0, 1.0]), 0)),
                      ("keep > n_cpm", (hd, rl, 2, bk, 3)), ("keep < 0", (hd, rl, 2, bk, -1)),
                      ("wrong length", (hd, rl, 2, bk[:5], 0))):
        try:
            K.block_from_key(*args)
            ok[f"T5 refuses {bad}"] = False
        except ValueError:
            ok[f"T5 refuses {bad}"] = True
    try:
        libs = np.full(120, 20000.0)
        K.rc_hybrid_block_k(C["profile"], libs, 0.1, np.random.default_rng(1), spread=0.25, overdispersion=4.0,
                            control=C["control"], key=C["key"], jitter_seed=1, call_set_first=True, call_set=None,
                            emit_counts=EMIT, block_key=rk)
        ok["T5 refuses a key without a call set"] = False
    except ValueError:
        ok["T5 refuses a key without a call set"] = True
    # T6: prescribed_mask
    ok["T6 prescribed_mask unset == production"] = np.array_equal(K.prescribed_mask(C["key"], C["control"], CS),
                                                                P.prescribed_mask(C["key"], C["control"], CS))
    ok["T6 prescribed_mask keyed == head | block_from_key"] = np.array_equal(
        K.prescribed_mask(C["key"], C["control"], CS, rk, 0), head | K.block_from_key(head, ref_lin, CS.n_cpm, rk, 0))
    for k, v in ok.items():
        print(("PASS " if v else "FAIL ") + k)
    n_pass = sum(bool(v) for v in ok.values())
    print(f"{n_pass}/{len(ok)} PASS")
    sys.exit(0 if n_pass == len(ok) else 1)


if __name__ == "__main__":
    main()
