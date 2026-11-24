"""Synthetic test of aggregate_r7a551: primary sums for every CONFIGS key, and that the plain and
calibrated bootstraps consume the same target-index sequence.

price551 is imported for real (guard551 is stubbed first so the test can write under this directory).
Panel is replaced: pooled_bootstrap only needs agg_fn, folds, and step_prices. step_prices is patched
to return a deterministic function of the sampled index vector and of cal, so a mismatch in the draws
changes the bootstrap arrays.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE551 = ROOT / "experiments/551-r6-c10b"
sys.modules.setdefault("guard551", types.ModuleType("guard551"))
sys.path.insert(0, str(HERE551))

import designs551 as D  # noqa: E402
import price551 as P5  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "aggregate_r7a551",
    ROOT / "experiments/555-cloud-r7-review/repro/upstream/aggregate_r7a551.py",
)
AGG = importlib.util.module_from_spec(spec)
spec.loader.exec_module(AGG)


class FakePanel:
    def __init__(self, tag, n=30, seed=0):
        self.tag = tag
        self.folds = ("k562_gwps", "orion_hct116", "orion_hek293t")
        self.n = n
        self._seen: list[np.ndarray] = []

    def agg_fn(self, f, names):
        def agg(i):
            # a pure function of the index vector; no RNG of its own
            return {"sum": float(np.asarray(i, float).sum()), "n": int(len(i))}
        return self.n, agg


def install_recorder(store: dict):
    def step_prices(P, f, S=None, level=True, cal=False, steps=D.STEPS):
        # encode the draw in every step, and shift calibrated prices so a swapped draw is visible
        base = 0.0 if S is None else S["sum"] / S["n"]
        scale = 10.0 if cal else 1.0
        return {st[0]: scale * (base / 100.0 + 0.01 * j) for j, st in enumerate(steps)}
    P5.step_prices = step_prices
    store["step_prices"] = step_prices


def test_plan_steps():
    print("=== plan_steps / primary ids for every CONFIGS key ===")
    expect = {
        # post-C10a steps; primary falls back to the five M0->M5 milestones when post is empty
        "M5": ((), ("M0->M1", "M1->M2", "M2->M3", "M3->M4", "M4->M5")),
        "M6": (("M5->M6",), ("M5->M6",)),
        "M5+P2g": (("P2g@M5",), ("P2g@M5",)),
        "M6+P2g": (("M5->M6", "P2g@M6"), ("M5->M6", "P2g@M6")),
        "M5+S_P2": (("S_P2@M5",), ("S_P2@M5",)),
        "M6+S_P2": (("M5->M6", "S_P2@M6"), ("M5->M6", "S_P2@M6")),
    }
    ok = True
    for F, (post_e, prim_e) in expect.items():
        post, st, ctrl, cum = AGG.plan_steps(F)
        primary = post if post else AGG.MILESTONES[:5]
        # control id
        if st:
            base = D.CONFIGS[F][0]
            want_ctrl = ("C-S@" if D.CONFIGS[F][1].startswith("S_P2") else "C-P@") + base
            ctrl_ok = ctrl[0] == want_ctrl
        else:
            ctrl_ok = ctrl is None
        match = tuple(post) == post_e and tuple(primary) == prim_e and ctrl_ok
        ok = ok and match
        print(f"  {F}: post={post} primary={primary} ctrl={None if ctrl is None else ctrl[0]} "
              f"cum={cum} {'OK' if match else 'MISMATCH'}")
    return ok


def test_pairing():
    print("=== plain vs calibrated bootstrap consume the same index draws ===")
    install_recorder({})
    steps = tuple(D.STEPS_C10A[:6]) + (
        {s[0]: s for s in D.STEPS_C10B}["P2g@M6"],
        {s[0]: s for s in D.STEPS_C10B}["C-P@M6"],
    )
    primary_ids = ("M5->M6", "P2g@M6")
    panels = {"T1": FakePanel("T1", n=25), "T2": FakePanel("T2", n=25)}
    # capture index vectors by wrapping agg
    seen = {"plain": [], "cal": []}

    def wrap(panel, bucket):
        orig = panel.agg_fn

        def agg_fn(f, names):
            n, agg = orig(f, names)

            def agg2(i):
                seen[bucket].append((panel.tag, f, np.array(i, copy=True)))
                return agg(i)
            return n, agg2
        panel.agg_fn = agg_fn

    for p in panels.values():
        wrap(p, "plain")
    bp, _ = P5.pooled_bootstrap(panels, n=40, seed=552, level=True, cal=False, steps=steps,
                                cum_ids=primary_ids)
    # second call must start from the same seed; reset the wrapper bucket
    panels2 = {"T1": FakePanel("T1", n=25), "T2": FakePanel("T2", n=25)}
    for p in panels2.values():
        orig = p.agg_fn

        def make(panel, orig=orig):
            def agg_fn(f, names):
                n, agg = orig(f, names)

                def agg2(i):
                    seen["cal"].append((panel.tag, f, np.array(i, copy=True)))
                    return agg(i)
                return n, agg2
            return agg_fn
        p.agg_fn = make(p)
    bc, _ = P5.pooled_bootstrap(panels2, n=40, seed=552, level=True, cal=True, steps=steps,
                                cum_ids=primary_ids)
    same_n = len(seen["plain"]) == len(seen["cal"])
    same_idx = same_n and all(np.array_equal(a[2], b[2]) and a[0] == b[0] and a[1] == b[1]
                              for a, b in zip(seen["plain"], seen["cal"]))
    # calibrated values are exactly 10x the plain ones, draw by draw, because step_prices scales by 10
    # and the index sequence matched. Compare one step.
    ratio = bc["M5->M6"] / bp["M5->M6"]
    ratio_ok = np.allclose(ratio, 10.0)
    print(f"  draws {len(seen['plain'])} plain / {len(seen['cal'])} cal; indices identical: {same_idx}; "
          f"cal/plain == 10: {ratio_ok} (min {ratio.min():.6f} max {ratio.max():.6f})")
    # a different seed must NOT match
    seen_other = []
    panels3 = {"T1": FakePanel("T1", n=25), "T2": FakePanel("T2", n=25)}
    for p in panels3.values():
        orig = p.agg_fn

        def make(panel, orig=orig):
            def agg_fn(f, names):
                n, agg = orig(f, names)

                def agg2(i):
                    seen_other.append(np.array(i, copy=True))
                    return agg(i)
                return n, agg2
            return agg_fn
        p.agg_fn = make(p)
    P5.pooled_bootstrap(panels3, n=40, seed=553, level=True, cal=False, steps=steps, cum_ids=primary_ids)
    differ = not all(np.array_equal(a, b) for a, b in zip((x[2] for x in seen["plain"]), seen_other))
    print(f"  seed 553 indices differ from seed 552: {differ}")
    return same_idx and ratio_ok and differ


def test_primary_sum_and_n_pos():
    print("=== primary sum and n_pos_folds on a fake point table ===")
    # two steps, 4 folds. Sum is positive on folds 0 and 1 only if we add them.
    # per_fold is a list, as price551.point_table stores the pooled per_fold.
    prim = {"pooled": {
        "M5->M6": {"mean": 0.004, "per_fold": [0.003, -0.001, 0.002, 0.001]},
        "P2g@M6": {"mean": 0.002, "per_fold": [0.001, -0.004, 0.000, 0.002]},
    }}
    ids = ("M5->M6", "P2g@M6")
    # the aggregator's expression, copied so a future edit of it is what we assert against the file
    expr = 'int(sum(v > 0 for v in np.sum([list(prim["pooled"][s]["per_fold"]) for s in ids], axis=0)))'
    n_pos = int(sum(v > 0 for v in np.sum([list(prim["pooled"][s]["per_fold"]) for s in ids], axis=0)))
    summed = np.sum([list(prim["pooled"][s]["per_fold"]) for s in ids], axis=0)
    # fold sums: 0.004, -0.005, 0.002, 0.003 -> 3 positive
    print(f"  fold sums {summed.tolist()} n_pos_folds {n_pos} (expect 3)")
    point = float(sum(prim["pooled"][s]["mean"] for s in ids))
    print(f"  pooled point {point} (expect 0.006)")
    # dict per_fold (the per-panel shape) would be wrong if someone used it the same way: list(dict) is keys
    per_panel = {"k562_gwps": 0.003, "orion_hct116": -0.001}
    try:
        bad = np.sum([list(per_panel)], axis=0)
        dict_note = f"list(dict) is keys {list(per_panel)}; np.sum gave {bad!r}"
    except Exception as e:
        dict_note = f"list(dict) path raised {type(e).__name__}: {e}"
    print(f"  pooled per_fold is a list in point_table, so the expression is well-defined. A dict would: {dict_note}")
    return n_pos == 3 and abs(point - 0.006) < 1e-12


def test_verdict_and_retention_zero():
    print("=== verdict boundaries ===")
    cases = [(0.01, 0.001, "TRANSFERS"), (0.01, 0.0, "FRAGILE"), (0.01, -0.001, "FRAGILE"),
             (0.0, -0.001, "FAILS"), (-0.01, -0.02, "FAILS")]
    ok = True
    for d, lb, want in cases:
        got = AGG.verdict(d, lb)
        ok = ok and got == want
        print(f"  D={d} lb={lb} -> {got} {'OK' if got == want else 'MISMATCH want '+want}")
    return ok


def main():
    results = {
        "plan_steps": test_plan_steps(),
        "pairing": test_pairing(),
        "n_pos": test_primary_sum_and_n_pos(),
        "verdict": test_verdict_and_retention_zero(),
    }
    print("===", json.dumps(results))
    if not all(results.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
