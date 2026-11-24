"""Pure tests for 342: rounding rules, the per-block adapter, the AST seam, the unchanged
floor path through production `emit_counts` and `rc_hybrid_block`, and injection-site
specificity. Synthetic fixtures only; no repository data and no .npz/.npy artifacts.

    uv run python experiments/342-hamilton-current-family/test342.py     (or pytest)
"""
from __future__ import annotations

import json
import os
import sys
import traceback
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "2")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(HERE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import rounding342 as RD                       # noqa: E402
import seam342 as SM                           # noqa: E402
from vcc2026 import donor_transfer as production  # noqa: E402
from vcc2026 import rc_emit                    # noqa: E402


# ----------------------------------------------------------------------------- helpers
def quasi_fixture(seed=42, genes=50, cells=25):
    rng = np.random.default_rng(seed)
    profile = rng.gamma(0.7, 1.0, genes)
    profile[:5] = 0.0
    profile /= profile.sum()
    libs = rng.integers(3000, 6000, cells).astype(np.float64)
    spread = np.where(rng.random(genes) < 0.5, 0.25, 1.0)
    return profile, libs, spread


def production_L(libs):
    return int(min(np.median(libs), 1_000_000))


def rng_state(rng):
    return json.dumps(rng.bit_generator.state, sort_keys=True, default=int)


def same_csr(a, b):
    return all(x.dtype == y.dtype and np.array_equal(x, y) for x, y in zip(a, b))


def to_dense(vals, cols, ptr, n, genes):
    block = np.zeros((n, genes), dtype=np.int64)
    for row in range(n):
        lo, hi = int(ptr[row]), int(ptr[row + 1])
        block[row, cols[lo:hi]] = vals[lo:hi].astype(np.int64)
    return block


def control_fixture(genes=60, controls=300, cells=20, tail=10, seed=3):
    rng = np.random.default_rng(seed)
    # Ten tail genes far below the scorer's 5-CPM gate in the controls (they are the
    # sub-gate slack the per-cell weld repairs with), fifty gated genes well above it.
    means = np.concatenate([rng.uniform(5.0, 200.0, genes - tail), np.full(tail, 0.002)])
    counts = rng.poisson(means[None, :] * rng.uniform(0.7, 1.3, (controls, 1))).astype(np.float64)
    control = rc_emit.control_block_from_counts(counts)
    assert int(control.keep.sum()) == genes - tail, "fixture tail genes must sit below 5 CPM"
    effect = rng.normal(0.0, 0.3, genes)
    effect[genes - tail:] = 0.0
    profile = means * np.exp(effect)
    # Each tail gene carries 2.5% of the gated mass (about 20% of the block in total), so
    # every cell holds hundreds of sub-gate counts for the weld to move.
    profile[genes - tail:] = profile[:genes - tail].sum() * 0.025
    profile /= profile.sum()
    libs = np.full(cells, 4000.0)
    spread = np.where(rng.random(genes) < 0.5, 0.25, 1.0)
    spread[genes - tail:] = 0.25
    key = np.where(control.keep, effect, 0.0)
    return control, profile, libs, spread, key


def rc_block(control, profile, libs, spread, key, emit_counts, seed=7):
    return rc_emit.rc_hybrid_block(
        profile, libs, 0.1, np.random.default_rng(seed), spread=spread, overdispersion=4.0,
        control=control, key=key, jitter_seed=1000, call_set_first=True,
        call_set=rc_emit.CallSet(5, 10, 1.9), emit_counts=emit_counts)


# ----------------------------------------------------------------------------- rules
def test_known_fractional_totals():
    cases = [
        ([1.2, 2.7, 3.1, 0.4], [1, 3, 3, 0]),          # sum 7.4 -> 7, one increment to .7
        ([0.5, 0.5, 0.5, 0.5], [1, 1, 0, 0]),          # sum 2.0 -> two, ties to lower index
        ([2.25, 0.75, 1.5, 3.5], [2, 1, 2, 3]),        # .75 first, then the .5 tie at index 2
        ([1.5, 0.5, 2.5], [2, 0, 2]),                  # round(4.5) = 4 (half to even) -> one
    ]
    for x, want in cases:
        x = np.array(x, dtype=np.float64)
        totals, chosen = RD.hamilton_totals(x)
        assert np.array_equal(totals, want), (x, totals)
        assert chosen.size == RD.residual_count(x)
        assert np.array_equal(RD.floor_totals(x), np.floor(x).astype(np.int64))


def test_zero_and_integer_entries():
    x = np.array([0.0, 2.0, 0.0, 3.0])
    assert RD.residual_count(x) == 0
    for rule in RD.RULES:
        totals = RD.BlockRounder(rule, 0)(x)
        assert np.array_equal(totals, [0, 2, 0, 3]), rule
    x = np.array([0.0, 1.6, 0.0, 2.4, 0.0])
    totals, chosen = RD.hamilton_totals(x)
    assert np.array_equal(totals, [0, 2, 0, 2, 0]) and chosen.tolist() == [1]
    for index in range(20):
        totals = RD.BlockRounder("random_up", index)(x)
        assert totals[[0, 2, 4]].tolist() == [0, 0, 0], "zero coordinates never move"
        assert int(totals.sum()) == 4


def test_stable_ties():
    x = np.full(7, 0.3)
    totals, chosen = RD.hamilton_totals(x)
    assert chosen.tolist() == [0, 1] and totals.tolist() == [1, 1, 0, 0, 0, 0, 0]
    x = np.array([0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.9])
    totals, chosen = RD.hamilton_totals(x)
    assert chosen.tolist() == [0, 1, 6], chosen
    a, _ = RD.hamilton_totals(x)
    b, _ = RD.hamilton_totals(x.copy())
    assert np.array_equal(a, b), "deterministic"


def test_mass_conservation_and_zero_one_differences():
    rng = np.random.default_rng(2024)
    for trial in range(20):
        x = rng.gamma(2.0, 30.0, 400)
        x[rng.random(400) < 0.1] = 0.0          # exact zeros
        x[rng.random(400) < 0.1] = 3.0          # exact integers
        k = RD.residual_count(x)
        frac = x - np.floor(x)
        for rule, index in (("hamilton", trial), ("random_up", trial)):
            totals = RD.BlockRounder(rule, index)(x)
            delta = totals - np.floor(x).astype(np.int64)
            assert int(totals.sum()) == int(round(float(x.sum()))), rule
            assert set(np.unique(delta).tolist()) <= {0, 1}, rule
            assert int(delta.sum()) == k, rule
            assert not np.any(delta[frac == 0]), rule
        summary = RD.summarize(x, *RD.hamilton_totals(x))
        assert summary["mass_conserved"] and summary["zero_one_differences"]
        assert summary["moved"] == k == summary["residual_k"]


def test_residual_count_bounds_and_definition():
    rng = np.random.default_rng(1)
    for _ in range(50):
        x = rng.gamma(1.5, 20.0, 300) * (rng.random(300) < 0.8)
        k = RD.residual_count(x)
        frac = x - np.floor(x)
        assert 0 <= k <= int((frac > 0).sum())
        assert k == int(round(float(x.sum()))) - int(np.floor(x).astype(np.int64).sum())
    for bad in (np.array([[1.0, 2.0]]), np.array([1.0, np.nan]), np.array([-1.0, 2.0])):
        try:
            RD.residual_count(bad)
            raise AssertionError("invalid input accepted")
        except ValueError:
            pass


def test_random_control_mass_and_independent_generator():
    rng = np.random.default_rng(0)
    x = rng.uniform(0.0, 50.0, 200)
    x[:20] = 0.0
    k = RD.residual_count(x)
    assert k > 20
    a = RD.BlockRounder("random_up", 3)(x)
    b = RD.BlockRounder("random_up", 3)(x)
    assert np.array_equal(a, b), "the dedicated generator is deterministic in target index"
    expected, chosen = RD.random_up_totals(x, RD.control_generator(3))
    assert np.array_equal(a, expected)
    assert chosen.size == k == np.unique(chosen).size, "without replacement"
    assert np.all(x[chosen] - np.floor(x[chosen]) > 0)
    c = RD.BlockRounder("random_up", 4)(x)
    assert not np.array_equal(a, c), "another target index draws another selection"
    delta = a - np.floor(x).astype(np.int64)
    assert set(np.unique(delta).tolist()) <= {0, 1} and int(delta.sum()) == k
    emission = np.random.default_rng(99)
    before = rng_state(emission)
    RD.BlockRounder("random_up", 5)(x)
    assert rng_state(emission) == before, "the emission stream is untouched"
    assert SM.SEAMED == "totals = __vcc342_round_totals__(n * L * profile_p)", \
        "the seam hands the rule ONE argument: the intended vector, never a generator"
    for rule in ("floor", "hamilton"):
        assert RD.BlockRounder(rule, 0).control_rng is None


def test_block_rounder_restriction_and_call_count():
    x = np.array([1.2, 2.7, 3.1, 0.4, 0.0])
    r = RD.BlockRounder("hamilton", 0)
    assert np.array_equal(r(x), [1, 3, 3, 0, 0])
    assert np.array_equal(r(np.array([0.0, 2.7, 0.0, 0.4, 0.0])), [0, 3, 0, 0, 0])
    assert r.calls == 2
    for bad in (np.array([1.2, 2.7, 3.1, 0.5, 0.0]), np.array([1.2, 2.7, 3.1])):
        try:
            r(bad)
            raise AssertionError("non-restriction accepted")
        except ValueError:
            pass
    f = RD.BlockRounder("floor", 0)
    assert np.array_equal(f(x), np.floor(x).astype(np.int64))
    assert np.array_equal(f(np.where(x > 1, 0.0, x)), [0, 0, 0, 0, 0])
    assert f.calls == 2
    try:
        RD.BlockRounder("ceil", 0)
        raise AssertionError("unknown rule accepted")
    except ValueError:
        pass


def test_integrity_treats_floor_deficit_as_intended():
    from run342 import rounding_invariants
    x = np.array([1.2, 2.7, 3.1, 0.4, 0.0])
    for rule in RD.RULES:
        r = RD.BlockRounder(rule, 0)
        r(x)
        r(np.where(x > 1, 0.0, x))
        record = r.summary()
        assert rounding_invariants(rule, [record])["all_pass"], rule
        assert record["mass_conserved"] == (rule != "floor")
        bad = {**record, "panel_total": record["panel_total"]+1}
        assert not rounding_invariants(rule, [bad])["all_pass"], rule
    assert not rounding_invariants("floor", [])["all_pass"]


# ----------------------------------------------------------------------------- seam
def test_unchanged_floor_path_through_emit_counts():
    profile, libs, spread = quasi_fixture()
    seam = SM.Seam()
    rounder = RD.BlockRounder("floor", 0)
    seamed = seam.emit_counts_with(rounder)
    sub = np.arange(len(profile)) % 3 == 0
    for seed, p in ((11, profile), (12, np.where(sub, 0.0, profile))):
        rng_a, rng_b = np.random.default_rng(seed), np.random.default_rng(seed)
        a = production.emit_counts(p, libs, 0.1, rng_a, mode="quasi", spread=spread,
                                   overdispersion=4.0)
        b = seamed(p, libs, 0.1, rng_b, mode="quasi", spread=spread, overdispersion=4.0)
        assert same_csr(a, b), "seamed floor differs from production"
        assert rng_state(rng_a) == rng_state(rng_b), "emission stream consumed identically"
    assert rounder.calls == 2
    x = len(libs) * production_L(libs) * profile
    assert np.array_equal(rounder.x, x), "the hook received production's own n * L * profile"


def test_hamilton_totals_flow_through_the_machinery():
    profile, libs, spread = quasi_fixture(seed=7)
    seam = SM.Seam()
    n, genes = len(libs), len(profile)
    x = n * production_L(libs) * profile
    expected, chosen = RD.hamilton_totals(x)
    assert chosen.size > 0
    rounder = RD.BlockRounder("hamilton", 4)
    vals, cols, ptr = seam.emit_counts_with(rounder)(
        profile, libs, 0.1, np.random.default_rng(3), mode="quasi", spread=spread,
        overdispersion=4.0)
    block = to_dense(vals, cols, ptr, n, genes)
    assert np.array_equal(block.sum(0), expected), "Hamilton totals were not delivered"
    assert np.array_equal(rounder.totals, expected) and rounder.calls == 1
    fv, fc, fp = production.emit_counts(profile, libs, 0.1, np.random.default_rng(3),
                                        mode="quasi", spread=spread, overdispersion=4.0)
    floor_block = to_dense(fv, fc, fp, n, genes)
    delta = block.sum(0) - floor_block.sum(0)
    assert set(np.unique(delta).tolist()) <= {0, 1} and int(delta.sum()) == chosen.size


def test_injection_site_specificity():
    seam = SM.Seam()
    att = seam.attestation
    assert att["replacements"] == 1 and att["n_line_differences"] == 1
    assert att["line_differences"] == [[SM.ORIGINAL, SM.SEAMED]]
    assert att["floor_sites_before"] == 4 and att["untouched_floor_sites_after"] == 3
    assert att["branch"] == "mode == 'quasi'" and att["function"] == "emit_counts"
    assert att["production_sha256"] == SM.attestation()["production_sha256"]

    class Trap:
        calls = 0

        def __call__(self, _x):
            Trap.calls += 1
            raise AssertionError("hook reached outside the quasi branch")

    profile, libs, spread = quasi_fixture(seed=9)
    seamed = seam.emit_counts_with(Trap())
    tight = np.arange(len(profile)) % 2 == 0
    for seed, kwargs in ((1, dict(mode="even")),
                         (2, dict(mode="mixed", tight_mask=tight)),
                         (5, dict(mode="nb"))):
        a = production.emit_counts(profile, libs, 0.1, np.random.default_rng(seed), **kwargs)
        b = seamed(profile, libs, 0.1, np.random.default_rng(seed), **kwargs)
        assert same_csr(a, b), kwargs
    assert Trap.calls == 0
    try:
        seam._base(profile, libs, 0.1, np.random.default_rng(0), mode="quasi", spread=spread,
                   overdispersion=4.0)
        raise AssertionError("an unbound seam must refuse to emit")
    except RuntimeError:
        pass
    assert SM.HOOK not in production.emit_counts.__code__.co_names, "production untouched"
    assert production.emit_counts.__code__ is not seamed.__code__
    assert production.emit_counts.__globals__.get(SM.HOOK) is None, "no monkeypatch"


def test_rc_hybrid_end_to_end():
    control, profile, libs, spread, key = control_fixture()
    seam = SM.Seam()
    ref, ref_diag = rc_block(control, profile, libs, spread, key, production.emit_counts)
    floor = RD.BlockRounder("floor", 0)
    got, got_diag = rc_block(control, profile, libs, spread, key, seam.emit_counts_with(floor))
    assert np.array_equal(ref, got), "seamed floor block differs from production"
    assert floor.calls == 2
    assert {k: float(v) for k, v in ref_diag.items()} == {k: float(v) for k, v in got_diag.items()}
    x = len(libs) * production_L(libs) * profile
    assert np.array_equal(floor.x, x)
    assert RD.hamilton_choice(x).size > 0
    for rule in ("hamilton", "random_up"):
        rounder = RD.BlockRounder(rule, 0)
        block, _ = rc_block(control, profile, libs, spread, key, seam.emit_counts_with(rounder))
        assert rounder.calls == 2, rule
        assert np.array_equal(block.sum(0), rounder.totals), rule
        assert np.array_equal(rounder.x, x), rule
        delta = rounder.totals - np.floor(x).astype(np.int64)
        assert set(np.unique(delta).tolist()) <= {0, 1} and int(delta.sum()) == RD.residual_count(x)
        assert np.array_equal(block.sum(0), ref.sum(0) + delta), "new totals through placement"
        assert block.min() >= 0
        cell_totals = block.sum(1)
        assert cell_totals.max() - cell_totals.min() <= 1, "per-cell weld held"
        assert int(cell_totals.sum()) == int(rounder.totals.sum())


TESTS = [v for k, v in list(globals().items()) if k.startswith("test_") and callable(v)]


def run_all() -> None:
    failures = []
    for test in TESTS:
        try:
            test()
            print(f"PASS {test.__name__}", flush=True)
        except Exception:  # noqa: BLE001 - report every failure, then raise
            failures.append(test.__name__)
            print(f"FAIL {test.__name__}\n{traceback.format_exc()}", flush=True)
    if failures:
        raise AssertionError(f"342 pure tests failed: {failures}")


if __name__ == "__main__":
    run_all()
