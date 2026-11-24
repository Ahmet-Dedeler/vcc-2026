#!/usr/bin/env python3
"""630: regression test for the 605 cloud-review blocker (negative-gain donors), plus the fail-safe guard.

The bug (pre-fix module, kept byte-for-byte at work/src_rank_up_prefix605.py, sha 8e1f05a6...): moving one count of
gene j from ballast row b to carrier a changes the gene's mean CPM by cw[a] - cw[b], cw = 1 / (n * library / 1e6). That is
positive only when the donor's library is LARGER than the carrier's. The pre-fix code clamped the ask's denominator at
1e-300 and let ZS.fill take from every donor, so with shallow donors the mean could FALL and an UP lfc flip DOWN while the
library, per-gene total and |dAUC| invariants all held. Production blocks cannot trigger it (ballast ~995k vs carriers
~20k counts on all 85 H1 targets; see RESULTS), but the module must be correct for any block.

  N1  the review's case (8 carriers at 22,000, 32 ballast at 18,000, every donor negative-gain): the PRE-fix module
      lowers the mean and flips the sign (the test can see the bug); the fixed module moves nothing, mean and sign
      unchanged, n_no_pos_donor counted
  N2  mixed donors (half deep, half shallow): the fixed module takes ONLY from deep donors, the mean rises, the sign holds,
      libraries and per-gene totals are exact
  N3  production-like (every donor deep): fixed == pre-fix bit for bit (the fix is a no-op where it cannot bind)
  N4  the guard fails SAFE: if calls still flip after GUARD_ROUNDS, apply returns the input block unchanged with
      ur_guard_failed = 1 instead of raising
  N5  any invariant assert inside the pass fails SAFE too (ur_invariant_failed = 1), same key set ([7772dc]'s review)
  N6  every path returns the same diagnostic keys (DIAG_KEYS), so the builder's per-key meta aggregation is uniform
Synthetic, < 1 s, < 200 MB.   uv run python experiments/630-rank-neutral-up/test_neggain630.py   (exit 0 = all PASS)
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS: list[tuple[str, bool, str]] = []


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


NEW = _load("rank_up_new", HERE / "src_rank_up.py")
OLD = _load("rank_up_old", HERE / "work/src_rank_up_prefix605.py")


def rec(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}", flush=True)


def make_case(bal_libs, n_car=8, car_lib=22_000, n_ctrl=400, seed=0):
    """Block (n_car + len(bal_libs)) x 5 genes: kept genes 0 (tested), 1, 2; sub-gate genes 3, 4 (library filler).
    Gene 0: carriers hold 2 counts (CPM ~91), one control just above the carriers, ballast holds a pile above the
    control max. ref_lin is set so the base lfc is +0.00144 (the review's case)."""
    rng = np.random.default_rng(seed)
    nb = len(bal_libs)
    n = n_car + nb
    X = np.zeros((n, 5), dtype=np.int64)
    X[:n_car, 0] = 2
    X[:n_car, 1] = X[:n_car, 2] = 1000
    X[:n_car, 3] = 10_000
    X[:n_car, 4] = car_lib - X[:n_car, :4].sum(1)
    for i, L in enumerate(bal_libs):
        r = n_car + i
        X[r, 0] = max(int(round(400 * L / 18_000)), 3)      # same CPM pile (~22k) at every depth
        X[r, 1] = X[r, 2] = int(L * 0.05)
        X[r, 3] = int(L * 0.4)
        X[r, 4] = L - X[r, :4].sum()
    assert (X >= 0).all()
    s = X.sum(1) / 1e6
    A = (X[:, :3] / s[:, None]).mean(0)
    ref = A.copy()
    ref[0] = A[0] / 2 ** 0.00144
    vcar = np.log1p(2 / s[0])
    cs0 = np.sort(np.r_[rng.uniform(0.0, vcar - 0.1, n_ctrl - 1), vcar + 0.05])
    cs = np.stack([cs0, np.sort(rng.uniform(0, 12, n_ctrl)), np.sort(rng.uniform(0, 12, n_ctrl))])
    cc = SimpleNamespace(keep=np.array([True, True, True, False, False]), ref_lin=ref, cs=cs)
    l_base = np.log2((A + 1e-9) / (ref + 1e-9))
    sel = np.array([True, False, False])
    target = l_base + np.array([1.0, 0.0, 0.0])
    return X, n_car, cc, l_base, target, sel


def mean_lfc(X, cc):
    s = X.sum(1) / 1e6
    A = (X[:, :3] / s[:, None]).mean(0)
    return A, np.log2((A + 1e-9) / (cc.ref_lin + 1e-9))


def run(mod, case, **kw):
    X, n_car, cc, l_base, target, sel = case
    return mod.up_pass2(X, n_car, cc, l_base, target, sel, **kw)


def main() -> int:
    for cross in (False, True):
        tag = f"cross={cross}"
        # N1: every donor negative-gain
        case = make_case([18_000] * 32)
        X0, _, cc, l_base, *_ = case
        A0, L0 = mean_lfc(X0, cc)
        Xo, so = run(OLD, case, cross=cross)
        Ao, Lo = mean_lfc(Xo, cc)
        rec(f"N1 {tag} pre-fix module shows the bug (mean falls, UP flips DOWN)",
            Ao[0] < A0[0] and np.sign(Lo[0]) != np.sign(L0[0]),
            f"lfc {L0[0]:+.5f} -> {Lo[0]:+.5f}, counts {so['counts']}")
        Xn, sn = run(NEW, case, cross=cross)
        An, Ln = mean_lfc(Xn, cc)
        rec(f"N1 {tag} fixed module: nothing moves, mean and sign unchanged",
            np.array_equal(Xn, X0) and sn["n_no_pos_donor"] == 1 and sn["n_moved"] == 0,
            f"lfc {Ln[0]:+.5f}, n_no_pos_donor {sn['n_no_pos_donor']}, n_neg_donor_rows {sn['n_neg_donor_rows']}")
        # N2: mixed donors
        libs = [18_000] * 16 + [995_000] * 16
        case = make_case(libs)
        X0, n_car, cc, *_ = case
        A0, L0 = mean_lfc(X0, cc)
        Xn, sn = run(NEW, case, cross=cross)
        An, Ln = mean_lfc(Xn, cc)
        changed = np.flatnonzero(np.any(Xn != X0, axis=1))
        donors = changed[changed >= n_car]
        shallow_touched = [int(r) for r in donors if X0[r].sum() < X0[:n_car].sum(1).max()]
        rec(f"N2 {tag} fixed module: only deep donors give, mean rises, sign holds, sums exact",
            sn["n_moved"] == 1 and not shallow_touched and An[0] > A0[0] and np.sign(Ln[0]) == np.sign(L0[0])
            and np.array_equal(Xn.sum(1), X0.sum(1)) and np.array_equal(Xn.sum(0), X0.sum(0)),
            f"lfc {L0[0]:+.5f} -> {Ln[0]:+.5f}, donors {donors.size}, shallow touched {len(shallow_touched)}")
        Xo, so = run(OLD, case, cross=cross)
        Ao, Lo = mean_lfc(Xo, cc)
        rec(f"N2 {tag} (info) pre-fix on the mixed case", True, f"lfc {L0[0]:+.5f} -> {Lo[0]:+.5f}")
        # N3: production-like
        case = make_case([995_000] * 4, n_car=396)
        Xo, so = run(OLD, case, cross=cross)
        Xn, sn = run(NEW, case, cross=cross)
        A0, L0 = mean_lfc(case[0], case[2])
        An, Ln = mean_lfc(Xn, case[2])
        rec(f"N3 {tag} all-deep donors: fixed == pre-fix bit for bit, mean rises",
            np.array_equal(Xn, Xo) and sn["n_neg_donor_rows"] == 0 and sn["n_moved"] == so["n_moved"] == 1
            and An[0] > A0[0], f"lfc {L0[0]:+.5f} -> {Ln[0]:+.5f}, counts {sn['counts']}")
    # N4: the guard fails safe (control flow only: up_pass2 and de_stats stubbed so calls always flip)
    X0 = np.arange(12, dtype=np.int64).reshape(4, 3)
    calls = {"n": 0}

    def fake_de_stats(X, cc, tie):
        calls["n"] += 1
        q = np.array([0.01, 0.5, 0.5]) if calls["n"] == 1 else np.array([0.5, 0.5, 0.5])   # base call, then always flipped
        return np.ones(3), q, np.array([0.1, 0.1, 0.1])

    def fake_up_pass2(X, n_car, cc, l0, t, sel, **kw):
        st = dict(n_sel=1, n_moved=1, n_reverted=0, n_no_candidate=0, counts=5, need=1.0, got=0.5, n_no_pos_donor=0,
                  n_neg_donor_rows=0, n_mean_revert=0, max_dAUC=1e-5, median_dAUC=1e-5, frac_got=0.5, capU=10)
        return X + 1, st

    saved = NEW.de_stats, NEW.up_pass2
    NEW.de_stats, NEW.up_pass2 = fake_de_stats, fake_up_pass2
    try:
        cc = SimpleNamespace(ref_lin=np.ones(3))
        Xg, dg = NEW.apply(X0, 2, cc, "b_k100", tie=None, lib_med=20_000.0)
        ok = np.array_equal(Xg, X0) and dg.get("ur_guard_failed") == 1.0 and dg.get("ur_n_moved") == 0.0
        rec("N4 guard exhausted -> input block returned, ur_guard_failed 1, no raise", ok,
            f"guard_rounds {dg.get('ur_guard_rounds')}, keys {len(dg)}")
    except AssertionError as e:
        rec("N4 guard exhausted -> input block returned, ur_guard_failed 1, no raise", False, f"raised {e}")
    finally:
        NEW.de_stats, NEW.up_pass2 = saved
    # N5: any invariant assert inside the pass fails SAFE with the same key set as a normal call ([7772dc]'s review)
    def boom(*a, **k):
        raise AssertionError("a library moved")

    calls["n"] = 0
    saved = NEW.de_stats, NEW.up_pass2
    NEW.de_stats, NEW.up_pass2 = fake_de_stats, boom
    try:
        Xi, di = NEW.apply(X0, 2, SimpleNamespace(ref_lin=np.ones(3)), "b_k100", tie=None, lib_med=20_000.0)
        ok = np.array_equal(Xi, X0) and di.get("ur_invariant_failed") == 1.0 and set(di) == set(dg)
        rec("N5 invariant assert -> input block returned, ur_invariant_failed 1, same keys as the guard path", ok,
            f"keys {len(di)} vs {len(dg)}")
    except AssertionError as e:
        rec("N5 invariant assert -> input block returned, ur_invariant_failed 1, same keys as the guard path", False,
            f"raised {e}")
    finally:
        NEW.de_stats, NEW.up_pass2 = saved
    # N6: DIAG_KEYS (the failure path's keys) == the guard path's keys == a normal call's keys
    rec("N6 DIAG_KEYS == the guard path's key set", {f"ur_{k}" for k in NEW.DIAG_KEYS} == set(dg),
        f"{sorted({f'ur_{k}' for k in NEW.DIAG_KEYS} ^ set(dg))}")
    # N7: a real apply() call (real de_stats, scorer tie sums) returns exactly DIAG_KEYS, invariant_failed 0
    X, n_car, cc, *_ = make_case([995_000] * 4, n_car=396)
    try:
        Xr, dr = NEW.apply(X, n_car, cc, "b_k100", tie=NEW.scorer_tie(cc), lib_med=20_000.0)
        rec("N7 real apply: keys == DIAG_KEYS, invariant_failed 0, guard_failed 0",
            set(dr) == {f"ur_{k}" for k in NEW.DIAG_KEYS} and dr["ur_invariant_failed"] == 0.0 and dr["ur_guard_failed"] == 0.0,
            f"moved {dr.get('ur_n_moved')}, sym diff {sorted(set(dr) ^ {f'ur_{k}' for k in NEW.DIAG_KEYS})}")
    except Exception as e:  # noqa: BLE001
        rec("N7 real apply: keys == DIAG_KEYS, invariant_failed 0, guard_failed 0", False, f"{type(e).__name__}: {e}")
    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} PASS")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
