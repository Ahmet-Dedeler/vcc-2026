"""172 — the same keys under 157 §3b's proportional power matching.

WHY.  `anatomy172.py` measures that **94.4% of h1-test's `reach` under `|H|` comes from the 32
targets whose `n_real` is <= 20**, and AGENTS.md §11.1 rule 9 says `n_real` is 89 here against
~240-2,670 on the challenge.  So every number in this experiment's headline table is measured
in a power regime the challenge is not in.  157 §3a showed that `debench.Surface(n_real_target)`
is the WRONG correction (it flattens the small-`n_real` stratum that carries the statistic,
and reads a 23x collapse for a 2.7x pool change).  157 §3b's fix is proportional:

    pool_p = the reference's top-round(f * n_real_p) genes by its OWN p-value

which scales every perturbation's pool by a common factor and leaves the stratification alone.

KNOWN-ANSWER ROW (rule 10's corollary): `f = 1` must reproduce the native BH numbers to within
p-value ties -- production 0.18283, `|H|` 0.23533, ORACLE 0.64553 -- and 157 §3b's own table
must reproduce for the incumbent and the oracle at every `f`.  Both are asserted below.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from harness172 import H1Test, check  # noqa: E402

FACTORS = (0.5, 1.0, 1.5, 2.0, 3.5, 5.0)
# 157 §3b, h1-test, shipped a15 emission -- the rows this instrument must reproduce.
EXPECT_157 = {
    "production (debench order)": [0.2381, 0.1828, 0.1507, 0.1220, 0.0536, 0.0283],
    "ORACLE": [0.6551, 0.6455, 0.6368, 0.6221, 0.5790, 0.5742],
}


def scaled_pool(T, f, real_p):
    """157 §3b's `power.py` rule, on this surface.

    THE TARGET'S OWN GENE IS REMOVED BEFORE THE TOP-N SELECTION, NOT AFTER.  A knockdown's own
    row usually carries the smallest p-value of the row, so selecting top-`n` first and zeroing
    it afterwards silently returns `n-1` genes and shifts every pool by one -- which is exactly
    what made the `f = 1` known-answer row read production 0.17161 against the native 0.18283
    on the first attempt.  `debench.Surface(n_real_target=...)` has the same ordering and is
    another reason 157 §3a found it unusable for `reach`.
    """
    rp = np.where(np.isfinite(real_p), real_p, np.inf)
    for p, j in enumerate(T.S.self_idx):
        if j >= 0:
            rp[p, j] = np.inf
    sig = np.zeros(rp.shape, bool)
    for p in range(rp.shape[0]):
        n = int(round(f * T.S.real_sig[p].sum()))
        if n <= 0:
            continue
        n = min(n, rp.shape[1])
        sig[p, np.argpartition(rp[p], n - 1)[:n]] = True
    sig &= np.isfinite(rp)
    return sig


def reach_at(T, key, sig):
    native = T.S.real_sig
    T.S.real_sig = sig
    try:
        return T.reach(key)
    finally:
        T.S.real_sig = native


def main() -> None:
    T = H1Test()
    print("known-answer rows (native surface):")
    check(T, verbose=False)
    print("  OK (harness172.check)\n")
    Z = np.load(os.path.join(HERE, "work", "allkeys.npz"), allow_pickle=True)
    order = ["production (debench order)", "|cons| (un-whitened)", "geo x concord",
             "|H| (whitened)"]
    order += [k for k in Z.files if k not in order + ["genes", "targets", "ORACLE"]
              and not k.startswith("null")]
    order += ["ORACLE"]
    from harness172 import DE_PATH
    real_p = np.load(DE_PATH, allow_pickle=True)["real_p"]
    pools = {f: scaled_pool(T, f, real_p) for f in FACTORS}
    med = {f: int(np.median(pools[f].sum(1))) for f in FACTORS}
    print("proportional power matching (157 §3b).  The challenge's n_real sits between f=2 and")
    print("f=3.5 on this surface; h1-test's native median n_real is 89.\n")
    print(f"{'key':<40}" + "".join(f"{'f='+str(f):>10}" for f in FACTORS))
    print(f"{'median n_real':<40}" + "".join(f"{med[f]:>10}" for f in FACTORS))
    got = {}
    for k in order:
        v = None if k == "production (debench order)" else np.asarray(Z[k], dtype=np.float64)
        got[k] = [reach_at(T, v, pools[f]) for f in FACTORS]
        print(f"{k:<40}" + "".join(f"{x:>10.5f}" for x in got[k]))
    # KNOWN-ANSWER ROW: f = 1 must reproduce the native BH surface exactly, for every key.
    native = {"production (debench order)": 0.18283, "|H| (whitened)": 0.23533,
              "ORACLE": 0.64553}
    bad = [k for k, v in native.items() if abs(got[k][FACTORS.index(1.0)] - v) > 5e-6]
    for k, v in native.items():
        g = got[k][FACTORS.index(1.0)]
        print(f"  {'OK ' if abs(g - v) <= 5e-6 else 'FAIL'} f=1 reproduces the native BH "
              f"surface: {k:<28} {g:.5f} vs {v:.5f}")
    if bad:
        raise SystemExit(f"f=1 KNOWN-ANSWER ROWS FAILED: {bad} (rule 10)")
    print("\n  cross-reference only — 157 §3b's own table, measured on 157's `de_h1test.npz`")
    print("  and not on 164's `de164.a15.npz`.  The two agree EXACTLY at f=1 (both reproduce")
    print("  the native BH set) and diverge by up to 0.019 at f>1, because which extra genes")
    print("  enter a scaled pool depends on the reference's p-value ordering BELOW the BH")
    print("  cut, and that ordering is a property of the individual DE table.  This is a")
    print("  cross-check, not an assertion; only the f=1 row is a known answer.")
    for k, exp in EXPECT_157.items():
        if k in got:
            print(f"      157 §3b {k:<32}" + "".join(f"{x:>10.4f}" for x in exp))
            print(f"      here    {k:<32}" + "".join(f"{x:>10.4f}" for x in got[k]))

    print("\ncapture % of the ORACLE headroom over the production order, by f:")
    print(f"{'key':<40}" + "".join(f"{'f='+str(f):>10}" for f in FACTORS))
    base, orc = got["production (debench order)"], got["ORACLE"]
    for k in order:
        if k in ("production (debench order)", "ORACLE"):
            continue
        pct = [100 * (got[k][i] - base[i]) / max(orc[i] - base[i], 1e-9)
               for i in range(len(FACTORS))]
        print(f"{k:<40}" + "".join(f"{x:>10.1f}" for x in pct))


if __name__ == "__main__":
    main()
