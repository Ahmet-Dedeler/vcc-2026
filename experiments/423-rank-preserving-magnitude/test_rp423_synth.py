#!/usr/bin/env python3
"""Synthetic invariant tests for rp423 with K > 1 ballast rows (the production 410b/417a bulk has K 1-16; the cap15
H1 bulk only exercises K = 1). Random blocks with mass-holder genes (a few carriers far above every control, the rest
zeros tied with control zeros) and dense genes; checks totals, libraries, non-negativity, no sign flips, untouched
genes' rank sums, and exact rank sums everywhere for zs/ex.  usage: uv run python test_rp423_synth.py"""
from __future__ import annotations

import sys
from types import SimpleNamespace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rp423  # noqa: E402


def make(seed, n_car=120, Kb=3, G=60, m=600):
    rng = np.random.default_rng(seed)
    keep = np.zeros(G, dtype=bool); keep[:40] = True
    # controls: libraries 3k..50k, poisson counts from gene rates; ~half the kept genes zero-inflated
    rate = np.where(np.arange(G) < 20, rng.uniform(1e-4, 5e-3, G), rng.uniform(2e-5, 3e-4, G))
    Lc = rng.integers(3000, 50000, m)
    C = rng.poisson(rate[None, :] * Lc[:, None])
    Lc = C.sum(1)
    V = np.log1p(C / (Lc / 1e6)[:, None])[:, keep].T
    V.sort(axis=1)
    cs = np.ascontiguousarray(V)
    cc = SimpleNamespace(keep=keep, cs=cs, tie=np.zeros(cs.shape[0]), ref_lin=(C / (Lc / 1e6)[:, None]).mean(0)[keep])
    # carriers: three library classes; dense genes poisson; zero-inflated genes: 3-6 mass holders far above the max
    lib = np.repeat([16000, 20000, 25000], n_car // 3 + 1)[:n_car]
    X = rng.poisson(rate[None, :] * lib[:, None] * rng.uniform(0.5, 1.5, G)[None, :])
    for g in range(20, 40):
        X[:, g] = 0
        hold = rng.choice(n_car, rng.integers(3, 7), replace=False)
        X[hold, g] = rng.integers(40, 120, hold.size)
    X[:, 40:] += rng.poisson(3.0, (n_car, G - 40))            # sub-gate reservoir in every carrier
    B = rng.poisson(rate * 900_000, (Kb, G)) + rng.poisson(200.0, (Kb, G)) * (~keep)
    X = np.vstack([X, B]).astype(np.int64)
    return X, n_car, cc


def lfc_of(X, cc):
    s = X.sum(1) / 1e6
    A = (X[:, cc.keep] / s[:, None]).mean(0)
    return np.log2((A + 1e-9) / (cc.ref_lin + 1e-9))


def R1(X, cc):
    s = X.sum(1) / 1e6
    out = np.zeros(cc.cs.shape[0])
    for j, col in enumerate(np.flatnonzero(cc.keep)):
        out[j] = rp423.ranks_of(X[:, col], s, cc.cs[j]).sum()
    return out


def main():
    fails = 0
    moved = {}
    for seed in range(6):
        for Kb in (1, 3, 7):
            X, n_car, cc = make(seed, Kb=Kb)
            l = lfc_of(X, cc)
            r0 = R1(X, cc)
            for v in ("all_k160", "dn_k250", "all_k160_zs", "all_k160_ex", "dn_k250_ex", "all_k200u120_ex", "all_k100_ex"):
                try:
                    X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
                except AssertionError as e:
                    print(f"[FAIL] seed {seed} K {Kb} {v}: assertion {e}"); fails += 1; continue
                l2 = lfc_of(X2, cc); r1 = R1(X2, cc); t = dg["_touched"]
                fin = np.isfinite(l) & np.isfinite(l2)
                ok = {"totals": np.array_equal(X2.sum(0), X.sum(0)), "libs": np.array_equal(X2.sum(1), X.sum(1)),
                      "nonneg": X2.min() >= 0, "signs": not np.any(np.sign(l[fin]) != np.sign(l2[fin])),
                      "untouched R1": np.array_equal(r1[~t], r0[~t])}
                if v.endswith(("_zs", "_ex")):
                    ok["R1 exact"] = np.array_equal(r1, r0)
                if v == "all_k100_ex":
                    ok["no-op"] = np.array_equal(X2, X)
                bad = [k for k, b in ok.items() if not b]
                fails += bool(bad)
                if bad:
                    print(f"[FAIL] seed {seed} K {Kb} {v}: {bad}")
                moved[v] = moved.get(v, 0) + int(dg["_touched"].sum())
    print("genes moved per variant over all synthetic blocks:", moved)
    print(f"{fails} failed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
