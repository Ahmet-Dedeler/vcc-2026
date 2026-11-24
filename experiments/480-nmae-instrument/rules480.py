"""480 step 3a: in-sample upper bound of class-based magnitude rules on H1 (our real construction).

Rule shape per class C (oracle-free classes only): lhat = m_C * base + s * c_C * min(1/sqrt(cnt), 1) * [cnt >= 0.3]
  shipped = (m 1, c 0.1) everywhere (= the coarse MDE's analytic form; 434's KA row).
The objective is separable across classes (nmae is additive over cells within a target), so the in-sample optimum is a
per-class grid search. The gain vs shipped is an UPPER BOUND for any rule of this shape on this class variable.
Also reports the analytic shipped row vs the emitted coarse arm (0.94270) as the known answer.

    uv run python rules480.py     # -> work/rules480.json
"""
from __future__ import annotations

import json

import numpy as np

from common480 import HERE, nmae
from decomp480 import contributions, h1_inputs, tss_classes

M_GRID = np.linspace(0, 2.5, 21)
C_GRID = np.linspace(0, 0.5, 26)


def class_vars(T, Lc, qc, cnt, cons, targets, genes, s):
    G = np.isfinite(T)
    q1, q2 = np.quantile(np.broadcast_to(cnt, T.shape)[G], [1 / 3, 2 / 3])
    expr = np.broadcast_to(np.where(cnt < q1, 0, np.where(cnt < q2, 1, 2)).astype(np.int8), T.shape)
    cis = tss_classes(targets, genes)
    cis2 = np.where(cis <= 1, 0, 1).astype(np.int8)
    strength = (np.abs(cons) > 0.2).sum(1)                          # donor-consensus strength (oracle-free)
    st = np.digitize(strength, np.quantile(strength, [1 / 3, 2 / 3]))
    mag = np.digitize(np.abs(cons), np.quantile(np.abs(cons)[G], [1 / 3, 2 / 3])).astype(np.int8)
    agree = np.where(cons == 0, 2, np.where(np.sign(cons) == s, 0, 1)).astype(np.int8)
    called = np.where(qc < 0.05, 0, 1).astype(np.int8)
    direc = np.where(s > 0, 0, 1).astype(np.int8)
    return {
        "none (one class)": (np.zeros(T.shape, np.int8), 1),
        "expression tertile": (expr, 3),
        "cis <100 kb vs rest": (cis2, 2),
        "donor-consensus strength tertile (target)": (np.broadcast_to(st[:, None], T.shape).astype(np.int8), 3),
        "predicted direction": (direc, 2),
        "called by us": (called, 2),
        "donor |cons| tertile": (mag, 3),
        "donor agrees with our sign": (agree, 3),
        "direction x agreement": ((direc * 3 + agree).astype(np.int8), 6),
        "direction x expression": ((direc * 3 + expr).astype(np.int8), 6),
        "called x direction x agreement": ((called * 6 + direc * 3 + agree).astype(np.int8), 12),
    }


def fit_classes(T, base, s, add1, cls, n, shape="m+c"):
    """Per-class in-sample optimum of (m, c). Returns (m_C, c_C, contribution arrays)."""
    G = np.isfinite(T)
    Tz = np.where(G, T, 0.0)
    D = np.abs(Tz).sum(1)
    D[D == 0] = np.nan
    ms = M_GRID if shape == "m+c" else np.array([1.0])
    best = np.full(n, np.inf)
    bm, bc = np.ones(n), np.full(n, 0.1)
    for m in ms:
        for c in C_GRID:
            v = contributions(np.abs(m * base + s * c * add1 - Tz), D, G, cls, n)
            better = v < best
            best[better], bm[better], bc[better] = v[better], m, c
    return bm, bc, best


def apply(base, s, add1, cls, bm, bc):
    return bm[cls] * base + s * bc[cls] * add1


def main():
    h, T, cnt, genes, targets, Lc, qc, L0, cons = h1_inputs()
    s = np.sign(Lc)
    s = np.where(s == 0, np.sign(cons), s)
    s = np.where(s == 0, 1.0, s)
    add1 = (np.minimum(1 / np.sqrt(np.maximum(cnt, 1e-3)), 1.0) * (cnt >= 0.3))[None, :]
    shipped = nmae(L0 + s * 0.1 * add1, T)
    print(f"KA: analytic shipped MDE on base {shipped:.5f} (434's row: 0.95237 - 0.00982 = 0.94255); "
          f"emitted coarse {nmae(Lc, T):.5f}")
    assert abs(shipped - 0.94255) < 2e-4
    out = {"shipped_analytic": shipped, "emitted_coarse": nmae(Lc, T)}
    for name, (cls, n) in class_vars(T, Lc, qc, cnt, cons, targets, genes, s).items():
        row = {}
        for shape in ["c", "m+c"]:
            bm, bc, _ = fit_classes(T, L0, s, add1, cls, n, shape)
            v = nmae(apply(L0, s, add1, cls, bm, bc), T)
            row[shape] = dict(nmae=v, gain=v - shipped, m=bm.tolist(), c=bc.tolist())
        out[name] = row
        print(f"{name:42s} c-only {row['c']['gain']:+.4f} (c {np.round(row['c']['c'], 2).tolist()}) | "
              f"m+c {row['m+c']['gain']:+.4f} (m {np.round(row['m+c']['m'], 2).tolist()}, "
              f"c {np.round(row['m+c']['c'], 2).tolist()})")
    json.dump(out, open(HERE / "work/rules480.json", "w"), indent=1)


if __name__ == "__main__":
    main()
