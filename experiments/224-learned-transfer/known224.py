"""224 Phase 0c — the KNOWN-ANSWER rows.  Nothing downstream is read until these pass.

Reproduces, from `data/donors/*/delta.npy` and nothing else:

  A. 159 §5.2's consensus-of-k table on the six-line intersection (1,592 targets x 5,206
     genes), the `1/cos^2(k) = A + B/k` fit and its k -> infinity limit;
  B. 159 §5.2's same-line positive control (K562-gwps predicted from K562-essential:
     pds 0.9468, cos 0.2627) and its cross-line comparator (pds 0.7839, cos 0.0897);
  C. kappa = cos_inf / sqrt(rho) = 0.591 and tau = kappa^2 = 0.349 (159 §3.3 route 3).

AGENTS.md §11.1 rule "every screen must carry a configuration whose answer you already
know, and abort if it does not reproduce".

    uv run python known224.py
"""
from __future__ import annotations

import itertools
import json
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

from substrate224 import LINES, ROOT, WORK, load_effects  # noqa: E402

SEED = 20260902                       # 159's seed, so the panel subsample is the same draw
PANEL_N = 300
SIX = ["orion_hct116", "orion_hek293t", "k562_gwps", "k562_essential", "rpe1", "nadig_jurkat"]
SAME_LINE = {"k562_gwps": "K562", "k562_essential": "K562"}

# 159 RESULTS.md §5.2 — the numbers this must reproduce.
EXPECT = {
    "orion_hct116": dict(cos=[0.1065, 0.1382, 0.1579, 0.1718, 0.1822], inf=0.253,
                         pds=[0.7861, 0.8207, 0.8375, 0.8492, 0.8559]),
    "k562_gwps": dict(cos=[0.1132, 0.1492, 0.1728, 0.1900], inf=0.303),
    "nadig_jurkat": dict(cos=[0.1232, 0.1602, 0.1841, 0.2013, 0.2146]),
}
EXPECT_CTRL = dict(same_pds=0.9468, same_cos=0.2627, cross_pds=0.7839, cross_cos=0.0897)
TOL = 0.0015


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def pds(P, T):
    A, B = unit(P.astype(np.float64)), unit(T.astype(np.float64))
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S)
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


def endpoints(P, T):
    A, B = unit(P.astype(np.float64)), unit(T.astype(np.float64))
    return dict(pds=pds(P, T), cos=float((A * B).sum(1).mean()))


def fit_inv_sq(ks, coss):
    """1/cos^2(k) = A + B/k  ->  cos_inf = 1/sqrt(A).  159 §5.1, R^2 = 0.99998."""
    x = 1.0 / np.asarray(ks, dtype=float)
    y = 1.0 / np.asarray(coss, dtype=float) ** 2
    M = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    pred = M @ coef
    ss = 1.0 - ((y - pred) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-300)
    return float(1.0 / np.sqrt(max(coef[0], 1e-12))), float(coef[0]), float(coef[1]), float(ss)


def main() -> None:
    rng = np.random.default_rng(SEED)
    raw = {}
    for n in SIX:
        p = ROOT / LINES[n]["path"]
        raw[n] = (
            [x.strip() for x in (p / "perts.txt").read_text().splitlines() if x.strip()],
            [x.strip() for x in (p / "genes.txt").read_text().splitlines() if x.strip()])
    ts = sorted(set.intersection(*[set(v[0]) for v in raw.values()]))
    gs = sorted(set.intersection(*[set(v[1]) for v in raw.values()]) - set(ts))
    print(f"six-line intersection: {len(ts)} targets x {len(gs)} genes "
          f"(159 reported 1592 x 5206)")
    assert len(ts) == 1592 and len(gs) == 5206, "substrate drifted from 159's"

    E = {n: load_effects(LINES[n]["path"], ts, gs) for n in SIX}
    C = {n: v - v.mean(0, keepdims=True) for n, v in E.items()}
    sel = np.sort(rng.choice(len(ts), min(PANEL_N, len(ts)), replace=False))

    out, fails = {"n_targets": len(ts), "n_genes": len(gs), "seed": SEED}, []

    print(f"\nA. CONSENSUS-OF-k  (panel of {len(sel)} targets, equal weight on unit rows)")
    print(f"{'reference':16s}{'k':>3}{'pds':>9}{'cos':>10}   {'159 cos':>9} {'d':>9}")
    cons = {}
    for ref in SIX:
        others = [k for k in SIX if SAME_LINE.get(k, k) != SAME_LINE.get(ref, ref)]
        T = C[ref][sel]
        cons[ref] = []
        for k in range(1, len(others) + 1):
            combos = list(itertools.combinations(others, k))
            if len(combos) > 10:
                combos = [combos[i] for i in rng.choice(len(combos), 10, replace=False)]
            vals = [endpoints(np.mean([unit(C[o])[sel] for o in combo], axis=0), T)
                    for combo in combos]
            e = {m: float(np.mean([v[m] for v in vals])) for m in vals[0]}
            cons[ref].append(dict(k=k, **e))
            exp = EXPECT.get(ref, {}).get("cos", [])
            ec = exp[k - 1] if k - 1 < len(exp) else None
            d = "" if ec is None else f"{e['cos'] - ec:+9.4f}"
            print(f"{ref:16s}{k:>3}{e['pds']:9.4f}{e['cos']:10.5f}   "
                  f"{'' if ec is None else f'{ec:9.4f}'} {d}")
            if ec is not None and abs(e["cos"] - ec) > TOL:
                fails.append(f"consensus cos {ref} k={k}: {e['cos']:.5f} vs {ec}")
    out["consensus"] = cons

    print(f"\n   1/cos^2 = A + B/k  extrapolation")
    for ref in SIX:
        ks = [r["k"] for r in cons[ref]]
        cs = [r["cos"] for r in cons[ref]]
        cinf, A, B, r2 = fit_inv_sq(ks, cs)
        exp = EXPECT.get(ref, {}).get("inf")
        d = "" if exp is None else f"   159: {exp:.3f}  d {cinf - exp:+.4f}"
        print(f"   {ref:16s} cos_inf {cinf:.4f}   R^2 {r2:.6f}{d}")
        out.setdefault("kinf", {})[ref] = dict(cos_inf=cinf, A=A, B=B, r2=r2)
        if exp is not None and abs(cinf - exp) > 0.006:
            fails.append(f"cos_inf {ref}: {cinf:.4f} vs {exp}")

    print("\nB. SAME-LINE POSITIVE CONTROL")
    T = C["k562_gwps"][sel]
    same = endpoints(unit(C["k562_essential"])[sel], T)
    cross = endpoints(unit(C["orion_hct116"])[sel], T)
    print(f"   K562-gwps <- K562-essential (same LINE): pds {same['pds']:.4f}  "
          f"cos {same['cos']:.5f}   159: 0.9468 / 0.2627")
    print(f"   K562-gwps <- HCT116        (cross-line): pds {cross['pds']:.4f}  "
          f"cos {cross['cos']:.5f}   159: 0.7839 / 0.0897")
    out["same_line_control"] = dict(same=same, cross=cross)
    for k, v, e in (("same_pds", same["pds"], EXPECT_CTRL["same_pds"]),
                    ("same_cos", same["cos"], EXPECT_CTRL["same_cos"]),
                    ("cross_pds", cross["pds"], EXPECT_CTRL["cross_pds"]),
                    ("cross_cos", cross["cos"], EXPECT_CTRL["cross_cos"])):
        if abs(v - e) > TOL:
            fails.append(f"control {k}: {v:.5f} vs {e}")

    print("\nC. kappa / tau  (159 §3.3 route 3)")
    rho = same["cos"]                       # two libraries of one line, true cosine 1
    cinf = out["kinf"]["k562_gwps"]["cos_inf"]
    kappa = cinf / np.sqrt(rho)
    print(f"   rho = same-line cos = {rho:.4f} ;  cos_inf(K562) = {cinf:.4f}")
    print(f"   kappa = {kappa:.4f}  (159: 0.591)   tau = kappa^2 = {kappa ** 2:.4f}  "
          f"(159: 0.349)")
    out["kappa"] = float(kappa)
    out["tau"] = float(kappa ** 2)
    if abs(kappa - 0.591) > 0.02:
        fails.append(f"kappa: {kappa:.4f} vs 0.591")

    out["fails"] = fails
    (WORK / "known224.json").write_text(json.dumps(out, indent=1))
    print("\n" + ("KNOWN-ANSWER ROWS PASS" if not fails else
                  "KNOWN-ANSWER FAILURES:\n  " + "\n  ".join(fails)))
    if fails:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
