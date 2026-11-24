#!/usr/bin/env python3
"""198 — adversarial check on 193's premise.

193's MANDATE says production whitens with `(C + lam I)^(-gamma/2)` and calls
`C + lam I` "the crudest possible covariance estimator: it shrinks every
eigenvalue by the same additive constant".

`scripts/make_submission_pm_cd4.py:493-494` ships WHITEN_GAMMA=0.65 and
WHITEN_LAM=**0.0**, and its own comment says lam "has sat at its 0.0 default in
every build this project has ever shipped". So there IS no ridge. The shipped
map is a pure power law with a RELATIVE FLOOR at ev.max()*1e-12
(`src/vcc2026/whiten.py:65`).

Measure what that map actually does to the spectrum. Read-only.
"""
from __future__ import annotations
import hashlib
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(exist_ok=True)

GAMMA, LAM, FLOOR_REL = 0.65, 0.0, 1e-12
OPS = [
    ("union_pruned (PRODUCTION)", "data/donors/whitener_union_pruned.npz"),
    ("union", "data/donors/whitener_union.npz"),
    ("union_full", "data/donors/whitener_union_full.npz"),
]

report = {}
for name, rel in OPS:
    p = ROOT / rel
    if not p.exists():
        print(f"MISSING {rel}")
        continue
    sha = hashlib.sha256(p.read_bytes()).hexdigest()      # rule 11: sha, not a name
    ev = np.load(p)["evals"].astype(np.float64)
    ev = np.sort(ev)[::-1]
    G = ev.size
    evmax = ev.max()
    floor = evmax * FLOOR_REL

    # whiten.py:65 applies the floor at FIT time; re-apply so we see the shipped map
    evf = np.maximum(ev, floor)
    w = (evf + LAM) ** (-GAMMA / 2.0)

    n_at_floor = int((ev <= floor).sum())
    # how much of C's TRACE lives in directions the map amplifies above 1.0?
    amp = w > 1.0
    print(f"\n=== {name} ===")
    print(f"  path {rel}")
    print(f"  sha256 {sha}")
    print(f"  G = {G:,}   max_ev = {evmax:.6f}   floor = evmax*1e-12 = {floor:.3e}")
    print(f"  eigenvalues AT the relative floor: {n_at_floor:,} "
          f"({100 * n_at_floor / G:.2f}% of the spectrum)")
    print(f"  applied weight  w = ev^(-{GAMMA / 2}):")
    print(f"     w(ev_max)  = {w[0]:.6g}")
    print(f"     w(median)  = {np.median(w):.6g}")
    print(f"     w(ev_min)  = {w[-1]:.6g}")
    print(f"  >>> CONDITION NUMBER OF THE APPLIED MAP = w_max/w_min = "
          f"{w.max() / w.min():,.0f}x")
    print(f"  directions AMPLIFIED (w>1): {int(amp.sum()):,} "
          f"({100 * amp.mean():.1f}%), carrying {100 * ev[amp].sum() / ev.sum():.4f}% "
          f"of trace(C)")
    # the tail that a ridge would kill
    for q in (1e-6, 1e-8, 1e-10, 1e-12):
        thr = evmax * q
        k = int((ev <= thr).sum())
        if k:
            wq = (max(thr, floor)) ** (-GAMMA / 2.0)
            print(f"     ev <= evmax*{q:.0e}: {k:6,d} dirs "
                  f"({100 * k / G:5.2f}%), amplified >= {wq / w[0]:,.0f}x the top")

    report[name] = dict(
        path=rel, sha256=sha, G=G, max_ev=float(evmax),
        n_at_floor=n_at_floor, frac_at_floor=n_at_floor / G,
        w_top=float(w[0]), w_min=float(w.min()), w_max=float(w.max()),
        condition_number=float(w.max() / w.min()),
        frac_amplified=float(amp.mean()),
        trace_share_amplified=float(ev[amp].sum() / ev.sum()),
    )

print("\n" + "=" * 74)
print("VERDICT for 193: the shipped operator carries NO ridge (lam = 0.0).")
print("Its eigenvalue map is a pure power law ev^(-0.325) with a RELATIVE FLOOR")
print("at 1e-12, so it is already a nonlinear shrinkage, not 'the crudest")
print("possible covariance estimator'. The identity row must be lam=0, and the")
print("untested lever is the FLOOR / RIDGE, not a replacement for one.")
json.dump(report, open(OUT / "whitener198.json", "w"), indent=1)
