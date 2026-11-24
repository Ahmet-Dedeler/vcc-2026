"""481 known answer: re-derive the prior gate's design strength and its q-quantile threshold from a stack dump.

The builder (468, VCC_DE_PRIOR_GATE) computes, per context and target, strength = count(|gated design row| > tau) with
tau 0.2, the gate being the context's render gate (pooled control CPM > min_cpm 10 and b_ctrl > 0), and the design row
being the carrier matrix (after the same-line blend; no targets are blended on these bundles except efr's D, 9 rows).
The dump holds the H1 x2 carriers (the stack-dump build), not the level builds' H1 x1 carriers, so this is an
approximation; it is checked against the builder's own A/B/C numbers in 469dr's meta (threshold 2, weak 86/83/83,
strength quartiles [0,1,3,7,max]).

usage: uv run python strength481.py <dump.npz> <ctrl tag prefix, e.g. def> [--q 0.33] [--tau 0.2] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def render_gate(tag: str, genes: np.ndarray, min_cpm: float = 10.0) -> np.ndarray:
    z = np.load(HERE / "work" / f"ctrl_{tag}.npz", allow_pickle=False)
    assert np.array_equal(z["genes"].astype(str), genes), f"gene axis differs for {tag}"
    tot = z["totals"]
    cpm = 1e6 * tot / tot.sum()
    b = np.log1p(5e4 * tot / tot.sum())
    return (cpm > min_cpm) & (b > 0)


def strength(carriers: np.ndarray, gate: np.ndarray, tau: float) -> np.ndarray:
    return (np.abs(np.where(gate[None, :], carriers, 0.0)) > tau).sum(1).astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("prefix")
    ap.add_argument("--q", type=float, default=0.33)
    ap.add_argument("--tau", type=float, default=0.2)
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    z = np.load(a.dump, allow_pickle=False)
    genes = z["genes"].astype(str)
    car = np.asarray(z["carriers"], dtype=np.float64)
    out = {}
    for c in ("D", "E", "F"):
        s = strength(car, render_gate(f"{a.prefix}_{c}", genes), a.tau)
        thr = float(np.quantile(s, a.q))
        weak = int((s < thr).sum())
        hist = {int(k): int((s == k).sum()) for k in range(0, 6)}
        out[c] = dict(threshold=thr, weak=weak, le1=int((s <= 1).sum()), eq0=int((s == 0).sum()),
                      q=[float(x) for x in np.quantile(s, [0, .25, .5, .75, 1])], hist0_5=hist,
                      frac_le1=float((s <= 1).mean()))
        print(f"  {c}: q{a.q:g} threshold {thr:g}  weak(<thr) {weak}  strength<=1 {out[c]['le1']} ({out[c]['frac_le1']:.1%})  "
              f"==0 {out[c]['eq0']}  quartiles {out[c]['q']}  hist0-5 {hist}")
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
