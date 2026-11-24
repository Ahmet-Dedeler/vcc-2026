#!/usr/bin/env python3
"""621 — where the switch's gain over production sits on h1: by delivery class (head = top-400 |key|, CPM block,
complement = uncalled), non-essential targets, identical gate pairs. Uses 521's state (435 production emitter)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
from inst636 import ESS, Inst, clean4_lolo  # noqa: E402
import lib521 as L  # noqa: E402


def main() -> None:
    I = Inst(); D = I.S["h1"]; T, G = D["T"], D["G"]
    Zb = np.load(HERE / "work/pred-h1.npz", allow_pickle=True)["Zbal"].astype(float)
    C4 = np.broadcast_to(clean4_lolo(I, "h1")[0], D["lfc"].shape).astype(float)
    S = L.load_state(0)
    ti = {t: i for i, t in enumerate(S["targets"])}; gi = {g: j for j, g in enumerate(S["genes"])}
    rt = np.array([ti[t] for t in T]); cg = np.array([gi.get(g, -1) for g in G]); ok = cg >= 0

    def remap(A, fill):
        out = np.full(D["lfc"].shape, fill, dtype=A.dtype)
        out[:, ok] = A[rt][:, cg[ok]]
        return out
    OL = remap(S["ol"].astype(float), np.nan)
    head = remap(S["head"], False); block = remap(S["block"] & ~S["head"], False); C = remap(S["C"], False)
    noness = np.array([t not in ESS for t in T])[:, None]
    m = D["gate"] & np.isfinite(C4) & np.isfinite(Zb) & np.isfinite(OL) & (OL != 0) & noness
    truth = np.sign(D["lfc"]); prod = np.sign(OL); c4 = np.where(C4 > 0.5, -1.0, 1.0); cons = np.sign(Zb)
    out = {}
    for c in (1, 2):
        sw = np.where(np.abs(Zb) >= c, cons, c4)
        for name, cm in (("head", head), ("block", block), ("complement", ~C)):
            mm = m & cm
            n = int(mm.sum())
            pr = float((prod[mm] == truth[mm]).mean()); s_ = float((sw[mm] == truth[mm]).mean())
            # switch applied ONLY where it overrides production with a confident consensus
            ov = mm & (np.abs(Zb) >= c) & (cons != prod)
            out[f"c{c}_{name}"] = dict(n=n, prod=pr, switch=s_, delta=s_ - pr, n_override=int(ov.sum()),
                                       override_acc_cons=float((cons[ov] == truth[ov]).mean()) if ov.any() else None)
            print(f"c{c} {name:10s} n={n:6d} prod {pr:.4f} switch {s_:.4f} Δ {s_ - pr:+.4f} | where confident consensus "
                  f"overrides prod: n={int(ov.sum())}, consensus right {out[f'c{c}_{name}']['override_acc_cons']}")
    (HERE / "work/byclass621.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
