#!/usr/bin/env python3
"""621 — confidence switch: per-target donor consensus sign where |Zbal| >= c, CLEAN4 (LOLO) elsewhere.

Scored on each surface's gate pairs that CLEAN4 covers (so every rule is compared on identical pairs), by |Zbal| bin
and as a switch rule over c. c is chosen on one blob-free surface (h1 or k562) and read on the other.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import Inst, clean4_lolo  # noqa: E402

CS = [0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 1e9]
BINS = [0, 1, 2, 3, 4, 6, 10, 20, 1e9]


def main() -> None:
    I = Inst()
    out = {}
    for s in I.S:
        D = I.S[s]
        z = np.load(HERE / f"work/pred-{s}.npz", allow_pickle=True)
        Zb = z["Zbal"].astype(float)
        C4 = np.broadcast_to(clean4_lolo(I, s)[0], D["lfc"].shape)
        m = D["gate"] & np.isfinite(C4) & np.isfinite(Zb)
        truth_dn = D["lfc"][m] < 0
        c4_dn = C4[m] > 0.5
        cons_dn = Zb[m] < 0
        a = np.abs(Zb[m])
        rows = {"n": int(m.sum()), "clean4_acc": float((c4_dn == truth_dn).mean()),
                "cons_acc": float((cons_dn == truth_dn).mean()), "bins": [], "switch": []}
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            k = (a >= lo) & (a < hi)
            if k.sum() == 0:
                continue
            rows["bins"].append({"absZ": [lo, hi], "n": int(k.sum()), "share": float(k.mean()),
                                 "clean4": float((c4_dn[k] == truth_dn[k]).mean()),
                                 "cons": float((cons_dn[k] == truth_dn[k]).mean())})
        for c in CS:
            pick = np.where(a >= c, cons_dn, c4_dn)
            rows["switch"].append({"c": c, "acc": float((pick == truth_dn).mean()), "share_cons": float((a >= c).mean())})
        out[s] = rows
        print(f"\n== {s}: gate pairs (CLEAN4-covered) {rows['n']}; CLEAN4 {rows['clean4_acc']:.4f}  consensus {rows['cons_acc']:.4f}")
        for b in rows["bins"]:
            print(f"   |Z| [{b['absZ'][0]:>4}, {b['absZ'][1]:>4}) n {b['n']:6d} ({b['share']:.3f})  CLEAN4 {b['clean4']:.4f}  consensus {b['cons']:.4f}  Δ {b['cons'] - b['clean4']:+.4f}")
        print("   switch c: " + "  ".join(f"{r['c']:g}:{r['acc']:.4f}({r['share_cons']:.2f})" for r in rows["switch"]))
    (HERE / "work/switch621.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
