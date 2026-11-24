#!/usr/bin/env python3
"""537 — the called-gene magnitude scale s (§6.1) on the block + CLEAN4 hook construction (CPH: run537 `lcap-kept+c4h`, 3 seeds),
analytic as in scale537.py (calls, signs, order held; only nmae moves): does s stack with CPH as it does on 446bj2z's base?
Parts: head (the 400 design genes among the calls) and all calls; s 1.25 / 1.5; 521's regimes. Writes work/scale537_cph.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(HERE))
import lib521 as L  # noqa: E402
import price537 as P7  # noqa: E402
from table521 import regime_pool  # noqa: E402

INV = 1.64787


def main():
    res = {}
    for seed in (0, 1, 2):
        S = L.load_state(seed)
        A = P7.load("lcap-kept+c4h", seed, True)
        call = (A["q"] < L.ALPHA) & ~S["own"]
        head = S["head"] & call
        for reg in ("1.0", "0.5", "0.35", "cmp"):
            pool = regime_pool(S, reg)[0]
            base = np.nanmean(L.nmae_rows(A["lfc"], pool, S["rl"]))
            for part, M in (("head", head), ("all", call)):
                for s in (1.25, 1.5):
                    v = np.nanmean(L.nmae_rows(np.where(M, A["lfc"] * s, A["lfc"]), pool, S["rl"]))
                    res.setdefault(f"{reg}|{part}|{s}", []).append(-(v - base) * INV / 6)
    out = {k: float(np.mean(v)) for k, v in res.items()}
    for reg in ("1.0", "0.5", "0.35", "cmp"):
        print(f"f {reg:4s} " + "  ".join(f"{p} x{s}: {out[f'{reg}|{p}|{s}']:+.5f}" for p in ("head", "all") for s in (1.25, 1.5)))
    (HERE / "work/scale537_cph.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
