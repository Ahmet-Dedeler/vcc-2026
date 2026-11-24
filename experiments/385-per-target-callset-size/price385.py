#!/usr/bin/env python3
"""385 step 3b: emit+score per-target-K arms. Usage: uv run python price385.py <arm>
Arms: Kc1 | Kc2 | Kshuf (from work/k385.json). Constant-K baseline = work/metrics-cap15.json.
One arm per process; run one at a time. ~142 s / ~5 GiB each.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run385 as H  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402

ARMS = {"Kc1": "K_c1", "Kc2": "K_c2", "Kshuf": "K_shuf"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", choices=list(ARMS))
    args = ap.parse_args()
    H.WORK.mkdir(parents=True, exist_ok=True)
    H.verify_with_known_drift()
    spec = json.loads((H.WORK / "k385.json").read_text())
    K = np.array(spec[ARMS[args.arm]], dtype=int)
    ktargets = spec["targets"]
    assert len(K) == 85 and (K >= 400).all() and (K <= 4474).all()

    def callset_fn(i, t):
        assert ktargets[i] == t, f"target order drift: {ktargets[i]} != {t}"
        k = int(K[i])
        head = min(400, k)
        return rc_emit.CallSet(head, k - head, 1.9)

    b = H.builder385()
    eff, genes, _ = H.native_cap15(b)
    out = H.emit_arm(args.arm, eff, genes, callset_fn=callset_fn)
    out["K_rule"] = {"key": ARMS[args.arm], "K": K.tolist(),
                     "n_floor400": int((K == 400).sum()), "n_cap4474": int((K == 4474).sum())}
    H.R.write_json(H.WORK / f"metrics-{args.arm}.json", out)
    print(args.arm, "K summary: min", int(K.min()), "med", float(np.median(K)),
          "mean", float(K.mean()), "max", int(K.max()))


if __name__ == "__main__":
    main()
