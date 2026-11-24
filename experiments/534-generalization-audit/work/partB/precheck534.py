"""534 Part B -- light pre-step for a panel: bind (frames from the bundle's controls, asserted == A/B/C's for a synthetic
bundle) and the ntclevel claims (Cp46 from the bundle's NTC cells). Caches bctrl-*.npy / claims.json in the panel's
work dir, so the heavy runs do not redo them. For efp, the claims must equal 438's stored efp claims (KA-CLAIMS).

    uv run python precheck534.py <panel>
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402

import lib534 as L  # noqa: E402


def main(panel):
    ctx = L.bind(panel)
    nc = L.ntc_claims(ctx)
    print(f"{panel}: contexts {ctx.contexts}; rules {ctx.rules}; claims " + " ".join(f"{s}={nc[s]['claim']:.6f}" for s in L.SLOTS))
    if panel == "efp":
        ref = json.loads((L.ROOT / "experiments/438-october-stack-repricing/work/efp/claims.json").read_text())
        d = max(abs(nc[s]["claim"] - ref[s]["claim"]) for s in L.SLOTS)
        print(f"KA-CLAIMS vs 438's efp claims.json: |d| {d:.1e}")
        assert d == 0.0


if __name__ == "__main__":
    main(sys.argv[1])
