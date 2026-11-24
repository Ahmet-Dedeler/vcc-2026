"""208 Phase 0 gate 1 — reproduce 196's CHALLENGE-side size tilt for both rounding rules.

196 measured `floor`'s delivered-vs-intended tilt at **-0.01234 log2** on challenge A/B/C
(decile 1 minus decile 10 of `log2(delivered/intended)` over the SCORER-gated pool, median over
blocks) and `ham`'s at **+0.00039**.  That number is the whole mechanism behind the `nmae` half
of the rider's price and the ratio it implies (challenge / h1-test = 0.913) is how this
experiment transfers the h1-test `nmae` gain, so it is re-measured here rather than quoted.

    uv run python tilt208.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
import recon196  # noqa: E402
from emit208 import tilt_stats  # noqa: E402
from ham196 import hamilton  # noqa: E402

KA = {"floor": -0.01234, "ham": +0.00039}   # 196 section 4a, challenge A/B/C


def main() -> None:
    acc = {r: [] for r in ("floor", "ham")}
    per_ctx = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        gate_scorer = D["mean_cpm"] > 5.0        # the SCORER's DE universe on this context
        sel = np.linspace(0, D["sums"].shape[0] - 1, 40).astype(int)
        rows = {r: [] for r in acc}
        for i in sel:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=400, seed=1000 + int(i))
            want = {"floor": R.totals(1.0), "ham": hamilton(R.x)}
            for r in acc:
                s = tilt_stats(want[r], R.x, gate_scorer, D["mean_cpm"])
                rows[r].append(s)
                acc[r].append(s)
        per_ctx[ctx] = {r: {k: float(np.median([x[k] for x in rows[r]]))
                            for k in rows[r][0]} for r in acc}
        print(f"ctx {ctx}: " + "  ".join(
            f"{r} tilt {per_ctx[ctx][r]['tilt_log2']:+.5f} "
            f"mean_rel {per_ctx[ctx][r]['mean_rel']:+.5f}" for r in acc), flush=True)

    print(f"\n{'rule':<8}{'mean rel':>11}{'mean |rel|':>12}{'common log2':>13}"
          f"{'dec1 log2':>11}{'dec10 log2':>12}{'TILT log2':>11}{'196 expects':>13}"
          f"{'':>6}")
    out = {}
    for r in acc:
        med = {k: float(np.median([x[k] for x in acc[r]])) for k in acc[r][0]}
        out[r] = med
        ok = abs(med["tilt_log2"] - KA[r]) <= 5e-4
        print(f"{r:<8}{med['mean_rel']:>+11.5f}{med['mean_abs_rel']:>12.5f}"
              f"{med['common_log2']:>+13.6f}{med['d1_log2']:>+11.5f}"
              f"{med['d10_log2']:>+12.5f}{med['tilt_log2']:>+11.5f}{KA[r]:>+13.5f}"
              f"{'  PASS' if ok else '  FAIL':>6}")
        if not ok:
            raise SystemExit(f"ABORT (rule 10): challenge tilt for {r} does not reproduce")
    h1 = -0.01352                                  # 196 / emit208.log, same statistic
    print(f"\n  challenge / h1-test tilt ratio = {out['floor']['tilt_log2']/h1:.3f}  "
          f"-- the factor this experiment transfers the h1-test `nmae` gain by")
    json.dump(dict(per_ctx=per_ctx, pooled=out, ka=KA,
                   ratio_to_h1=out["floor"]["tilt_log2"] / h1),
              open(WORK / "tilt208.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
