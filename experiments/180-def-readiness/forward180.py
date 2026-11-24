"""180 — what the knee rule is worth on a D/E/F context that lands off the knee.

The A/B/C gain is small because A/B/C are already close to the knee (rho 0.999 / 0.940 /
0.878 under the global 0.982, which was FITTED on A/B/C board points). A new panel has no
such guarantee. This prices the rule as a function of where a context lands, using 096's
identity and A/B/C's own measured terms:

    raw mse_c - 1 = ( s^2 A_c - 2 s X + Kc_c - min(s^2 B_c, claim_c) ) / sumD

with `claim_c` set so that `rho = B_c / claim_c` takes the value on each row, and `s*` the
knee `sqrt(claim_c / B_c) = 1 / sqrt(rho)`.

WHAT IS CARRIED AND WHAT IS NOT. `A_c`, `B_c`, `Kc_c` are A/B/C's measured values; `X` and
`sumD` are the A/B/C board-anchored solve. A D/E/F context will have its OWN `X` and `sumD`
and neither is knowable before the round. So this is the SHAPE of the gain, sized on our
current alignment -- not a forecast of the final round's number. It is quoted that way.

    uv run python forward180.py
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SUMD = 0.229141
SPAN = 0.950162566
RHOS = (1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50, 1.10, 1.25)


def main() -> None:
    knee = json.loads((WORK / "knee180.json").read_text())
    X = knee["mse"]["X"]
    rows = []
    print(f"[fwd] A/B/C-anchored X {X:.6f}, sumD {SUMD}, mse span {SPAN}")
    print(f"[fwd] a context with context C's own A / B / Kc, placed at each rho:\n")
    A = knee["per_context"]["C"]["before"]["A"]
    B = knee["per_context"]["C"]["before"]["B"]
    Kc = knee["Kc"]["C"]
    print(f"      A {A:.6f}  B {B:.6f}  Kc {Kc:.6f}")
    print(f"\n  {'rho':>6}{'claim':>10}{'knee s*':>10}{'alpha*':>10}"
          f"{'d raw mse':>12}{'d score_mse':>13}{'d overall':>12}")
    for rho in RHOS:
        claim = B / rho
        s = (claim / B) ** 0.5
        num0 = A + Kc - min(B, claim) - 2 * X
        num1 = s * s * A + Kc - min(s * s * B, claim) - 2 * s * X
        d_raw = (num1 - num0) / SUMD
        d_score = -d_raw / SPAN
        rows.append(dict(rho=rho, claim=claim, s=s, alpha=0.982 * s,
                         d_raw_mse=d_raw, d_score_mse=d_score,
                         d_overall_one_context=d_score / 3 / 6))
        print(f"  {rho:>6.2f}{claim:>10.5f}{s:>10.5f}{0.982 * s:>10.5f}"
              f"{d_raw:>+12.5f}{d_score:>+13.5f}{d_score / 3 / 6:>+12.6f}")
    print("\n  `d overall` is ONE context of three on ONE member of six: d_score / 3 / 6.")
    print("  The curve is ASYMMETRIC and steeper on the high side. Below the knee the "
          "un-refunded part of the amplitude is\n  only `A - B` (0.0012 here, 0.5% of A), "
          "so being 30% short costs 0.030 of scaled mse; ABOVE the knee the whole "
          "`2sA`\n  is charged with no refund at all, so being 10% long costs 0.090. A "
          "D/E/F context that lands high is the expensive\n  case, and a global alpha "
          "fitted on A/B/C has no way to know which side it will land on.")
    (WORK / "forward180.json").write_text(json.dumps(
        {"X": X, "sumD": SUMD, "A": A, "B": B, "Kc": Kc, "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
