"""216 — the 1/k arithmetic, so the mandate's "a fourth line cannot explain x2.06" is a
number on disk rather than a sentence.

159's law as 207 §7.2 calibrated it, with the board->surface anchor at 197's a15 point:

    1 / cos_surface^2 = 66.5 + 1920 / k ,     cos_surface = cos_board / 2.078

Known-answer rows: it must return 207's published `k_eff` of 2.46 (CD4 off), 3.09 (shipped),
4.49 (f = 0.34), 3.84 (the best equal consensus of all four panel-covering public lines) and
10.41 (`jw_13288`'s floor).

    uv run python keff216.py
"""
from __future__ import annotations

import json
from pathlib import Path

WORK = Path(__file__).resolve().parent / "work"
ANCHOR, A_LAW, B_LAW = 2.078, 66.5, 1920.0
EXPECT = {"CD4 off (207 §7.2)": 2.46, "shipped f=0.10": 3.09, "f=0.34 (incumbent)": 4.49,
          "best 4-line equal consensus (197 §4.2)": 3.84, "jw_13288 floor (215)": 10.41}


def k_of(board_cos: float) -> float:
    cs = board_cos / ANCHOR
    return B_LAW / (1.0 / cs ** 2 - A_LAW)


def board_of(k: float) -> float:
    return ANCHOR * (A_LAW + B_LAW / k) ** -0.5


def main() -> None:
    price = json.loads((WORK / "price216.json").read_text())
    grid = {f"f={r['f']:g}": r["board_cos_k1"] for r in price["rows"]}
    rows = [
        ("CD4 off (207 §7.2)", 0.07138),
        ("shipped f=0.10", grid["f=0.1"]),
        ("f=0.20", grid["f=0.2"]),
        ("f=0.34 (incumbent)", grid["f=0.34"]),
        ("f=0.40 (unconstrained argmax)", grid["f=0.4"]),
        ("f=0.55", grid["f=0.55"]),
        ("f=0.80", grid["f=0.8"]),
        ("f=1.00 (cosine argmax)", grid["f=1"]),
        ("f=1.00, k_cos-calibrated",
         [r for r in price["rows"] if r["f"] == 1.0][0]["board_cos_kcal"]),
        ("best 4-line equal consensus (197 §4.2)", 0.08730),
        ("jw_13288 floor (215)", 0.13117),
        ("cqawesome floor (215)", 0.13211),
    ]
    out, gates = [], []
    print(f"{'configuration':<40}{'board cos':>11}{'k_eff':>9}")
    for name, b in rows:
        k = k_of(b)
        out.append(dict(name=name, board_cos=b, k_eff=k))
        print(f"{name:<40}{b:>11.6f}{k:>9.2f}")
        if name in EXPECT:
            ok = abs(round(k, 2) - EXPECT[name]) <= 0.01
            gates.append(dict(name=name, expected=EXPECT[name], measured=round(k, 2),
                              passed=bool(ok)))

    k0 = k_of(grid["f=0.1"])
    summary = dict(
        k_shipped=k0,
        k_plus_one_line=k0 + 1.0,
        board_cos_plus_one_line=board_of(k0 + 1.0),
        rel_plus_one_line=board_of(k0 + 1.0) / grid["f=0.1"],
        rel_sqrt_4_over_3=(4 / 3) ** 0.5,
        rel_dose_f034=grid["f=0.34"] / grid["f=0.1"],
        rel_dose_cos_argmax=grid["f=1"] / grid["f=0.1"],
        leader_move_rel=0.13211 / 0.06427,
        k_needed_for_leader=k_of(0.13211),
        lines_needed=k_of(0.13211) - k0)
    print()
    print(f"our k_eff = {k0:.3f}.  ONE more equal line -> k = {k0+1:.3f} -> board cosine "
          f"{summary['board_cos_plus_one_line']:.6f} = x{summary['rel_plus_one_line']:.4f}")
    print(f"159's sqrt(4/3) rule of thumb for 3 -> 4 lines:            "
          f"x{summary['rel_sqrt_4_over_3']:.4f}")
    print(f"the DOSE alone at f = 0.34:                                "
          f"x{summary['rel_dose_f034']:.4f}")
    print(f"the DOSE at its cosine argmax (f = 1.0, -17.2 sigma pds):  "
          f"x{summary['rel_dose_cos_argmax']:.4f}")
    print(f"the leaders' measured move (215 §3.2):                     "
          f"x{summary['leader_move_rel']:.4f}")
    print(f"k_eff needed to stand where they stand: {summary['k_needed_for_leader']:.2f} "
          f"against our {k0:.2f} -- {summary['lines_needed']:.1f} more effective lines")
    print(f"\nknown-answer gates: {sum(g['passed'] for g in gates)}/{len(gates)} pass")
    (WORK / "keff216.json").write_text(json.dumps(
        dict(rows=out, gates=gates, **summary), indent=2))


if __name__ == "__main__":
    main()
