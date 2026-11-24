"""218 Phase 2 pricing — every arm on all six members, under AGENTS.md section 8.4's rule.

Section 8.4 is explicit that pricing a DIRECTION arm on `mse` alone understates it ~2.5x and
that `nmae`/`fid`/`reach`/`jac` must never be carried at zero.  The cleanest calibration this
project has is the `cd4f034` board row itself (`210/BOARD-RESULT.md`): one arm, one variable,
six measured member deltas and a measured cosine move.  Dividing that row by its own cosine
gain gives a per-unit-cosine rate card:

    cos 0.072280 -> 0.090270   (197's one-sided inversion of raw `mse`), d = +0.017990

    member   d scaled overall     per +0.01 cosine
    mse       +0.006537            +0.003634
    reach     +0.004400            +0.002446
    nmae      +0.001762            +0.000980
    fid       +0.001524            +0.000847
    jac       +0.000179            +0.000100
    pds       -0.001345            -0.000748
    TOTAL     +0.013057            +0.007259

`mse`'s rate reproduces 197's independently derived +0.00365 overall per +0.01 cosine to
0.5%, which is the check that the card is not circular.

The PLOLO -> board transfer is section 8.4's measured **1.80x on the RELATIVE cosine gain**
(212 read +13.83% at `f = 0.34`; the board delivered +24.9%).  Both the uncalibrated (1.0x)
and calibrated (1.80x) columns are printed and the uncalibrated one is the headline, as in
207.

CAVEAT CARRIED FORWARD, not buried: the rate card rests on ONE board point.  Section 8.4's
own lesson is that a class constant estimated from one point should be treated as a LOWER
BOUND on magnitude, so the `pds` cost may be understated as easily as the gains.

    uv run python price218.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

BOARD_COS = 0.090270           # cd4f034's inverted board cosine (the current incumbent)
K_COS = 1.80                   # section 8.4's PLOLO -> board relative-cosine factor
RATE = dict(mse=+0.003634, reach=+0.002446, nmae=+0.000980, fid=+0.000847,
            jac=+0.000100, pds=-0.000748)          # overall per +0.01 ABSOLUTE cosine
SPAN_PDS = 0.452               # 162's board reference span for `pds`
SIGMA_PDS = 0.001494           # 162's leave-one-out sigma on raw PLOLO `pds`
GATE_SHAPE = 0.0019            # rule 13, shape-only class
GATE_AMPL = 0.0089             # rule 13, amplitude-changing class
SEED = 20260904


def price(rel_cos_gain: float, d_pds_plolo: float, k: float = 1.0) -> dict:
    """Six members from a relative PLOLO cosine gain and a PLOLO `pds` delta."""
    d_cos_abs = BOARD_COS * (K_COS * k) * rel_cos_gain if k else BOARD_COS * rel_cos_gain
    members = {m: RATE[m] * d_cos_abs / 0.01 for m in RATE}
    # `pds` twice: from the rate card (cosine-driven) and from PLOLO directly (dose-driven,
    # 207's route: raw delta / span / 6 metrics x 162's 1.011, then section 8.4's 2.0x).
    pds_plolo = 2.0 * (d_pds_plolo / SPAN_PDS / 6.0 * 1.011)
    members["pds_from_plolo"] = pds_plolo
    total_card = sum(members[m] for m in RATE)
    total_mixed = sum(members[m] for m in RATE if m != "pds") + pds_plolo
    return dict(d_cos_abs=d_cos_abs, **members, total_rate_card=total_card,
                total_pds_from_plolo=total_mixed,
                pds_sigma=abs(d_pds_plolo) / SIGMA_PDS)


def main() -> None:
    rows = json.loads((WORK / "compose218.json").read_text())
    folds = sorted({r["fold"] for r in rows})
    by = {(r["fold"], r["arm"]): r for r in rows}

    for cf in (0.34, 0.10):
        base = {f: by.get((f, f"cd4{cf:g}+none0")) for f in folds}
        if any(v is None for v in base.values()):
            continue
        arms = sorted({r["arm"] for r in rows if r["cd4_f"] == cf and r["second"] != "none"})
        print("=" * 116)
        print(f"SECOND SOURCE ON TOP OF cd4 f = {cf:g}   "
              f"(baseline 4-fold mean cos_agg "
              f"{np.mean([base[f]['cos_agg'] for f in folds]):.6f}, "
              f"pds {np.mean([base[f]['pds'] for f in folds]):.6f})")
        print("=" * 116)
        print(f"{'arm':34s} {'nf':>3s} {'rel dcos%':>10s} {'d pds':>10s} {'d q_w':>10s} "
              f"{'sigma':>6s} {'price k=1':>10s} {'price k=1.8':>11s} {'gate':>6s}")
        out = []
        for arm in arms:
            sel = [(f, by[(f, arm)]) for f in folds if (f, arm) in by]
            if len(sel) < 2:
                continue
            rel = float(np.mean([r["cos_agg"] / base[f]["cos_agg"] - 1.0 for f, r in sel]))
            dp = float(np.mean([r["pds"] - base[f]["pds"] for f, r in sel]))
            dq = float(np.mean([r["a15_sgn_w"] - base[f]["a15_sgn_w"] for f, r in sel]))
            p1 = price(rel, dp, k=0)          # k=0 -> no 1.80x, PLOLO carried straight across
            p2 = price(rel, dp, k=1)
            verdict = "PASS" if p1["total_pds_from_plolo"] > GATE_SHAPE else "fail"
            out.append(dict(arm=arm, n_folds=len(sel), rel_cos=rel, d_pds=dp, d_qw=dq,
                            price_k1=p1, price_k18=p2, gate=verdict,
                            per_fold={f: dict(cos_agg=r["cos_agg"], pds=r["pds"],
                                              rel=r["cos_agg"] / base[f]["cos_agg"] - 1.0,
                                              d_pds=r["pds"] - base[f]["pds"])
                                      for f, r in sel}))
            print(f"{arm:34s} {len(sel):3d} {100*rel:+10.3f} {dp:+10.6f} {dq:+10.6f} "
                  f"{p1['pds_sigma']:6.2f} {p1['total_pds_from_plolo']:+10.5f} "
                  f"{p2['total_pds_from_plolo']:+11.5f} {verdict:>6s}")
        (WORK / f"price218_cd4{cf:g}.json").write_text(json.dumps(out, indent=2))
        print()

    # ---------------------------------------------------------- known-answer: cd4f034 itself
    print("=" * 116)
    print("KNOWN-ANSWER on the PRICE: re-price the cd4f034 arm this machinery was "
          "calibrated on")
    print("=" * 116)
    r = json.loads((WORK / "resid218.json").read_text())
    b = {x["fold"]: x for x in r if x["arm"] == "f0"}
    a = {x["fold"]: x for x in r if x["arm"] == "cd4_gn0.34"}
    rel = float(np.mean([a[f]["cos_agg"] / b[f]["cos_agg"] - 1.0 for f in a]))
    dp = float(np.mean([a[f]["pds"] - b[f]["pds"] for f in a]))
    p = price(rel, dp, k=1)
    print(f"  PLOLO rel cosine gain vs CD4-off : {100*rel:+.3f}%   "
          f"(207 measured +18.0% at f=0.34 against the SHIPPED f=0.10 baseline)")
    print(f"  priced total (k=1.80)            : {p['total_pds_from_plolo']:+.5f}")
    print(f"  cd4f034 realised on the board    : +0.01306   "
          f"[against f=0.10, not f=0; this row is a scale check, not an identity]")


if __name__ == "__main__":
    main()
