#!/usr/bin/env python3
"""702 reader: the preregistered bins (PREREGISTRATION.md, amendments 1-2) for the PIE bulk arm's board entry vs ALLSW's.

Primary: Δ(pds + mse) in overall units = (scaled Δpds + scaled Δmse) / 6, against ALLSW's own entry (4WNmsNfdQhYx6bfib0dR).
Bins: TRANSFERS ≥ +0.0015; WEAK +0.0005 … +0.0015 (adopt only if both held-out reads gained at w* — they did, amendment 1 —
and the DE members did not fall by more than 0.0005 overall); NULL −0.0005 … +0.0005; HARM ≤ −0.0005.
Falsifiers: DE members moving by more than 0.0015 overall in total (a bulk-only arm moves them only through the re-roll),
reported as STOP. Prediction (amendment 1): primary +0.0008 … +0.0017; Δmse +0.0013 … +0.0031, Δpds −0.0005 … −0.0014.
usage: uv run python experiments/702-pie-bulk-arm/read702.py <702 entry json> [<ALLSW entry json>]
       (entry jsons as written by scripts/capture_board_entry.py into data/board-entries/)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data/board-entries/4WNmsNfdQhYx6bfib0dR.json"
SCALED = ("score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac")


def entry(p):
    d = json.loads(Path(p).read_text())
    return d.get("entry", d)


def main():
    arm = entry(sys.argv[1])
    base = entry(sys.argv[2] if len(sys.argv) > 2 else BASE)
    d = {k: arm[k] - base[k] for k in SCALED}
    primary = (d["score_pds"] + d["score_mse"]) / 6.0
    de = (d["score_nmae"] + d["score_fid"] + d["score_reach"] + d["score_jac"]) / 6.0
    total = arm["score_avg"] - base["score_avg"]
    if abs(de) > 0.0015:
        verdict = "STOP: DE members moved by more than the re-roll allows for a bulk-only arm (diagnose before reading)"
    elif primary >= 0.0015:
        verdict = "TRANSFERS: adopt for October (PIE dump hook on D/E/F at w 0.25, seen rule; def_day pass-through + rehearsal)"
    elif primary >= 0.0005:
        verdict = ("WEAK: adopt (both held-out reads gained at w*)" if de >= -0.0005
                   else "WEAK but the DE members fell by more than 0.0005 overall: do not adopt")
    elif primary > -0.0005:
        verdict = "NULL: not adopted; PIE stays a candidate only through 702c / per-context follow-ups"
    else:
        verdict = "HARM: not adopted; record the held-out instruments' over-read"
    out = {"arm_entry": arm.get("entry_id"), "base_entry": base.get("entry_id"), "arm_overall": arm["score_avg"],
           "base_overall": base["score_avg"], "delta_overall": total, "primary_pds_plus_mse": primary,
           "de_members_overall": de, "delta_scaled": d,
           "delta_raw": {"pds": arm["pds_cosine"] - base["pds_cosine"],
                         "mse": arm["expr_mse_unbiased_capped_norm"] - base["expr_mse_unbiased_capped_norm"]},
           "prediction": "primary +0.0008 … +0.0017 (Δmse +0.0013 … +0.0031, Δpds −0.0005 … −0.0014)", "verdict": verdict}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
