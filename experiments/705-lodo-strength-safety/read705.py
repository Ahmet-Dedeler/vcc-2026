#!/usr/bin/env python3
"""705 reader: the preregistered SAFETY bins (PREREGISTRATION.md) for the LODO strength rider's board entry vs ALLSW's.

Primary: Δ(pds + mse) in overall units = (scaled Δpds + scaled Δmse) / 6, against ALLSW's own entry (4WNmsNfdQhYx6bfib0dR).
The direction members reproduce exactly on the board, so this read is exact.
Bins: HARM (primary ≤ −0.0005) VETOES the rider for October; otherwise it is CLEARED (ships once [17adab]'s pass-through
is applied with its known answers and one rehearsal). Falsifier, read first: DE members moving by more than 0.0015 overall
in total -> STOP (a bulk-only arm moves them only through the re-roll, σ ≈ 0.00025).
Prediction: primary ≈ +0.0001 (−0.0003 … +0.0004).
usage: uv run python experiments/705-lodo-strength-safety/read705.py <705 entry json> [<ALLSW entry json>]
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
    if abs(de) > 0.0015:
        verdict = "STOP: DE members moved by more than the re-roll allows for a bulk-only arm (diagnose before reading)"
    elif primary <= -0.0005:
        verdict = "HARM: VETO, the rider stays off for October; record what the October-like instruments missed"
    else:
        verdict = ("CLEARED: lodo_strength theta 0.25 may ship for October once [17adab]'s pass-through is applied with its "
                   "known answers and one rehearsal (else it ships off)")
    out = {"arm_entry": arm.get("entry_id"), "base_entry": base.get("entry_id"), "arm_overall": arm["score_avg"],
           "base_overall": base["score_avg"], "delta_overall": arm["score_avg"] - base["score_avg"],
           "primary_pds_plus_mse": primary, "de_members_overall": de, "delta_scaled": d,
           "delta_raw": {"pds": arm["pds_cosine"] - base["pds_cosine"],
                         "mse": arm["expr_mse_unbiased_capped_norm"] - base["expr_mse_unbiased_capped_norm"]},
           "prediction": "primary ≈ +0.0001 (−0.0003 … +0.0004); DE re-roll only", "verdict": verdict}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
