#!/usr/bin/env python3
"""710 reader: the preregistered bins (PREREGISTRATION.md, "the coordinator's") for the H1-dose board entry vs ALLSW's.

Written by the coordinator [2d6f46] on 2026-10-16, before 710's build or read exists; the bins are copied verbatim from the
preregistration.
Primary: Δoverall (all six members) vs ALLSW's own entry (4WNmsNfdQhYx6bfib0dR). Reported beside it: Δ(pds + mse) in
overall units (exact) and the DE members in overall units (true effect plus re-roll, σ ≈ 0.00025).
Bins:
  HARM      Δoverall ≤ −0.0003           -> BLOCKS any 708 move that raises H1 above ×2 (708 AMENDMENT 2; decide708 R710=HARM)
  NULL      −0.0003 < Δoverall < +0.0003 -> no block; 708's rule decides on D/E/F's panel alone
  POSITIVE  Δoverall ≥ +0.0003           -> the instrument's H1 step is calibrated on Arc's panel; no block
Prediction: total ≈ +0.0002, in −0.0005 … +0.0008 (bulk +0.0001 … +0.0004; DE prior 0 ± 0.0005).
Falsifiers are build-side (dump KA, config diff, final rho band, zero-gene audit): check build_h1s710.log's gates first;
this reader does not see them.
usage: uv run python experiments/710-h1-dose-board/read710.py <710 entry json> [<ALLSW entry json>]
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
    primary = arm["score_avg"] - base["score_avg"]
    bulk = (d["score_pds"] + d["score_mse"]) / 6.0
    de = (d["score_nmae"] + d["score_fid"] + d["score_reach"] + d["score_jac"]) / 6.0
    if primary <= -0.0003:
        verdict = "HARM: R710=HARM, blocks any 708 move that raises H1 above x2 (decide708)"
    elif primary >= 0.0003:
        verdict = "POSITIVE: no block; the instrument's H1 step is calibrated on Arc's panel"
    else:
        verdict = "NULL: no block; 708's rule decides on D/E/F's panel alone"
    out = {"arm_entry": arm.get("entry_id"), "base_entry": base.get("entry_id"), "arm_overall": arm["score_avg"],
           "base_overall": base["score_avg"], "primary_delta_overall": primary, "bulk_pds_plus_mse_overall": bulk,
           "de_members_overall": de, "delta_scaled": d,
           "prediction": "total ≈ +0.0002 (−0.0005 … +0.0008)", "verdict": verdict}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
