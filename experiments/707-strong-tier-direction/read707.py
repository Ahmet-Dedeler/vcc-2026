#!/usr/bin/env python3
"""707 reader: the preregistered bins (PREREGISTRATION.md) for the strong-tier direction arm's board entry vs ALLSW's.

Primary: Δ(pds + mse) in overall units = (scaled Δpds + scaled Δmse) / 6 vs ALLSW (4WNmsNfdQhYx6bfib0dR); exact read.
Bins: TRANSFERS ≥ +0.0010 (adopt for October only with dumphook707 applied + one rehearsal by 10-20); WEAK +0.0005 … +0.0010
(not adopted); NULL −0.0005 … +0.0005 (not adopted); HARM ≤ −0.0005 (not adopted; record the over-read).
Falsifier, read first: DE members moving by more than 0.0015 overall in total -> STOP.
Prediction (best-of-N corrected): primary ≈ −0.0004 … +0.0001 (raw A/B/C folds +0.0003 … +0.0010).
usage: uv run python experiments/707-strong-tier-direction/read707.py <707 entry json> [<ALLSW entry json>]
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
    elif primary >= 0.0010:
        verdict = "TRANSFERS: adopt for October only with dumphook707 applied (DEF_DAY_DUMP_HOOK) + one rehearsal by 10-20"
    elif primary >= 0.0005:
        verdict = "WEAK: not adopted (below the late-switch bar); recorded"
    elif primary > -0.0005:
        verdict = "NULL: not adopted; the October-like gain was a panel artefact on Arc's panel"
    else:
        verdict = "HARM: not adopted; record what the October-like folds and h1-test over-read"
    out = {"arm_entry": arm.get("entry_id"), "base_entry": base.get("entry_id"), "arm_overall": arm["score_avg"],
           "base_overall": base["score_avg"], "delta_overall": arm["score_avg"] - base["score_avg"],
           "primary_pds_plus_mse": primary, "de_members_overall": de, "delta_scaled": d,
           "delta_raw": {"pds": arm["pds_cosine"] - base["pds_cosine"],
                         "mse": arm["expr_mse_unbiased_capped_norm"] - base["expr_mse_unbiased_capped_norm"]},
           "prediction": "primary ≈ −0.0004 … +0.0001 (best-of-N corrected)", "verdict": verdict}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
