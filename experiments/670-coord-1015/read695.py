#!/usr/bin/env python3
"""Readers for R-ABC2 (695) and R-ABC3 (696): the composition bins from their preregistrations (readiness [4592ee]).

Written by the coordinator [2d6f46] on 2026-10-16, before either build or read exists; the formula and bins are copied verbatim.
695 (R-ABC2): pred = R-ABC 0.2743318 (PooEbFClgTLMQoXiOCQP) + Σ over adopted switches of that switch's own board Δoverall vs ALLSW
  (4WNmsNfdQhYx6bfib0dR). Pass --s705 <entry json> if 705 was adopted (θ 0.25) and --s707 <entry json> if 707 was adopted (hook).
  D = R-ABC2 − pred: ABOVE ≥ +0.0005 (ships as decided; record the interaction), COMPOSES (−0.0005, +0.0005) (ships as decided),
  BELOW ≤ −0.0005 (revert the switches R-ABC2 added, latest adopted first: 707's hook, then 705's θ, until R-ABC3 / 10-21).
696 (R-ABC3, --rabc3): pred = R-ABC2's own board read + ku30's Δoverall vs 446bj2z (+0.0011785); the same ±0.0005 bins; and the bulk
  members (pds, mse raw) must equal R-ABC2's (a K-only change must leave the bulk identical): otherwise STOP.
697 (R-ABC2H, --rabc2h): pred = R-ABC2's own board read + 707's Δoverall vs ALLSW (sw707 AYbAiHRj, +0.0010330); ±0.0005 bins; BELOW ⇒ the
  6a route ships WITHOUT the swap (697 PREREGISTRATION).
698 (C1, --c1): pred = R-ABC2H's own board read + 0.00153 / 1.7232 (lab PREREG-C1-ABC d1f55f78); |D| ≤ 0.0005 AS PREDICTED, D > +0.0005 ABOVE,
  D ≤ −0.0005 HARM; a second HARM condition is C1 < R-ABC2H − 0.0005. Either HARM closes 811's 6a route.
699 (C0, --c0; added 2026-10-20 02:52Z, before C0 is built): pred = R-ABC2's own board read + 0.00269 / 1.1280 (lab PREREG-C0-ABC
  e603fe73); the same bins as C1; the second HARM condition is C0 < R-ABC2 − 0.0005 (SIXA-DECISION Amendment 1). HARM closes X2 on the
  swap-OFF branch only. Also printed (information): C0 − R-ABC2 beside C1 − R-ABC2H (811 with and without the swap), with --c1-entry.
673 (W1, --w1; added 2026-10-20 03:33Z, before W1 is built): pred = R-ABC2's own board read + 710's Δoverall vs ALLSW (h1s710
  entry, +0.0004824); the ±0.0005 bins; BELOW blocks any 708 H1 move above ×2 in the final (673 PREREGISTRATION 5c4f4bc5). Guard only.
674 (W2, --w2; added 2026-10-20 04:32Z, before W2 is built): pred = W1's own board read + 457dr's Δoverall vs 453dr (+0.0031484); the ±0.0005
  bins; a second BELOW condition is W2 < W1 − 0.0005. BELOW → X1 may not carry CD4 > ×1 (674 PREREGISTRATION 12562add, SIXA Amendment 3).
usage: uv run python experiments/670-coord-1015/read695.py <rabc2 entry json> [--s705 <705 entry json>] [--s707 <707 entry json>]
       uv run python experiments/670-coord-1015/read695.py --rabc3 <rabc3 entry json> <rabc2 entry json>
       uv run python experiments/670-coord-1015/read695.py --rabc2h <rabc2h entry json> <rabc2 entry json>
       uv run python experiments/670-coord-1015/read695.py --c1 <c1 entry json> <rabc2h entry json>
       uv run python experiments/670-coord-1015/read695.py --w1 <w1 entry json> <rabc2 entry json>
       uv run python experiments/670-coord-1015/read695.py --w2 <w2 entry json> <w1 entry json>
       uv run python experiments/670-coord-1015/read695.py --c0 <c0 entry json> <rabc2 entry json> [--c1-entry <c1 json> --rabc2h-entry <json>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BOARD = ROOT / "data/board-entries"
RABC, ALLSW, BASE446, KU30 = "PooEbFClgTLMQoXiOCQP", "4WNmsNfdQhYx6bfib0dR", "be90smZUDcIV7j5MgQ03", "HAbelEdgYYziqwfjTgdI"
SW707 = "AYbAiHRjuHdil2zCjlsk"
H1S710 = "UW4LsxTUkeLumDCBotfU"
CD4_457DR, BASE453DR = "1FfirZsGsAZvaXYTd39J", "q3VuTImh4FRGKXH2waCp"
C1_INC = 0.00153 / 1.7232   # lab PREREG-C1-ABC (d1f55f78), fixed before R-ABC2H's read
C0_INC = 0.00269 / 1.1280   # lab PREREG-C0-ABC (e603fe73), fixed before C0 was built (class-median rho of 707 and C1)


def entry(p):
    d = json.loads(Path(p).read_text())
    return d.get("entry", d)


def ov(e):
    return e["score_avg"]


def bin_of(d):
    if d >= 0.0005:
        return "ABOVE: ships as decided; record the interaction"
    if d <= -0.0005:
        return "BELOW: revert the switches added, latest adopted first, one at a time (coordinator's rule)"
    return "COMPOSES: ships as decided"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm")
    ap.add_argument("ref", nargs="?")
    ap.add_argument("--s705")
    ap.add_argument("--s707")
    ap.add_argument("--rabc3", action="store_true")
    ap.add_argument("--rabc2h", action="store_true")
    ap.add_argument("--c1", action="store_true")
    ap.add_argument("--c0", action="store_true")
    ap.add_argument("--w1", action="store_true")
    ap.add_argument("--w2", action="store_true")
    ap.add_argument("--c1-entry")
    ap.add_argument("--rabc2h-entry")
    a = ap.parse_args()
    if a.w2:
        arm, w1 = entry(a.arm), entry(a.ref)
        k = ov(entry(BOARD / f"{CD4_457DR}.json")) - ov(entry(BOARD / f"{BASE453DR}.json"))
        pred = ov(w1) + k
        d = ov(arm) - pred
        below_w1 = ov(arm) < ov(w1) - 0.0005
        v = ("BELOW: X1 may not carry CD4 > x1 (SIXA Amendment 3)" + (" (W2 below W1 by >= 0.0005)" if below_w1 else "")
             if d <= -0.0005 or below_w1 else
             "ABOVE: no block (guard only; A/B/C-positive by 438)" if d >= 0.0005 else
             "COMPOSES: no block (guard only; A/B/C-positive by 438)")
        out = {"read": "674 W2", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "w1_overall": ov(w1),
               "cd4_457dr_delta": k, "pred": pred, "D": d, "w2_minus_w1": ov(arm) - ov(w1), "verdict": v}
        print(json.dumps(out, indent=1)); return
    if a.w1:
        arm, r2 = entry(a.arm), entry(a.ref)
        k = ov(entry(BOARD / f"{H1S710}.json")) - ov(entry(BOARD / f"{ALLSW}.json"))
        pred = ov(r2) + k
        d = ov(arm) - pred
        v = ("BELOW: no 708 move may raise H1 above x2 in the final" if d <= -0.0005 else
             "ABOVE: no block (guard only; the weights stay as 708 decides)" if d >= 0.0005 else
             "COMPOSES: no block (guard only; the weights stay as 708 decides)")
        out = {"read": "673 W1", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "rabc2_overall": ov(r2),
               "h1s710_delta": k, "pred": pred, "D": d, "verdict": v}
        print(json.dumps(out, indent=1)); return
    if a.c0:
        arm, r2 = entry(a.arm), entry(a.ref)
        pred = ov(r2) + C0_INC
        d = ov(arm) - pred
        below_r2 = ov(arm) < ov(r2) - 0.0005
        if d <= -0.0005 or below_r2:
            v = "HARM: X2 on the swap-OFF branch is closed" + (" (C0 below R-ABC2 by >= 0.0005)" if below_r2 else "")
        elif d > 0.0005:
            v = "ABOVE: the swap-OFF X2 has its A/B/C read (no Arc-drawn HARM); record the interaction"
        else:
            v = "AS PREDICTED: the swap-OFF X2 has its A/B/C read (no Arc-drawn HARM)"
        out = {"read": "699 C0", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "rabc2_overall": ov(r2),
               "increment": C0_INC, "pred": pred, "D": d, "c0_minus_rabc2": ov(arm) - ov(r2), "verdict": v}
        if a.c1_entry and a.rabc2h_entry:
            c1, h = entry(a.c1_entry), entry(a.rabc2h_entry)
            out["info_c1_minus_rabc2h"] = ov(c1) - ov(h)
            out["info_c1_minus_c0"] = ov(c1) - ov(arm)
        print(json.dumps(out, indent=1)); return
    if a.rabc2h:
        arm, rabc2 = entry(a.arm), entry(a.ref)
        k = ov(entry(BOARD / f"{SW707}.json")) - ov(entry(BOARD / f"{ALLSW}.json"))
        pred = ov(rabc2) + k
        d = ov(arm) - pred
        v = bin_of(d)
        if d <= -0.0005:
            v = "BELOW: the 6a route ships WITHOUT the swap (the composition does not hold on A/B/C)"
        out = {"read": "697 R-ABC2H", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "rabc2_overall": ov(rabc2),
               "sw707_delta": k, "pred": pred, "D": d, "verdict": v}
        print(json.dumps(out, indent=1)); return
    if a.c1:
        arm, h = entry(a.arm), entry(a.ref)
        pred = ov(h) + C1_INC
        d = ov(arm) - pred
        below_h = ov(arm) < ov(h) - 0.0005
        if d <= -0.0005 or below_h:
            v = "HARM: 811's 6a route is closed" + (" (C1 below R-ABC2H by >= 0.0005)" if below_h else "")
        elif d > 0.0005:
            v = "ABOVE: 811's 6a route stays open (no Arc-drawn HARM); record the interaction"
        else:
            v = "AS PREDICTED: 811's 6a route stays open (no Arc-drawn HARM)"
        out = {"read": "698 C1", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "rabc2h_overall": ov(h),
               "increment": C1_INC, "pred": pred, "D": d, "c1_minus_rabc2h": ov(arm) - ov(h), "verdict": v}
        print(json.dumps(out, indent=1)); return
    if a.rabc3:
        arm, rabc2 = entry(a.arm), entry(a.ref)
        k = ov(entry(BOARD / f"{KU30}.json")) - ov(entry(BOARD / f"{BASE446}.json"))
        pred = ov(rabc2) + k
        d = ov(arm) - pred
        leak = {m: arm[f] - rabc2[f] for m, f in (("pds", "pds_cosine"), ("mse", "expr_mse_unbiased_capped_norm"))}
        verdict = (f"STOP: bulk moved vs R-ABC2 ({leak}); a K-only change must leave the bulk (pds, mse) identical"
                   if max(abs(v) for v in leak.values()) > 1e-6 else bin_of(d)
                   + ("; BELOW reverts K (rc_callset_cpm off) for R0" if d <= -0.0005 else ""))
        out = {"read": "696 R-ABC3", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "rabc2_overall": ov(rabc2),
               "ku30_delta": k, "pred": pred, "D": d, "bulk_delta_raw_vs_rabc2": leak, "verdict": verdict}
    else:
        arm = entry(a.arm)
        base_rabc, allsw = entry(BOARD / f"{RABC}.json"), entry(BOARD / f"{ALLSW}.json")
        adopted = {}
        for name, path in (("705 theta 0.25", a.s705), ("707 hook", a.s707)):
            if path:
                adopted[name] = ov(entry(path)) - ov(allsw)
        if not adopted:
            raise SystemExit("no adopted switch given: R-ABC2 would equal R-ABC (the driver stops 'nothing new'); nothing to read")
        pred = ov(base_rabc) + sum(adopted.values())
        d = ov(arm) - pred
        order = [n for n in ("707 hook", "705 theta 0.25") if n in adopted]
        verdict = bin_of(d) + (f"; revert order: {' then '.join(order)}" if d <= -0.0005 else "")
        out = {"read": "695 R-ABC2", "arm_entry": arm.get("entry_id"), "arm_overall": ov(arm), "rabc_overall": ov(base_rabc),
               "adopted_deltas_vs_allsw": adopted, "pred": pred, "D": d, "verdict": verdict}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
