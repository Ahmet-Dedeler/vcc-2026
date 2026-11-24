#!/usr/bin/env python3
"""10-05 board reader (coordinator, written 2026-10-04 ~20:15Z, BEFORE either 10-05 entry exists). Prints bins, writes nothing.

usage: uv run python experiments/663-idall3-arm/read1005.py [--a496 <entry id>] [--idall3 <entry id>]

It reads captured board entries (data/board-entries/<id>.json, written by scripts/capture_board_entry.py) and applies the
registered rules unchanged:

496a vs 446bj2 (YvVru3DAF89Jiyr34olR), per 496-uncalled-sign-hook/PREREGISTRATION-board.md:
  - falsifier first: pds and mse raw equal to the base's within 1e-6 (the prereg says "exactly"; DE-only arms move
    mse by O(1e-8) raw, so the project's DE-arm falsifier is 1e-6: 537 §5.2, 592, 653 R3; exact equality is reported);
  - bins on Δnmae raw: TRANSFERS <= -0.0010, PARTIAL in (-0.0010, -0.0003], FAILS > -0.0003;
  - leak check: |Δfid| or |Δjac| > 0.001 raw, or |Δreach| > 0.01 raw, means diagnose;
  - plus PRED-662-496a's instrument reading (I-662 h1-test -0.0022 [-0.0027, -0.0016], 618 -0.0014 [-0.0017, -0.0011]) and
    CPH's T = Δnmae / -0.001979 (662/PRED-662-CPH.md).
idall3 vs 498p2g (JWSghc0o7MN5D6o8Iadc), per 612-key-gated-hook/PREREGISTRATION-draft-IDALL3.md:
  - falsifier: pds or mse differ from the base's by > 1e-6 raw, or |Δjac| > 0.001 raw;
  - primary: Δ(fid + nmae) overall = (Δscore_fid + Δscore_nmae) / 6, arm - base, read against Δ_496a (the same statistic
    for 496a - 446bj2);
  - HARM if the primary is < -0.0002 against its own base; otherwise WIN if it is >= Δ_496a + 0.0002 AND Δreach raw lies
    inside ±0.0042, TIE if |primary - Δ_496a| < 0.0002, BELOW if it is < Δ_496a - 0.0002. A primary >= Δ_496a + 0.0002 with
    reach outside the band does not meet the WIN bin's condition and prints as "WIN-on-primary, reach OUT";
  - plus the registered predictions: 618 TIE +0.00007 [-0.00002, +0.00015] (PRED-661-idall3), I-662 +0.00020
    [+0.00006, +0.00028] (its I-662 ALTERNATIVE), h1-test WIN +0.0007 ... +0.0010 (the prereg).
Scaled members compare only between entries with the same anchor_version; the reader checks it.
"""
import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASE_446BJ2 = "YvVru3DAF89Jiyr34olR"
BASE_498P2G = "JWSghc0o7MN5D6o8Iadc"
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard", "nmae": "de_wilcoxon_lfc_nmae"}


def entry(eid):
    p = REPO / "data/board-entries" / f"{eid}.json"
    if not p.exists():
        sys.exit(f"{p} is missing: capture it first (scripts/capture_board_entry.py {eid})")
    e = json.loads(p.read_text())["entry"]
    if e.get("status") != "published" or e.get("score_avg") is None:
        sys.exit(f"{eid}: status {e.get('status')}, score_avg {e.get('score_avg')} -- not a scored entry")
    return e


def deltas(arm, base):
    if arm["anchor_version"] != base["anchor_version"]:
        print(f"  WARNING anchor_version differs ({arm['anchor_version']} vs {base['anchor_version']}): scaled members are not comparable")
    raw = {k: arm[v] - base[v] for k, v in RAW.items()}
    sc = {k: arm[f"score_{k}"] - base[f"score_{k}"] for k in RAW}
    sc["avg"] = arm["score_avg"] - base["score_avg"]
    print("  raw Δ:    " + "  ".join(f"{k} {raw[k]:+.6f}" for k in RAW))
    print("  scaled Δ: " + "  ".join(f"{k} {sc[k]:+.6f}" for k in RAW) + f"  | overall {sc['avg']:+.7f}"
          f" ({arm['score_avg']:.7f} vs {base['score_avg']:.7f})")
    return raw, sc


def read496a(eid):
    a, b = entry(eid), entry(BASE_446BJ2)
    print(f"496a {eid} vs 446bj2 {BASE_446BJ2}")
    raw, sc = deltas(a, b)
    # 496's prereg (09-27) says "exactly". Corrected 10-05 01:10Z, before idall3's read: DE-only arms move the scorer's
    # jackknife-corrected mse by O(1e-8) raw (446bja +2.4e-8, 458d/g -5.3e-8, 446bj2zcp -2.3e-8, K20 +4.1e-8), so the
    # project's DE-arm falsifier is 1e-6 raw (537 section 5.2, 592's 10-01 correction in 535/BOARD-RESULT-446bj2zcp, 653 R3).
    # Exact equality is still reported.
    bitid = raw["pds"] == 0 and raw["mse"] == 0
    exact = abs(raw["pds"]) <= 1e-6 and abs(raw["mse"]) <= 1e-6
    print(f"  falsifier (pds, mse raw within 1e-6; bit-identical: {bitid}): {'PASS' if exact else 'FAIL -- diagnose before reading anything else'}"
          f" (Δpds {raw['pds']:+.3e}, Δmse {raw['mse']:+.3e})")
    leak = abs(raw["fid"]) > 0.001 or abs(raw["jac"]) > 0.001 or abs(raw["reach"]) > 0.01
    print(f"  leak check (|Δfid|, |Δjac| <= 0.001, |Δreach| <= 0.01 raw): {'DIAGNOSE' if leak else 'clean'}")
    d = raw["nmae"]
    b496 = "TRANSFERS" if d <= -0.0010 else ("PARTIAL" if d <= -0.0003 else "FAILS")
    print(f"  496 bin (Δnmae raw {d:+.5f}): {b496}")
    inst = ("both I-662 panels hold" if -0.0017 <= d <= -0.0016 else
            "below -0.0027: both under-read (transfer above the delivery multipliers)" if d < -0.0027 else
            "below -0.0017: h1-test's level reads it, 618 under-reads" if d < -0.0017 else
            "above -0.0016: 618 reads it, h1-test over-reads")
    print(f"  PRED-662-496a: h1-test -0.0022 [-0.0027, -0.0016], 618 -0.0014 [-0.0017, -0.0011] -> {inst}")
    print(f"  CPH's T = Δnmae / -0.001979 = {d / -0.001979:.3f} (I-662 h1-test point gives 1.12)")
    prim = (sc["fid"] + sc["nmae"]) / 6
    print(f"  Δ_496a = Δ(fid + nmae) overall = {prim:+.6f}")
    return prim, exact and not leak


def readidall3(eid, d496):
    a, b = entry(eid), entry(BASE_498P2G)
    print(f"idall3 {eid} vs 498p2g {BASE_498P2G}")
    raw, sc = deltas(a, b)
    fals = abs(raw["pds"]) > 1e-6 or abs(raw["mse"]) > 1e-6 or abs(raw["jac"]) > 0.001
    print(f"  falsifier (pds, mse within 1e-6 raw; |Δjac| <= 0.001 raw): {'FAIL -- a delivery leak, stop and audit' if fals else 'PASS'}")
    prim = (sc["fid"] + sc["nmae"]) / 6
    inband = abs(raw["reach"]) <= 0.0042
    print(f"  primary Δ(fid + nmae) overall {prim:+.6f}; Δreach raw {raw['reach']:+.5f} ({'inside' if inband else 'OUTSIDE'} ±0.0042)")
    if d496 is None:
        print("  Δ_496a unknown (pass --a496): only HARM can be read")
        print(f"  {'HARM' if prim < -0.0002 else 'not HARM'}")
        return
    x = prim - d496
    if prim < -0.0002:
        b612 = "HARM"
    elif x >= 0.0002:
        b612 = "WIN" if inband else "WIN-on-primary, reach OUT (the WIN bin's reach condition fails)"
    elif x > -0.0002:
        b612 = "TIE"
    else:
        b612 = "BELOW"
    print(f"  primary - Δ_496a = {x:+.6f} -> {b612}")
    near = ("618 TIE band [-0.00002, +0.00015]" if -0.00002 <= x <= 0.00015 else None,
            "I-662 band [+0.00006, +0.00028]" if 0.00006 <= x <= 0.00028 else None,
            "h1-test WIN band [+0.0007, +0.0010]" if 0.0007 <= x <= 0.0010 else None)
    print(f"  registered predictions containing the read: {', '.join(n for n in near if n) or 'none'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a496")
    ap.add_argument("--idall3")
    a = ap.parse_args()
    d496 = None
    if a.a496:
        d496, ok = read496a(a.a496)
        if not ok:
            print("  (496a failed its falsifier or leak check: Δ_496a is reported but must not anchor idall3's bins)")
            d496 = None
    if a.idall3:
        readidall3(a.idall3, d496)


if __name__ == "__main__":
    main()
