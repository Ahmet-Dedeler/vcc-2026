#!/usr/bin/env python3
"""653 SUGGESTED v4 of 537's rule1007_draft537.py v3 (sha 8367c89d), from 653's review. Changes, each tested in selftest():
  R1  the 10-08 list excludes arms actually scored/taken (--scored=TAG,..., plus CPH/K20 when their ids are given), not the
      heads of the recomputed 10-07 lists (v3 dropped 459dr whenever K20's id was passed on 10-08)
  R2  a missing board-entry capture reads MISSING (treated as UNKNOWN, printed) instead of crashing the rule
  R3  every DE input must keep pds and mse within 1e-6 raw of its base (537 section 5.2); otherwise its bin is LEAK (off)
  R4  CPH SHIP also needs the overall delta vs 446bj2z > 0 (fid + jac alone ignored reach and nmae)
  R8  MZ / CPMZ (the Arc hook) only with --arc-hook (652 F1 recommends dropping it from the October path)
  R11 every compared pair must share anchor_version (scaled deltas are anchor-dependent)
Original v3 docstring follows.

537 — a DRAFT of the 10-07 / 10-08 slot rule for the coordinator (not armed; the coordinator fixes the real rule before
any 10-06 entry exists). It reads board entries that do not exist yet (like 498/rule1004.py) and prints candidate lists; the
coordinator's chain takes the first READY, not-yet-scored package of each list (novelty and gates are the chain's job).

Inputs (entry ids or `none`):
  cp   446bj2zcp (10-01 s2) vs 446bj2z            535's bins on d(fid + jac) overall: POS >= +0.0012, PARTIAL >= +0.0004,
                                                   NULL > -0.0004, else NEG; "rescue" when NULL/NEG with d jac raw >= +0.002
                                                   and d fid < 0 (the composition transferred, the unsteered signs cost fid)
  h    496a (10-05 s1) vs 446bj2                  496a's bins on d nmae raw: TRANSFERS <= -0.0010, PARTIAL <= -0.0003
  k16  446bj2zk16 (10-06 s2, or s1) vs 446bj2z    523's bins on d(fid + jac) overall: POS >= +0.0012, NEG <= -0.0005, else MID
  k20  446bj2zk20 if it was scored on 10-06, else none
  c    495 (10-03 s1) vs 446bj2                   495's bins: FAILS if overall < +0.0004 or d mse > -0.002 raw
  cph  446bj2zcph (10-07 s1) vs 446bj2z, for the 10-08 list only (v3): SHIP when d(fid + jac) overall >= +0.0012 (535's
       bin with the hook on: the block ships in October only then); none before 10-07's read
Arms (tags):
  CPH   446bj2zcph   = 446bj2z + VCC_RC_CALLSET_BLOCK=prop:<528 abc artifact> + VCC_DE_MDE_SIGN=clean4:0.05     (build after 10-01)
  CPMZ  446bj2zmzcp  = 446bj2z + VCC_RC_CALLSET_BLOCK=prop:<...> + VCC_DE_MDE_SIGN=file:<517 mz artifact>:...:0.05:0.05:le1
  CPHK20 446bj2zcphk20 = CPH + VCC_RC_CALLSET_CPM=1600 (the hook removes the block's smaller-K leak; RESULTS §3.2 P5)
  MZ    446bj2zmz    (READY, 527)
  K20   446bj2zk20   (READY, 523)
  ALLSW the all-switches production arm (built after 10-05's read; RESULTS §5.2)
  S1    518's strength-key arm (only if 495 FAILS)
  459dr (READY)
Why this order (RESULTS §3.2, §5): on H1 the CLEAN4 hook repairs the call-propensity block's fid cost (the hook's signs reach
the called, unsteered swapped-in genes) and the two are super-additive; so if the block reads >= PARTIAL, the combination
is the October DE form to read first. v2 (2026-09-30 19:25Z, RESULTS §3.2 P12): the hook SOURCE test (CPMZ vs CPH) moves
from 10-07 s2 to 10-08: without sign information the Arc hook leaks on the block (−0.0030 vs CPH at f 0.35) while CLEAN4's
costs ≈ 0, so a CPMZ read changes October's source only if it beats its r 0.37 prediction by ~2× — little D/E/F value for a
day-1 slot. 459dr (October's H1-dose rule) takes 10-07 s2 and heads the tail; S1 (A/B/C-only) goes last.
Without the hook the block does not stack with a smaller K, so the K dose partner (K20) is read only when the block is off;
with the hook, CPHK20 is the 10-08 arm when the board reads the small-n_real regime (k16 POS).
v3 (2026-10-01 01:15Z, after the 10-01 read: 446bj2zcp NULL with the rescue flag): the 10-08 list takes CPH's read. The
block (with the hook) ships only if CPH reads SHIP; then CPHK20 (K POS) / CPMZ (H on) follow; otherwise the block is off for
October and the K dose partner K20 (k16 >= MID) and MZ (H on) take 10-08, exactly as in the block-off branch of 10-07.
usage: uv run python rule1007_draft537.py <cp id|none> <h id|none> <k16 id|none> <k20 id|none> <c id|none> [<cph id|none>] [--selftest]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
Z, BJ2 = "be90smZUDcIV7j5MgQ03", "YvVru3DAF89Jiyr34olR"      # 446bj2z, 446bj2


MISSING = []


def load(i, fake=None):
    if i in (None, "none"):
        return None
    if fake is not None:
        x = fake.get(i)
        if x is None:
            MISSING.append(i)
        return x
    f = ROOT / f"data/board-entries/{i}.json"
    if not f.is_file():                      # R2: a missing capture must not crash the rule
        MISSING.append(i)
        print(f"WARNING: no capture data/board-entries/{i}.json (run scripts/capture_board_entry.py); read as UNKNOWN",
              file=sys.stderr)
        return None
    d = json.loads(f.read_text())
    return d.get("entry", d)


def pinned(a, b):
    """R3 / R11: a DE arm must keep pds and mse within 1e-6 raw of its base, on the same anchors."""
    if a.get("anchor_version") != b.get("anchor_version"):
        raise SystemExit(f"ANCHOR MISMATCH {a.get('anchor_version')} vs {b.get('anchor_version')}: recompute the bins from raw")
    return abs(a[RAW["pds"]] - b[RAW["pds"]]) <= 1e-6 and abs(a[RAW["mse"]] - b[RAW["mse"]]) <= 1e-6


def dd(a, b):
    return {m: a[RAW[m]] - b[RAW[m]] for m in RAW}, {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}, a["score_avg"] - b["score_avg"]


def bins(cp, h, k16, c, fake=None, cph=None):
    z, bj2 = load(Z, fake), load(BJ2, fake)
    out = {}
    a = load(cp, fake)
    if a is None:
        out["CP"] = "UNKNOWN"
    else:
        r, s, _ = dd(a, z)
        p = (s["fid"] + s["jac"]) / 6
        out["CP"] = "POS" if p >= 0.0012 else "PARTIAL" if p >= 0.0004 else "NULL" if p > -0.0004 else "NEG"
        if not pinned(a, z):
            out["CP"] = "LEAK"
        # the composition transferred (jac up) but the unsteered swapped-in genes' signs cost fid: the hook is the repair
        out["CP_rescue"] = bool(out["CP"] in ("NULL", "NEG") and r["jac"] >= 0.002 and r["fid"] < 0)
    a = load(h, fake)
    if a is None:
        out["H"] = "UNKNOWN"
    else:
        r, _, _ = dd(a, bj2)
        out["H"] = "TRANSFERS" if r["nmae"] <= -0.0010 else "PARTIAL" if r["nmae"] <= -0.0003 else "FAILS"
        if not pinned(a, bj2):
            out["H"] = "LEAK"
    a = load(k16, fake)
    if a is None:
        out["K"] = "UNKNOWN"
    else:
        _, s, _ = dd(a, z)
        p = (s["fid"] + s["jac"]) / 6
        out["K"] = "POS" if p >= 0.0012 else "NEG" if p <= -0.0005 else "MID"
        if not pinned(a, z):
            out["K"] = "LEAK"
    a = load(c, fake)
    if a is None:
        out["C"] = "UNKNOWN"
    else:
        r, _, dov = dd(a, bj2)
        out["C"] = "FAILS" if (dov < 0.0004 or r["mse"] > -0.002) else "OK"
    a = load(cph, fake)
    if a is None:
        out["CPH"] = "UNKNOWN"
    else:
        _, s, dov = dd(a, z)
        out["CPH"] = "SHIP" if ((s["fid"] + s["jac"]) / 6 >= 0.0012 and dov > 0) else "NO"   # R4: overall floor
        if not pinned(a, z):
            out["CPH"] = "LEAK"
    if MISSING:
        out["MISSING"] = ",".join(dict.fromkeys(MISSING))
    return out


def rule(B, k20_scored, scored=None, arc_hook=False):
    """Ordered candidate lists for 10-07 s1 / s2 and 10-08 (the chain takes the first READY, unscored package).
    On H1 the block + hook dominates K 1,600 / 2,000 in every regime and the block does not stack with a smaller K
    (RESULTS §3.2 P2/P3), so K20 is read only when the block is off."""
    cp_on = B["CP"] in ("POS", "PARTIAL") or B.get("CP_rescue", False)
    h_on = B["H"] in ("TRANSFERS", "PARTIAL") and arc_hook   # R8: h_on only gates the Arc-hook arms (MZ, CPMZ)
    # v3: a block that reached 10-07 only through the rescue flag is predicted NOT to ship (RESULTS §0 10-01 READ), so the
    # K dose partner is not held back by it: K20 runs beside CPH when k16 >= MID
    k_on = B["K"] in ("POS", "MID") and not k20_scored and B["CP"] not in ("POS", "PARTIAL")
    s1 = (["CPH"] if cp_on else []) + (["MZ"] if h_on and not cp_on else []) + (["K20"] if k_on else [])
    s2 = (["K20"] if k_on else []) + (["MZ"] if h_on and not cp_on else [])
    tail = ["459dr", "ALLSW"] + (["S1"] if B["C"] == "FAILS" else [])
    s1 = s1 + [x for x in tail if x not in s1]
    s2 = [x for x in s2 if x != s1[0]]
    s2 = s2 + [x for x in tail if x not in s2 and x != s1[0]]
    # 10-08 (its rule reads 10-06's k16 and, v3, 10-07's CPH): the block stays on only if CPH ships (or, before CPH is
    # read, if the block itself read >= PARTIAL); then the block + hook at K 2,000 when the board reads the small-n_real
    # regime and the hook source; otherwise the block-off arms
    block_on = B.get("CPH") == "SHIP" if B.get("CPH", "UNKNOWN") != "UNKNOWN" else B["CP"] in ("POS", "PARTIAL")
    k_off = B["K"] in ("POS", "MID") and not k20_scored and not block_on
    # R1: exclude what is KNOWN to be scored or taken. On 10-07 (CPH unread) the planned picks are the list heads. On 10-08
    # pass --scored=<tags scored on 10-07>; without it only CPH / K20 (their ids given) are excluded, never a recomputed head,
    # and the chain's novelty gate refuses any package that was in fact scored.
    if scored is not None:
        done = set(scored) | ({"K20"} if k20_scored else set())
    elif B.get("CPH", "UNKNOWN") == "UNKNOWN":
        done = {s1[0], s2[0]} | ({"K20"} if k20_scored else set())
    else:
        done = {"CPH"} | ({"K20"} if k20_scored else set())
    cand = ((["CPHK20"] if block_on and B["K"] == "POS" else []) + (["CPMZ"] if block_on and h_on else [])
            + (["K20"] if k_off else []) + (["MZ"] if h_on and not block_on else []) + tail)
    d8 = [x for x in dict.fromkeys(cand) if x not in done]
    return s1, s2, d8


def selftest():
    base = {m: 0.5 for m in RAW.values()} | {f"score_{m}": 0.0 for m in RAW} | {"score_avg": 0.27}

    def e(**kw):
        x = dict(base)
        for k, v in kw.items():
            x[k] = v
        return x
    fake = {Z: base, BJ2: base,
            "cpPOS": e(score_fid=0.006, score_jac=0.006), "cpNUL": e(score_fid=0.0), "hTR": e(de_wilcoxon_lfc_nmae=0.498),
            "hFA": e(de_wilcoxon_lfc_nmae=0.4999),
            "cpRES": e(score_fid=-0.012, score_jac=0.009, de_wilcoxon_sig_jaccard=0.503, de_wilcoxon_direction_fidelity_yield_raw=0.497), "kPOS": e(score_fid=0.01), "kNEG": e(score_fid=-0.01),
            "cOK": e(score_avg=0.273, expr_mse_unbiased_capped_norm=0.49), "cF": e(score_avg=0.2701)}
    cases = [(("cpPOS", "hTR", "kPOS", "cOK"), False, (["CPH", "459dr", "ALLSW"], ["459dr", "ALLSW"])),
             (("cpPOS", "hFA", "kNEG", "cOK"), False, (["CPH", "459dr", "ALLSW"], ["459dr", "ALLSW"])),
             (("cpNUL", "hTR", "kPOS", "cF"), False, (["MZ", "K20", "459dr", "ALLSW", "S1"], ["K20", "459dr", "ALLSW", "S1"])),
             (("cpNUL", "hFA", "kPOS", "cOK"), True, (["459dr", "ALLSW"], ["ALLSW"])),
             (("cpRES", "hFA", "kNEG", "cOK"), False, (["CPH", "459dr", "ALLSW"], ["459dr", "ALLSW"])),
             (("cpRES", "hFA", "kPOS", "cOK"), False, (["CPH", "K20", "459dr", "ALLSW"], ["K20", "459dr", "ALLSW"])),
             (("cpNUL", "hFA", "kPOS", "cOK"), False, (["K20", "459dr", "ALLSW"], ["459dr", "ALLSW"]))]
    # 10-08 (v2): CPHK20 when the block is on and k16 reads POS, then the hook-source test CPMZ when the hook transfers
    B = bins("cpPOS", "hTR", "kPOS", "cOK", fake)
    assert rule(B, False, arc_hook=True)[2][:2] == ["CPHK20", "CPMZ"], rule(B, False, arc_hook=True)[2]
    # 10-08 (v3): the rescue branch (CP NULL + rescue) with CPH read: SHIP keeps the block arms, NO switches to K20 / MZ
    fake["cphSHIP"] = e(score_fid=0.006, score_jac=0.004, score_avg=0.27 + 0.010 / 6)   # v4: overall consistent with the members
    fake["cphNO"] = e(score_fid=0.001, score_jac=0.001)
    B = bins("cpRES", "hTR", "kPOS", "cOK", fake, cph="cphSHIP")
    assert rule(B, False, arc_hook=True)[2][:2] == ["CPHK20", "CPMZ"], rule(B, False, arc_hook=True)[2]
    B = bins("cpRES", "hTR", "kPOS", "cOK", fake, cph="cphNO")
    s1, s2, d8 = rule(B, False, scored=["CPH", "K20"], arc_hook=True)     # v4: 10-08 states what 10-07 scored
    assert (s1[0], s2[0], d8[:2]) == ("CPH", "K20", ["MZ", "459dr"]), (s1, s2, d8)
    # R1: 10-08 with K20 scored on 10-07 and no hook: 459dr (never scored) heads the list, whichever way K20 is passed
    B = bins("cpRES", "hFA", "kPOS", "cOK", fake, cph="cphNO")
    assert rule(B, True)[2][:2] == ["459dr", "ALLSW"], rule(B, True)
    # without --scored and with k20 passed as none, K20 stays listed (the novelty gate skips it if 10-07 scored it)
    assert rule(B, False)[2][:3] == ["K20", "459dr", "ALLSW"], rule(B, False)
    assert rule(B, True, scored=["CPH", "K20"])[2][:2] == ["459dr", "ALLSW"]
    # ... and 459dr drops out once it was the one actually scored
    assert rule(B, True, scored=["CPH", "459dr"])[2][:1] == ["ALLSW"]
    # R8: without --arc-hook a transferring 496a adds no Arc-hook arm anywhere
    B = bins("cpRES", "hTR", "kPOS", "cOK", fake, cph="cphNO")
    s1, s2, d8 = rule(B, True)
    assert "MZ" not in s1 + s2 + d8 and "CPMZ" not in d8 and d8[:1] == ["459dr"], (s1, s2, d8)
    # R2: a missing capture reads UNKNOWN and is reported
    B = bins("cpRES", "hNOPE", "kPOS", "cOK", fake)
    assert B["H"] == "UNKNOWN" and "hNOPE" in B["MISSING"], B
    # R3: a bulk leak in k16 bins LEAK and is treated as off (no K20)
    fake["kLEAK"] = e(score_fid=0.01, pds_cosine=0.49)
    B = bins("cpRES", "hFA", "kLEAK", "cOK", fake)
    assert B["K"] == "LEAK" and "K20" not in rule(B, False)[1], B
    # R4: fid + jac over the bar but the overall delta <= 0 (reach / nmae paid) does not SHIP
    fake["cphREACH"] = e(score_fid=0.006, score_jac=0.002, score_reach=-0.04, score_avg=0.2695)
    B = bins("cpRES", "hFA", "kPOS", "cOK", fake, cph="cphREACH")
    assert B["CPH"] == "NO", B
    for (cp, h, k, c), k20s, want in cases:
        B = bins(cp, h, k, c, fake)
        s1, s2, _ = rule(B, k20s, arc_hook=True)          # v3's branches, with the Arc hook allowed as in v3
        assert (s1, s2) == want, (B, s1, s2, want)
    print(f"SELFTEST PASS ({len(cases)} synthetic branches + 10-08 checks + R1/R2/R3/R4/R8 regressions)")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    cp, h, k16, k20, c, cph = ([x for x in sys.argv[1:] if not x.startswith("--")] + ["none"] * 6)[:6]
    scored = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--scored=")), None)
    if cph not in (None, "none") and scored is None:
        print("WARNING: computing the 10-08 list without --scored=<tags scored on 10-07>; only CPH / K20 are excluded",
              file=sys.stderr)
    B = bins(cp, h, k16, c, cph=cph)
    s1, s2, d8 = rule(B, k20 not in (None, "none"), scored=scored, arc_hook="--arc-hook" in sys.argv)
    print(f"{' '.join(s1)}|{' '.join(s2)}|10-08: {' '.join(d8)}|" + ",".join(f"{k}={v}" for k, v in B.items()))


if __name__ == "__main__":
    main()
