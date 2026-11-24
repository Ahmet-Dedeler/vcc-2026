"""Synthetic captures for every FINAL-PHASE-FEEDBACK-PLAN section 2 row, including the edge cases the plan names.

  uv run python experiments/531-cloud-final-phase-tooling/test_readdef.py
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import readdef as R

RAW_KEYS = list(R.RAW.values())
SCALED_KEYS = list(R.SCALED.values())


def entry(eid, overall=0.0, raw=None, scaled=None, hide=()):
    raw = dict(raw or {})
    scaled = dict(scaled or {})
    e = {"entry_id": eid, "status": "published", "score_avg": overall}
    for m, k in R.RAW.items():
        e[k] = raw.get(m, 0.5)
    for m, k in R.SCALED.items():
        e[k] = scaled.get(m, 0.0)
    e["score_avg"] = overall
    for k in hide:
        e[k] = None
    return e


def wrap(e):
    return {"captured_at": "2026-10-23T00:00:00+00:00", "entry": e, "note": "synthetic"}


def shifted(base, eid, draw=0.0, sdraw=None, overall=None):
    """Copy base and add raw/scaled deltas. overall defaults to base overall + sum(scaled deltas)/6."""
    e = entry(eid)
    for m, k in R.RAW.items():
        e[k] = base[k] + (draw.get(m, 0.0) if isinstance(draw, dict) else 0.0)
    sd = sdraw or {}
    for m, k in R.SCALED.items():
        e[k] = base[k] + sd.get(m, 0.0)
    if overall is None:
        e["score_avg"] = base["score_avg"] + sum(sd.get(m, 0.0) for m in R.RAW) / 6.0
    else:
        e["score_avg"] = overall
    return e


def check(name, cond, detail=""):
    if not cond:
        raise SystemExit(f"FAIL {name} {detail}")
    print(f"PASS {name}")


def main() -> int:
    n = 0

    def t(name, row, arm, base, alt, code, adopt=None):
        nonlocal n
        v = R.verdict(row, arm, base, alt)
        n += 1
        if v["code"] != code or (adopt is not None and v["adopt"] is not adopt):
            raise SystemExit(
                f"FAIL {name}: code {v['code']!r} adopt {v['adopt']!r} bin {v['bin']!r}\n"
                f"ADOPTION: {v['adoption']}\nexpected code {code!r} adopt {adopt!r}"
            )
        # The rendered verdict carries the bin sentence the plan uses.
        text = R.render(v)
        if "BIN:" not in text or "ADOPTION:" not in text:
            raise SystemExit(f"FAIL {name}: render missing BIN/ADOPTION")
        print(f"PASS {name} [{v['code']}]")

    base = entry("R0", overall=0.20, scaled={m: 0.20 for m in R.RAW})

    # --- row 1: baseline, no adoption even if the pair differs ---
    t("r1 baseline", "1", shifted(base, "X", sdraw={"pds": 1.0}, overall=0.9), base, None, "baseline", None)

    # --- row 2: bulk, exact. both lose; better side; tie; one side; tiny exact win ---
    lose_a = shifted(base, "r98", sdraw={"pds": -0.006, "mse": -0.006})
    lose_b = shifted(base, "r102", sdraw={"pds": -0.003, "mse": -0.003})
    t("r2 both lose", "2", lose_a, base, lose_b, "base_stands", False)
    v = R.verdict("2", lose_a, base, lose_b)
    check("r2 both-lose words", "R0's ρ* stands" in v["adoption"] and "both sides lose" in v["bin"])

    win_a = shifted(base, "r98", sdraw={"pds": 0.012, "mse": 0.006})  # bulk (0.018)/6
    win_b = shifted(base, "r102", sdraw={"pds": 0.006, "mse": 0.0})
    t("r2 better side is first", "2", win_a, base, win_b, "adopt_arm", True)
    t("r2 better side is second", "2", win_b, base, win_a, "adopt_alt", True)
    t("r2 one side wins", "2", win_a, base, lose_b, "adopt_arm", True)
    t("r2 one side only, wins", "2", win_a, base, None, "adopt_arm", True)
    t("r2 one side only, loses", "2", lose_a, base, None, "base_stands", False)

    tie_a = shifted(base, "r98", sdraw={"pds": 0.006, "mse": 0.0})
    tie_b = shifted(base, "r102", sdraw={"mse": 0.006, "pds": 0.0})
    t("r2 tie both beat", "2", tie_a, base, tie_b, "tie", False)

    tiny = shifted(base, "tiny", sdraw={"pds": 1e-12})
    t("r2 exact tiny bulk win", "2", tiny, base, None, "adopt_arm", True)
    check("r2 tiny not rounded away in the bin", "beats the base" in R.verdict("2", tiny, base, None)["bin"])
    zero = shifted(base, "zero", sdraw={})
    t("r2 exact zero is not a win", "2", zero, base, None, "base_stands", False)

    # --- row 3: q curve ---
    q15 = shifted(base, "q15", sdraw={"pds": 0.030})  # overall +0.005
    q20_worse = shifted(base, "q20", sdraw={"pds": -0.030})
    q20_better = shifted(base, "q20", sdraw={"pds": 0.060})
    q20_mid = shifted(base, "q20", sdraw={"pds": 0.012})  # beats base, loses to q15
    t("r3 curve turns, higher loses to base", "3", q15, base, q20_worse, "adopt_arm", True)
    check("r3 turn words", "turned" in R.verdict("3", q15, base, q20_worse)["bin"])
    t("r3 curve turns, higher still beats base but loses to lower", "3", q15, base, q20_mid, "adopt_arm", True)
    t("r3 curve has not turned", "3", q15, base, q20_better, "adopt_alt", True)
    check("r3 not-turned words", "has not turned" in R.verdict("3", q15, base, q20_better)["bin"])
    t("r3 both lose", "3", q20_worse, base, shifted(base, "q20b", sdraw={"mse": -0.01}), "base_stands", False)
    t("r3 lower loses, higher wins", "3", q20_worse, base, q20_better, "adopt_alt", True)
    t("r3 one q only", "3", q15, base, None, "adopt_arm", True)
    t("r3 tie", "3", q15, base, shifted(base, "q20", sdraw={"pds": 0.030}), "tie", False)
    # exact overall equality with the base is not better
    t("r3 exact overall tie with base", "3", shifted(base, "q", overall=base["score_avg"]), base, None, "base_stands", False)

    # --- row 4: two independent switches, exact overall ---
    up = shifted(base, "cap", sdraw={"mse": 1e-9})
    down = shifted(base, "cd4", sdraw={"mse": -1e-9})
    t("r4 tiny exact win", "4", up, base, None, "adopt_arm", True)
    t("r4 loss", "4", down, base, None, "base_stands", False)
    t("r4 both win", "4", up, base, shifted(base, "cd4", sdraw={"pds": 0.01}), "adopt_both", True)
    t("r4 only second wins", "4", down, base, up, "adopt_alt", True)
    t("r4 neither", "4", down, base, shifted(base, "cd4", sdraw={"pds": -0.01}), "base_stands", False)

    # --- row 5: overall exact, DE veto beyond 0.009 when bulk rises ---
    # bulk +0.020, DE -0.008, overall +0.012 → adopt
    ok_w = shifted(base, "w0", sdraw={"pds": 0.060, "mse": 0.060, "fid": -0.012, "reach": -0.012, "jac": -0.012, "nmae": -0.012})
    # bulk (0.12)/6 = 0.02; de (-0.048)/6 = -0.008; overall = +0.012
    t("r5 adopt under the DE cap", "5", ok_w, base, None, "adopt", True)
    # DE cost exactly 0.009 is not beyond
    edge = shifted(base, "w1", sdraw={"pds": 0.060, "fid": -0.054})
    # bulk 0.01, de -0.009, overall +0.001
    d_edge = R.delta(edge, base)
    check("r5 fixture DE cost is exactly 0.009", abs(d_edge["de"] + 0.009) < 1e-15, str(d_edge["de"]))
    t("r5 DE cost of exactly 0.009 does not veto", "5", edge, base, None, "adopt", True)
    # just beyond
    beyond = shifted(base, "w1", sdraw={"pds": 0.120, "fid": -0.054, "nmae": -1e-9})
    d_beyond = R.delta(beyond, base)
    check("r5 fixture DE cost is beyond 0.009", d_beyond["de"] < -0.009 and d_beyond["bulk"] > 0 and d_beyond["overall"] > 0, str(d_beyond))
    t("r5 veto when DE cost is beyond 0.009", "5", beyond, base, None, "veto_de", False)
    check("r5 veto words", "NOT adopted" in R.verdict("5", beyond, base, None)["adoption"])
    # overall not better, even with a healthy DE
    t("r5 not overall better", "5", shifted(base, "w", sdraw={"fid": 0.01, "pds": -0.03}), base, None, "not_better", False)
    # two w: the vetoed one loses, the clean one is adopted
    t("r5 two w, only the clean one", "5", beyond, base, ok_w, "adopt_alt", True)
    # both clear: the better overall
    ok_w2 = shifted(base, "w0", sdraw={"pds": 0.030, "mse": 0.030})
    t("r5 two w, better overall", "5", ok_w, base, ok_w2, "adopt_arm", True)

    # --- row 6: prior-only. 0.009 does not apply. threshold on fid+nmae+reach overall ---
    # fnr = (0.003+0.003+0.003)/6 = 0.0015 exactly, overall contribution of all six can be under 0.009
    hit = shifted(base, "prog", sdraw={"fid": 0.003, "nmae": 0.003, "reach": 0.003})
    d_hit = R.delta(hit, base)
    check("r6 fixture fnr exactly 0.0015 and overall under 0.009",
          abs(d_hit["fnr"] - 0.0015) < 1e-15 and d_hit["overall"] < 0.009, str(d_hit))
    t("r6 exact threshold adopts under the noise floor", "6", hit, base, None, "adopt", True)
    check("r6 says prior-only", "Prior-only" in R.verdict("6", hit, base, None)["adoption"])
    miss = shifted(base, "prog", sdraw={"fid": 0.003, "nmae": 0.003, "reach": 0.002999})
    t("r6 just under 0.0015", "6", miss, base, None, "below_bin", False)
    # large overall from pds does not adopt when the three-member sum misses
    bulk_only = shifted(base, "prog", sdraw={"pds": 0.20})
    check("r6 bulk-only overall is large", R.delta(bulk_only, base)["overall"] > 0.009)
    t("r6 large overall does not adopt when fid+nmae+reach miss", "6", bulk_only, base, None, "below_bin", False)

    # --- row 7: read468 boundaries. prior-only ---
    def g(eid, overall_d, raw_d=None):
        e = shifted(base, eid, draw=raw_d or {}, sdraw={}, overall=base["score_avg"] + overall_d)
        return e

    t("r7 transfers", "7", g("gate", 0.0015), base, None, "adopt_transfers", True)
    t("r7 just under transfers is weak", "7", g("gate", 0.0014999), base, None, "adopt_weak", True)
    t("r7 weak boundary 0.0005", "7", g("gate", 0.0005), base, None, "adopt_weak", True)
    t("r7 just under 0.0005 is NULL and still adopts", "7", g("gate", 0.0004999), base, None, "adopt_null", True)
    t("r7 zero is NULL", "7", g("gate", 0.0), base, None, "adopt_null", True)
    t("r7 just above harm is NULL", "7", g("gate", -0.0004999), base, None, "adopt_null", True)
    t("r7 harm boundary", "7", g("gate", -0.0005), base, None, "harm", False)
    t("r7 harm", "7", g("gate", -0.001), base, None, "harm", False)
    check("r7 harm closes", "closes the gate" in R.verdict("7", g("gate", -0.0005), base, None)["adoption"])
    check("r7 null keeps", "NULL read keeps" in R.verdict("7", g("gate", 0.0), base, None)["adoption"])
    # under 0.009 and still a real read
    check("r7 weak is under 0.009 and adopts",
          R.verdict("7", g("gate", 0.0006), base, None)["adopt"] is True)
    # falsifiers: exactly 1e-6 is ok; just over is a bug. A positive overall must not override a bug.
    t("r7 pds at exactly 1e-6 is not a bug", "7", g("gate", 0.002, {"pds": 1e-6, "mse": -1e-6}), base, None, "adopt_transfers", True)
    t("r7 pds just over 1e-6 is a bug", "7", g("gate", 0.002, {"pds": 1.1e-6}), base, None, "bug", False)
    t("r7 jac at exactly 0.0010 is not a bug", "7", g("gate", 0.002, {"jac": 0.0010}), base, None, "adopt_transfers", True)
    t("r7 jac just over 0.0010 is a bug", "7", g("gate", 0.002, {"jac": 0.0010001}), base, None, "bug", False)
    # mechanism words
    mech = g("gate", 0.002, {"reach": 0.003, "fid": -0.0001})
    check("r7 mechanism transfers", "H1 mechanism transfers" in R.verdict("7", mech, base, None)["bin"])
    mech_bad = g("gate", -0.001, {"reach": -0.001, "fid": -0.001})
    check("r7 mechanism removes signal", "only removes signal" in R.verdict("7", mech_bad, base, None)["bin"])

    # --- row 8-9: nmae <= -0.002, pinned pds/mse, overall > 0.009. NOT prior-only ---
    def scale(eid, nmae, overall_d, pds=0.0, mse=0.0):
        return shifted(base, eid, draw={"nmae": nmae, "pds": pds, "mse": mse}, sdraw={"nmae": overall_d * 6}, overall=base["score_avg"] + overall_d)

    good = scale("s125", -0.002, 0.0090001)
    t("r8-9 exact nmae and just over the noise", "8-9", good, base, None, "adopt", True)
    t("r8 alias", "8", good, base, None, "adopt", True)
    t("r9 alias", "9", good, base, None, "adopt", True)
    # 0.20 + 0.009 is a hair over 0.009 once subtracted back. Store the threshold itself.
    exact_noise = scale("s", -0.003, 0.0)
    exact_noise["score_avg"] = 0.009
    base_zero = dict(base)
    base_zero["score_avg"] = 0.0
    t("r8-9 overall exactly 0.009 does not clear", "8-9", exact_noise, base_zero, None, "below_noise", False)
    check("r8-9 below-noise asks for a second seed",
          "two seeds" in R.verdict("8-9", scale("s", -0.003, 0.005), base, None)["adoption"])
    t("r8-9 nmae just worse than -0.002", "8-9", scale("s", -0.001999, 0.02), base, None, "nmae_miss", False)
    t("r8-9 pds not pinned", "8-9", scale("s", -0.01, 0.02, pds=1.1e-6), base, None, "not_pinned", False)
    t("r8-9 pds at exactly 1e-6 is pinned", "8-9", scale("s", -0.01, 0.02, pds=1e-6, mse=-1e-6), base, None, "adopt", True)
    # a small overall with a good nmae must NOT adopt (the prior-only exception does not apply)
    t("r8-9 small overall is not prior-only", "8-9", scale("s", -0.01, 0.002), base, None, "below_noise", False)
    # two scales: both clear, the more negative nmae wins
    s_lo = scale("s125", -0.004, 0.02)
    s_hi = scale("s16", -0.003, 0.03)
    t("r8-9 both clear, more negative nmae", "8-9", s_lo, base, s_hi, "adopt_arm", True)
    t("r8-9 only the second clears", "8-9", scale("s125", -0.001, 0.02), base, s_hi, "adopt_alt", True)
    t("r8-9 neither clears", "8-9", scale("s125", -0.001, 0.02), base, scale("s16", -0.01, 0.001), "none", False)
    same = scale("s16", -0.004, 0.02)
    t("r8-9 tie of nmae", "8-9", s_lo, base, same, "tie", False)

    # --- row 10-12: combination must be scored; exact overall ---
    t("r10 adopt when better", "10-12", shifted(base, "combo", sdraw={"pds": 1e-12}), base, None, "adopt", True)
    t("r10 keep base when worse", "10", shifted(base, "combo", sdraw={"pds": -0.01}), base, None, "scored_not_better", False)
    t("r10 exact tie keeps the base", "11", shifted(base, "combo"), base, None, "scored_not_better", False)
    t("r12 alias is the combination read", "12", shifted(base, "combo", sdraw={"mse": 0.01}), base, None, "adopt", True)
    check("r10 mentions rule 2 and the row-2 re-check",
          "rule 2" in R.verdict("10-12", shifted(base, "combo", sdraw={"pds": 0.01}), base, None)["adoption"]
          and "row-2" in R.verdict("10-12", shifted(base, "combo", sdraw={"pds": 0.01}), base, None)["adoption"])

    # --- row 13-14: reserve, not a read ---
    t("r13 reserve", "13-14", shifted(base, "X", sdraw={"pds": 1}), base, None, "reserve", None)
    t("r14 alias", "14", base, base, None, "reserve", None)

    # --- members hidden: every reading row, and it must not adopt ---
    hidden = entry("hid", hide=("score_avg", "score_pds"))
    for row in ("2", "3", "4", "5", "6", "7", "8-9", "10-12"):
        t(f"hidden {row}", row, hidden, base, None, "hidden", None)
    # alt hidden on a two-sided row blocks the both-lose / better-side sentence
    t("r2 alt hidden", "2", win_a, base, hidden, "hidden", None)
    # a hidden base blocks even a strong arm
    t("r6 base hidden", "6", hit, hidden, None, "hidden", None)
    check("hidden adoption names the dark plan",
          "dark plan" in R.verdict("6", hidden, base, None)["adoption"])
    # rows that do not read members still return their own code
    t("r1 ignores hidden", "1", hidden, base, None, "baseline", None)
    t("r13 ignores hidden", "13", hidden, hidden, None, "reserve", None)

    # --- loader: wrapper, bare entry, entry id, missing file ---
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        wrapped = td / "wrapped.json"
        bare = td / "bare.json"
        wrapped.write_text(json.dumps(wrap(hit)))
        bare.write_text(json.dumps(base))
        loaded = R.load_entry(str(wrapped))
        check("loader wrapper", loaded["entry_id"] == "prog" and loaded["score_avg"] == hit["score_avg"])
        v = R.verdict("6", R.load_entry(str(wrapped)), R.load_entry(str(bare)), None)
        check("loader through verdict", v["code"] == "adopt")
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parent / "readdef.py"), "6", str(wrapped), str(bare)],
            capture_output=True, text=True,
        )
        check("cli exit 0", proc.returncode == 0, proc.stderr)
        check("cli prints bin", "BIN:" in proc.stdout and "ADOPT the programme source" in proc.stdout, proc.stdout)
        proc_h = subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parent / "readdef.py"), "2", str(wrapped), str(td / "nope.json")],
            capture_output=True, text=True,
        )
        # nope.json does not exist: usage-style SystemExit from load_entry
        check("cli missing file", proc_h.returncode == 1, proc_h.stderr)
        hid_path = td / "hid.json"
        hid_path.write_text(json.dumps(wrap(hidden)))
        proc_hid = subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parent / "readdef.py"), "8-9", str(hid_path), str(bare)],
            capture_output=True, text=True,
        )
        check("cli hidden exit 2", proc_hid.returncode == 2 and "MEMBERS HIDDEN" in proc_hid.stdout, proc_hid.stdout)

    real = R.load_entry("be90smZUDcIV7j5MgQ03")
    check("entry id resolves a real capture", abs(real["score_avg"] - 0.2728259609816559) < 1e-15, str(real.get("score_avg")))
    # two real captures still produce a row-2 verdict rather than crashing
    other = R.load_entry("iUziqE1gewfcPKzYUDaB")
    v = R.verdict("2", real, other, None)
    check("real pair yields a bulk verdict", v["code"] in ("adopt_arm", "base_stands") and "mse+pds" in v["bin"], v["bin"])

    try:
        R.verdict("0", base, base, None)
        raise SystemExit("FAIL unknown row should exit")
    except SystemExit as e:
        check("unknown row", "unknown row" in str(e))

    print(f"\n{n} row-verdicts passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
