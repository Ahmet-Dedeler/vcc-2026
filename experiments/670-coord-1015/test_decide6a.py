#!/usr/bin/env python3
"""Tests for decide6a.py (SIXA-DECISION + Amendments 1-5) on synthetic inputs, plus one end-to-end run on the 10-16 rehearsal's
real decision708 file and the real board entries.   uv run python experiments/670-coord-1015/test_decide6a.py"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import decide6a as D  # noqa: E402

OK = dict(rabc2h="COMPOSES: ships", c1="ABOVE: open", c0="ABOVE: open", w1="NOT READ", w2="NOT READ")


def rec(point, lo, folds, n=3, lt3=0.0001, per=None):
    per = per if per is not None else [point] * n
    return dict(point=point, ci95=[lo, point + 0.002], folds_pos=folds, n_folds=n, leave_top3_out=lt3,
                per_fold=[dict(fold=f"F{i}", overall=v) for i, v in enumerate(per)])


def agg(c1b1=True, c0b0=True, c1c0=None):
    r = {}
    for cand, base, good in (("C1", "B1", c1b1), ("C0", "B0", c0b0)):
        r[f"{cand}|{base}|DEFo|pds"] = rec(0.002, 0.0005, 3) if good else rec(0.0005, -0.0004, 2, lt3=-0.0001)
        r[f"{cand}|{base}|DEF|pds"] = rec(0.0015, 0.0002, 3) if good else rec(-0.0002, -0.001, 1)
    if c1c0 is not None:
        lo = c1c0[2] if len(c1c0) > 2 else -0.001
        r["C1|C0|DEFo|pds"] = rec(c1c0[0], lo, 0, per=c1c0[1])
    return r


def d708(weights="h2c1", rider=True, vetoed=False, swap=False):
    return dict(weights=weights, rider=rider, vetoed=vetoed, swap=swap)


def run(d, bins=OK, res=None, lz="PASS", clean="YES", g1="NOTBUILT", g2="NOTBUILT"):
    return D.decide(d, bins, res if res is not None else agg(), lz, clean, g1, g2)


CASES = []


def case(f):
    CASES.append(f)
    return f


@case
def t_swap_off_no_move():
    o = run(d708())
    assert not o["X1"]["eligible"] and o["X2"]["construction"] == "C0" and o["X2"]["eligible"], o
    assert o["X2"]["env"] == {"DEF_DAY_DUMP_HOOK": D.HOOK_C0} and o["ship"] == "PENDING: build X2", o


@case
def t_swap_on_c1_stays():
    o = run(d708(swap=True), res=agg(c1c0=(0.0009, [0.001, 0.0008, 0.0009], 0.0002)))
    assert o["X1"]["env"] == {"DEF_DAY_DUMP_HOOK": D.HOOK707} and o["X2"]["construction"] == "C1", o
    assert o["build_order"] == ["X1", "X2"] and o["X2"]["lz_increment"] == "C1 - B1", o


@case
def t_swap_on_c1c0_fails_ships_c0():
    o = run(d708(swap=True), res=agg(c1c0=(-0.0003, [-0.001, 0.0002, -0.0001])))
    assert o["X2"]["construction"] == "C0" and o["X2"]["lz_increment"] == "C0 - B0", o


@case
def t_swap_on_c1c0_leans_positive_but_bound_below_zero():
    """Amendment 4b: the LABC / ABCp pattern (point > 0, 2/3 folds, ci95[0] < 0) ships C0."""
    o = run(d708(swap=True), res=agg(c1c0=(0.00023, [0.0006, 0.0003, -0.0002], -0.0004)))
    assert o["X2"]["construction"] == "C0", o


@case
def t_swap_on_c1c0_two_of_three_but_point_neg():
    o = run(d708(swap=True), res=agg(c1c0=(-0.0001, [0.0001, 0.0002, -0.0007])))
    assert o["X2"]["construction"] == "C0", o


@case
def t_swap_on_c1c0_fails_c0_fails_keeps_c1():
    o = run(d708(swap=True), res=agg(c0b0=False, c1c0=(-0.0003, [-0.001, 0.0002, -0.0001])))
    assert o["X2"]["construction"] == "C1", o


@case
def t_swap_on_c1b1_fails_no_x2():
    o = run(d708(swap=True), res=agg(c1b1=False, c1c0=(0.001, [0.001] * 3)))
    assert not o["X2"]["eligible"] and o["X2"]["construction"] is None and o["X1"]["eligible"], o


@case
def t_swap_on_missing_c1c0_waits():
    o = run(d708(swap=True), res=agg())
    assert not o["X2"]["eligible"] and any("C1:C0 missing" in w for w in o["X2"]["why"]), o


@case
def t_swap_off_c0_harm_closes():
    o = run(d708(), bins=dict(OK, c0="HARM: closed"))
    assert not o["X2"]["eligible"], o


@case
def t_lz_fail_and_none():
    assert not run(d708(), lz="FAIL")["X2"]["eligible"]
    o = run(d708(), lz="NONE")
    assert not o["X2"]["eligible"] and any("not read yet" in w for w in o["X2"]["why"]), o


@case
def t_lab_not_clean():
    assert not run(d708(), clean="NO")["X2"]["eligible"]


@case
def t_weights_move_w2_below_drops_cd4():
    o = run(d708(weights="h2.5c1.5"), bins=dict(OK, w1="COMPOSES: ok", w2="BELOW: blocks"))
    assert o["X1"]["env"] == {"DEF_DAY_H1_STACK": "x2.5"} and not o["nested"], o
    assert o["X2"]["construction"] == "C0" and o["X2"]["env"] == {"DEF_DAY_DUMP_HOOK": D.HOOK_C0}, o


@case
def t_weights_move_w1_below_backstop():
    o = run(d708(weights="h3c1"), bins=dict(OK, w1="BELOW: blocks"))
    assert not o["X1"]["eligible"] and o["nested"], o


@case
def t_rider_veto_non_nested():
    o = run(d708(rider=False, vetoed=True))
    assert "DEF_DAY_LODO_STRENGTH" in o["X1"]["env"] and not o["nested"] and o["X2"]["construction"] == "C0", o


@case
def t_ship_order():
    assert run(d708(swap=True), res=agg(c1c0=(0.001, [0.001] * 3, 0.0003)), g1="PASS", g2="PASS")["ship"] == "X2"
    assert run(d708(swap=True), res=agg(c1c0=(0.001, [0.001] * 3, 0.0003)), g1="PASS", g2="FAIL")["ship"] == "X1"
    assert run(d708(swap=True), res=agg(c1c0=(0.001, [0.001] * 3, 0.0003)), g1="FAIL", g2="FAIL")["ship"] == "R0 (stays final)"
    assert run(d708(swap=True), res=agg(c1c0=(0.001, [0.001] * 3, 0.0003)), g1="PASS", g2="NOTBUILT")["ship"] == "PENDING: build X2"
    assert run(d708(swap=True), res=agg(c1b1=False), g1="NOTBUILT")["ship"] == "PENDING: build X1"
    assert run(d708(), lz="FAIL")["ship"] == "R0 (stays final)"


@case
def t_rabc2h_below_turns_swap_off():
    o = run(d708(swap=True), bins=dict(OK, rabc2h="BELOW: x"))
    assert o["swap_branch"] == "OFF" and not o["X1"]["eligible"] and o["X2"]["construction"] == "C0", o


@case
def t_end_to_end_rehearsal_file():
    """The 10-16 rehearsal's real decision708 + the real board entries (W1/W2 not yet read) + a synthetic agg json."""
    with tempfile.TemporaryDirectory() as t:
        p = Path(t) / "agg.json"
        p.write_text(json.dumps(dict(results=agg())))
        f = HERE.parents[1] / "experiments/708-def-panel-repricing/work/decision708-d708abcdef.txt"
        out = subprocess.run([sys.executable, str(HERE / "decide6a.py"), "--decision708", str(f), "--agg", str(p),
                              "--lz", "NONE", "--lab-clean", "YES", "--w1-entry", "", "--w2-entry", ""],
                             capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
        o = json.loads(out.stdout)
        assert o["swap_branch"] == "OFF" and o["weights"] == "h2c1" and not o["X1"]["eligible"], o
        assert o["board"]["c1"].startswith("ABOVE") and o["board"]["c0"].startswith("ABOVE"), o["board"]
        assert o["board"]["rabc2h"].startswith("COMPOSES"), o["board"]
        assert o["X2"]["construction"] == "C0" and not o["X2"]["eligible"], o   # LZ not read yet


if __name__ == "__main__":
    bad = 0
    for f in CASES:
        try:
            f(); print(f"PASS {f.__name__}")
        except AssertionError as e:
            bad += 1; print(f"FAIL {f.__name__}: {str(e)[:400]}")
    print(f"{len(CASES) - bad}/{len(CASES)} PASS")
    sys.exit(1 if bad else 0)
