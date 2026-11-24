"""804 — the old instruments' record on the known answers: what the instrument used at the time predicted BEFORE the board
read (quoted from each row's PREREGISTRATION / BOARD-RESULT / ledger text), graded by PROTOCOL.md's hit rule.
    uv run python experiments/804-known-answers/old_record804.py   -> work/old_record.json + a table
"""
import json
from pathlib import Path
W = Path(__file__).resolve().parent / "work"
B = json.loads((W / "board_pairs.json").read_text())["pairs"]
# row -> (instrument, predicted overall (lo, hi) or point, would it have shipped, source)
OLD = {
    "495-vs-446bj2": ("PLOLO A/B/C folds (420/518 route)", (0.0019, 0.0029), True, "BOARD-LEDGER 10-03: predicted +0.0019 … +0.0029"),
    "479p2g-vs-446bj2": ("PLOLO + October-like folds (479 prereg)", (0.0008, 0.0012), True, "479 PREREGISTRATION: +0.0012 / +0.0008"),
    "498p2g-vs-495": ("additivity of 479p2g on cap 6 (540c) / R6 C10b calibrated", (0.0006, 0.0013), True, "C10b calibrated +0.00093 [+0.00061, +0.00128]"),
    "479sp2-vs-446bj2": ("PLOLO + October-like folds (479 prereg)", (0.0008, 0.0013), True, "479 PREREGISTRATION: +0.0013 / +0.0008"),
    "702-vs-ALLSW": ("K562 GWPS seen + h1-test seen (702 amendment 1)", (0.0008, 0.0017), True, "702 RESULTS: predicted +0.0008 … +0.0017"),
    "446bj2z-vs-446bj2": ("H1 real emitter (456)", (0.00036, 0.0008), True, "ledger: nmae −0.0013 … −0.0029 raw predicted"),
    "496a-vs-446bj2": ("618 I-662 / h1-test", (0.00038, 0.0006), True, "ledger: 618 −0.0014, h1-test −0.0022 raw nmae"),
    "DN-vs-630arm": ("H1 real emitter (631 prereg)", (0.0005, 0.0013), True, "ledger: −0.0018 … −0.0048 raw nmae predicted"),
    "630arm-vs-446bj2z": ("H1 real emitter (630)", (0.0014, 0.0030), True, "ledger: H1 centre −0.005 … −0.011 raw nmae"),
    "446bj2zk20-vs-446bj2z": ("610 instrument card (f 0.35)", (0.0044, 0.0044), True, "ledger: 610 card predicted +0.0044"),
    "446bj2zk16-vs-446bj2z": ("I-662 (f 0.61)", (-0.0035, -0.0035), False, "ledger: I-662 Δfid −0.0116 predicted"),
    "458d-vs-446b": ("H1 programme on H1 (462l / 467)", (0.0067, 0.0110), True, "FINDINGS §1: Arc programme +0.0067 (ungated) on H1"),
    "471c-vs-446b": ("bulk instrument (471 prereg)", (0.0001, 0.0005), True, "ledger: mse −0.0002 raw vs predicted −0.0005 … −0.0030"),
}
rows, hits, graded = [], 0, 0
for k, (inst, rng, ship, src) in OLD.items():
    bd = B[k]["d_overall"]
    if abs(bd) < 3e-4:
        verdict = "not graded"
    else:
        graded += 1
        ok = ship if bd > 0 else (not ship)
        hits += ok
        verdict = "HIT" if ok else "MISS"
    lo, hi = rng
    rows.append(dict(row=k, board=bd, instrument=inst, predicted=[lo, hi], shipped=ship, verdict=verdict,
                     over_read=(None if bd == 0 or (lo + hi) / 2 * bd <= 0 else (lo + hi) / 2 / bd), source=src))
    print(f"{k:24s} board {bd:+.5f}  old {lo:+.4f}..{hi:+.4f} ship={str(ship):5s} {verdict:10s} "
          f"rho {'' if rows[-1]['over_read'] is None else format(rows[-1]['over_read'], '.2f')}  [{inst}]")
print(f"old instruments: {hits}/{graded} hits")
(W / "old_record.json").write_text(json.dumps(dict(rows=rows, hits=hits, graded=graded), indent=1))
