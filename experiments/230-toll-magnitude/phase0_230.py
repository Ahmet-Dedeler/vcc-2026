"""230 Phase 0 — the mandatory known-answer rows.

Six rows, every one of them a number some other experiment already recorded.  If any of
them fails to reproduce, the instrument is wrong and nothing downstream may be believed
(AGENTS.md rule 10's corollary).

  KA1  the scaling map            a8 raw reach 0.1890 -> scaled 0.124349   (BOARD-LEDGER head)
  KA2  the toll's calibration     mean(a21, a22) reach == -0.014569        (AGENTS.md 11.1 r13)
  KA3  214's renderer factor      0.002069 / 0.001124 == 1.84x             (214 section 4)
  KA4  208's rider-A composite    +0.001197, and 219's correction +0.002203/+0.003909
  KA5  the per-gene panel-total digest on the frozen challenge artifact
       (A, ABCD1) == bb6befc01caf7493                                      (191/202/208/210)
  KA6  a21 and a22 moved board `pds` and `mse` by EXACTLY zero             (162 section 3)

    uv run python phase0_230.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
ROOT = HERE.parent.parent

# ---------------------------------------------------------------------------------------
# The board record.  Every number here is transcribed from a *BOARD-RESULT / LOG entry and
# is a scored measurement, not a model output.  Sources are in the `src` field.
# ---------------------------------------------------------------------------------------
BOARD = {
    # scaled six-member deltas against the arm's own incumbent
    "a8":   dict(pds=0.0, mse=0.0, nmae=-0.000038, fid=+0.006393, reach=+0.021398,
                 jac=-0.003439, overall=+0.0040523, src="162 RESULTS s3 table"),
    "a10":  dict(pds=0.0, mse=0.0, nmae=+0.000037, fid=+0.000797, reach=-0.000754,
                 jac=-0.001140, overall=-0.0001766, src="162 RESULTS s3 table"),
    "a21":  dict(pds=-0.0, mse=-0.0, nmae=-0.064214, fid=-0.029464, reach=-0.018482,
                 jac=-0.016582, overall=-0.0214572, src="162 RESULTS s3 / LOG 09-01"),
    "a22":  dict(pds=-0.0, mse=-0.0, nmae=+0.001520, fid=-0.004447, reach=-0.010656,
                 jac=-0.001523, overall=-0.0025176, src="162 RESULTS s3 / LOG 09-01"),
    "a23":  dict(reach=-0.007933, overall=-0.000727, src="LEDGER (0.107435-0.115368); 162 s5"),
    "rc_hybrid":  dict(pds=0.0, mse=0.0, nmae=+0.00005118, fid=+0.000639,
                       reach=+0.054940, jac=+0.000366, overall=+0.0093325,
                       src="BOARD-LEDGER s3b, vs a15"),
    "rc_callset": dict(pds=0.0, mse=0.0, nmae=+0.00000025, fid=+0.000490,
                       reach=+0.022590, jac=+0.027588, overall=+0.0084448,
                       src="BOARD-LEDGER s3b, vs rc_hybrid"),
}

SCALE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532,
         "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
REACH_B = 0.0790          # the `reach` anchor, BOARD-LEDGER head


def sha16(a: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(a, dtype=np.int64).tobytes()).hexdigest()[:16]


def row(ok: bool, name: str, got, want, note: str = "") -> dict:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<46} got {got!s:<26} want {want!s:<24}{note}")
    return dict(name=name, ok=bool(ok), got=str(got), want=str(want), note=note)


def main() -> None:
    out: list[dict] = []
    print("=== 230 PHASE 0 — known-answer rows ===\n")

    # -- KA1: the scaling map ------------------------------------------------------------
    got = (0.1890 - REACH_B) * SCALE["reach"]
    out.append(row(abs(got - 0.124349) < 2e-4, "KA1 scaling map, a8 raw reach -> scaled",
                   f"{got:.6f}", "0.124349", " (BOARD-LEDGER's own check)"))

    # -- KA2: the toll IS mean(a21, a22) -------------------------------------------------
    toll = 0.5 * (BOARD["a21"]["reach"] + BOARD["a22"]["reach"])
    out.append(row(abs(toll - (-0.014569)) < 1e-6, "KA2 toll == mean(a21, a22) reach",
                   f"{toll:.6f}", "-0.014569", " (AGENTS.md 11.1 rule 13 / 162 s5)"))
    out.append(row(abs(toll / 6.0 - (-0.002428)) < 1e-6, "KA2b toll on overall (toll/6)",
                   f"{toll / 6.0:.6f}", "-0.002428", ""))

    # -- KA3: 214's renderer factor ------------------------------------------------------
    f = 0.002069 / 0.001124
    out.append(row(abs(f - 1.84) < 0.01, "KA3 214 shipped-renderer factor",
                   f"{f:.4f}", "1.84", " (Hamilton PLOLO d_pds, f=0.10)"))

    # -- KA4: 208's rider-A composite and 219's correction -------------------------------
    # Read 208's OWN measured legs off its JSON rather than re-deriving them (rule 11:
    # a path, not a description).  219's correction is a multiplicative 1.8407x on the
    # `pds` leg only; every other leg is unchanged.
    p208 = json.loads((ROOT / "experiments/208-free-riders/work/price208.json").read_text())
    ham = p208["rows"]["ham_tier"]
    c208 = p208["gate"]["ham_tier"][
        "none (196's reading; 169 measured +0.0006 for profile arms)"]["priced"]
    RFAC = 0.002069 / 0.001124                          # 214 s4, shipped vs a23 renderer
    c219 = c208 * RFAC
    six208, six219 = ham["sum"], ham["sum"] + ham["pds"] * (RFAC - 1.0)
    out.append(row(abs(c208 - 0.001197) < 1e-6, "KA4 208 rule-13 composite (208 renderer)",
                   f"{c208:.6f}", "+0.001197", " = pds_overall*1.011 + mse leg"))
    out.append(row(abs(c219 - 0.002203) < 1e-5, "KA4b 219 composite (shipped renderer)",
                   f"{c219:.6f}", "+0.002203", " CLEARS the +0.0019 floor at 1.16x"))
    out.append(row(abs(six208 - 0.002899) < 1e-5, "KA4c 208 six-member sum, as measured",
                   f"{six208:.6f}", "+0.002899", " (208 s6, `reach` leg at 1.290)"))
    out.append(row(abs(six219 - 0.003909) < 1e-5, "KA4d 219 six-member, corrected",
                   f"{six219:.6f}", "+0.003909", " (only the `pds` leg is rescaled)"))
    out.append(row(abs((c219 - 0.00243) - (-0.000227)) < 1e-5,
                   "KA4e composite with the -0.00243 toll charged",
                   f"{c219 - 0.00243:.6f}", "-0.000227", " FAILS -> the arm dies on the toll"))

    # -- KA5: the per-gene panel-total digest on the frozen challenge artifact -----------
    z = np.load(ROOT / "experiments/196-emitted-depth/work/chal-totals-a15.npz",
                allow_pickle=True)
    blocks = np.asarray(z["blocks"])
    totals = np.asarray(z["totals"])
    i = int(np.flatnonzero(blocks == "A|ABCD1")[0])
    d = sha16(totals[i])
    out.append(row(d == "bb6befc01caf7493", "KA5 digest (A, ABCD1) on the frozen totals",
                   d, "bb6befc01caf7493", f"  {totals.shape[0]} blocks x {totals.shape[1]} genes"))

    # -- KA6: a21/a22 moved pds and mse by EXACTLY zero -----------------------------------
    #   This is the row the whole experiment turns on, so it is a known-answer row.
    z6 = all(BOARD[a][m] == 0.0 for a in ("a21", "a22", "rc_hybrid", "rc_callset")
             for m in ("pds", "mse"))
    out.append(row(z6, "KA6 a21/a22/rc_* board d_pds == d_mse == 0 exactly",
                   "0.0 x 8", "0.0 x 8", " (LOG 09-01: 'to 16 significant digits')"))

    ok = all(r["ok"] for r in out)
    (WORK / "phase0_230.json").write_text(json.dumps(
        dict(rows=out, all_pass=ok, board=BOARD, scale=SCALE), indent=2))
    print(f"\n=== {'ALL PASS' if ok else 'FAILURES PRESENT — STOP'} ===")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
