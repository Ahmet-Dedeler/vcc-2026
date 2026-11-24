#!/usr/bin/env python3
"""708 decisions: apply the PREREGISTRATION (+ AMENDMENT 1) and 707 AMENDMENT 1 to the day's reads.

usage: decide708.py <tag> <dayof438 printout> <Part A design or ""> <R705: CLEARED|HARM|NONE> <R707: TRANSFERS|WEAK|NULL|HARM|NONE>
                    [<run690.log>]
Reads work/score708-<tag>.json and readiness's run690.log (default experiments/690-def-readiness-1015/work/run690.log).
A WEIGHTS move needs the day path's byte proof: the last R-KA verdict PASS, its "R-KA tree:" line quoting the CURRENT
scripts/def_day.sh sha256, and, for any dose other than CD4 x1.5 at H1 x2 (h2c1.5, 530's knob), "690i APPLIED" (readiness's
DEF_DAY_H1_STACK / widened CD4 grid). Otherwise the move is BLOCKED until 690i gets its own all-unset def_day KA.
Prints one line per switch:
  WEIGHTS  : Part A's rule result (from dayof438's RECOMMENDATION) and AMENDMENT 1's confirmation (A:<design> mean >= 0)
  RIDER    (AMENDMENT 3): 705 CLEARED -> ON (it ships in R0); 708 VETOES it iff the D/E/F-panel mean of r <= -0.0005
  SWAP707  (AMENDMENT 2, relaxed): ON iff 707's A/B/C read is not HARM AND, on the D/E/F panel, mean >= +0.0005, >= 2/3 folds
             positive and leave-top-3-out > 0; judged on its MARGINAL arm s|r when the rider is ON, else on s (AMENDMENT 3)
  Any replacement entry is gated by OCT22-CHECKLIST.md section C step 6a (245b8b1c).
Env R710 (HARM|NOT_HARM|NONE; 710's A/B/C board read of H1 x2.5): HARM blocks any WEIGHTS move that raises H1 above x2.
AMENDMENT 4 (robustness): every ACTION (WEIGHTS confirmation, RIDER veto, SWAP707 ON) needs its rule to pass on BOTH the
all-target arms ("arms") and the floor-kept arms ("arms_floor", score708 --floor-json). A rule passing in one version only
is "floor-fragile": no action. If arms_floor is missing, every action is blocked (and the FLOOR line says so).
A switch that is ON still needs its day-path pass-through, known answer and rehearsal before a final entry uses it.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def rka_proof(log_path: Path, des: str):
    """(ok, quote): the last R-KA verdict + tree line, checked against the live def_day.sh and the dose's 690i need."""
    if not log_path.exists():
        return False, f"no {log_path}"
    L = log_path.read_text().splitlines()
    vi = [i for i, l in enumerate(L) if " R-KA PASS:" in l or " R-KA FAIL:" in l]
    if not vi:
        return False, "no R-KA verdict in run690.log"
    v = L[vi[-1]]
    tree = next((l for l in L[vi[-1] + 1:] if " R-KA tree:" in l), None)
    if tree is None:
        return False, f"no R-KA tree line after: {v[:120]}"
    m = re.search(r"def_day\.sh ([0-9a-f]{64})", tree)
    live = hashlib.sha256((ROOT / "scripts/def_day.sh").read_bytes()).hexdigest()
    ok = " R-KA PASS:" in v and m is not None and m.group(1) == live
    needs690i = des not in ("", "h2c1.5")
    if needs690i:
        ok = ok and "690i APPLIED" in tree
    why = ("" if ok else " [" + ("verdict FAIL; " if " R-KA FAIL:" in v else "")
           + ("def_day.sh moved since the proof; " if (m is None or m.group(1) != live) else "")
           + ("690i not applied in the proven tree; " if needs690i and "690i APPLIED" not in tree else "") + "]")
    return ok, tree.split(" R-KA tree:", 1)[-1].strip() + why


AMEND5_TS = "2026-10-16T12:51Z"     # 708 AMENDMENT 5: the rehearsal on file must postdate it


def rehearsal_on_file(wdir: Path):
    """(ok, why): 708 AMENDMENT 5's precondition for any 708 action on a real bundle: kaA708's Part A PASS, rehearse708's
    REHEARSAL PASS after AMENDMENT 5, and rehearse708b's completed day path (day708 rc 0) after that pass."""
    import os
    if os.environ.get("R708_REHEARSAL") == "1":
        return True, "rehearsal run (synthetic bundle): the rehearsal-on-file check does not apply"
    ka = wdir / "kaA708-r708abcdef.json"
    if not ka.exists() or json.loads(ka.read_text()).get("verdict") != "PASS":
        return False, f"no Part A KA PASS in {ka.name}"
    def lines(name):
        f = wdir / name
        return f.read_text().splitlines() if f.exists() else []
    passes = [l for l in lines("rehearse708.log") if l.endswith(" REHEARSAL PASS") and l[:17] >= AMEND5_TS[:17]]
    if not passes:
        return False, f"no REHEARSAL PASS after AMENDMENT 5 ({AMEND5_TS}) in rehearse708.log"
    t_pass = passes[-1][:20]
    done = [l for l in lines("rehearse708b.log") if l[:20] > t_pass and " day708 rc 0" in l]
    if not done:
        return False, "no completed day path (day708 rc 0) after that pass in rehearse708b.log"
    return True, f"Part A KA PASS (Amendment 5); REHEARSAL PASS {t_pass}; day path rc 0 {done[-1][:20]}"


def main():
    tag, dayof, des, r705, r707 = sys.argv[1:6]
    rlog = Path(sys.argv[6]) if len(sys.argv) > 6 else ROOT / "experiments/690-def-readiness-1015/work/run690.log"
    R = json.loads((HERE / "work" / f"score708-{tag}.json").read_text())
    S, SF = R["arms"], R.get("arms_floor")          # SF: AMENDMENT 4's floor-kept version (None: every action blocked)
    import os
    reh_ok, reh_why = rehearsal_on_file(Path(os.environ.get("R708_REHEARSAL_DIR", HERE / "work")))
    if not reh_ok:                                  # AMENDMENT 5: 708 is OFF without the rehearsal on file
        SF = None
    rec = [l for l in Path(dayof).read_text().splitlines() if l.startswith("RECOMMENDATION")]
    rec = rec[-1] if rec else "RECOMMENDATION missing"
    lines = [f"REHEARSAL {'ON FILE' if reh_ok else 'MISSING -> 708 OFF: every action BLOCKED (AMENDMENT 5)'}: {reh_why}"]
    if not reh_ok:
        lines.append("FLOOR    (not consulted: 708 is OFF without the rehearsal)")
    elif SF is None:
        lines.append("FLOOR    no floor-kept version in score708's json: every 708 action is BLOCKED (AMENDMENT 4)")
    else:
        fl = R.get("floor", {})
        lines.append(f"FLOOR    F {fl.get('threshold_cpm', float('nan')):.4f} CPM from the panel's own controls; floor-kept "
                     + " / ".join(f"{k} {v}" for k, v in fl.get("kept", {}).items()) + "; an action needs both versions")
    moved = False
    if des:
        a = S.get(f"A:{des}")
        af = SF.get(f"A:{des}") if SF is not None else None
        conf_all = a is not None and a["mean"] >= 0.0
        conf_floor = af is not None and af["mean"] >= 0.0
        conf = conf_all and conf_floor
        pok, quote = rka_proof(rlog, des)
        import os
        mh = re.match(r"h([0-9.]+)c([0-9.]+)$", des)
        h1_up = bool(mh) and float(mh.group(1)) > 2.0
        r710 = os.environ.get("R710", "NONE")
        b710 = h1_up and r710 == "HARM"
        moved = conf and pok and not b710
        if h1_up:
            lines.append(f"710      A/B/C board read of H1 x2.5: {r710} -> {'BLOCKS this H1 move' if b710 else 'does not block'}")
        why = ("MOVE" if moved else ("BLOCKED by 710's HARM" if (conf and pok) else ("BLOCKED by the day-path proof" if conf
               else ("BLOCKED: floor-fragile (AMENDMENT 4), stay h2c1" if (conf_all and not conf_floor) else "BLOCKED: stay h2c1"))))
        lines.append((f"WEIGHTS  move to {des} ({rec.split(':', 1)[-1].strip()[:60]}); confirmation under the production key: "
                      f"mean {a['mean']:+.5f} (floor-kept {af['mean']:+.5f}) -> " if (a and af) else
                      (f"WEIGHTS  move to {des} ({rec.split(':', 1)[-1].strip()[:60]}); confirmation under the production key: "
                       f"mean {a['mean']:+.5f} (floor-kept: none) -> " if a else f"WEIGHTS  move to {des} but no A:{des} score -> "))
                     + why)
        lines.append(f"R-KA     {quote}")
    else:
        lines.append("WEIGHTS  stay H1 x2 / CD4 x1 (no mover passed 438's rule)")
    r, s, rs = S["r"], S["s"], S["rs"]

    def relaxed(a):
        """708 AMENDMENT 2 (b)-(d) on the D/E/F panel's folds."""
        return a["mean"] >= 0.0005 and a["n_pos"] >= 2 and a["lt3_mean"] > 0

    veto_all = r["mean"] <= -0.0005
    veto_floor = SF is not None and SF["r"]["mean"] <= -0.0005
    rider = r705 == "CLEARED" and not (veto_all and veto_floor)
    rf = f"{SF['r']['mean']:+.5f}" if SF is not None else "none"
    lines.append(f"RIDER    705={r705}: r mean {r['mean']:+.5f} ({r['n_pos']}/3), floor-kept {rf}; veto iff mean <= -0.0005 "
                 f"in both -> " + ("ON" + (" (veto floor-fragile or floor missing: no veto, AMENDMENT 4)"
                                          if (r705 == "CLEARED" and veto_all) else "")
                                   if rider else ("VETOED (OFF in the replacement entry)" if r705 == "CLEARED" else "OFF")))
    swn = "s|r" if (rider and "s|r" in S) else "s"
    sw = S[swn]
    swf = SF.get(swn) if SF is not None else None
    sw_all = relaxed(sw)
    sw_floor = swf is not None and relaxed(swf)
    swap = r707 in ("TRANSFERS", "WEAK", "NULL") and sw_all and sw_floor
    ftxt = (f"floor-kept mean {swf['mean']:+.5f} ({swf['n_pos']}/3) lt3 {swf['lt3_mean']:+.5f} relaxed {sw_floor}"
            if swf is not None else "floor-kept: none")
    lines.append(f"SWAP707  707={r707}: judged on {swn}: mean {sw['mean']:+.5f} ({sw['n_pos']}/3) lt3 {sw['lt3_mean']:+.5f} "
                 f"relaxed {sw_all}; {ftxt} -> {'ON' if swap else 'OFF'}"
                 + (" (floor-fragile, AMENDMENT 4)" if (r707 in ("TRANSFERS", "WEAK", "NULL") and sw_all and not sw_floor) else "")
                 + f"   (s alone {s['mean']:+.5f}; rs {rs['mean']:+.5f})")
    lines.append(f"FINAL    weights {des if moved else 'h2c1'} | rider {'ON' if rider else 'OFF'} | "
                 f"swap707 {'ON' if swap else 'OFF'}")
    lines.append("PROTOCOL any replacement entry: experiments/690-def-readiness-1015/OCT22-CHECKLIST.md section C step 6a (245b8b1c)")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
