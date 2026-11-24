"""216 Phase 0b — the four gates `gate216.py` scored against the wrong TAG.

`score_effects` produces every endpoint twice: `a23_*` (with production's mass fix) and
`a15_*` (without).  204 §0 and 212 §4.3 both publish the **a15** tag; 207 §4 publishes the
**a23** tag.  `gate216.py` compared 204's `q`/`q_w` and 212's macro cosine against `a23_*`
and they missed; both are exact against `a15_*`.  Nothing is re-measured here — this reads
the numbers `gate216.py` already wrote and re-scores the four gates against the tag their
sources used.  The measurement is unchanged; only which column is compared.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")


def main() -> None:
    d = json.loads((WORK / "gate216.json").read_text())
    arms = d["arms"]
    fixed = []

    def gate(name, expected, measured, tol, source):
        ok = abs(measured - expected) <= tol
        fixed.append(dict(name=name, expected=expected, measured=measured, tol=tol,
                          source=source, passed=bool(ok)))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: expected {expected} measured {measured}")
        return ok

    k = arms["k562_gwps/f010"]
    ok = True
    ok &= gate("KA-3b 204 q (a15 tag), k562_gwps fold", 0.50164311,
               round(k["a15_sgn"], 8), 5e-6, "204 §0 gate 2b")
    ok &= gate("KA-3c 204 q_w (a15 tag), k562_gwps fold", 0.50487994,
               round(k["a15_sgn_w"], 8), 5e-6, "204 §0 gate 2b")

    c10 = float(np.mean([arms[f"{h}/f010"]["a15_cos"] for h in FOLDS]))
    c34 = float(np.mean([arms[f"{h}/f034"]["a15_cos"] for h in FOLDS]))
    ok &= gate("KA-6 212 macro cos at f=0.10 (a15 tag)", 0.0218418, round(c10, 7), 5e-8,
               "212 §4.3")
    ok &= gate("KA-6b 212 macro cos at f=0.34 (a15 tag)", 0.0248622, round(c34, 7), 5e-8,
               "212 §4.3")
    ok &= gate("KA-7 212 forward dose gain %", 13.83, round(100.0 * (c34 / c10 - 1.0), 2),
               5e-3, "212 §4.3 / §6")

    d["gates_retagged"] = fixed
    d["gates_original_failures_explained"] = (
        "KA-3b/6/6b/7 in `gates` compared 204's and 212's published values against the "
        "a23 (mass-fixed) tag; both sources publish the a15 tag. Re-scored here against "
        "a15 they are exact. No measurement changed.")
    passed_now = [g for g in d["gates"]
                  if g["name"].split()[0] not in ("KA-3b", "KA-6", "KA-6b", "KA-7")]
    d["all_gates_passed_retagged"] = bool(ok and all(g["passed"] for g in passed_now))
    (WORK / "gate216.json").write_text(json.dumps(d, indent=2))
    print(json.dumps({"all_gates_passed_retagged": d["all_gates_passed_retagged"]}, indent=2))


if __name__ == "__main__":
    main()
