#!/usr/bin/env python3
"""KA630 audit -- P2 of the coordinator's "10-06 APPLY DAY v2" ([4d0207], with [d4cd5f]; 2026-10-02).

A run_ka build never runs def_day.sh, so 630dd's FAILSAFE line is never written for it. This tool runs 630dd's OWN report
snippet on the build's rc_audit dump. The snippet is extracted from `patch630dd.C_NEW` at run time, so it is byte-for-byte
the code def_day runs and cannot drift. The tool prints the snippet's lines exactly as def_day logs them ("630dd: ...").

Exit 0 ONLY if BOTH hold (KA PASS):
- the line `630dd: uprn FAILSAFE_TOTAL 0` is present. n > 0 means a fail-safe fired; n = -1 means the dump is missing or
  unreadable, or has no ur_ rows; a missing line also fails;
- every expected context (`--contexts`, default A,B,C) is in the audit and clears 630's prereg capacity floor: mean
  ur_frac_got >= 0.127 AND mean ur_n_moved >= 500 per target ([d4cd5f]: a pre-slot condition; reference 446bj2z v3 blocks
  A 0.284 / 1,550, B 0.295 / 970, C 0.325 / 1,985, about 2x margin). The floor is the snippet's own `capacity_floor_met`
  from its `uprn <ctx>: {...}` lines of this run, the same thresholds def_day reports. A miss prints `FLOOR FAIL <ctx>`.
Exit 1 otherwise.

The KA build must write the dump: run it with VCC_RC_AUDIT_DUMP=<path>.jsonl. The builder APPENDS, so remove a stale file
first. It also needs VCC_DE_UPRN set (otherwise there are no ur_ rows, and the audit reports -1).

usage: uv run python experiments/630-rank-neutral-up/ka630_audit.py <build>.rc_audit.jsonl [--out-dir DIR]
       (DIR, default: next to the dump, receives uprn_audit.json and, if any context is effectively off,
        uprn_effectively_off.txt, exactly as def_day writes them into its work dir)
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import patch630dd as P  # noqa: E402

START = "<<'PY' 2>&1 | while IFS= read -r _l; do log \"630dd: $_l\"; done || true\n"


def report_snippet() -> str:
    """The python heredoc def_day runs for 630dd's report, cut from the patch's own C_NEW text."""
    c = P.C_NEW
    i = c.index(START) + len(START)
    j = c.index("\nPY\n", i)
    return c[i:j + 1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dump", type=Path, help="the build's rc_audit jsonl (VCC_RC_AUDIT_DUMP)")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--contexts", default="A,B,C", help="contexts that must be present and clear the capacity floor")
    a = ap.parse_args()
    out = a.out_dir or a.dump.parent
    out.mkdir(parents=True, exist_ok=True)
    # def_day runs the snippet with the system python3 (stdlib only), stdout+stderr merged, each line logged "630dd: <line>"
    r = subprocess.run(["python3", "-", str(a.dump), str(out / "uprn_audit.json"), str(out / "uprn_effectively_off.txt")],
                       input=report_snippet(), capture_output=True, text=True)
    lines = [f"630dd: {ln}" for ln in (r.stdout + r.stderr).splitlines()]
    for ln in lines:
        print(ln)
    tot = [m.group(1) for ln in lines if (m := re.fullmatch(r"630dd: uprn FAILSAFE_TOTAL (-?\d+)", ln))]
    ok = len(tot) == 1 and tot[0] == "0"
    # capacity floor per expected context, parsed from the snippet's `uprn <c>: {json}` lines of THIS run (never from a
    # uprn_audit.json that an earlier run may have left in out_dir)
    rep = {m.group(1): json.loads(m.group(2)) for ln in lines if (m := re.fullmatch(r"630dd: uprn (\S+): (\{.*\})", ln))}
    for c in [x for x in a.contexts.split(",") if x]:
        v = rep.get(c)
        if v is None:
            print(f"FLOOR FAIL {c} (absent from the audit)")
            ok = False
            continue
        print(f"floor {c}: mean_frac_got {v['mean_frac_got']:.4f} (>= 0.127) mean_n_moved {v['mean_n_moved']:.1f} (>= 500) "
              f"targets {v['targets']} -> {'ok' if v['capacity_floor_met'] else 'MISS'}")
        if not v["capacity_floor_met"]:
            print(f"FLOOR FAIL {c}")
            ok = False
    print(f"KA630 audit: {'PASS' if ok else 'FAIL'} (FAILSAFE_TOTAL {tot[0] if tot else 'line missing'}; contexts "
          f"{a.contexts}; module {P.UPRN_MODULE_SHA256[:12]}; dump {a.dump})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
