#!/usr/bin/env python3
"""690k — the corrected R-KA verdict (10-16). R-KA (run690, 09:12:07Z) logged FAIL because it asked the novelty checker whether
its D/E/F-labelled package is 524dr2's board-scored one. That can never match: the board scored relabel453's A/B/C clone,
and novelty hashes `context`. This script judges the same run against the right stage and, only if every link holds, appends
a verdict + tree line to run690.log in the shape decide708.py reads (" R-KA PASS:" then " R-KA tree: def_day.sh <sha> …").
The chain, every link read from a file:
  1. SUBMITTED_LEDGER: entry dVHLws6KdMxudq0Yn5ap (524dr2, board 0.2730889) has h5ad_sha256 H (its meta's output.sha256)
  2. submissions/524dr2-october-recipe-abc.h5ad.meta.json: output.sha256 == H, relabel453.source_output_sha256 == S
  3. 658's 524dr2_final.h5ad.meta.json (the def_day output it was relabelled from): output.sha256 == S
  4. R-KA's dr_ka/r0def_final.h5ad.meta.json: output.sha256 == S        (byte-identical def_day output)
  5. def_metadiff (content-aware) R-KA final vs 524dr2 final: no real diff
  6. scripts/def_day.sh's sha now == the sha R-KA ran on (run690.log's 09:12:07Z tree line); 690i applied in it
usage: uv run python experiments/690-def-readiness-1015/rka690k.py [--write]
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
W = ROOT / "experiments/690-def-readiness-1015/work"
LOG = W / "run690.log"
ENTRY = "dVHLws6KdMxudq0Yn5ap"
META_PKG = ROOT / "submissions/524dr2-october-recipe-abc.h5ad.meta.json"
RUN_524 = ROOT / "experiments/658-prereg-1007-builds/work/dr2/run"
RUN_KA = W / "dr_ka"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    out, ok = [], True

    def chk(name, cond, detail):
        nonlocal ok
        ok &= bool(cond)
        out.append(f"  {'ok  ' if cond else 'FAIL'} {name}: {detail}")

    led = json.loads((ROOT / "submissions/SUBMITTED_LEDGER.json").read_text())
    def walk(o):
        if isinstance(o, dict):
            if o.get("entry_id") == ENTRY:
                yield o
            for v in o.values():
                yield from walk(v)
        elif isinstance(o, list):
            for v in o:
                yield from walk(v)
    rec = next(walk(led), None)
    H = (rec or {}).get("h5ad_sha256")
    chk("1 ledger entry", bool(H), f"{ENTRY} h5ad_sha256 {str(H)[:16]}")
    m = json.loads(META_PKG.read_text())
    S = (m.get("relabel453") or {}).get("source_output_sha256")
    chk("2 board package meta", (m.get("output") or {}).get("sha256") == H and bool(S),
        f"output.sha256 {str((m.get('output') or {}).get('sha256'))[:16]} == ledger; relabel453.source_output_sha256 {str(S)[:16]}")
    m524 = json.loads((RUN_524 / "524dr2_final.h5ad.meta.json").read_text())
    chk("3 524dr2's def_day output", (m524.get("output") or {}).get("sha256") == S, f"{str((m524.get('output') or {}).get('sha256'))[:16]}")
    mka = json.loads((RUN_KA / "r0def_final.h5ad.meta.json").read_text())
    ka_out = (mka.get("output") or {}).get("sha256")
    chk("4 R-KA's def_day output", ka_out == S, f"{str(ka_out)[:16]} (byte-identical iff == {str(S)[:16]})")
    md = subprocess.run([sys.executable, str(ROOT / "scripts/def_metadiff.py"), str(RUN_KA / "r0def_final.h5ad.meta.json"), str(RUN_KA),
                         str(RUN_524 / "524dr2_final.h5ad.meta.json"), str(RUN_524), "--expect", ""], capture_output=True, text=True)
    chk("5 def_metadiff R-KA vs 524dr2 finals", md.returncode == 0, md.stdout.strip().splitlines()[-1][:200] if md.stdout.strip() else md.stderr[-200:])
    lg = LOG.read_text().splitlines()
    tl = [ln for ln in lg if " R-KA tree: def_day.sh " in ln]
    ran = re.search(r"def_day\.sh ([0-9a-f]{64})", tl[-1]).group(1) if tl else None
    live = sha(ROOT / "scripts/def_day.sh")
    applied = "690i: the doses come from 438's grid" in (ROOT / "scripts/def_day.sh").read_text()
    chk("6 def_day.sh unchanged since R-KA", ran == live and applied, f"ran {str(ran)[:16]}, live {live[:16]}, 690i {'APPLIED' if applied else 'NOT applied'}")
    print("\n".join(out))
    print("ALL LINKS", "PASS" if ok else "FAIL")
    if "--write" in sys.argv:
        if not ok:
            print("not writing (a link failed)"); return 1
        ts = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        bsha = sha(ROOT / "scripts/make_submission_pm_cd4.py")[:16]
        with open(LOG, "a") as f:
            f.write(f"{ts} R-KA PASS: (corrected verdict, 690k) the Oct 22 chain rebuilt 524dr2's def_day output byte for byte: "
                    f"r0def_final output.sha256 {S[:16]} == 524dr2_final's == relabel453.source_output_sha256 of the board-scored "
                    f"package (entry {ENTRY}, 0.2730889; ledger h5ad {H[:16]}); def_metadiff diff []. The 09:12:07Z FAIL compared the "
                    f"D/E/F package with the A/B/C relabel (novelty hashes context).\n")
            f.write(f"{ts} R-KA tree: def_day.sh {live} helpers {sha(ROOT / 'scripts/def_day_helpers.py')} builder {bsha} "
                    f"at {head} (R-KA ran at cf0feef3, same def_day.sh); 690i APPLIED\n")
        print("appended the corrected verdict + tree line to", LOG)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
