#!/usr/bin/env python3
"""690h tests: the report block cut from a patched scratch def_day.sh, run in a bash harness with a stub log.
  H1 --check/apply/second apply refused; bash -n
  H2 on the DN build's real dump (A/B/C, experiments/671-dnmag-arm/work/dn631.rc_audit.jsonl): three PASS contexts,
     `LOW_CONTEXTS none`, no dnmag_low.txt
  H3 context B's dn_frac scaled to 0.05: `LOW_CONTEXTS B`, dnmag_low.txt = B, the warning line fires
  H4 dn_moved < 200 in C: LOW C          H5 a missing dump: `LOW_CONTEXTS -1`, rc 0 (never fatal)
  H6 knob unset: nothing runs
usage: uv run python experiments/693-dnmag-speed/test690h.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DUMP = ROOT / "experiments/671-dnmag-arm/work/dn631.rc_audit.jsonl"
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def main() -> int:
    td = Path(tempfile.mkdtemp(prefix="t690h_"))
    try:
        (td / "scripts").mkdir()
        shutil.copy2(ROOT / "scripts/def_day.sh", td / "scripts/def_day.sh")
        run = lambda *a: subprocess.run([sys.executable, str(HERE / "patch690h.py"), "--root", str(td), *a],
                                        capture_output=True, text=True)
        rec(run("--check").returncode == 0, "H1 --check on the live def_day.sh")
        rec(run().returncode == 0, "H1 apply")
        rec(run().returncode == 1, "H1 second apply refused")
        rec(subprocess.run(["bash", "-n", str(td / "scripts/def_day.sh")]).returncode == 0, "H1 bash -n")
        dd = (td / "scripts/def_day.sh").read_text()
        i = dd.index('if [ -n "${VCC_DE_DNMAG:-}" ]; then   # 690h')
        j = dd.index("PY\nfi\n", i) + len("PY\nfi\n")
        blk = dd[i:j]
        wline = next(ln for ln in dd.splitlines() if ln.startswith('[ -s "$W/dnmag_low.txt" ]'))
        rows = [json.loads(l) for l in open(DUMP)]

        def case(name: str, rows_, knob: bool = True) -> tuple[subprocess.CompletedProcess, Path]:
            w = td / name
            w.mkdir()
            if rows_ is not None:
                (w / "fin.rc_audit.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows_))
            env = {k: v for k, v in os.environ.items() if not k.startswith(("VCC_", "DEF_DAY_"))}
            env.update(W=str(w), FINAL=str(w / "fin"))
            if knob:
                env["VCC_DE_DNMAG"] = "d_m200"
            # the warning line is `[ -s f ] && log ...`, rc 1 when no LOW file exists (harmless in def_day.sh, no set -e)
            sh = 'log(){ echo "LOG: $*"; }\n' + blk + wline + "\ntrue\n"
            return subprocess.run(["bash", "-c", sh], capture_output=True, text=True, env=env), w

        r, w = case("real", rows)
        rec(r.returncode == 0 and "LOW_CONTEXTS none" in r.stdout and r.stdout.count('"low": false') == 3
            and not (w / "dnmag_low.txt").exists() and (w / "dnmag_audit.json").exists(),
            "H2 the real A/B/C dump: 3 contexts pass, LOW_CONTEXTS none")
        rb = [dict(x, dn_frac=0.05) if x["context"] == "B" else x for x in rows]
        r, w = case("lowfrac", rb)
        rec(r.returncode == 0 and "LOW_CONTEXTS B" in r.stdout and (w / "dnmag_low.txt").read_text() == "B\n"
            and "WARNING DNMAG under-delivered" in r.stdout, "H3 B at dn_frac 0.05: LOW B + warning")
        rc_ = [dict(x, dn_moved=150) if x["context"] == "C" else x for x in rows]
        r, w = case("lowmoved", rc_)
        rec(r.returncode == 0 and "LOW_CONTEXTS C" in r.stdout, "H4 C at 150 genes/target: LOW C")
        r, w = case("missing", None)
        rec(r.returncode == 0 and "LOW_CONTEXTS -1" in r.stdout, "H5 missing dump: -1, never fatal")
        r, w = case("unset", rows, knob=False)
        rec(r.returncode == 0 and "dnmag" not in r.stdout, "H6 knob unset: silent")
    finally:
        shutil.rmtree(td, ignore_errors=True)
    n = sum(ok for ok, _ in RESULTS)
    print(f"\n{n}/{len(RESULTS)} PASS")
    return 0 if n == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
