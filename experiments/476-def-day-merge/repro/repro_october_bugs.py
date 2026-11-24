#!/usr/bin/env python3
"""Reproductions for the October day-of review (experiment 485).

Run from the repo root:
  uv run python experiments/485-cloud-october-review/repro/repro_october_bugs.py

Prints one block per finding. Exit 0 means every asserted reproduction behaved
as this file claims. Does not touch production scripts or data files.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))

os.chdir(ROOT)


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def run(cmd, **kw):
    return subprocess.run(cmd, text=True, capture_output=True, **kw)


def write_ntc(path: Path, rhos: dict[str, float], scales: dict[str, float]) -> None:
    per = {
        c: {"rho": rhos[c], "scale": scales[c], "B": 1.0, "claim": 1.0}
        for c in rhos
    }
    path.write_text(json.dumps({"per_context": per}))


def main() -> int:
    failures = []

    # ------------------------------------------------------------------
    section("F1 setup: bash is GNU 5.x here; production is bash 3.2. Emulation notes inline.")
    ver = run(["bash", "--version"])
    print(ver.stdout.splitlines()[0])
    print("Production traps from AGENTS.md §6.2 were checked by reading def_day.sh")
    print("(no associative arrays, no '$' before '\\|' in grep). See F-negatives.")

    # ------------------------------------------------------------------
    section("F1: level-build exit code is ignored; a pre-existing h5ad is success")
    # Mirrors scripts/def_day.sh build() lines 113-116 and the driver's
    # scripts/def_pipeline.py:696-697 (refuse to overwrite, then the shell
    # only tests [ -s h5ad ]).
    td = Path(tempfile.mkdtemp(prefix="f1-level-"))
    try:
        h5 = td / "arm-lvl1.h5ad"
        log = td / "driver-lvl1.log"
        # Previous attempt wrote a COMPLETE artifact at q=1.75, factors 1.2.
        h5.write_text("BUILT q=1.75 factors=D=1.2000,E=1.2000,F=1.2000\n")
        stub = td / "driver_stub.sh"
        stub.write_text(
            "#!/bin/bash\n"
            f"echo 'REFUSING: {h5} exists; the builder refuses to overwrite' > \"$1\"\n"
            "exit 1\n"
        )
        stub.chmod(0o755)
        script = td / "retry.sh"
        # The requested factors are the NEW run's (lower q). The file on disk is the old build.
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
fact='D=1.0500,E=1.0500,F=1.0500'
"{stub}" "{log}"
rc=$?
if [ -s "{h5}" ]; then
  echo "SHELL_ACCEPTS rc=$rc file=$(cat "{h5}" | tr '\\n' ' ') requested=$fact"
else
  echo "SHELL_DIES rc=$rc"
fi
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script)])
        print(out.stdout.strip())
        print("driver log:", log.read_text().strip())
        if "SHELL_ACCEPTS rc=1" not in out.stdout or "q=1.75" not in out.stdout:
            failures.append("F1 level stale accept")
        if "requested=D=1.0500" not in out.stdout:
            failures.append("F1 did not show the mis-labelled factor")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F1b: a build the arm guard refuses is still shipped (file exists, rc ignored)")
    td = Path(tempfile.mkdtemp(prefix="f1-arms-"))
    try:
        h5 = td / "arm-lvl0.h5ad"
        h5.write_bytes(b"complete-h5ad")
        script = td / "accept.sh"
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
# driver wrote the h5ad, then assert_build_carries_arms raised SystemExit
( echo 'REFUSING THE BUILD: artifact does not carry the configuration' > "{td}/driver.log"; exit 1 )
if [ -s "{h5}" ]; then echo SHIP; else echo DIE; fi
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script)])
        print(out.stdout.strip(), "driver:", (td / "driver.log").read_text().strip())
        if out.stdout.strip() != "SHIP":
            failures.append("F1b arm-guard bypass")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F1c: dump stage always crashes the arm assert (no meta); shell keeps any non-empty npz")
    from def_pipeline import assert_build_carries_arms  # noqa: E402

    missing = Path(tempfile.mkdtemp(prefix="f1-dump-")) / "no-such.h5ad.meta.json"
    try:
        assert_build_carries_arms(missing, {"VCC_H1_UNION_G": "0.56"})
        print("UNEXPECTED: assert returned")
        failures.append("F1c assert did not raise")
    except FileNotFoundError as e:
        print(f"assert_build_carries_arms on a dump (no meta) raises {type(e).__name__}: {e}")
    print(
        "Builder path: make_submission_pm_cd4.py writes the npz then `raise SystemExit(0)` "
        "BEFORE write_submission (the only meta writer). r.run() therefore returns 0, "
        "then assert_build_carries_arms runs on a path that was never created."
    )
    # Stale npz: new dump stage fails BEFORE savez; old npz remains; shell only tests -s.
    td = Path(tempfile.mkdtemp(prefix="f1-dumpstale-"))
    try:
        dump = td / "dump_h1x2_cd4x1.npz"
        dump.write_bytes(b"OLD_DUMP_LODO_OFF")
        script = td / "dump.sh"
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
# memory gate / early driver failure: npz not touched, exit 1
( echo 'free+inactive memory 10.0 GiB < 23 GiB; refusing' > "{td}/driver-dump.log"; exit 1 )
if [ -s "{dump}" ]; then echo "DUMP_KEPT bytes=$(wc -c < "{dump}") content=$(cat "{dump}")"; else echo DUMP_DIE; fi
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script)])
        print(out.stdout.strip())
        if "OLD_DUMP_LODO_OFF" not in out.stdout:
            failures.append("F1c stale dump")
    finally:
        shutil.rmtree(td)
        missing.parent.rmdir()

    # ------------------------------------------------------------------
    section("F2: out-of-bound next factor dies inside the loop; the safe-build fallback never runs")
    # rho* and the first-step rule: f = 1 * ntclevel scale. scale > 1.4 when un-levelled rho < ~0.50.
    import scipy.stats

    q = 0.31 / (0.31 + 1.78)
    rho_star = math.exp(2.0 * 0.01 * float(scipy.stats.norm.ppf(q)))
    print(f"rho* = {rho_star:.6f} (shell fallback hardcodes 0.9793; |delta|={abs(rho_star-0.9793):.3e})")
    # Known-answer from 452's def_day.log line 29: extrapolate (1.1948, 0.9609)-(1.0000, 0.6861) -> 1.2069
    td = Path(tempfile.mkdtemp(prefix="f2-relevel-"))
    try:
        j0 = td / "lvl0.ntclevel.json"
        j1 = td / "lvl1.ntclevel.json"
        write_ntc(j0, {"D": 0.6861, "E": 0.6206, "F": 0.5911}, {"D": 1.1948, "E": 1.2562, "F": 1.2871})
        write_ntc(j1, {"D": 0.9609, "E": 0.9566, "F": 1.0568}, {"D": 1.0096, "E": 1.0118, "F": 0.9627})
        ka = run([
            "uv", "run", "python", "experiments/448-def-rehearsal-q175/relevel_iter.py",
            f"{j0}:D=1,E=1,F=1",
            f"{j1}:D=1.1948,E=1.2562,F=1.2871",
        ])
        print("452 known-answer (want D next 1.2069):")
        print(ka.stdout)
        if "next f 1.2069" not in ka.stdout:
            failures.append("F2 known-answer did not reproduce 452's 1.2069")
        # Flat exponent: rho barely moves, secant explodes past 1.4.
        jf = td / "flat0.ntclevel.json"
        jg = td / "flat1.ntclevel.json"
        write_ntc(jf, {"D": 0.80, "E": 0.80, "F": 0.80}, {"D": 1.1, "E": 1.1, "F": 1.1})
        write_ntc(jg, {"D": 0.801, "E": 0.801, "F": 0.801}, {"D": 1.0, "E": 1.0, "F": 1.0})
        flat = run([
            "uv", "run", "python", "experiments/448-def-rehearsal-q175/relevel_iter.py",
            f"{jf}:D=1.0000,E=1.0000,F=1.0000",
            f"{jg}:D=1.1500,E=1.1500,F=1.1500",
        ])
        print("flat-rho secant:")
        print(flat.stdout or flat.stderr)
        # First step from a low un-levelled rho. scale stored as ntclevel would: sqrt(rho*/rho) roughly,
        # here set explicitly to 1.55 so the printed factor is 1.5500.
        jlow = td / "low.ntclevel.json"
        write_ntc(jlow, {"D": 0.45, "E": 0.70, "F": 0.70}, {"D": 1.55, "E": 1.20, "F": 1.20})
        low = run([
            "uv", "run", "python", "experiments/448-def-rehearsal-q175/relevel_iter.py",
            f"{jlow}:D=1,E=1,F=1",
        ])
        print("low-rho first step:")
        print(low.stdout)
        script = td / "loop.sh"
        script.write_text(
            r"""#!/bin/bash
set -uo pipefail
die(){ echo "STOP: $*"; exit 1; }
# fallback that WOULD ship lvl0, copied from def_day.sh:146-163, simplified
ship_fallback(){ echo "FALLBACK_SHIPPED lvl0"; }
OUT="$1"
NEXT="$(echo "$OUT" | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | cut -d= -f2-)"
[ -n "$NEXT" ] || die "relevel_iter printed no factors"
uv run python scripts/def_day_helpers.py inbound "$NEXT" \
  || die "next factors $NEXT leave the builder bound [0.7, 1.4]: lower q and rerun"
ship_fallback
"""
        )
        script.chmod(0o755)
        loop = run(["bash", str(script), low.stdout])
        print("shell on that proposal:", loop.stdout.strip(), "rc", loop.returncode)
        if "FALLBACK_SHIPPED" in loop.stdout:
            failures.append("F2 fallback ran despite out-of-bound factors")
        if loop.returncode == 0 or "STOP:" not in loop.stdout:
            failures.append("F2 shell did not die")
        # Identical-f division by zero (retry of the same factors, rho moved).
        jz = td / "z0.ntclevel.json"
        jy = td / "z1.ntclevel.json"
        write_ntc(jz, {"D": 0.90, "E": 0.90, "F": 0.90}, {"D": 1.0, "E": 1.0, "F": 1.0})
        write_ntc(jy, {"D": 0.91, "E": 0.91, "F": 0.91}, {"D": 1.0, "E": 1.0, "F": 1.0})
        zdiv = run([
            "uv", "run", "python", "experiments/448-def-rehearsal-q175/relevel_iter.py",
            f"{jz}:D=1.2000,E=1.2000,F=1.2000",
            f"{jy}:D=1.2000,E=1.2000,F=1.2000",
        ])
        print("same-f different-rho rc", zdiv.returncode)
        err = (zdiv.stderr or "") + (zdiv.stdout or "")
        print(err.strip().splitlines()[-3:])
        if zdiv.returncode == 0 or "VCC_BULK_CTX_SCALE=" in (zdiv.stdout or ""):
            failures.append("F2 same-f did not crash")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F3: lvlC is accepted without reading the new build's rho")
    td = Path(tempfile.mkdtemp(prefix="f3-lvlc-"))
    try:
        script = td / "lvlc.sh"
        script.write_text(
            r"""#!/bin/bash
set -uo pipefail
# Condensed from def_day.sh:131-140. build() deletes the h5ad only when past_cliff
# (rho>1). A combination that lands at rho 0.999 is kept and labelled CONVERGED.
build(){
  # simulated measurement: rho D=0.999 E=0.999 F=0.999, all <= 1, all outside the 1% band of 0.9793
  echo "measured rho 0.999 0.999 0.999" > "$1.meas"
  echo present > "$1.h5ad"
}
FINAL=""
NEXT="D=1.3000,E=1.3000,F=1.3000"
# no historical build carries NEXT
build "$2/lvlC"
if [ -s "$2/lvlC.h5ad" ]; then FINAL="$2/lvlC"; fi
if [ -n "$FINAL" ]; then echo "CONVERGED_LABEL ${FINAL##*/} $(cat "$2/lvlC.meas")"; else echo "NOT_ACCEPTED"; fi
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script), "x", str(td)])
        print(out.stdout.strip())
        if "CONVERGED_LABEL" not in out.stdout or "0.999" not in out.stdout:
            failures.append("F3 lvlC gate")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F4: past_cliff exit 1 is both 'safe' and 'crashed'")
    td = Path(tempfile.mkdtemp(prefix="f4-cliff-"))
    try:
        good = td / "ok.ntclevel.json"
        bad = td / "bad.ntclevel.json"
        missing = td / "nope.ntclevel.json"
        write_ntc(good, {"D": 0.97, "E": 0.97, "F": 0.97}, {"D": 1, "E": 1, "F": 1})
        write_ntc(bad, {"D": 1.05, "E": 0.97, "F": 0.97}, {"D": 1, "E": 1, "F": 1})
        # crash case: no per_context
        missing.write_text("{}\n")
        def cliff(p):
            r = run(["uv", "run", "python", "scripts/def_day_helpers.py", "past_cliff", str(p)])
            return r.returncode
        c_ok, c_bad, c_crash = cliff(good), cliff(bad), cliff(missing)
        print(f"exit codes: safe={c_ok} past={c_bad} crashed={c_crash}")
        script = td / "if.sh"
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
decide(){{
  if uv run python scripts/def_day_helpers.py past_cliff "$1"; then echo "DELETE $2"; else echo "KEEP $2"; fi
}}
decide "{good}" safe
decide "{bad}" past
decide "{missing}" crashed
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script)])
        print(out.stdout.strip())
        if "KEEP crashed" not in out.stdout or "DELETE past" not in out.stdout or "KEEP safe" not in out.stdout:
            failures.append("F4 past_cliff polarity")
        if c_ok != 1 or c_bad != 0 or c_crash != 1:
            failures.append(f"F4 exit codes {c_ok, c_bad, c_crash}")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F5: zero-gene audit and prep dry-run do not stop the script")
    td = Path(tempfile.mkdtemp(prefix="f5-audit-"))
    try:
        # audit-shaped program: summary line, then detail, then a json sidecar.
        fake = td / "audit_fake.py"
        fake.write_text(
            "import pathlib, sys\n"
            "print('blocks 900; gated genes with ZERO block total (target excluded): 24 in 10 blocks; lfc < -5: 24')\n"
            "print('detail ' * 400)\n"
            "pathlib.Path(sys.argv[1]).write_text('json-written')\n"
            "raise SystemExit(0)\n"
        )
        side = td / "side.json"
        script = td / "audit.sh"
        # Exact pipeline from def_day.sh:173 and :181-182.
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
uv run python "{fake}" "{side}" 2>&1 | head -1 | tee "{td}/log.txt"
echo "pipeline_rc=${{PIPESTATUS[0]}}"
if [ -f "{side}" ]; then echo "JSON_WRITTEN"; else echo "JSON_MISSING"; fi
# prep failure is a warning, then DONE
( echo 'prep failed: Unknown context' > "{td}/prep.log"; exit 1 ) \
  && echo "prep ok" || echo "WARN prep dry run failed"
echo DONE
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script)])
        print(out.stdout.strip())
        if "DONE" not in out.stdout or "WARN prep" not in out.stdout:
            failures.append("F5 script stopped")
        if "24 in 10 blocks" not in out.stdout:
            failures.append("F5 summary missing")
        # The zero count is visible to a human and invisible to the script: DONE still printed.
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F6: nctx failure is 'same-line off', not a stop")
    td = Path(tempfile.mkdtemp(prefix="f6-nctx-"))
    try:
        bad = td / "sameline.npz"
        bad.write_bytes(b"not a zip")
        r = run(["uv", "run", "python", "scripts/def_day_helpers.py", "nctx", str(bad)])
        print("nctx rc", r.returncode, "stdout", repr(r.stdout), "stderr_tail", (r.stderr or "")[-200:])
        script = td / "sl.sh"
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
NCTX="$(uv run python scripts/def_day_helpers.py nctx "{bad}" 2>/dev/null)"
if [ "$NCTX" -gt 0 ]; then echo "MODE_ON nctx=$NCTX"; else echo "MODE_OFF nctx=$(printf %q "$NCTX")"; fi
echo "continued"
"""
        )
        script.chmod(0o755)
        out = run(["bash", str(script)])
        print(out.stdout.strip())
        if "MODE_OFF" not in out.stdout or "continued" not in out.stdout:
            failures.append("F6 nctx fail-open")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F7: EXIT trap does not fail the script when restore leaves the wrong builder")
    td = Path(tempfile.mkdtemp(prefix="f7-trap-"))
    try:
        script = td / "trap.sh"
        script.write_text(
            r"""#!/bin/bash
set -uo pipefail
BSHA0=aaaaaaaaaaaaaaaa
restore(){
  echo "restore rc=$1" >> "$2"
  # last command succeeds, hashes differ, no comparison
  echo "builder after restore: bbbbbbbbbbbbbbbb (want $BSHA0)"
}
trap 'restore fail '"$td"'/restore.txt' EXIT
exit 0
"""
        )
        # The trap above doesn't capture the incoming $? because we pass a literal.
        # Faithful version: capture $? first, then a failing restore, then a successful log, no hash check.
        script.write_text(
            f"""#!/bin/bash
set -uo pipefail
BSHA0=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
restore(){{
  local why=$?
  echo "incoming_exit=$why" >> "{td}/trace"
  false
  echo "restore_command_rc=$?" >> "{td}/trace"
  echo "builder after restore: bbbb (want ${{BSHA0:0:4}})"
}}
trap restore EXIT
exit 0
"""
        )
        out = run(["bash", str(script)])
        print("script rc", out.returncode)
        print("stdout:", out.stdout.strip())
        print("trace:", (td / "trace").read_text().strip())
        if out.returncode != 0:
            failures.append("F7 trap unexpectedly failed the script")
        if "bbbb" not in out.stdout or "aaaa" not in out.stdout:
            failures.append("F7 log line missing")
    finally:
        shutil.rmtree(td)

    # ------------------------------------------------------------------
    section("F8: CD4 gate has no pass-through; preflight does not require the names D,E,F")
    text = (ROOT / "scripts/def_day.sh").read_text()
    drv = [ln for ln in text.splitlines() if "def_pipeline.py" in ln and "DRV=" in ln or "min-free-gib" in ln]
    print("DRV lines:")
    for ln in text.splitlines():
        if "def_pipeline.py" in ln or "cd4-thresholds" in ln or "FACT=" in ln or "contexts D,E,F" in ln:
            print(" ", ln.strip())
    if "--cd4-thresholds" in text:
        failures.append("F8 unexpectedly found a cd4-thresholds pass-through")
    pre = (ROOT / "scripts/def_preflight.py").read_text()
    if 'len(contexts) == 3' not in pre:
        failures.append("F8 preflight context-count check moved")
    if '"D"' in pre and "contexts ==" in pre:
        print("preflight mentions D by name (inspect)")
    else:
        print("preflight gates on len(contexts)==3, not on the labels D,E,F.")
    print("extract_cd4.py defaults: --min-union-targets 290, --min-targets 275 (raises ValueError).")
    print("def_day.sh never passes --cd4-thresholds relaxed, so a low-coverage panel stops at setup.")

    # ------------------------------------------------------------------
    section("F-negatives")
    # Environment isolation of the dump subshell.
    script = r"""#!/bin/bash
set -uo pipefail
unset VCC_H1_UNION_G VCC_PERT_SCALE VCC_BULK_DUMP
( export VCC_H1_UNION_G=0.56 VCC_BULK_DUMP=/tmp/d VCC_PERT_SCALE=lodo_lo0.5_hi3_def
  echo "inside G=$VCC_H1_UNION_G scale=$VCC_PERT_SCALE" )
echo "outside G=${VCC_H1_UNION_G-<unset>} scale=${VCC_PERT_SCALE-<unset>} dump=${VCC_BULK_DUMP-<unset>}"
"""
    out = run(["bash", "-c", script])
    print(out.stdout.strip())
    if "outside G=<unset>" not in out.stdout:
        failures.append("negative: dump subshell leaked")
    # bash 3.2 empty-array note, demonstrated on THIS bash (5.2 does not throw).
    empty = run(["bash", "-c", "set -u; a=(); echo \"n=${#a[@]}\"; echo \"${a[@]}\"; echo ok"])
    print("empty array under set -u on bash 5.2 rc", empty.returncode, "stdout", repr(empty.stdout))
    print("bash 3.2 (not installed here) errors on \"${a[@]}\" when a is empty and set -u.")
    print("def_day.sh only expands POINTS after build() appends lvl0, so the 3.2 error is not reached.")
    # grep patterns
    day = (ROOT / "scripts/def_day.sh").read_text()
    if "\\|" in day:
        print("WARNING: def_day.sh contains \\| which BSD grep mishandles before $")
        failures.append("negative grep")
    else:
        print("def_day.sh has no \\| grep alternation (the 2026-09-25 BSD grep bug is not in this file).")
    # .4f round-trip vs same() tolerance
    bad_rt = []
    for i in range(7000, 14001):
        s = f"{i/10000:.4f}"
        if f"{float(s):.4f}" != s:
            bad_rt.append(s)
            break
    print(f".4f round-trip over [0.7000, 1.4000]: {'FAIL '+bad_rt[0] if bad_rt else 'all hold'}")
    if bad_rt:
        failures.append("round-trip")
    # EXIT status preserved when the trap's last command succeeds (bash 5.2).
    trap = run(["bash", "-c", "set -uo pipefail; trap 'echo trapped' EXIT; exit 7"])
    print(f"exit 7 with a successful EXIT trap: rc={trap.returncode} (5.2 preserves it; trap does not mask die)")
    if trap.returncode != 7:
        failures.append("trap masked exit on this bash")

    # ------------------------------------------------------------------
    section("SUMMARY")
    if failures:
        print("FAILURES:", failures)
        return 1
    print("all reproductions matched the claims")
    return 0


if __name__ == "__main__":
    sys.exit(main())
