#!/usr/bin/env python3
"""659 tests: chain1008.sh in DRY mode on synthetic 10-07 slot logs, tail drivers, gate logs and packages (work/test/), with
the novelty command stubbed (CH1008_NOVELTY, honoured only with DRY=1). The rule itself runs for real on today's captures.
usage: uv run python experiments/659-chain1008/test659.py"""
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = HERE / "work/test"
OK, BAD = [], []


def check(name, cond, out=""):
    (OK if cond else BAD).append(name)
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        print(out)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tail_arm(tag, gate_tag, ready=True, mismatch=False, no_gate_sha=False):
    """a fake driver log + gate log + package for one tail arm; returns its tail1008.tsv row"""
    d = T / tag
    d.mkdir(parents=True, exist_ok=True)
    pkg = d / f"{tag}.vcc"
    pkg.write_bytes(os.urandom(64))
    s = sha(pkg)
    logged = ("0" * 64) if mismatch else s
    drv, gate = d / "driver.log", d / "gate.log"
    lines = ["2026-10-07T10:00:00Z start", f"2026-10-07T12:00:00Z package {pkg} sha256 {logged}"]
    lines.append(f"2026-10-07T12:00:01Z READY {tag} {logged}" if ready else "2026-10-07T12:00:01Z NOT READY bulk is NOT 498p2g's")
    drv.write_text("\n".join(lines) + "\n")
    gate.write_text(f"[{gate_tag}] PASS: exactly ['de_mag']\n[{gate_tag}] CLIFF OK {{}}\nnovelty rc=0\n"
                    + ("" if no_gate_sha else f"package sha256: {logged}\n"))
    return f"{tag}\t-\t{drv}\t{gate}\t{gate_tag}"


def run(spec_rows, s1007, hold=()):
    if T.exists():
        for x in T.glob("HOLD-*"):
            x.unlink()
    (T / "spec.tsv").write_text("# tag\trule_tag\tdriver_log\tgate_log\tgate_tag\n" + "\n".join(spec_rows) + "\n")
    for h in hold:
        (T / h).write_text("held\n")
    log = T / "dry-chain1008.log"
    if log.exists():
        log.unlink()
    env = dict(os.environ, DRY="1", CH1008_W=str(T), CH1008_S1007=str(s1007), CH1008_TAIL=str(T / "spec.tsv"),
               CH1008_NOVELTY="true")
    r = subprocess.run(["bash", str(HERE / "chain1008.sh")], cwd=ROOT, env=env, capture_output=True, text=True)
    return r.returncode, log.read_text() if log.exists() else r.stderr


def main():
    if T.exists():
        shutil.rmtree(T)
    T.mkdir(parents=True)
    s1007 = T / "s1007"
    s1007.mkdir()
    # 10-07: slot 1 submitted and captured CPH; slot 2's K20 refused its gates (not scored)
    (s1007 / "submit1007s1.log").write_text(
        "2026-10-06T20:00:05Z ARMED slot 1: 446bj2zcph pkg=x name='n' (pid 1)\n2026-10-07T00:20:00Z SUBMITTED ok\n"
        "2026-10-07T00:21:00Z CAPTURED own entry AAAAAAAAAAAAAAAAAAAA -> data/board-entries/AAAAAAAAAAAAAAAAAAAA.json\n2026-10-07T00:21:01Z END\n")
    (s1007 / "submit1007s2.log").write_text(
        "2026-10-06T20:00:06Z ARMED slot 2: 446bj2zk20 pkg=y name='m' (pid 2)\n2026-10-06T20:00:07Z NOT SUBMITTING: cliff check not OK\n")
    rows = [tail_arm("630arm", "630arm"), tail_arm("524dr2", "524dr", ready=False), tail_arm("allsw658", "allsw658")]
    rc, out = run(rows, s1007)
    check("T1 scored = the 10-07 slot that submitted (CPH), not the refused one", rc == 0 and "--scored=CPH" in out
          and "scored on 10-07 (build tags): 446bj2zcph" in out, out)
    check("T2 candidates = the rule's 10-08 list then the tail, deduplicated", "candidates (in order): 459dr allsw658 630arm 524dr2" in out, out)
    check("T3 plan: s1 459dr (rule list first), s2 allsw658 (READY driver, sha matches)",
          "would arm 10-08: s1 459dr" in out and "s2 allsw658" in out, out)
    # slot 1 held: slot 2's pick goes as slot 1
    rc, out = run(rows, s1007, hold=("HOLD-1008-s1",))
    check("T4 HOLD-1008-s1: slot 1 dropped, slot 2's pick goes as slot 1",
          "HOLD-1008-s1 present" in out and "slot 2's pick 459dr goes as slot 1" in out and "would arm 10-08: s1 459dr" in out, out)
    rc, out = run(rows, s1007, hold=("HOLD-1008-s1", "HOLD-1008-s2"))
    check("T5 both held: nothing armed", "no eligible 10-08 arm" in out and out.rstrip().endswith("END"), out)
    # a READY line whose sha does not match the package; a gate log without the sha line
    rows2 = [tail_arm("630arm", "630arm", mismatch=True), tail_arm("524dr2", "524dr", no_gate_sha=True),
             tail_arm("allsw658", "allsw658", ready=False)]
    rc, out = run(rows2, s1007)      # slot 1 takes 459dr, so slot 2 walks all three tail arms
    check("T6 sha mismatch / missing gate sha / NOT READY driver are all 'not ready'",
          "630arm: package sha != the driver's READY sha" in out and "524dr2: gate log" in out
          and "allsw658: driver log" in out, out)
    # gate tag: 524dr2 reads its gate log's "[524dr]" lines; with the wrong gate tag it is not ready
    rows3 = [tail_arm("630arm", "630arm", ready=False), tail_arm("524dr2", "524dr"), tail_arm("allsw658", "allsw658", ready=False)]
    rc, out = run(rows3, s1007)
    check("T7 524dr2 READY through its gate tag [524dr] (slot 2)", "would arm 10-08: s1 459dr" in out and "s2 524dr2" in out, out)
    bad = rows3[1].rsplit("\t", 1)[0] + "\t524dr2"
    rc, out = run([rows3[0], bad, rows3[2]], s1007)
    check("T8 the wrong gate tag makes 524dr2 not ready", "524dr2: no '[524dr2] PASS'" in out, out)
    # nothing scored on 10-07 (no logs): --scored= empty, the rule still gets the flag
    empty = T / "s1007-empty"
    empty.mkdir()
    rc, out = run(rows, empty)
    check("T9 no 10-07 logs: --scored= (empty) is still passed", "scored on 10-07 (build tags): none; --scored=" in out
          and "rule: " in out, out)
    # a TBD row reads not ready and never fails the chain
    rc, out = run([f"630arm\t630arm\tTBD\tTBD\t630arm", rows[1], rows[2]], s1007)
    check("T10 a TBD tail row reads not ready", rc == 0 and "630arm: driver log (TBD)" in out, out)
    # QUEUEFIRST: the tail (630arm, 524dr2, allsw658) before the rule's fallbacks (459dr, ALLSW)
    rows = [tail_arm("630arm", "630arm"), tail_arm("524dr2", "524dr", ready=False), tail_arm("allsw658", "allsw658")]  # fresh fixtures
    (T / "QUEUEFIRST-1008").write_text("on\n")
    rc, out = run(rows, s1007)
    (T / "QUEUEFIRST-1008").unlink()
    check("T12 QUEUEFIRST-1008: candidates 630arm 524dr2 allsw658 459dr; plan s1 630arm, s2 allsw658",
          "candidates (in order): 630arm 524dr2 allsw658 459dr" in out and "would arm 10-08: s1 630arm" in out
          and "s2 allsw658" in out, out)
    # bash -n on both scripts
    for f in ("chain1008.sh", "arm1008.sh"):
        r = subprocess.run(["bash", "-n", str(HERE / f)], capture_output=True, text=True)
        check(f"T11 bash -n {f}", r.returncode == 0, r.stderr)
    print(f"\n{len(OK)} passed, {len(BAD)} failed")
    (HERE / "work/test659.out").write_text(f"{len(OK)} passed, {len(BAD)} failed\n" + "\n".join("PASS " + x for x in OK)
                                           + "".join("\nFAIL " + x for x in BAD) + "\n")
    return 0 if not BAD else 1


if __name__ == "__main__":
    sys.exit(main())
