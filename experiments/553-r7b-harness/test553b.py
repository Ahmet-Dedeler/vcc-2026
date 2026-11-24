#!/usr/bin/env python3
"""553 — tests of the 2026-10-02 fixes from 555's review (light, seconds; no R7b truth is opened).

  B1 r7_ledger_section: a sha256 counts only inside the ledger's `## Reserve R7` section (a `###` sub-heading does not end it;
     the next `## ` does); a mention in another section does not count.
  B2 a SEALED access refuses an --owner-go claim whose sha256 is not in that section, before plan.json is written; with the
     sha256 recorded there it passes the check and stops at the next gate (a stand-in surface: fx's files, sealed=True,
     check_r7b_outputs stubbed; arms that do not exist, so the run ends at "no stats" before any phase).
  B3 a malformed --pairs entry is refused before plan.json (it used to fail in `report`, after the truth phases).
  B4 band553: an effect missing on a control half is RECORDED with band_total None and the reason (it used to be dropped
     silently); every other effect's row is unchanged.
usage: uv run python experiments/553-r7b-harness/test553b.py
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import access553 as A  # noqa: E402
import lib553 as S  # noqa: E402


def run_main(argv: list[str]) -> str:
    """access553.main() in-process; returns the SystemExit message ('' if it returned)."""
    old = sys.argv
    sys.argv = ["access553.py", *argv]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            A.main()
        return ""
    except SystemExit as e:
        return str(e)
    finally:
        sys.argv = old


def b1() -> None:
    sha = "ab" * 32
    led = ("# Ledger\n\n## Reserve R6\n- C10b: do not arm with " + "cd" * 32 + "\n\n## Reserve R7 — DECLARED\n- declared\n"
           "### R7 claim\n- recorded: sha256 `" + sha + "`\n\n## Something else\n- " + "cd" * 32 + "\n")
    sec = A.r7_ledger_section(led)
    assert sha in sec and ("cd" * 32) not in sec, sec
    assert A.r7_ledger_section("## Reserve R6\n- x\n") == ""
    tail = "x\n## Reserve R7\n- " + sha + "\n### sub\n- y\n"        # R7 last in the file
    assert A.r7_ledger_section(tail).endswith("- y\n")
    print("B1 PASS: the R7 section is bounded by level-2 headings only; other sections never count")


def b2(tmp: Path) -> None:
    claim = tmp / "claim.md"
    claim.write_text("a test claim, not R7's\n")
    sha = hashlib.sha256(claim.read_bytes()).hexdigest()
    fx = S.surface_spec("fx")
    orig_spec, orig_led, orig_chk = S.surface_spec, A.LEDGER, S.check_r7b_outputs
    S.surface_spec = lambda name: dict(fx, name=name, sealed=True)
    S.check_r7b_outputs = lambda check_truth_sha=True: {"stub": True}
    try:
        for name, text, want in (
                ("absent", "## Reserve R7\n- nothing recorded\n", "not in the ledger"),
                ("elsewhere", "## Reserve R6\n- " + sha + "\n## Reserve R7\n- nothing\n", "not in the ledger"),
                ("recorded", "## Reserve R7\n### claim\n- sha256 `" + sha + "`\n", "no stats for")):
            led = tmp / f"ledger-{name}.md"
            led.write_text(text)
            A.LEDGER = led
            out = tmp / f"out-{name}"
            msg = run_main(["fx", "--arms", "nonexistent-arm", "--seeds", "0", "--owner-go", str(claim), "--out", str(out)])
            assert want in msg, (name, msg)
            assert not (out / "plan.json").exists(), (name, "plan.json written")
            print(f"B2 {name:9s} -> refused: {msg[:110]}")
    finally:
        S.surface_spec, A.LEDGER, S.check_r7b_outputs = orig_spec, orig_led, orig_chk
    print("B2 PASS: a sealed access arms only with the claim's sha256 in the ledger's R7 section")


def b3(tmp: Path) -> None:
    for bad in ("armA", "armA:armA", "armA:nonexistent", "armA:armB:cpm"):
        out = tmp / "out-pairs"
        msg = run_main(["fx", "--arms", "armA,armB", "--seeds", "0", "--pairs", bad, "--out", str(out)])
        assert "--pairs" in msg and not (out / "plan.json").exists(), (bad, msg)
    print("B3 PASS: malformed --pairs refused before plan.json (4 cases)")


def b4(tmp: Path) -> None:
    """band553 on synthetic access jsons under work/ (tag -t553b; removed afterwards)."""
    W = HERE / "work"
    tag = "-t553b"
    base_rec = json.loads((HERE / "work/access553_fx/access553.json").read_text())
    dirs = []
    try:
        def put(name, rec):
            d = W / f"access553_{name}{tag}"
            d.mkdir(parents=True, exist_ok=False)
            dirs.append(d)
            (d / "access553.json").write_text(json.dumps(rec))
        put("fx", base_rec)
        for k in (0, 1):
            for h in "AB":
                rec = json.loads(json.dumps(base_rec))
                if (k, h) == (1, "B"):
                    rec["arms"].pop("armB")                 # armB is missing on one half
                put(f"fxc{k}{h}", rec)
        r = subprocess.run([sys.executable, str(HERE / "band553.py"), "fx", "--splits", "0,1", "--tag", tag],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-800:]
        res = json.loads((W / f"band553_fx{tag}.json").read_text())
        assert res["armB vs cpm"]["band_total"] is None and "1B" in res["armB vs cpm"]["no_band"], res["armB vs cpm"]
        assert res["armA vs cpm"]["band_total"] is not None and "no_band" not in res["armA vs cpm"], res["armA vs cpm"]
        print("B4 PASS: an effect missing on a half is recorded (band_total None, the half named); the others keep their band")
    finally:
        for d in dirs:
            shutil.rmtree(d)
        (W / f"band553_fx{tag}.json").unlink(missing_ok=True)


def main() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="test553b-"))
    try:
        b1()
        b2(tmp)
        b3(tmp)
        b4(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("ALL PASS")


if __name__ == "__main__":
    main()
