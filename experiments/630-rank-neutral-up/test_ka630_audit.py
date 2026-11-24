#!/usr/bin/env python3
"""Tests for ka630_audit.py (P2; [4d0207]). Light, no data.
  A1  the snippet the tool runs appears VERBATIM in the composed def_day.sh (606's chain, scratch tree): same code
  A2  exit 0 + `630dd: uprn FAILSAFE_TOTAL 0` on a clean dump
  A3  exit 1 + FAILSAFE_TOTAL n > 0 when a guard / invariant fail-safe fired
  A4  exit 1 + FAILSAFE_TOTAL -1 for a dump with no ur_ rows, an unreadable dump, a missing dump
  A5  exit 1 + `FLOOR FAIL <ctx>` when a context misses the capacity floor (frac_got < 0.127 or moved < 500) or is absent;
      exit 0 exactly at the floor (0.127 / 500: the prereg's ">="); a passing uprn_audit.json left in out_dir is never read
usage: uv run python experiments/630-rank-neutral-up/test_ka630_audit.py"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
sys.path.insert(0, str(ROOT / "experiments/606-applyday-composition"))
import ka630_audit as K  # noqa: E402
import test598 as T598  # noqa: E402
import compose606 as C606  # noqa: E402

RES: list[tuple[bool, str]] = []


def rec(ok, name):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def row(c, t, **kw):
    base = dict(context=c, target=t, ur_n_sel=10.0, ur_frac_got=0.3, ur_n_moved=900.0, ur_guard_failed=0.0,
                ur_invariant_failed=0.0, ur_n_no_pos_donor=0.0)
    return json.dumps({**base, **kw})


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="tka630-") as tds:
        td = Path(tds)
        t = td / "tree"
        T598.canon_tree(t)
        (t / "data/identity_refs").mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", t / "data/identity_refs/hek293_nonflex_ref642.npz")
        (t / "src/vcc2026").mkdir(parents=True, exist_ok=True)
        for f in (ROOT / "src/vcc2026").glob("*.py"):
            shutil.copy2(f, t / "src/vcc2026" / f.name)
        for name in C606.DEFAULT_ORDER:
            c = subprocess.run([sys.executable, str(ROOT / C606.PATCHES[name][0]), "--root", str(t), "--check"], capture_output=True, text=True)
            co = (c.stdout + c.stderr).lower()
            if c.returncode != 0 and "already" in co and "partially" not in co:
                continue
            r = subprocess.run([sys.executable, str(ROOT / C606.PATCHES[name][0]), "--root", str(t)], capture_output=True, text=True)
            assert r.returncode == 0, name
        day = (t / "scripts/def_day.sh").read_text()
        snip = K.report_snippet()
        rec(day.count(snip) == 1 and "FAILSAFE_TOTAL" in snip, f"A1 the audited snippet is verbatim in the composed def_day.sh ({len(snip)} chars)")
        cases = {
            "clean": ("\n".join([row("A", "T1"), row("B", "T2"), row("C", "T3")]) + "\n", 0, "630dd: uprn FAILSAFE_TOTAL 0"),
            "floor_miss": ("\n".join([row("A", "T1"), row("B", "T2"), row("C", "T3", ur_frac_got=0.10)]) + "\n", 1, "FLOOR FAIL C"),
            "moved_miss": ("\n".join([row("A", "T1"), row("B", "T2", ur_n_moved=400.0), row("C", "T3")]) + "\n", 1, "FLOOR FAIL B"),
            "ctx_absent": ("\n".join([row("A", "T1"), row("B", "T2")]) + "\n", 1, "FLOOR FAIL C (absent from the audit)"),
            # the dump is missing, so the snippet writes no uprn_audit.json; the passing one an earlier run left must not be read
            "stale_json": (None, 1, "FLOOR FAIL A (absent from the audit)"),
            "floor_edge": ("\n".join([row("A", "T1"), row("B", "T2"), row("C", "T3", ur_frac_got=0.127, ur_n_moved=500.0)]) + "\n",
                           0, "KA630 audit: PASS (FAILSAFE_TOTAL 0; contexts A,B,C; module " + K.P.UPRN_MODULE_SHA256[:12] + "; dump "),
            "fired": ("\n".join([row("A", "T1"), row("B", "T2", ur_guard_failed=1.0), row("C", "T3", ur_invariant_failed=1.0)]) + "\n",
                      1, "630dd: uprn FAILSAFE_TOTAL 2"),
            "no_ur": (json.dumps(dict(context="A", target="T1", zs_x=1.0)) + "\n", 1, "630dd: uprn FAILSAFE_TOTAL -1"),
            "bad": ("{not json\n", 1, "630dd: uprn FAILSAFE_TOTAL -1"),
            "missing": (None, 1, "630dd: uprn FAILSAFE_TOTAL -1"),
        }
        for name, (text, want_rc, want_line) in cases.items():
            dump = td / f"{name}.jsonl"
            if text is not None:
                dump.write_text(text)
            if name == "stale_json":   # a passing uprn_audit.json left by an earlier run must not be read
                (td / f"o_{name}").mkdir()
                (td / f"o_{name}" / "uprn_audit.json").write_text(json.dumps({c: dict(targets=1, mean_frac_got=0.3, mean_n_moved=900.0,
                    guard_failed=0, invariant_failed=0, no_pos_donor_genes=0, capacity_floor_met=True) for c in "ABC"}))
            r = subprocess.run([sys.executable, str(HERE / "ka630_audit.py"), str(dump), "--out-dir", str(td / f"o_{name}")],
                               capture_output=True, text=True)
            lab = {"clean": "A2", "fired": "A3", "floor_miss": "A5", "moved_miss": "A5", "ctx_absent": "A5", "stale_json": "A5", "floor_edge": "A5"}.get(name, "A4")
            hit = any(ln.startswith(want_line) for ln in r.stdout.splitlines()) if name == "floor_edge" else want_line in r.stdout.splitlines()
            rec(r.returncode == want_rc and hit, f"{lab} {name}: exit {r.returncode} (want {want_rc}), '{want_line[:60]}' present")
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS" if n == len(RES) else f"\n{len(RES) - n} FAILED of {len(RES)}")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
