#!/usr/bin/env python3
"""630 composition tests, on temp trees only (the live repo is never written). Light: no builds, no scorer, < 1 GB.

  T0  canonical POST-bundle tree = 598's own construction (test598.canon_tree: git HEAD's copies + 538 -> 554 -> 590)
  T1  patch630 on the canonical tree WITHOUT 598/612: --check and apply both refuse (order 598 -> 612 -> 630), nothing written
  T2  598 then 612 applied (their own patch scripts), then patch630 --check passes, apply passes, a second apply is refused,
      the patched builder compiles (py_compile), the patched text == patch630.patched_text(post-612 text), and
      src/vcc2026/rank_up.py is installed with the pinned sha
  T3  the diff post-612 -> post-630 is exactly the patch's inserted blocks (every replaced anchor line reappears verbatim,
      no other builder line moves)
  T4  the knob block (exec'd with stubs): unset -> None; b_k100 / s_k075 / b_k300 / b_k050 accepted; b_k040, b_k310, x_k100,
      b_100, '' (whitespace) handled; a set knob without VCC_DE_MAG refuses
  T5  the builder's pinned sha == the installed module's sha == experiments/630's src_rank_up.py

  uv run python experiments/630-rank-neutral-up/test630.py      (exit 0 = all PASS)
"""
from __future__ import annotations

import ast
import difflib
import hashlib
import importlib.util
import os
import py_compile
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/630-rank-neutral-up"
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
import test598 as T598  # noqa: E402  (canon_tree, run, PY, copy_tree: read-only reuse of 598's tree construction)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P630 = _load("patch630", X / "patch_builder630.py")
RESULTS: list[tuple[str, bool, str]] = []
B = "scripts/make_submission_pm_cd4.py"
M = "src/vcc2026/rank_up.py"


def rec(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}", flush=True)


def knob_block(text):
    """The module-level `if DE_UPRN_ENV:` statement plus the three assignments before it."""
    tree = ast.parse(text)
    src = []
    for node in tree.body:
        s = ast.get_source_segment(text, node) or ""
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") in ("DE_UPRN_ENV", "DE_UPRN", "DE_UPRN_MODULE_PATH",
                                                                         "DE_UPRN_MODULE_SHA256") for t in node.targets):
            src.append(s)
        if isinstance(node, ast.If) and s.startswith("if DE_UPRN_ENV:"):
            src.append(s)
    return "\n".join(src)


def run_knob(block_src, env_val, de_mag):
    ns = {"os": os, "ROOT": Path("/tmp/x"), "DE_MAG": de_mag}
    old = os.environ.get("VCC_DE_UPRN")
    try:
        if env_val is None:
            os.environ.pop("VCC_DE_UPRN", None)
        else:
            os.environ["VCC_DE_UPRN"] = env_val
        exec(compile(block_src, "<knob>", "exec"), ns)
        return ("ok", ns["DE_UPRN"])
    except SystemExit as e:
        return ("exit", str(e))
    finally:
        if old is None:
            os.environ.pop("VCC_DE_UPRN", None)
        else:
            os.environ["VCC_DE_UPRN"] = old


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="test630-") as td:
        td = Path(td)
        canon = td / "canon"
        T598.canon_tree(canon)
        rec("T0 canonical post-bundle tree", (canon / B).exists(), f"bundle {T598.BUNDLE_STATE}")
        # T1: refuses without 598/612
        before = (canon / B).read_bytes()
        c = T598.run([T598.PY, X / "patch_builder630.py", "--root", canon, "--check"])
        a = T598.run([T598.PY, X / "patch_builder630.py", "--root", canon])
        rec("T1 refuses without 612 (check)", c.returncode != 0 and "612" in (c.stdout + c.stderr), (c.stdout + c.stderr)[-160:])
        rec("T1 refuses without 612 (apply), nothing written",
            a.returncode != 0 and (canon / B).read_bytes() == before and not (canon / M).exists())
        # T2: 598 -> 612 -> 630
        t2 = td / "t2"
        T598.copy_tree(canon, t2)
        r598 = T598.run([T598.PY, ROOT / "experiments/598-tercilec-knob/patch598.py", "--root", t2])
        rec("T2 598 applies", r598.returncode == 0, (r598.stdout + r598.stderr)[-160:])
        r612 = T598.run([T598.PY, ROOT / "experiments/612-key-gated-hook/patch612.py", "--root", t2])
        rec("T2 612 applies after 598", r612.returncode == 0, (r612.stdout + r612.stderr)[-160:])
        post612 = (t2 / B).read_text()
        c = T598.run([T598.PY, X / "patch_builder630.py", "--root", t2, "--check"])
        rec("T2 630 --check passes after 612", c.returncode == 0, (c.stdout + c.stderr)[-200:])
        rec("T2 --check wrote nothing", (t2 / B).read_text() == post612 and not (t2 / M).exists())
        a = T598.run([T598.PY, X / "patch_builder630.py", "--root", t2])
        rec("T2 630 applies", a.returncode == 0, (a.stdout + a.stderr)[-200:])
        a2 = T598.run([T598.PY, X / "patch_builder630.py", "--root", t2])
        rec("T2 second apply refused", a2.returncode != 0 and "already" in (a2.stdout + a2.stderr))
        post630 = (t2 / B).read_text()
        try:
            py_compile.compile(str(t2 / B), doraise=True)
            ok = True
        except py_compile.PyCompileError as e:
            ok, _ = False, e
        rec("T2 patched builder compiles", ok)
        rec("T2 patched text == patched_text(post-612)", post630 == P630.patched_text(post612))
        msha = hashlib.sha256((t2 / M).read_bytes()).hexdigest() if (t2 / M).exists() else None
        rec("T2 module installed with the pinned sha", msha == P630.MODULE_SHA, f"{msha}")
        # T3: diff is exactly the inserted blocks
        a_lines, b_lines = post612.splitlines(), post630.splitlines()
        sm = difflib.SequenceMatcher(a=a_lines, b=b_lines, autojunk=False)
        bad = [op for op in sm.get_opcodes() if op[0] in ("replace", "delete")]
        inserted = sum(j2 - j1 for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag == "insert")
        expected = sum(new.count("\n") - old.count("\n") for old, new in P630.EDITS)
        rec("T3 only insertions (no builder line replaced or deleted)", not bad, f"{bad[:2]}")
        rec("T3 inserted line count == the patch's", inserted == expected, f"{inserted} vs {expected}")
        # T4: the knob
        kb = knob_block(post630)
        cases = [(None, "BASE", ("ok", None)), ("b_k100", "BASE", ("ok", "b_k100")), ("s_k075", "BASE", ("ok", "s_k075")),
                 ("b_k300", "BASE", ("ok", "b_k300")), ("b_k050", "BASE", ("ok", "b_k050")), ("   ", "BASE", ("ok", None)),
                 ("b_k100_q20r3p10", "BASE", ("ok", "b_k100_q20r3p10")), ("s_k075_q40r3p10", "BASE", ("ok", "s_k075_q40r3p10"))]
        for val, dm, want in cases:
            got = run_knob(kb, val, dm)
            rec(f"T4 knob {val!r} -> {want[1]!r}", got == want, f"{got}")
        for val in ("b_k040", "b_k310", "x_k100", "b_100", "b_k1000", "b_k100_q20r3", "b_k100_q20",
                    "b_k100_q20r0p10", "b_k100_q0r3p10", "b_k100_q20r3p0", "b_k100_q020r3p10"):   # 607: no zero (or 0-led) suffix
            got = run_knob(kb, val, "dn_k250_zs")
            rec(f"T4 knob {val!r} refused", got[0] == "exit", f"{got}")
        got = run_knob(kb, "b_k100", None)
        rec("T4 knob set without VCC_DE_MAG refused", got[0] == "exit" and "VCC_DE_MAG" in got[1], f"{got}")
        # T5: sha chain
        pinned = [ln for ln in post630.splitlines() if ln.startswith("DE_UPRN_MODULE_SHA256 = ")]
        src_sha = hashlib.sha256((X / "src_rank_up.py").read_bytes()).hexdigest()
        rec("T5 pinned sha == installed module == 630/src_rank_up.py",
            len(pinned) == 1 and src_sha in pinned[0] and msha == src_sha, f"{pinned} {src_sha[:16]}")
    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} PASS")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
