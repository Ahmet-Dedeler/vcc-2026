#!/usr/bin/env python3
"""631 composition tests (adapted from 630's), on temp trees only (the live repo is never written). < 1 GB, no builds.

  T0  canonical POST-bundle tree = 598's construction (test598.canon_tree)
  T1  patch631 on 598 -> 612 WITHOUT 630: --check and apply both refuse (order ... -> 630 -> 631), nothing written
  T2  598 -> 612 -> 630 applied, then patch631 --check passes, apply passes, a second apply is refused, the patched builder
      compiles, the patched text == patch631.patched_text(post-630 text), src/vcc2026/dn_mag.py installed with the pin
  T3  the diff post-630 -> post-631 is insertions only, exactly the patch's line count
  T4  the knob: unset -> None; d_m200 / d_m125 / d_m400 / d_m250_q3 accepted; d_m100, d_m410, d_m200_q0, d_m200_q03,
      x_m200, d_200 refused; a set knob without VCC_DE_MAG refuses
  T5  the builder's pinned sha == the installed module == 631/src_dn_mag631.py
  T6  the DOWN block sits AFTER 630's UP call in the per-block code and at the zs level of indentation
  uv run python experiments/631-nmae-headroom/test631.py      (exit 0 = all PASS)
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
X = ROOT / "experiments/631-nmae-headroom"
X630 = ROOT / "experiments/630-rank-neutral-up"
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
import test598 as T598  # noqa: E402  (canon_tree, run, PY, copy_tree: read-only reuse of 598's tree construction)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P631 = _load("patch631", X / "patch_builder631.py")
RESULTS: list[tuple[str, bool, str]] = []
B = "scripts/make_submission_pm_cd4.py"
M = "src/vcc2026/dn_mag.py"


def rec(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}", flush=True)


def knob_block(text):
    """The module-level `if DE_UPRN_ENV:` statement plus the three assignments before it."""
    tree = ast.parse(text)
    src = []
    for node in tree.body:
        s = ast.get_source_segment(text, node) or ""
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") in ("DE_DNMAG_ENV", "DE_DNMAG", "DE_DNMAG_MODULE_PATH",
                                                                         "DE_DNMAG_MODULE_SHA256") for t in node.targets):
            src.append(s)
        if isinstance(node, ast.If) and s.startswith("if DE_DNMAG_ENV:"):
            src.append(s)
    return "\n".join(src)


def run_knob(block_src, env_val, de_mag):
    ns = {"os": os, "ROOT": Path("/tmp/x"), "DE_MAG": de_mag}
    old = os.environ.get("VCC_DE_DNMAG")
    try:
        if env_val is None:
            os.environ.pop("VCC_DE_DNMAG", None)
        else:
            os.environ["VCC_DE_DNMAG"] = env_val
        exec(compile(block_src, "<knob>", "exec"), ns)
        return ("ok", ns["DE_DNMAG"])
    except SystemExit as e:
        return ("exit", str(e))
    finally:
        if old is None:
            os.environ.pop("VCC_DE_DNMAG", None)
        else:
            os.environ["VCC_DE_DNMAG"] = old


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="test631-") as td:
        td = Path(td)
        canon = td / "canon"
        T598.canon_tree(canon)
        rec("T0 canonical post-bundle tree", (canon / B).exists(), f"bundle {T598.BUNDLE_STATE}")
        t1 = td / "t1"
        T598.copy_tree(canon, t1)
        T598.run([T598.PY, ROOT / "experiments/598-tercilec-knob/patch598.py", "--root", t1])
        T598.run([T598.PY, ROOT / "experiments/612-key-gated-hook/patch612.py", "--root", t1])
        before = (t1 / B).read_bytes()
        c = T598.run([T598.PY, X / "patch_builder631.py", "--root", t1, "--check"])
        a = T598.run([T598.PY, X / "patch_builder631.py", "--root", t1])
        rec("T1 refuses without 630 (check)", c.returncode != 0 and "630" in (c.stdout + c.stderr), (c.stdout + c.stderr)[-160:])
        rec("T1 refuses without 630 (apply), nothing written",
            a.returncode != 0 and (t1 / B).read_bytes() == before and not (t1 / M).exists())
        t2 = td / "t2"
        T598.copy_tree(canon, t2)
        for script in ("experiments/598-tercilec-knob/patch598.py", "experiments/612-key-gated-hook/patch612.py",
                       "experiments/630-rank-neutral-up/patch_builder630.py"):
            r = T598.run([T598.PY, ROOT / script, "--root", t2])
            rec(f"T2 {Path(script).name} applies", r.returncode == 0, (r.stdout + r.stderr)[-160:])
        post630 = (t2 / B).read_text()
        c = T598.run([T598.PY, X / "patch_builder631.py", "--root", t2, "--check"])
        rec("T2 631 --check passes after 630", c.returncode == 0, (c.stdout + c.stderr)[-200:])
        rec("T2 --check wrote nothing", (t2 / B).read_text() == post630 and not (t2 / M).exists())
        a = T598.run([T598.PY, X / "patch_builder631.py", "--root", t2])
        rec("T2 631 applies", a.returncode == 0, (a.stdout + a.stderr)[-200:])
        a2 = T598.run([T598.PY, X / "patch_builder631.py", "--root", t2])
        rec("T2 second apply refused", a2.returncode != 0 and "already" in (a2.stdout + a2.stderr))
        post631 = (t2 / B).read_text()
        try:
            py_compile.compile(str(t2 / B), doraise=True)
            ok = True
        except py_compile.PyCompileError:
            ok = False
        rec("T2 patched builder compiles", ok)
        rec("T2 patched text == patched_text(post-630)", post631 == P631.patched_text(post630))
        msha = hashlib.sha256((t2 / M).read_bytes()).hexdigest() if (t2 / M).exists() else None
        rec("T2 module installed with the pinned sha", msha == P631.MODULE_SHA, f"{msha}")
        a_lines, b_lines = post630.splitlines(), post631.splitlines()
        sm = difflib.SequenceMatcher(a=a_lines, b=b_lines, autojunk=False)
        bad = [op for op in sm.get_opcodes() if op[0] in ("replace", "delete")]
        inserted = sum(j2 - j1 for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag == "insert")
        expected = sum(new.count("\n") - old.count("\n") for old, new in P631.EDITS)
        rec("T3 only insertions", not bad, f"{bad[:2]}")
        rec("T3 inserted line count == the patch's", inserted == expected, f"{inserted} vs {expected}")
        kb = knob_block(post631)
        for val, want in ((None, None), ("d_m200", "d_m200"), ("d_m125", "d_m125"), ("d_m400", "d_m400"),
                          ("d_m250_q3", "d_m250_q3"), ("   ", None)):
            got = run_knob(kb, val, "dn_k250_zs")
            rec(f"T4 knob {val!r} -> {want!r}", got == ("ok", want), f"{got}")
        for val in ("d_m100", "d_m410", "d_m200_q0", "d_m200_q03", "x_m200", "d_200", "d_m2000"):
            got = run_knob(kb, val, "dn_k250_zs")
            rec(f"T4 knob {val!r} refused", got[0] == "exit", f"{got}")
        got = run_knob(kb, "d_m200", None)
        rec("T4 knob set without VCC_DE_MAG refused", got[0] == "exit" and "VCC_DE_MAG" in got[1], f"{got}")
        pinned = [ln for ln in post631.splitlines() if ln.startswith("DE_DNMAG_MODULE_SHA256 = ")]
        src_sha = hashlib.sha256((X / "src_dn_mag631.py").read_bytes()).hexdigest()
        rec("T5 pinned sha == installed module == 631/src_dn_mag631.py",
            len(pinned) == 1 and src_sha in pinned[0] and msha == src_sha, f"{pinned} {src_sha[:16]}")
        i_up = post631.index("block, _ur_diag = _rank_up.apply(")
        i_dn = post631.index("block, _dn_diag = _dn_mag.apply(")
        line_if = [ln for ln in post631.splitlines() if ln.strip() == "if DE_DNMAG is not None:" and ln.startswith(" " * 28)
                   and not ln.startswith(" " * 29)]
        line_up = [ln for ln in post631.splitlines() if ln.strip() == "if DE_UPRN is not None:" and ln.startswith(" " * 28)
                   and not ln.startswith(" " * 29)]
        rec("T6 DOWN pass after 630's UP call, at the same (zs-level) indentation", i_dn > i_up and len(line_if) == 1
            and len(line_up) == 1, f"{i_up} {i_dn} {len(line_if)} {len(line_up)}")
    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} PASS")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
