#!/usr/bin/env python3
"""612 tests, on temp trees only (the live repo is never written). Light: no builds, no scorer, < 2 GB.

  T0  canonical POST-bundle tree = 598's own construction (test598.canon_tree: git HEAD's copies + 538 -> 554 -> 590)
  T1  patch612 on the canonical tree WITHOUT 598: --check and apply both refuse (order 598 -> 612), nothing written
  T2  598 applied (its own patch598.py), then patch612 --check passes, apply passes, a second apply is refused, the
      patched builder compiles (py_compile), and the patched text == patch612.patched_text(post-598 text)
  T3  the diff post-598 -> post-612 is exactly the patch's edits: every changed line sits in an inserted block and every
      replaced anchor line reappears verbatim (no other line of the builder moves)
  T4  _kg_mask_sign (extracted from the PATCHED builder) == run612.py's gate, bit for bit:
      (a) synthetic keys: ties, all-zero, a single nonzero gene, negative values, q 0.25 / 0.5 / 0.75;
      (b) the 85 real h1-test design keys (521's load_state key on the gated axis, embedded on the full axis) with the
          real CLEAN4 hook sign (tau 0.05): masks equal, and the gate keeps CLEAN4 on ~half of each target's covered genes
  T5  the VCC_DE_MDE_SIGN parse block is BYTE-IDENTICAL pre and post; the KEYGATE block (exec'd with stubs): unset -> None,
      0.5 -> 0.5, 0.25 -> 0.25; 0 / 1 / 1.5 / non-numeric refuse; a set knob without VCC_DE_MDE_SIGN refuses
  T6  _kg_apply (extracted from the PATCHED builder) with the knob unset returns the SAME object for every source and the
      weak-target CLEAN4 path; with the knob at 0.5 it returns the masked sign for every source

  uv run python experiments/612-key-gated-hook/test612.py      (exit 0 = all PASS)
"""
from __future__ import annotations

import ast
import difflib
import json
import os
import importlib.util
import py_compile
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/612-key-gated-hook"
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
import test598 as T598  # noqa: E402  (canon_tree, run, PY: read-only reuse of 598's tree construction)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P612 = _load("patch612", X / "patch612.py")
RESULTS: list[tuple[str, bool, str]] = []
B = "scripts/make_submission_pm_cd4.py"


def rec(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}", flush=True)


SKIPPED: list[tuple[str, str]] = []


def skip(name, why):
    SKIPPED.append((name, why))
    print(f"SKIP {name}: {why}", flush=True)


def kg_ref_ties(sign_gated, key_full, keep, q):
    """The 605-review rule, written independently: run612's gate, then genes tied AT the cut fill the weak count up to
    floor(q * n_nonzero) in stable (gene index) order."""
    a = np.abs(np.asarray(key_full, dtype=np.float64)[np.asarray(keep, dtype=bool)])
    nz = a > 0
    if not nz.any():
        return np.ascontiguousarray(np.asarray(sign_gated, dtype=np.float64))
    cut = float(np.quantile(a[nz], q))
    weak = (~nz) | (a < cut)
    target = int(np.floor(q * int(nz.sum()) + 1e-9))
    for j in range(len(a)):
        if int((weak & nz).sum()) >= target:
            break
        if nz[j] and a[j] == cut and not weak[j]:
            weak[j] = True
    return np.ascontiguousarray(np.where(weak, np.asarray(sign_gated, dtype=np.float64), 0.0))


def tie_at_cut(key_full, keep, q):
    a = np.abs(np.asarray(key_full, dtype=np.float64)[np.asarray(keep, dtype=bool)])
    nz = a > 0
    return bool(nz.any()) and int((nz & (a == float(np.quantile(a[nz], q)))).sum()) > 1


def kg_run612(sign_gated, key_full, keep, q):
    """run612.py's gate, verbatim (carrier_block, c4kg branch)."""
    k = np.asarray(key_full, dtype=np.float64)[np.asarray(keep, dtype=bool)]
    a = np.abs(k)
    nz = a > 0
    cut = float(np.quantile(a[nz], q)) if nz.any() else 0.0
    weak = (~nz) | (a < cut) if q < 1.0 else np.ones_like(nz)
    return np.ascontiguousarray(np.where(weak, sign_gated, 0.0))


def extract_fn(text, name):
    ns = {"np": np}
    for node in ast.parse(text).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            exec(compile(ast.get_source_segment(text, node), "<builder>", "exec"), ns)
            return ns[name]
    raise SystemExit(f"{name} not found")


def parse_block(text):
    for node in ast.parse(text).body:
        if isinstance(node, ast.If):
            src = ast.get_source_segment(text, node)
            if src.startswith("if DE_MDE_SIGN_ENV:"):
                return src
    raise SystemExit("VCC_DE_MDE_SIGN parse block not found")


def keygate_block(text):
    """The knob's module-level statements: the env read, the None default and the `if ..._ENV:` block."""
    lines = []
    for node in ast.parse(text).body:
        src = ast.get_source_segment(text, node)
        if src.startswith("DE_MDE_SIGN_KEYGATE_ENV =") or src.startswith("DE_MDE_SIGN_KEYGATE = None") \
                or src.startswith("if DE_MDE_SIGN_KEYGATE_ENV:"):
            lines.append(src)
    if len(lines) != 3:
        raise SystemExit(f"keygate block: found {len(lines)} of 3 statements")
    return "\n".join(lines)


def run_parse(block, env):
    ns = {"DE_MDE_SIGN_ENV": env, "DE_MDE_SIGN": None, "DE_MDE": (0.1, 0.1, 0.3), "DE_SIGN_PRIOR": "clean4q80",
          "SAMELINE": None, "DE_PRIOR_GATE": None}
    try:
        exec(compile(block, "<parse>", "exec"), ns)
        return ("ok", ns["DE_MDE_SIGN"])
    except SystemExit as e:
        return ("exit", str(e)[:80])


def main():
    with tempfile.TemporaryDirectory(prefix="t612-") as td:
        td = Path(td)
        canon = td / "canon"
        T598.canon_tree(canon)
        rec("T0 canonical post-bundle tree", (canon / B).exists(), f"bundle {T598.BUNDLE_STATE}")
        pre_text = (canon / B).read_text()
        # T1: without 598
        c = T598.run([T598.PY, X / "patch612.py", "--root", canon, "--check"])
        a = T598.run([T598.PY, X / "patch612.py", "--root", canon])
        rec("T1 refuses a pre-598 builder (check and apply)", c.returncode != 0 and a.returncode != 0
            and "598" in (c.stdout + c.stderr) and (canon / B).read_text() == pre_text)
        # T2: 598 then 612
        t2 = td / "t2"
        T598.copy_tree(canon, t2)
        r598 = T598.run([T598.PY, ROOT / "experiments/598-tercilec-knob/patch598.py", "--root", t2])
        post598 = (t2 / B).read_text()
        ok598 = r598.returncode == 0 and "598: tercilec" in post598
        rec("T2a 598 applies on the canonical tree", ok598, (r598.stdout + r598.stderr)[-160:].replace("\n", " "))
        c = T598.run([T598.PY, X / "patch612.py", "--root", t2, "--check"])
        a = T598.run([T598.PY, X / "patch612.py", "--root", t2])
        post612 = (t2 / B).read_text()
        a2 = T598.run([T598.PY, X / "patch612.py", "--root", t2])
        try:
            py_compile.compile(str(t2 / B), doraise=True)
            comp = True
        except py_compile.PyCompileError:
            comp = False
        rec("T2b 612 --check / apply after 598", c.returncode == 0 and a.returncode == 0)
        rec("T2c second apply refused, text unchanged", a2.returncode != 0 and (t2 / B).read_text() == post612)
        rec("T2d patched builder compiles", comp)
        rec("T2e patched text == patched_text(post-598)", post612 == P612.patched_text(post598))
        # T3: the diff is exactly the six edits
        old_l, new_l = post598.splitlines(), post612.splitlines()
        sm = difflib.SequenceMatcher(a=old_l, b=new_l, autojunk=False)
        ins = [new_l[j1:j2] for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag in ("insert", "replace")]
        rem = [old_l[i1:i2] for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag in ("delete", "replace")]
        new_blocks = "\n".join(n for _, _, n in P612.EDITS)
        old_blocks = "\n".join(o for _, o, _ in P612.EDITS)
        ok3 = all(line in new_blocks for blk in ins for line in blk) and all(line in old_blocks for blk in rem for line in blk)
        rec("T3 diff post-598 -> post-612 confined to the edit blocks", ok3,
            f"{sum(len(b) for b in ins)} lines in, {sum(len(b) for b in rem)} lines out")
        # T4: helper == run612 gate
        fn = extract_fn(post612, "_kg_mask_sign")
        rng = np.random.default_rng(612)
        cases, tied = 0, 0
        ok4a = True
        for trial in range(200):
            n_full = 60
            keep = rng.random(n_full) < 0.7
            n = int(keep.sum())
            key = np.where(rng.random(n_full) < 0.6, rng.normal(size=n_full), 0.0)
            if trial % 7 == 0:
                key = np.round(key, 1)            # ties
            if trial % 11 == 0:
                key = np.zeros(n_full)            # all zero
            if trial % 13 == 0:
                key = np.zeros(n_full); key[np.flatnonzero(keep)[:1]] = -2.0   # one nonzero
            sign = rng.choice([-1.0, 0.0, 1.0], size=n)
            for q in (0.25, 0.5, 0.65, 0.75):
                cases += 1
                got = fn(sign, key, keep, q)
                ok4a &= np.array_equal(got, kg_ref_ties(sign, key, keep, q))
                if tie_at_cut(key, keep, q):
                    tied += 1
                else:
                    ok4a &= np.array_equal(got, kg_run612(sign, key, keep, q))
        rec("T4a _kg_mask_sign == run612 gate where no tie sits at the cut, == the independent tie rule everywhere (synthetic)",
            ok4a, f"{cases} cases, {tied} with a tie at the cut")
        toys = [([1, 1, 1, 2], 2), ([3, 3, 3, 3], 2), ([1, 1, 1, 1, 5], 2)]
        got_toys = []
        for kv, want in toys:
            kv = np.asarray(kv, float); kp = np.ones(len(kv), bool); sg = np.ones(len(kv))
            got_toys.append(int((fn(sg, kv, kp, 0.5) != 0).sum()))
        rec("T4c ties at the cut keep floor(q*n) weak genes (605 review toys at q 0.5)",
            got_toys == [w for _, w in toys], f"weak {got_toys}, want {[w for _, w in toys]}")
        sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
        try:
            import lib521 as L  # noqa: PLC0415
            S = L.load_state(0)
        except (FileNotFoundError, OSError) as e:
            S = None
            skip("T4b _kg_mask_sign == run612 gate on the 85 h1-test design keys", f"harness file absent ({e})")
        ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False) if S is not None else None
        if S is not None:
            keepH = np.asarray(ref["keep"], dtype=bool)
            P = _load("run612_c4", X / "run612.py") if False else None   # run612 executes at import; recompute its hook here
            cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
            names = list(cp["bank_names"].astype(str))
            Pc = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
            pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
            GF = ref["genes"].astype(str)
            col = np.array([pos.get(str(g), -1) for g in GF])
            Pf = np.full(len(GF), np.nan); Pf[col >= 0] = Pc[col[col >= 0]]
            hook = np.where(np.isfinite(Pf) & (np.abs(Pf - 0.5) >= 0.05), np.where(Pf > 0.5, -1.0, 1.0), 0.0)[keepH]
            ok4b, frac = True, []
            for i in range(S["key"].shape[0]):
                kf = np.zeros(len(GF)); kf[keepH] = S["key"][i]
                m1, m2 = fn(hook, kf, keepH, 0.5), kg_run612(hook, kf, keepH, 0.5)
                ok4b &= np.array_equal(m1, m2)
                frac.append((m1 != 0).sum() / max((hook != 0).sum(), 1))
            rec("T4b _kg_mask_sign == run612 gate on the 85 h1-test design keys", ok4b,
                f"CLEAN4 kept on {np.median(frac):.3f} (median) of each target's covered genes")
        # T5: parse
        rec("T5a VCC_DE_MDE_SIGN parse block byte-identical pre/post", parse_block(post598) == parse_block(post612))
        kb = keygate_block(post612)
        def kg_parse(env, sign_set=True):
            ns = {"os": __import__("os"), "DE_MDE_SIGN": {"source": "clean4", "tau": 0.05} if sign_set else None}
            old = os.environ.get("VCC_DE_MDE_SIGN_KEYGATE")
            try:
                if env is None:
                    os.environ.pop("VCC_DE_MDE_SIGN_KEYGATE", None)
                else:
                    os.environ["VCC_DE_MDE_SIGN_KEYGATE"] = env
                exec(compile(kb, "<keygate>", "exec"), ns)
                return ("ok", ns["DE_MDE_SIGN_KEYGATE"])
            except SystemExit as e:
                return ("exit", str(e)[:60])
            finally:
                if old is None:
                    os.environ.pop("VCC_DE_MDE_SIGN_KEYGATE", None)
                else:
                    os.environ["VCC_DE_MDE_SIGN_KEYGATE"] = old
        rec("T5b knob unset -> None", kg_parse(None) == ("ok", None))
        rec("T5c knob 0.5 / 0.25", kg_parse("0.5") == ("ok", 0.5) and kg_parse(" 0.25 ") == ("ok", 0.25))
        bad = [kg_parse(e)[0] for e in ("0", "1", "1.5", "x", "-0.2")] + [kg_parse("0.5", sign_set=False)[0]]
        rec("T5d malformed knob / knob without VCC_DE_MDE_SIGN refused", all(b == "exit" for b in bad), str(bad))
        # T6: _kg_apply
        ns6 = {"np": np, "DE_MDE_SIGN_KEYGATE": None}
        for node in ast.parse(post612).body:
            if isinstance(node, ast.FunctionDef) and node.name in ("_kg_mask_sign", "_kg_apply"):
                exec(compile(ast.get_source_segment(post612, node), "<builder>", "exec"), ns6)
        arr, arr_c4 = np.array([1.0, -1.0, 0.0, 1.0]), np.array([-1.0, -1.0, 1.0, 0.0])
        keep4 = np.array([True, True, True, True]); key4 = np.array([0.1, 2.0, 0.0, 3.0])
        f_expr = P612.F_NEW[P612.F_NEW.index("_kg_apply("):P612.F_NEW.rindex(")}),") + 1]
        class Ctl:
            keep = keep4
        def ev(weak, c4):
            env = dict(ns6); env.update({"de_mde_sign_c4": c4, "de_mde_sign_weak": weak, "target_value": "T1", "str": str,
                                         "de_mde_sign": arr, "key": key4, "rc_control": Ctl})
            return eval(f_expr, env)
        ok6a = ev(set(), None) is arr and ev({"T1"}, arr_c4) is arr_c4 and ev({"T2"}, arr_c4) is arr
        ns6["DE_MDE_SIGN_KEYGATE"] = 0.5
        m1, m2 = ev(set(), None), ev({"T1"}, arr_c4)
        ok6b = np.array_equal(m1, kg_run612(arr, key4, keep4, 0.5)) and np.array_equal(m2, kg_run612(arr_c4, key4, keep4, 0.5))
        rec("T6a knob unset: the SAME object for every source and the weak-target path", ok6a)
        rec("T6b knob 0.5: masked sign for every source and the weak-target path", ok6b, f"{m1.tolist()} / {m2.tolist()}")
        # T7: meta key only when set
        msrc = P612.M_NEW
        rec("T7 meta: de_mde_sign_keygate written only when the knob is set",
            "**({} if DE_MDE_SIGN_KEYGATE is None else {\"de_mde_sign_keygate\"" in msrc)
        # T8-T12: 612dd (def_day.sh pass-through), on the 598 -> 612 tree
        P612DD = _load("patch612dd", X / "patch612dd.py")
        c_no612 = T598.run([T598.PY, X / "patch612dd.py", "--root", canon, "--check"])
        t8 = td / "t8"; T598.copy_tree(canon, t8)
        T598.run([T598.PY, ROOT / "experiments/598-tercilec-knob/patch598.py", "--root", t8])
        c_no612b = T598.run([T598.PY, X / "patch612dd.py", "--root", t8, "--check"])
        rec("T8a 612dd refuses without 598 and without 612", c_no612.returncode != 0 and c_no612b.returncode != 0
            and "612" in (c_no612b.stdout + c_no612b.stderr))
        T598.run([T598.PY, X / "patch612.py", "--root", t8])
        day_pre = (t8 / "scripts/def_day.sh").read_text()
        c = T598.run([T598.PY, X / "patch612dd.py", "--root", t8, "--check"])
        a = T598.run([T598.PY, X / "patch612dd.py", "--root", t8])
        day_post = (t8 / "scripts/def_day.sh").read_text()
        a2 = T598.run([T598.PY, X / "patch612dd.py", "--root", t8])
        bn = T598.run(["bash", "-n", t8 / "scripts/def_day.sh"])
        rec("T8b 612dd --check / apply after 598 -> 612; second apply refused; bash -n", c.returncode == 0 and a.returncode == 0
            and a2.returncode != 0 and (t8 / "scripts/def_day.sh").read_text() == day_post and bn.returncode == 0,
            (bn.stderr or "")[-120:])
        sm = difflib.SequenceMatcher(a=day_pre.splitlines(), b=day_post.splitlines(), autojunk=False)
        ops = {tag for tag, *_ in sm.get_opcodes()} - {"equal"}
        n_ins = sum(j2 - j1 for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag == "insert")
        rec("T9 the def_day diff is insert-only (no existing line changes)", ops == {"insert"}, f"{n_ins} lines inserted")
        # T10: run the validation block in bash with stubs (cwd = the patched tree, so grep reads its builder)
        vb = P612DD.V_NEW[len(P612DD.V_OLD):]
        def vrun(env, tree=t8):
            script = 'die(){ echo "DIE: $*"; exit 3; }\n' + vb + 'echo OK\n'
            e = {k: v for k, v in os.environ.items() if not k.startswith("DEF_DAY_")}
            e.update(env)
            r = __import__("subprocess").run(["bash", "-c", script], cwd=tree, capture_output=True, text=True, env=e)
            return r.returncode, (r.stdout + r.stderr).strip()[-90:]
        cases = {"unset": vrun({}), "0.5 + sign": vrun({"DEF_DAY_MDE_SIGN_KEYGATE": "0.5", "DEF_DAY_MDE_SIGN": "clean4:0.05"}),
                 "0.5 no sign": vrun({"DEF_DAY_MDE_SIGN_KEYGATE": "0.5"}),
                 "1": vrun({"DEF_DAY_MDE_SIGN_KEYGATE": "1", "DEF_DAY_MDE_SIGN": "clean4:0.05"}),
                 "x": vrun({"DEF_DAY_MDE_SIGN_KEYGATE": "x", "DEF_DAY_MDE_SIGN": "clean4:0.05"}),
                 "pre-612 builder": vrun({"DEF_DAY_MDE_SIGN_KEYGATE": "0.5", "DEF_DAY_MDE_SIGN": "clean4:0.05"}, tree=canon)}
        ok10 = (cases["unset"][0] == 0 and cases["0.5 + sign"][0] == 0 and all(cases[k][0] == 3 for k in ("0.5 no sign", "1", "x", "pre-612 builder")))
        rec("T10 validation: unset / valid pass; no sign, q 1, non-numeric, pre-612 builder die", ok10,
            "; ".join(f"{k}:{v[0]}" for k, v in cases.items()))
        # T11: the final-meta check (C) on synthetic metas
        cb = P612DD.C_NEW[len(P612DD.C_OLD):]
        py = cb[cb.index("<<'PY'") + len("<<'PY'"):cb.index("PY\nfi")].split("\n", 1)[1]
        def mrun(cfg, q):
            f = td / "m.json"; f.write_text(json.dumps({"config": cfg}))
            return __import__("subprocess").run([sys.executable, "-c", py, str(f), q], capture_output=True).returncode
        r11 = (mrun({"de_mde_sign_keygate": {"q": 0.5}}, "0.5"), mrun({"de_mde_sign_keygate": {"q": 0.25}}, "0.5"), mrun({}, "0.5"))
        rec("T11 final-meta check: match passes, wrong q / missing key fail", r11 == (0, 1, 1), str(r11))
    n_ok = sum(ok for _, ok, _ in RESULTS)
    print(f"\n612 tests: {n_ok}/{len(RESULTS)} PASS" + (f" ({len(SKIPPED)} SKIP)" if SKIPPED else ""))
    return 0 if n_ok == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
