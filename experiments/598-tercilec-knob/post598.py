#!/usr/bin/env python3
"""598 post-apply check: run after `patch598.py --root <repo>` on the live tree and BEFORE committing it. Read-only on the
live tree (temp files only). Light (~0.3 GB).

1. The six files equal a canonical application: git HEAD's copies (+ any 10-05 bundle step not yet on HEAD, as test598
   builds it) + 598, rebuilt in a temp tree. HEAD must still be the pre-598 commit.
2. On the real A/B/C dump (479p2g's 00ccfcba), the LIVE builder's functions against the canonical pre-598 builder's:
   tercile and stercile q bytes and audits identical; tercilec's audit + key record == work/ka598c_expected.json.
3. bash -n on the live def_day.sh; py_compile on the five live Python files.

  uv run python experiments/598-tercilec-knob/post598.py        (exit 0 = ALL PASS)
"""
from __future__ import annotations

import hashlib
import json
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import patch598  # noqa: E402
import test598 as T  # noqa: E402

ROOT = T.ROOT


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    fails = []
    with tempfile.TemporaryDirectory(prefix="post598-") as td:
        canon = Path(td) / "canon"
        T.canon_tree(canon)
        print(f"canonical tree: HEAD {subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip()}"
              f" + bundle steps {T.BUNDLE_STATE}")
        try:
            want = patch598.patched_texts(canon)
        except SystemExit as e:
            print(f"[1] cannot rebuild the canonical application: {e} (is 598 already committed on HEAD?)")
            return 1
        for rel in patch598.FILES:
            got = (ROOT / rel).read_text()
            ok = got == want[rel]
            print(f"[1] {'OK  ' if ok else 'DIFF'} {rel} live {sha(ROOT / rel)[:16]} canonical "
                  f"{hashlib.sha256(want[rel].encode()).hexdigest()[:16]}")
            if not ok:
                fails.append(f"sha {rel}")
        old = T.module_ns((canon / patch598.BUILDER).read_text())
        new = T.module_ns((ROOT / patch598.BUILDER).read_text())
    with np.load(T.DUMP, allow_pickle=False) as z:
        bulk, car = np.asarray(z["bulk"], np.float64), np.asarray(z["carriers"], np.float64)
        perts, genes = [str(x) for x in z["perts"]], [str(x) for x in z["genes"]]
        env = json.loads(str(z["env"]))
    fo, fn = old["sharpen_q_by_tercile"], new["sharpen_q_by_tercile"]
    for mode in ("tercile", "stercile"):
        ok = T.same_q(fo, fn, bulk, car, genes, perts, mode=mode)
        print(f"[2] {mode}: live == pre-598 (q bytes, audit): {'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"{mode} identity")
    carried, key = new["dump_key_carried"](env, perts, bulk, car, genes)
    _, a = fn(bulk, car, genes, perts, {"mode": "tercilec", "q": [2.5, 2.0, 1.75]}, carried=carried)
    got = {**a, **key}
    exp = json.loads((HERE / "work/ka598c_expected.json").read_text())
    bad = [k for k in exp if (got.get(k) != exp[k] if k not in ("edges", "key_ratio_c", "key_ratio_maxdev")
                              else not np.allclose(np.asarray(got.get(k), float), np.asarray(exp[k], float), rtol=1e-12, atol=1e-15))]
    print(f"[2] tercilec audit == ka598c_expected.json: {'PASS' if not bad and set(got) == set(exp) else 'FAIL ' + str(bad)}"
          f" (counts {a['counts']}, moved {sorted(a['moved_vs_tercile'])})")
    if bad or set(got) != set(exp):
        fails.append("tercilec audit")
    bn = subprocess.run(["/bin/bash", "-n", str(ROOT / patch598.DAY)], capture_output=True, text=True)
    okc = bn.returncode == 0
    with tempfile.TemporaryDirectory() as td:
        for rel in patch598.FILES:
            if rel.suffix == ".py":
                try:
                    py_compile.compile(str(ROOT / rel), cfile=str(Path(td) / (rel.name + "c")), doraise=True)
                except py_compile.PyCompileError as e:
                    okc = False
                    print(e)
    print(f"[3] bash -n def_day.sh, py_compile x5: {'PASS' if okc else 'FAIL ' + bn.stderr.strip()}")
    if not okc:
        fails.append("syntax")
    print("ALL PASS" if not fails else f"FAILED: {fails}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
