#!/usr/bin/env python3
"""537 — the builder knob `VCC_EMISSION_SEED` (NOT applied): re-roll ONLY the count-emission stream of a build, so the board
can read the re-roll noise of one design directly (RESULTS §2.4 / §5: the re-roll pair).

Today `make_rngs(SEED)` returns (default_rng(SEED), default_rng(SEED + EMISSION_SEED_OFFSET)) = (library picks, emission)
and the meta records `"emission_seed": SEED + EMISSION_SEED_OFFSET` (1,000,003). The patch adds

    EMISSION_SEED_ENV = os.environ.get("VCC_EMISSION_SEED", "").strip()
    EMISSION_SEED = int(EMISSION_SEED_ENV) if EMISSION_SEED_ENV else SEED + EMISSION_SEED_OFFSET

uses EMISSION_SEED for the emission stream in `make_rngs` when it is called with the module SEED (the one production call,
`make_rngs(SEED)`), and records it in the meta's `emission_seed`. Unset: the same streams and the same meta, byte for byte.
Set (e.g. 1000004): the library picks are unchanged, every count draw is new, and the config diff vs the base is exactly
`emission_seed`.

usage: uv run python experiments/537-arm-pipeline/patch_builder537_seed.py <builder.py> [--check]
  --check: patch a temporary copy and run the self-test (never touches the given file).
The coordinator applies it (after a build logs "builder after restore"), then runs KA-1 (446bj2z's env, knob unset ->
ec3a5ae6 byte for byte) before any arm.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

A1 = 'SEED = 0\nEMISSION_SEED_OFFSET = 1_000_003\n'
B1 = ('SEED = 0\nEMISSION_SEED_OFFSET = 1_000_003\n'
      '# `VCC_EMISSION_SEED` (537, 2026-09-30): re-roll ONLY the count-emission stream (library picks unchanged), to read the\n'
      '# board\'s re-roll noise of one design. Unset = SEED + EMISSION_SEED_OFFSET, the stream every build has used.\n'
      'EMISSION_SEED_ENV = os.environ.get("VCC_EMISSION_SEED", "").strip()\n'
      'if EMISSION_SEED_ENV and not EMISSION_SEED_ENV.isdigit():\n'
      '    raise SystemExit(f"VCC_EMISSION_SEED must be a non-negative integer when set: {EMISSION_SEED_ENV!r}")\n'
      'EMISSION_SEED = int(EMISSION_SEED_ENV) if EMISSION_SEED_ENV else SEED + EMISSION_SEED_OFFSET\n')
A2 = '        np.random.default_rng(seed + EMISSION_SEED_OFFSET),\n'
B2 = '        np.random.default_rng(EMISSION_SEED if seed == SEED else seed + EMISSION_SEED_OFFSET),\n'
A3 = '            "emission_seed": SEED + EMISSION_SEED_OFFSET,\n'
B3 = '            "emission_seed": EMISSION_SEED,\n'


def patch_text(src: str) -> str:
    for a, b in ((A1, B1), (A2, B2), (A3, B3)):
        n = src.count(a)
        if n != 1:
            raise SystemExit(f"anchor found {n} times (need exactly 1): {a!r}")
        src = src.replace(a, b)
    return src


SELFTEST = r'''
import importlib.util, os, sys, numpy as np
def load(path, env):
    for k in ("VCC_EMISSION_SEED",):
        os.environ.pop(k, None)
    os.environ.update(env)
    src = open(path).read()
    # execute only the constants and make_rngs (the module's full import is heavy); both come before any use
    i0 = src.index("SEED = 0\n"); i1 = src.index("GIB = 1024**3")
    j0 = src.index("def make_rngs("); j1 = src.index("def assert_frozen_contract")
    ns = {"os": os, "np": np}
    exec(src[i0:i1], ns); exec(src[j0:j1], ns)
    return ns
orig, new = sys.argv[1], sys.argv[2]
o = load(orig, {}); n0 = load(new, {}); n1 = load(new, {"VCC_EMISSION_SEED": "1000004"})
po, eo = o["make_rngs"](o["SEED"]); pn, en = n0["make_rngs"](n0["SEED"]); pn1, en1 = n1["make_rngs"](n1["SEED"])
a = [x.random(8) for x in (po, eo, pn, en, pn1, en1)]
assert np.array_equal(a[0], a[2]) and np.array_equal(a[1], a[3]), "unset knob changed a stream"
assert np.array_equal(a[0], a[4]) and not np.array_equal(a[1], a[5]), "set knob must change only the emission stream"
assert n0["EMISSION_SEED"] == 1000003 and n1["EMISSION_SEED"] == 1000004
try:
    load(new, {"VCC_EMISSION_SEED": "x1"}); raise AssertionError("a malformed seed was accepted")
except SystemExit:
    pass
print("SELFTEST PASS: unset == HEAD streams; set changes only the emission stream; malformed refused")
'''


def main():
    path = Path(sys.argv[1])
    check = "--check" in sys.argv
    src = path.read_text()
    new = patch_text(src)
    if not check:
        path.write_text(new)
        print(f"patched {path}: {hashlib.sha256(src.encode()).hexdigest()[:16]} -> {hashlib.sha256(new.encode()).hexdigest()[:16]}")
        return
    with tempfile.TemporaryDirectory() as td:
        o, p = Path(td) / "orig.py", Path(td) / "patched.py"
        shutil.copy(path, o)
        p.write_text(new)
        subprocess.run([sys.executable, "-m", "py_compile", str(p)], check=True)
        d = subprocess.run(["diff", str(o), str(p)], capture_output=True, text=True).stdout
        hunks = [ln for ln in d.splitlines() if ln and ln[0].isdigit()]
        assert len(hunks) == 3, d
        (Path(td) / "t.py").write_text(SELFTEST)
        r = subprocess.run([sys.executable, str(Path(td) / "t.py"), str(o), str(p)], capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr.strip())
        assert r.returncode == 0
        print(f"--check OK on a copy: 3 hunks, compiles; {hashlib.sha256(src.encode()).hexdigest()[:16]} -> "
              f"{hashlib.sha256(new.encode()).hexdigest()[:16]}")


if __name__ == "__main__":
    main()
