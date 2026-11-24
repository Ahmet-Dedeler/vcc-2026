#!/usr/bin/env python3
"""616: run an existing instrument script unchanged, inside a process-local sandbox.

  * READ redirects: an explicit map {original file -> replacement} (the clean Orion gate tables), plus "sandbox first":
    a file the sandboxed run itself wrote earlier (its mirror exists) is read back from the mirror. So a chain of
    instruments (decomp639 -> block613's KA, ds612b -> ds612c) sees its own upstream outputs, not the committed ones.
  * WRITE redirects: every write to a path under the repo but outside experiments/616-clean-orion-folds goes to
    work/sandbox/<tag>/<same relative path>. No other experiment's file is ever touched.
Hooks: numpy.load / save / savez / savez_compressed, io.open (= builtins.open, used by pathlib's read_text /
write_text and by pandas). Every redirect is logged to work/sandbox/<tag>/redirects.json.

    uv run python sandbox616.py <tag> <clean|identity> <script.py> [script args ...]

`identity` maps nothing (writes are still sandboxed). It is the known-answer mode: each instrument's sandboxed output
must equal its committed output bit for bit (ka616.py). `clean` maps the two Orion gate tables to work/gate-*_clean.npz.
"""
from __future__ import annotations

import builtins
import io
import json
import os
import runpy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TAG, MODE, SCRIPT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.dont_write_bytecode = True                      # no .pyc written into other experiments' __pycache__
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
assert MODE in ("clean", "identity", "map"), MODE
SANDBOX = HERE / "work" / "sandbox" / TAG
SANDBOX.mkdir(parents=True, exist_ok=True)
G480 = ROOT / "experiments/480-nmae-instrument/work"
MAP: dict[str, str] = {}
if MODE == "clean" and not (HERE / "work" / "RELEASE-S2").exists():
    # coordinator 19:06Z: no clean-mode RECOMPUTATION until the list of Orion-dependent readings exists and is re-ruled;
    # the clean gate tables (S1) are still built. Created by hand when S2 is released.
    print(f"[sandbox616] {TAG} clean: HELD (work/RELEASE-S2 absent); nothing run", flush=True)
    sys.exit(3)
if MODE == "map":
    # 618: an explicit read map from a json file {original repo path: replacement path} (both relative to the repo root or
    # absolute), e.g. h1-test's reference / key / base stats / gate table -> the h1-trainval surface's
    _m = json.loads(Path(os.environ["SANDBOX616_MAP"]).read_text())
    for a, b in _m.items():
        a_, b_ = (str(Path(x) if Path(x).is_absolute() else ROOT / x) for x in (a, b))
        assert os.path.exists(b_), b_
        MAP[a_] = b_
if MODE == "clean":
    for name in ("orion_hek293t", "orion_hct116"):
        rep = HERE / f"work/gate-{name}_clean.npz"
        assert rep.exists(), rep
        MAP[str(G480 / f"gate-{name}.npz")] = str(rep)
LOG: list[tuple[str, str, str]] = []


def _abs(p) -> str | None:
    if isinstance(p, (str, os.PathLike)):
        return os.path.abspath(os.fspath(p))
    return None


def _mirror(s: str) -> str | None:
    if s.startswith(str(ROOT) + os.sep) and not s.startswith(str(HERE) + os.sep):
        return str(SANDBOX / os.path.relpath(s, ROOT))
    return None


TRACE = os.environ.get("SANDBOX616_TRACE") == "1"   # log EVERY repo-file read (used to scope 618's 613d re-check)


def _read(p):
    s = _abs(p)
    if s is None:
        return p
    if TRACE and s.startswith(str(ROOT) + os.sep):
        LOG.append(("read-any", s, ""))
    if s in MAP:
        LOG.append(("read-map", s, MAP[s]))
        return MAP[s]
    m = _mirror(s)
    if m and os.path.exists(m):
        LOG.append(("read-sandbox", s, m))
        return m
    return p


def _write(p):
    s = _abs(p)
    if s is None:
        return p
    m = _mirror(s)
    if m:
        os.makedirs(os.path.dirname(m), exist_ok=True)
        LOG.append(("write", s, m))
        return m
    return p


_open = io.open


def open_(file, mode="r", *a, **k):
    w = any(c in mode for c in "wax+")
    return _open(_write(file) if w else _read(file), mode, *a, **k)


io.open = open_
builtins.open = open_
_load, _save, _savez, _savezc = np.load, np.save, np.savez, np.savez_compressed
np.load = lambda file, *a, **k: _load(_read(file), *a, **k)
np.save = lambda file, *a, **k: _save(_write(file), *a, **k)
np.savez = lambda file, *a, **k: _savez(_write(file), *a, **k)
np.savez_compressed = lambda file, *a, **k: _savezc(_write(file), *a, **k)

script = os.path.abspath(SCRIPT)
sys.argv = [script] + sys.argv[4:]
sys.path.insert(0, os.path.dirname(script))
rc = 0
try:
    runpy.run_path(script, run_name="__main__")
except SystemExit as e:
    rc = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
finally:
    io.open = _open
    builtins.open = _open
    hits = {"tag": TAG, "mode": MODE, "script": os.path.relpath(script, ROOT), "map": MAP, "rc": rc,
            "redirects": sorted(set(LOG))}
    with _open(SANDBOX / "redirects.json", "w") as f:
        json.dump(hits, f, indent=1)
    n_map = sum(1 for x in set(LOG) if x[0] == "read-map")
    print(f"[sandbox616] {TAG} {MODE}: rc {rc}; mapped reads {n_map}; sandbox reads "
          f"{sum(1 for x in set(LOG) if x[0] == 'read-sandbox')}; writes {sum(1 for x in set(LOG) if x[0] == 'write')}",
          flush=True)
sys.exit(rc)
