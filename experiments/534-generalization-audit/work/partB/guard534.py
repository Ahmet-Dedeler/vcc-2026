"""534 Part B -- write guard. Imported FIRST by every Part B script (before numpy, h5py or any foreign experiment module).

The mandate's write scope is `experiments/534-generalization-audit/` only. 438's instrument imports ~15 modules from other
experiments, several of which cache files under their own `work/` dirs when a cache is missing (lib438.bind writes
`bctrl-*.npy`, lib438.bank_stats writes `bankstats/*.npz`, plolo writes `bctrl-challenge-*.npy` and `fold-*.npz`,
sens415 writes `full-*.npz`, lib420 writes `pert_scale_415pp_lo0.5_hi3.npz` ...). Rather than trust an audit of every
code path, this makes any such write RAISE: every Python-level write under the repo root that is not under
`experiments/534-generalization-audit/` is refused with PermissionError. Paths outside the repo (temp dirs, caches)
are untouched. Child processes are not covered: their output paths are set explicitly by the Part B scripts.

Covered: builtins.open / io.open in any write mode (numpy save/savez/memmap, zipfile, pathlib write_text/bytes, pandas,
json dumps to files), os.open with write flags, os.mkdir/makedirs of a NEW directory, os.remove/unlink/rmdir/rename/
replace, shutil.rmtree, h5py.File in a write mode.
"""
from __future__ import annotations

import builtins
import io
import os
import shutil
import sys

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

_HERE = os.path.dirname(os.path.realpath(__file__))
ALLOW = os.path.realpath(os.path.join(_HERE, "..", ".."))          # experiments/534-generalization-audit
ROOT = os.path.realpath(os.path.join(ALLOW, "..", ".."))            # the repo
assert os.path.basename(ALLOW) == "534-generalization-audit", ALLOW
assert os.path.isfile(os.path.join(ROOT, "AGENTS.md")), ROOT
REFUSED: list[str] = []


def _inside(p: str, d: str) -> bool:
    return p == d or p.startswith(d + os.sep)


def check(path, what: str) -> None:
    if isinstance(path, int):
        return
    try:
        p = os.path.realpath(os.fspath(path))
    except TypeError:
        return
    if _inside(p, ROOT) and not _inside(p, ALLOW):
        REFUSED.append(f"{what}: {p}")
        raise PermissionError(f"guard534: refused {what} outside experiments/534-generalization-audit: {p}")


_open = builtins.open


def _guarded_open(file, mode="r", *a, **k):
    if isinstance(mode, str) and any(c in mode for c in "wax+"):
        check(file, f"open(mode={mode!r})")
    return _open(file, mode, *a, **k)


builtins.open = _guarded_open
io.open = _guarded_open

_WFLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC
_os_open = os.open


def _guarded_os_open(path, flags, *a, **k):
    if flags & _WFLAGS:
        check(path, "os.open(write)")
    return _os_open(path, flags, *a, **k)


os.open = _guarded_os_open


def _wrap(name, what, only_new=False):
    orig = getattr(os, name)

    def f(path, *a, **k):
        if not (only_new and os.path.exists(path)):
            check(path, what)
        return orig(path, *a, **k)
    setattr(os, name, f)


_wrap("mkdir", "mkdir", only_new=True)
for _n in ("remove", "unlink", "rmdir"):
    _wrap(_n, _n)
for _n in ("rename", "replace"):
    _orig = getattr(os, _n)

    def _mv(src, dst, *a, _o=_orig, _w=_n, **k):
        check(src, _w)
        check(dst, _w)
        return _o(src, dst, *a, **k)
    setattr(os, _n, _mv)
_rmtree = shutil.rmtree


def _guarded_rmtree(path, *a, **k):
    check(path, "rmtree")
    return _rmtree(path, *a, **k)


shutil.rmtree = _guarded_rmtree

try:
    import h5py as _h5py

    _h5init = _h5py.File.__init__

    def _h5_guarded(self, name, mode="r", *a, **k):
        if mode not in ("r",):
            check(name, f"h5py.File(mode={mode!r})")
        return _h5init(self, name, mode, *a, **k)
    _h5py.File.__init__ = _h5_guarded
except ImportError:                                               # pragma: no cover
    pass
