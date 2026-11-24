"""801 -- write guard + data-closure recorder. A COPY of experiments/534-generalization-audit/work/partB/guard534.py (git blob
545769bbdf07a09f4fe1ae336e5c8acf4ae7de30) with ALLOW pointed at experiments/801-lab-folds and one addition at the end:
`sys.addaudithook` records every file under the repo that the process opens (reads included) into CLOSURE; set the env
GUARD801_CLOSURE=<path> to have the sorted list written there at exit (that list is what the cloud runs stage).
guard534's docstring follows.

534 Part B -- write guard. Imported FIRST by every Part B script (before numpy, h5py or any foreign experiment module).
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
ALLOW = os.path.realpath(_HERE)                                     # experiments/801-lab-folds
ROOT = os.path.realpath(os.path.join(ALLOW, "..", ".."))            # the repo
assert os.path.basename(ALLOW) == "801-lab-folds", ALLOW
assert os.path.isfile(os.path.join(ROOT, "pyproject.toml")), ROOT          # 801: the cloud tarball carries no AGENTS.md
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
        raise PermissionError(f"guard801: refused {what} outside experiments/801-lab-folds: {p}")


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


if os.environ.get("GUARD801_ALLOW_MKDIR") != "1":           # 801: on EC2 only (cloud801.sh sets it): foreign modules
    _wrap("mkdir", "mkdir", only_new=True)                    # mkdir their work dirs at import; an empty dir is harmless
                                                              # there, and every FILE write outside 801 still raises
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

# ---------------------------------------------------------------------------------------------- 801: the data closure
import atexit as _atexit  # noqa: E402

CLOSURE: set[str] = set()


def _audit(event, args):
    if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        try:
            p = os.path.realpath(os.fsdecode(args[0]))
        except (TypeError, ValueError):
            return
        if _inside(p, ROOT) and os.path.isfile(p):
            CLOSURE.add(os.path.relpath(p, ROOT))


sys.addaudithook(_audit)


def _dump_closure():
    out = os.environ.get("GUARD801_CLOSURE")
    if out:
        with _open(out, "w") as fh:
            fh.write("\n".join(sorted(CLOSURE)) + "\n")


_atexit.register(_dump_closure)

try:                                     # h5py opens through the HDF5 C library, which the audit hook never sees
    import h5py as _h5py2

    _h5init2 = _h5py2.File.__init__

    def _h5_record(self, name, mode="r", *a, **k):
        try:
            p = os.path.realpath(os.fspath(name))
            if _inside(p, ROOT) and os.path.isfile(p):
                CLOSURE.add(os.path.relpath(p, ROOT))
        except TypeError:
            pass
        return _h5init2(self, name, mode, *a, **k)
    _h5py2.File.__init__ = _h5_record
except ImportError:                                               # pragma: no cover
    pass
