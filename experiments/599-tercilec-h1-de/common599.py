"""599 — shared pieces: the H1 key (410's construction, cached), the coverage-thinned carried set, the energy-preserving bulk
at a clip, and the per-target q of each arm from 598's patched sharpen_q_by_tercile.

Modes (one variable between t and tc: the tier rule; s175 anchors both on the shipped scalar render):
  ka   : 410b's ps_bulk exactly (carried as fitted, clip 0.3 / 3), no sharpening -> must reproduce 435's base run
  s175 : thinned carried, clip 0.5 / 6 (October-like), scalar q 1.75 on every target (sharpen_pt off)
  t    : thinned, clip 0.5 / 6, `tercile:2.5,2.0,1.75`  (the DECIDED October value)
  tc   : thinned, clip 0.5 / 6, `tercilec:2.5,2.0,1.75` (598)
Thinning: a seeded random subset of the fitted carried targets loses its key (w = 1), so that 43 of the 85 targets are
non-carried (efp: ~153-160 of ~300); the median is re-taken over the remaining carried targets (key410.weights).
"""
from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
KEYFILE = WORK / "key599.npz"
BUILDER598 = ROOT / "experiments/598-tercilec-knob/work/builder.post598.py"
BUILDER598_SHA = "5cf0e0124452e1b158770ab28a8d1fdb681d3dfeedc5b8d926237832a03a7a68"
THIN_SEED, N_NONCARRIED = 599, 43
Q3 = [2.5, 2.0, 1.75]
MODES = {"ka": dict(thin=False, lo=0.3, hi=3.0, q=None),
         "s175": dict(thin=True, lo=0.5, hi=6.0, q="scalar:1.75"),
         "t": dict(thin=True, lo=0.5, hi=6.0, q="tercile"),
         "tc": dict(thin=True, lo=0.5, hi=6.0, q="tercilec")}
for _p in ("experiments/410-pertscale-revisit", "experiments/162-retrodicting-validator"):
    sys.path.insert(0, str(ROOT / _p))


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def tercile_fn():
    """598's patched sharpen_q_by_tercile (the applied text, 5cf0e012), extracted by ast; self-contained (numpy only)."""
    assert sha256(BUILDER598) == BUILDER598_SHA, "598's patched builder text moved"
    text = BUILDER598.read_text()
    ns: dict = {"np": np}
    for node in ast.parse(text).body:
        if isinstance(node, ast.FunctionDef) and node.name == "sharpen_q_by_tercile":
            exec(compile(ast.get_source_segment(text, node), "<builder598>", "exec"), ns)
    return ns["sharpen_q_by_tercile"]


def fit_key(targets) -> dict:
    """run496.ps_bulk's key, verbatim (key410.shared_energy_key on the H1 operator, run496.DONORS); cached in KEYFILE."""
    if KEYFILE.exists():
        z = np.load(KEYFILE, allow_pickle=False)
        assert list(map(str, z["targets"])) == list(map(str, targets)), "cached key is for other targets"
        return {k: z[k] for k in z.files}
    import key410  # noqa: PLC0415
    sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
    import run496 as M  # noqa: PLC0415  (H1Op and DONORS; importing it has no side effects beyond its WORK dir)
    op = M.H1Op()
    key, carried, nd = key410.shared_energy_key(list(map(str, targets)), M.DONORS, op, op.cols, cols=None)
    del op
    WORK.mkdir(parents=True, exist_ok=True)
    np.savez(KEYFILE, targets=np.asarray(list(map(str, targets))), key=key, carried=carried, n_donors=nd)
    return {"targets": np.asarray(list(map(str, targets))), "key": key, "carried": carried, "n_donors": nd}


def thinned(carried: np.ndarray) -> np.ndarray:
    """Drop a seeded random subset of the carried targets so that N_NONCARRIED of them end non-carried."""
    carried = np.asarray(carried, dtype=bool)
    need = N_NONCARRIED - int((~carried).sum())
    assert need >= 0, f"already {int((~carried).sum())} non-carried"
    drop = np.random.default_rng(THIN_SEED).choice(np.flatnonzero(carried), need, replace=False)
    out = carried.copy()
    out[drop] = False
    return out


def mode_bulk(mode: str, eff, genes, targets, kz: dict):
    """(eff_bulk, w, carried, c) for a mode: key410.weights at the mode's clip, energy-preserving over non-panel genes
    (run496.ps_bulk's arithmetic verbatim)."""
    import key410  # noqa: PLC0415
    spec = MODES[mode]
    carried = thinned(kz["carried"]) if spec["thin"] else np.asarray(kz["carried"], dtype=bool)
    w, med = key410.weights(np.asarray(kz["key"], dtype=np.float64), carried, spec["lo"], spec["hi"])
    panel = set(map(str, targets))
    nonp = np.array([str(g) not in panel for g in genes])
    e2 = np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp])
    c = float(np.sqrt(e2.sum() / (w * w * e2).sum()))
    return eff * (w * c)[:, None], w, carried, c


def mode_q(mode: str, eff_bulk, eff, genes, targets, carried):
    """Per-target q (None = no sharpening) and the tier audit."""
    q = MODES[mode]["q"]
    if q is None:
        return None, {"mode": "none"}
    if q.startswith("scalar:"):
        v = float(q.split(":", 1)[1])
        return np.full(len(targets), v), {"mode": q}
    fn = tercile_fn()
    kw = {"carried": carried} if q == "tercilec" else {}
    qt, audit = fn(eff_bulk, eff, list(map(str, genes)), list(map(str, targets)), {"mode": q, "q": list(Q3)}, **kw)
    return np.asarray(qt, dtype=np.float64), audit


def tiers_of(qt) -> np.ndarray:
    """q -> tier index (0 LOW 2.5, 1 MID 2.0, 2 HIGH 1.75)."""
    return np.array([Q3.index(float(x)) for x in qt])


def dump_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) + "\n")
