#!/usr/bin/env python3
"""661: 625's octile-flip driver (flip625.py) on 618's h1-trainval surface. A minimal copy; the differences are:
  * SURF = tv : evaluate323's reference dir -> 618/work/ref_tv, run400.native_cap15 -> 618's faithful tv key
    (hash-checked against its compose.json), CLASSES3 from ref_tv (asserted equal to H1's: the same controls), the four
    CLEAN4 banks (asserted);
    SURF = h1 : no swaps, i.e. flip625 exactly. This is the KA: OCT none must equal 612's c4kg stats bit for bit.
  * REF (the gated axis / gene names) = the surface's reference.npz (618 recorded keep_equals_h1 = true);
  * outputs -> 661/work/<SURF>_oct<OCT>/ (stats, metrics; dumps unused), never into 625's dir.
Everything else is flip625's: base = 612's c4kg (CLEAN4 q80 block prior, the hook on each target's weak-key half),
octile k = the CLEAN4 sign flipped in BOTH the block prior (run400.prior_signs_full) and the hook source.
usage: uv run python flip618.py <tv|h1> <none|0..7> base <seed> --bulk ps <MDE> --zs dn_k250_zs [--zero-repair -5]
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SURF = sys.argv.pop(1)
assert SURF in ("tv", "h1"), SURF
OCT = sys.argv.pop(1)
assert OCT == "none" or OCT in {str(k) for k in range(8)}, OCT
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
os.environ["RUN426_EMITTER"] = "prod"
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; WORK / HERE redirected below)
from vcc2026 import rc_emit_k as KPROD  # noqa: E402

EMITTER_SHA = hashlib.sha256(Path(KPROD.__file__).read_bytes()).hexdigest()
assert EMITTER_SHA.startswith("0f56209d57e8385a"), f"emitter is not the committed production 0f56209d: {EMITTER_SHA}"
assert M.K is KPROD
ODIR = HERE / "work" / f"{SURF}_oct{OCT}"         # 661: this experiment's own dir
WORK = ODIR
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = ODIR                                     # run496 dumps to M.HERE / "work/blocks" (21dd4c's trap: NOT WORK)
W618 = ROOT / "experiments/618-h1-trainval-surface/work"
REF = np.load((W618 / "ref_tv/reference.npz") if SURF == "tv" else
              (ROOT / "experiments/323-current-family-validation/work/reference.npz"), allow_pickle=False)
KEEPH = np.asarray(REF["keep"], dtype=bool)
GF = REF["genes"].astype(str)
V8 = np.load(ROOT / "experiments/625-module-readout/work/modules625_v8.npz", allow_pickle=True)
OCTMAP = dict(zip(V8["genes"].astype(str), V8["octile"].astype(int)))
OCT_FULL = np.array([OCTMAP.get(str(g), -1) for g in GF])
FLIP_FULL = (OCT_FULL == int(OCT)) if OCT != "none" else np.zeros(len(GF), bool)


def _clean4_P_full():
    """612/run537's CLEAN4 P on the reference axis (mean of the four banks; NaN where uncovered)."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in GF])
    out = np.full(len(GF), np.nan)
    out[col >= 0] = P[col[col >= 0]]
    return out


P_FULL = _clean4_P_full()
HOOK_FULL = np.where(np.isfinite(P_FULL) & (np.abs(P_FULL - 0.5) >= 0.05), np.where(P_FULL > 0.5, -1.0, 1.0), 0.0)
HOOK_FULL = np.where(FLIP_FULL, -HOOK_FULL, HOOK_FULL)                 # octile k: the hook's sign flipped too
HOOK = np.ascontiguousarray(HOOK_FULL[KEEPH])

_prior = M.H.prior_signs_full


def prior_signs_full(genes_full):
    out, cov = _prior(genes_full)
    if OCT != "none":
        g = np.asarray([str(x) for x in genes_full])
        flip = np.array([OCTMAP.get(x, -1) == int(OCT) for x in g])
        out = np.where(flip, -out, out)
    return out, cov


M.H.prior_signs_full = prior_signs_full
if SURF == "tv":
    # 618's run618e swaps, verbatim in substance
    H, R, E = M.H, M.R, M.E
    _wj = R.write_json
    def _nan_to_null(v):
        if isinstance(v, float) and v != v:
            return None
        if isinstance(v, dict):
            return {k: _nan_to_null(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [_nan_to_null(x) for x in v]
        return v
    R.write_json = lambda path, value: _wj(path, _nan_to_null(value))
    REFDIR = W618 / "ref_tv"
    _load, _score = E.load_reference, E.score
    def _with_ref(fn):
        def wrapped(*a, **k):
            old = R.WORK
            R.WORK = REFDIR
            try:
                return fn(*a, **k)
            finally:
                R.WORK = old
        return wrapped
    E.load_reference = _with_ref(_load)
    E.score = _with_ref(_score)
    TVK = W618 / "tv" / "effects.npz"
    def native_cap15_tv(builder_module):
        with np.load(TVK, allow_pickle=False) as z:
            comp = json.loads((TVK.parent / "compose.json").read_text())
            if builder_module.array_sha256(z["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
                raise ValueError("tv cap15 final effect hash drift vs compose.json")
            return z["cap15"].copy(), z["genes"].copy(), z["targets"].copy()
    _H1C3 = tuple(H.CLASSES3)
    H.native_cap15 = native_cap15_tv
    r0, _, _ = E.load_reference()
    with np.load(TVK, allow_pickle=False) as z:
        assert list(map(str, r0["targets"])) == list(z["targets"].astype(str)), "tv key / reference target order"
    assert tuple(H.PRIOR_BANKS) == ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
    L_ = np.asarray(r0["libraries"], float)
    H.CLASSES3 = tuple((1 / 3, float(np.quantile(L_, q))) for q in (0.25, 0.5, 0.75))
    assert H.CLASSES3 == _H1C3, "tv control classes differ from H1's"
REC = {"surf": SURF, "oct": OCT, "emitter_sha256": EMITTER_SHA, "n_flip_full": int(FLIP_FULL.sum()),
       "n_flip_gated": int(FLIP_FULL[KEEPH].sum()), "hook_n": int((HOOK != 0).sum()),
       "hook_sha256": hashlib.sha256(HOOK.tobytes()).hexdigest()}
_carrier = M.carrier_block


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    """612's c4kg gate, verbatim: the hook's sign on the target's weak-key half (|key| == 0 or below the median nonzero |key|)."""
    k = np.asarray(key, dtype=np.float64)[np.asarray(emit_c.keep, dtype=bool)]
    assert k.shape == HOOK.shape, f"gated key {k.shape} vs hook {HOOK.shape}"
    a = np.abs(k)
    nz = a > 0
    cut = float(np.quantile(a[nz], 0.5)) if nz.any() else 0.0
    weak = (~nz) | (a < cut)
    M.MDE_ADD_SIGN = np.ascontiguousarray(np.where(weak, HOOK, 0.0))
    return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)


M.carrier_block = carrier_block

if __name__ == "__main__":
    (WORK / "flip625.json").write_text(json.dumps(REC, indent=1))
    print(f"[625d] {json.dumps(REC)}", flush=True)
    M.main()
