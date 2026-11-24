#!/usr/bin/env python3
"""618 stage C: run496's driver (production emitter) on the h1-trainval surface, with the FAITHFUL tv key.

Modes (all but h1ka swap ONLY what is listed; built from run614's reference wrapper and run612's hook injection):
  h1ka      no swaps at all: must reproduce 435's H1 stats (ka618e.py; proves the wrapper inert)
  tv        evaluate323.load_reference / score -> work/ref_tv; run400.native_cap15 -> work/tv/effects.npz cap15 (hash-checked
            against work/tv/compose.json, as run400 checks 323's); CLASSES3 from the reference's control libraries
            (h1-test's own controls, so they must EQUAL H1's constants: asserted); the four CLEAN4 banks (asserted)
  tvcp      + 528's call-propensity block (k_emit528v3, block_steer "kept") with H1's lcap row (keys528.npz; the same
            control set, so the same key), as run614's k562cp mode wires it
  tvc4h     + the plain CLEAN4 hook (run612's c4kgall: the four banks' mean P_down at tau 0.05) on uncalled genes
  tvidall3  + idall3 (run612's idall: idsign612v3's sign, ungated)
  tvidkg3   + idkg3 (run612's idkg at q 0.5: idsign612v3's sign on each target's weak-key half)
Every mode writes into work/<name>/ only.
    uv run python run618e.py <mode> <name> <run496 args ...>
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODE = sys.argv.pop(1)
NAME = sys.argv.pop(1)
assert MODE in ("h1ka", "tv", "tvcp", "tvc4h", "tvidall3", "tvidkg3"), MODE
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
os.environ["RUN426_EMITTER"] = "prod"
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; WORK redirected below)

H, R, E = M.H, M.R, M.E


def _nan_to_null(v):
    if isinstance(v, float) and v != v:
        return None
    if isinstance(v, dict):
        return {k: _nan_to_null(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_nan_to_null(x) for x in v]
    return v


WORK = HERE / "work" / NAME
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = WORK
REC = {"mode": MODE, "name": NAME, "argv": sys.argv[1:]}
H1_CLASSES3 = tuple(H.CLASSES3)
if MODE != "h1ka":
    _wj = R.write_json
    R.write_json = lambda path, value: _wj(path, _nan_to_null(value))
    REFDIR = HERE / "work" / "ref_tv"
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
    TVK = HERE / "work" / "tv" / "effects.npz"

    def native_cap15_tv(builder_module):
        """run400.native_cap15 with the tv artifacts: the same integrity check against compose.json."""
        with np.load(TVK, allow_pickle=False) as z:
            comp = json.loads((TVK.parent / "compose.json").read_text())
            if builder_module.array_sha256(z["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
                raise ValueError("tv cap15 final effect hash drift vs compose.json")
            return z["cap15"].copy(), z["genes"].copy(), z["targets"].copy()

    H.native_cap15 = native_cap15_tv
    r0, _, _ = E.load_reference()
    with np.load(TVK, allow_pickle=False) as z:
        kt = list(z["targets"].astype(str))
    if list(map(str, r0["targets"])) != kt:
        raise SystemExit("tv key targets differ from the reference's (order matters: run496 indexes eff by r['targets'])")
    C4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
    assert tuple(H.PRIOR_BANKS) == C4 and tuple(M.CLEAN4) == C4, "CLEAN4 is not the four banks"
    L = np.asarray(r0["libraries"], float)
    H.CLASSES3 = tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75))
    assert H.CLASSES3 == H1_CLASSES3, f"tv control classes {H.CLASSES3} differ from H1's {H1_CLASSES3}"
    REC.update(refdir=str(REFDIR), key=str(TVK), key_sha256=R.sha(TVK), prior_banks=list(H.PRIOR_BANKS),
               classes3=[c[1] for c in H.CLASSES3], n_targets=len(kt))
    with np.load(REFDIR / "reference.npz", allow_pickle=False) as zr:
        KEEP = np.asarray(zr["keep"], bool)
        GF = zr["genes"].astype(str)
    with np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False) as zh:
        same_keep = bool(np.array_equal(np.asarray(zh["keep"], bool), KEEP) and np.array_equal(zh["genes"].astype(str), GF))
    REC["keep_equals_h1"] = same_keep
else:
    r0, _, _ = E.load_reference()
    L = np.asarray(r0["libraries"], float)
    assert tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75)) == tuple(H.CLASSES3), "H1 classes rule"
    REC.update(classes3=[c[1] for c in H.CLASSES3])

if MODE in ("tvc4h", "tvidall3", "tvidkg3"):
    def _clean4_P_full():
        """run612._clean4_P_full, copied: the four CLEAN4 banks' mean P_down on the reference axis, NaN where uncovered."""
        cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
        names = list(cp["bank_names"].astype(str))
        P = cp["P_banks"].astype(float)[[names.index(b) for b in C4]].mean(0)
        pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
        col = np.array([pos.get(str(g), -1) for g in GF])
        out = np.full(len(GF), np.nan)
        out[col >= 0] = P[col[col >= 0]]
        return out

    P4 = _clean4_P_full()
    HOOK_C4 = np.ascontiguousarray(np.where(np.isfinite(P4) & (np.abs(P4 - 0.5) >= 0.05),
                                            np.where(P4 > 0.5, -1.0, 1.0), 0.0)[KEEP])
    if MODE in ("tvidall3", "tvidkg3"):
        with np.load(ROOT / "experiments/612-key-gated-hook/work/idsign612v3.npz", allow_pickle=False) as zi:
            assert list(zi["genes"].astype(str)) == list(GF), "idsign612v3 is not on the reference gene axis"
            WEAK_SRC = np.ascontiguousarray(np.asarray(zi["sign"], np.float64)[KEEP])
    KIND = {"tvc4h": "c4kgall", "tvidall3": "idall", "tvidkg3": "idkg"}[MODE]
    _carrier = M.carrier_block

    def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
        # as run612: the design key on the full axis, gated by control.keep; mde_add_sign lives on the gated axis
        k = np.asarray(key, dtype=np.float64)[np.asarray(emit_c.keep, dtype=bool)]
        assert k.shape == HOOK_C4.shape, f"gated key {k.shape} vs hook {HOOK_C4.shape}"
        if KIND == "c4kgall":
            sign = HOOK_C4
        elif KIND == "idall":
            sign = WEAK_SRC
        else:
            a = np.abs(k)
            nz = a > 0
            cut = float(np.quantile(a[nz], 0.5)) if nz.any() else 0.0
            sign = np.ascontiguousarray(np.where((~nz) | (a < cut), WEAK_SRC, 0.0))
        M.MDE_ADD_SIGN = np.ascontiguousarray(sign)
        return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)

    M.carrier_block = carrier_block
    REC.update(hook_kind=KIND, hook_c4_n=int((HOOK_C4 != 0).sum()))

if MODE == "tvcp":
    sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
    import k_emit528v3 as K528  # noqa: E402
    with np.load(ROOT / "experiments/528-callprop-block/work/keys528.npz", allow_pickle=False) as zk:
        lc = np.asarray(zk["lcap"], np.float64)
    assert np.all(lc == lc[0]), "keys528 lcap rows are not one tiled control-based key"
    assert lc.shape[1] == int(KEEP.sum()), "lcap is not on this reference's gated axis"
    BLOCK_KEY = np.ascontiguousarray(np.tile(lc[0], (len(kt), 1)))
    CUR = [None]

    class KWrap:
        LAST_PRIOR_REC = K528.LAST_PRIOR_REC

        @staticmethod
        def rc_hybrid_block_k(*a, **kw):
            kw["block_key"] = BLOCK_KEY[CUR[0]]
            kw["block_keep"] = 0
            kw["block_steer"] = "kept"
            return K528.rc_hybrid_block_k(*a, **kw)

    M.K = KWrap
    _carrier_cp = M.carrier_block

    def carrier_block_cp(profile, libs, emit_rng, spread, key, i, emit_c):
        CUR[0] = i
        assert BLOCK_KEY.shape[1] == int(np.asarray(emit_c.keep).sum()), "key is not on the gated axis"
        return _carrier_cp(profile, libs, emit_rng, spread, key, i, emit_c)

    M.carrier_block = carrier_block_cp
    REC.update(cp_key="keys528.npz lcap row 0 (H1 controls = tv controls)", cp_emitter=str(Path(K528.__file__).name),
               block_steer="kept")
(WORK / "run618e.json").write_text(json.dumps(REC, indent=1))
print(f"[618e] {json.dumps(REC)}", flush=True)

if __name__ == "__main__":
    M.main()
