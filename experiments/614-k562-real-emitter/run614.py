#!/usr/bin/env python3
"""614 stage 5: run496's driver (production emitter) on a harness reference other than H1's.
Swaps, all made BEFORE run496.main() and only in the non-KA modes:
  - evaluate323.load_reference / score read the harness reference dir (wrapped; run323.WORK itself is never repointed,
    because the builder's production inputs live there);
  - run400.native_cap15 -> the harness key (k562: work/key614.npz; its targets must equal the reference's, in order);
  - run400.PRIOR_BANKS and run496.CLEAN4 -> CLEAN4 minus both K562 banks (k562 mode);
  - run400.CLASSES3 -> the 25/50/75 % quantiles of the reference's control libraries (H1's constants are exactly its own).
Every mode sends run496's WORK / HERE (dumps) to work/<name>/.
Mode `h1ka` swaps nothing else, and must reproduce 435's H1 stats (the wrapper's known answer; ka614.py).
    uv run python run614.py <k562|h1ka> <name> <run496 args ...>
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
assert MODE in ("k562", "h1ka", "k562cp", "h1t11k", "k562c4h", "k562idall3", "k562f", "k562fcp"), MODE
# k562f / k562fcp (TRUSTGATE faithful-key re-run, 3d851ae7): identical to k562 / k562cp except the key, which is
# fkey614.py's cap15 (work/fkey/effects.npz, hash-checked against its compose.json as run400 checks 323's)
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


if MODE != "h1ka":
    # K562 targets can lack a key run496 summarises (39 % have no `mde_lambda`), so its metrics median is NaN and
    # run323.write_json (allow_nan=False) raises after the run. Write NaN as null; stats and raw members are untouched.
    _wj = R.write_json
    R.write_json = lambda path, value: _wj(path, _nan_to_null(value))
WORK = HERE / "work" / NAME
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = WORK
REC = {"mode": MODE, "name": NAME, "argv": sys.argv[1:]}
if MODE in ("k562", "k562cp", "k562c4h", "k562idall3", "k562f", "k562fcp"):
    REFDIR = HERE / "work" / "ref_k562s"      # amendment A1: the truth with the H1-shaped sub-gate reservoir (sub614.py)
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
    if MODE in ("k562f", "k562fcp"):
        _fk = HERE / "work/fkey/effects.npz"
        with np.load(_fk, allow_pickle=False) as z:
            _comp = json.loads((_fk.parent / "compose.json").read_text())
            import hashlib as _hl
            # the builder's array_sha256 (canonical <f8 bytes), checked against compose.json as run400 checks 323's
            if _hl.sha256(np.ascontiguousarray(z["cap15"], dtype="<f8").tobytes()).hexdigest() != _comp["arms"]["cap15"]["final_sha256"]:
                raise SystemExit("fkey cap15 hash drift vs compose.json")
            KEY = (np.asarray(z["cap15"], np.float64), z["genes"].astype(str), z["targets"].astype(str))
        KEYFILE = _fk
    else:
        with np.load(HERE / "work/key614.npz", allow_pickle=False) as z:
            KEY = (np.asarray(z["eff"], np.float64), z["genes"].astype(str), z["targets"].astype(str))
        KEYFILE = HERE / "work/key614.npz"
    r0, _, _ = E.load_reference()
    if list(map(str, r0["targets"])) != list(KEY[2]):
        raise SystemExit("key614 targets differ from the reference's (order matters: run496 indexes eff by r['targets'])")
    H.native_cap15 = lambda b: (KEY[0].copy(), KEY[1].copy(), KEY[2].copy())
    H.PRIOR_BANKS = ("nadig_jurkat", "nadig_hepg2")
    M.CLEAN4 = H.PRIOR_BANKS
    L = np.asarray(r0["libraries"], float)
    H.CLASSES3 = tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75))
    REC.update(refdir=str(REFDIR), key_file=str(KEYFILE), key_sha256=R.sha(KEYFILE), prior_banks=list(H.PRIOR_BANKS),
               classes3=[c[1] for c in H.CLASSES3], n_targets=int(len(KEY[2])))
elif MODE == "h1t11k":
    # descriptive control: H1 thinned to K562's depth (thin614.py); H1's own key and all four CLEAN4 banks, so the only
    # changes vs H1@20k are the reference loader and the library classes, both from the thinned truth
    REFDIR = HERE / "work" / "ref_h1t11k"
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
    r0, _, _ = E.load_reference()
    L = np.asarray(r0["libraries"], float)
    H.CLASSES3 = tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75))
    REC.update(refdir=str(REFDIR), classes3=[c[1] for c in H.CLASSES3])
else:
    r0, _, _ = E.load_reference()
    L = np.asarray(r0["libraries"], float)
    assert tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75)) == tuple(H.CLASSES3), "H1 classes rule"
    REC.update(classes3=[c[1] for c in H.CLASSES3])
if MODE in ("k562c4h", "k562idall3"):
    # hook forms (PRED-614-levers.md): the MDE add's direction on uncalled genes from idsign614.npz, ungated, as run612
    # injects it (M.MDE_ADD_SIGN set before each target's carrier block; the production emitter reads it)
    SRCN = {"k562c4h": "c4h", "k562idall3": "idall3"}[MODE]
    with np.load(HERE / "work/idsign614.npz", allow_pickle=False) as zs_:
        assert list(zs_["genes"].astype(str)) == list(r0["genes"].astype(str)), "idsign614 is not on the reference axis"
        # the gated (DE-universe) axis: reference.npz's `keep` (= the emitter's control.keep; ref614 asserts they are
        # equal). E.load_reference() does not return it (the first hook run died here, 17:55Z, KeyError 'keep').
        with np.load(REFDIR / "reference.npz", allow_pickle=False) as zr_:
            KEEP_REF = np.asarray(zr_["keep"], bool)
        HOOK_SIGN = np.ascontiguousarray(np.asarray(zs_[SRCN], np.float64)[KEEP_REF])
    _carrier_h = M.carrier_block

    def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
        assert HOOK_SIGN.shape == (int(np.asarray(emit_c.keep).sum()),), "hook sign is not on the gated axis"
        M.MDE_ADD_SIGN = HOOK_SIGN
        return _carrier_h(profile, libs, emit_rng, spread, key, i, emit_c)

    M.carrier_block = carrier_block
    REC.update(hook_source=SRCN, hook_sha256=R.sha(HERE / "work/idsign614.npz"), hook_n_sign=int((HOOK_SIGN != 0).sum()))
if MODE in ("k562cp", "k562fcp"):
    # G3: 528's call-propensity block (lcap, block_steer "kept") on the K562 harness: run528's KWrap / carrier_block
    # pattern verbatim, with 528's final emitter (k_emit528v3) and the K562 key (key614cp.py; lcap KA'd on H1)
    sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
    import k_emit528v3 as K528  # noqa: E402
    with np.load(HERE / "work/keys614cp.npz", allow_pickle=False) as zk:
        assert list(zk["targets"].astype(str)) == list(KEY[2]), "cp key targets differ from the reference order"
        BLOCK_KEY = np.ascontiguousarray(zk["lcap"], dtype=np.float64)
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
    _carrier = M.carrier_block

    def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
        CUR[0] = i
        assert BLOCK_KEY.shape[1] == int(np.asarray(emit_c.keep).sum()), "key is not on the gated axis"
        return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)

    M.carrier_block = carrier_block
    REC.update(cp_key_sha256=R.sha(HERE / "work/keys614cp.npz"), cp_emitter=str(Path(K528.__file__).name), block_steer="kept")
(WORK / "run614.json").write_text(json.dumps(REC, indent=1))
print(f"[614] {json.dumps(REC)}", flush=True)

if __name__ == "__main__":
    M.main()
