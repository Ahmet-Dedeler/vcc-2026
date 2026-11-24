#!/usr/bin/env python3
"""615 stage 5: run496's driver (production emitter) on the Orion raw HEK293T harness (a copy of 614's run614.py;
the hek modes differ from 614's k562 modes ONLY in the reference dir, the key file, and keeping all four CLEAN4 banks).
Swaps, all made BEFORE run496.main() and only in the non-KA modes:
  - evaluate323.load_reference / score read the harness reference dir (wrapped; run323.WORK itself is never repointed,
    because the builder's production inputs live there);
  - run400.native_cap15 -> the harness key (work/key615.npz: k562_gwps x s; its targets must equal the reference's);
  - run400.PRIOR_BANKS / run496.CLEAN4 are NOT changed: none of the four CLEAN4 banks is Orion-derived (asserted);
  - run400.CLASSES3 -> the 25/50/75 % quantiles of the reference's control libraries (H1's constants are exactly its own).
Every mode sends run496's WORK / HERE (dumps) to work/<name>/.
Mode `h1ka` swaps nothing else, and must reproduce 435's H1 stats (the wrapper's known answer; ka614.py).
    uv run python run615.py <hek|hekcp|hekc4h|hekidall3|h1ka> <name> <run496 args ...>
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
assert MODE in ("hek", "h1ka", "hekcp", "hekc4h", "hekidall3"), MODE
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
    # harness targets can lack a key run496 summarises (on K562, 39 % had no `mde_lambda`), so its metrics median is NaN and
    # run323.write_json (allow_nan=False) raises after the run. Write NaN as null; stats and raw members are untouched.
    _wj = R.write_json
    R.write_json = lambda path, value: _wj(path, _nan_to_null(value))
WORK = HERE / "work" / NAME
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = WORK
REC = {"mode": MODE, "name": NAME, "argv": sys.argv[1:]}
if MODE in ("hek", "hekcp", "hekc4h", "hekidall3"):
    REFDIR = HERE / "work" / "ref_hek"        # the P' truth (truth615.py); the full 18,533 axis has its own sub-gate genes
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
    with np.load(HERE / "work/key615.npz", allow_pickle=False) as z:
        KEY = (np.asarray(z["eff"], np.float64), z["genes"].astype(str), z["targets"].astype(str))
    r0, _, _ = E.load_reference()
    if list(map(str, r0["targets"])) != list(KEY[2]):
        raise SystemExit("key615 targets differ from the reference's (order matters: run496 indexes eff by r['targets'])")
    H.native_cap15 = lambda b: (KEY[0].copy(), KEY[1].copy(), KEY[2].copy())
    C4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
    assert tuple(H.PRIOR_BANKS) == C4 and tuple(M.CLEAN4) == C4, "CLEAN4 is not the four Orion-free banks"
    L = np.asarray(r0["libraries"], float)
    H.CLASSES3 = tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75))
    REC.update(refdir=str(REFDIR), key_sha256=R.sha(HERE / "work/key615.npz"), prior_banks=list(H.PRIOR_BANKS),
               classes3=[c[1] for c in H.CLASSES3], n_targets=int(len(KEY[2])))
else:
    r0, _, _ = E.load_reference()
    L = np.asarray(r0["libraries"], float)
    assert tuple((1 / 3, float(np.quantile(L, q))) for q in (0.25, 0.5, 0.75)) == tuple(H.CLASSES3), "H1 classes rule"
    REC.update(classes3=[c[1] for c in H.CLASSES3])
if MODE in ("hekc4h", "hekidall3"):
    # hook forms: the MDE add's direction on uncalled genes from idsign615.npz, ungated, as run612
    # injects it (M.MDE_ADD_SIGN set before each target's carrier block; the production emitter reads it)
    SRCN = {"hekc4h": "c4h", "hekidall3": "idall3"}[MODE]
    with np.load(HERE / "work/idsign615.npz", allow_pickle=False) as zs_:
        assert list(zs_["genes"].astype(str)) == list(r0["genes"].astype(str)), "idsign615 is not on the reference axis"
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
    REC.update(hook_source=SRCN, hook_sha256=R.sha(HERE / "work/idsign615.npz"), hook_n_sign=int((HOOK_SIGN != 0).sum()))
if MODE == "hekcp":
    # G3: 528's call-propensity block (lcap, block_steer "kept") on the HEK harness: run528's KWrap / carrier_block
    # pattern verbatim, with 528's final emitter (k_emit528v3) and the HEK key (key615cp.py; lcap KA'd on H1)
    sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
    import k_emit528v3 as K528  # noqa: E402
    with np.load(HERE / "work/keys615cp.npz", allow_pickle=False) as zk:
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
    REC.update(cp_key_sha256=R.sha(HERE / "work/keys615cp.npz"), cp_emitter=str(Path(K528.__file__).name), block_steer="kept")
(WORK / "run615.json").write_text(json.dumps(REC, indent=1))
print(f"[615] {json.dumps(REC)}", flush=True)

if __name__ == "__main__":
    M.main()
