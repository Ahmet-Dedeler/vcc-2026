#!/usr/bin/env python3
"""528 — REAL production-emitter runs on H1 with the CallSet's CPM block chosen by a call-propensity key.

435's harness (run496.py with RUN426_EMITTER=prod: the shipped DE construction on H1 — production rc_emit_k, CallSet(400,
2200, 1.9), CLEAN4 q80 block prior, coarse MDE; `--zs dn_k250_zs` adds the incumbent's zs pass on the same blocks) with ONE
variable: the emitter is 528's patched copy (k_emit528.py = patch_emitter528.py on production 1658ae2c) and the block key.
Nothing of 496's / 435's is written: outputs go to work/h1_<MODE>/.

MODE = cpm | <key>[-k<keep>]
  cpm      : block_key unset (the shipped block; KA: stats bit-equal to 435's production run, pre-zs and zs, same seed)
  <key>    : a key in work/keys528.npz (per target, gated axis): lcap (lc + 22 x, beta chosen on h1-trainval capture), shuf (its
             falsifier: x shuffled within expression bins), ltv / ptv (fitted on h1-trainval), l2f / p2f (517 2-fold on
             h1-test), rand (matched random control)
  -k<keep> : block_keep = the number of top-CPM block genes kept before the key fills the rest (default 0 = full swap;
             k1980 = the bottom 10 % swapped, k1650 = 25 %, k1100 = 50 %)
  -dec     : block_steer="cpm" (k_emit528v2.py): only the CALL set follows the key; the steered block (CLEAN4 signs, q80,
             kf, the MDE split) stays the CPM block -- a pure call-set composition change
  -kept    : block_steer="kept" (k_emit528v3.py): the call set follows the key; only its CPM part is steered (genes swapped
             in are called unsteered; genes swapped out are neither called nor steered)
  -n<cpm>  : interaction check with 523's size knob: CallSet(400, <cpm>, 1.9) instead of (400, 2200, 1.9)
  +mzg     : interaction check with 517's uncalled-gene hook as 527 builds it (446bj2zmz): the MDE add on complement genes
             takes meanz_ne's sign at tau 0.05 on 467's strong targets and CLEAN4's at 0.05 on the weak third (run517's
             `c4-mzg0.05`, its functions copied here; KA: `cpm+mzg` == 517's real c4-mzg0.05 run, same seed)
usage: uv run python experiments/528-callprop-block/run528.py <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 \
         --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")     # before any import that sets it (run400/run496)
os.environ["RUN426_EMITTER"] = "prod"      # run496 imports production rc_emit_k as K (tag `_prodem`); replaced below
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
MODE = sys.argv.pop(1)
m_ = re.fullmatch(r"(cpm|ltv|ptv|l2f|p2f|lcap|shuf|rand)(?:-k([0-9]+))?(-dec|-kept)?(?:-n([0-9]+))?(\+mzg)?", MODE)
assert m_, MODE
KEYNAME = m_.group(1)
KEEP = int(m_.group(2) or 0)
STEER = {"-dec": "cpm", "-kept": "kept"}.get(m_.group(3) or "", "key")
DEC = STEER != "key"
N_CPM = int(m_.group(4) or 2200)
HOOK = bool(m_.group(5))
assert not (DEC and KEYNAME == "cpm"), "-dec / -kept need a key"
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; its WORK is redirected below)
V3CHECK = os.environ.get("RUN528_EMITTER", "") == "v3"   # the final emitter (v3) on any mode, into h1_<MODE>-v3check
if STEER == "kept" or V3CHECK:
    import k_emit528v3 as K528  # noqa: E402  (+ block_steer "kept"; == v2 for key / cpm, test_emitter528v2.py)
elif STEER == "cpm":
    import k_emit528v2 as K528  # noqa: E402  (528's emitter with block_steer key / cpm; == k_emit528 for "key")
else:
    import k_emit528 as K528  # noqa: E402

WORK = HERE / "work" / (f"h1_{MODE}" + ("-v3check" if V3CHECK else ""))
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
BLOCK_KEY = None
if KEYNAME != "cpm":
    with np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False) as zr:
        REF_TARGETS = list(zr["targets"].astype(str))
    with np.load(HERE / "work/keys528.npz", allow_pickle=False) as z:
        assert list(z["targets"].astype(str)) == REF_TARGETS, "key targets are not in the reference order"
        BLOCK_KEY = np.ascontiguousarray(z[KEYNAME], dtype=np.float64)
    if KEYNAME == "rand" and KEEP == 0:
        raise SystemExit("rand needs its matched keep (keys528.json: info.rand.keep)")
CUR = [None]
REC = {"mode": MODE, "key": KEYNAME, "block_keep": KEEP, "block_steer": STEER, "n_cpm": N_CPM,
       "hook_mzg": HOOK,
       "key_sha256": None if BLOCK_KEY is None else hashlib.sha256(BLOCK_KEY.tobytes()).hexdigest(),
       "emitter_sha256": hashlib.sha256(Path(K528.__file__).read_bytes()).hexdigest()}


class KWrap:
    """run496 calls `K.rc_hybrid_block_k` and pops `K.LAST_PRIOR_REC`: both are 528's patched module, with the block key
    of the current target injected (nothing injected in `cpm` mode, so the verbatim code path runs)."""
    LAST_PRIOR_REC = K528.LAST_PRIOR_REC

    @staticmethod
    def rc_hybrid_block_k(*a, **kw):
        if BLOCK_KEY is not None:
            kw["block_key"] = BLOCK_KEY[CUR[0]]
            kw["block_keep"] = KEEP
            if DEC:
                kw["block_steer"] = STEER
        return K528.rc_hybrid_block_k(*a, **kw)


M.K = KWrap
_carrier = M.carrier_block
if N_CPM != 2200:
    M.H.CALLSET = (400, N_CPM, 1.9)
HK: dict[str, np.ndarray] = {}
if HOOK:
    # run517's `c4-mzg0.05` (copied, read-only sources): 467's strong set; Arc's P_down over h1_2025_full minus h1-test's
    # rows; 517's meanz_ne pseudo-P with Arc P_down's tau-.05 coverage on the gated axis; CLEAN4's P for the weak targets
    _REF = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    _H1T = set(_REF["targets"].astype(str))
    _KEEPH = np.asarray(_REF["keep"], dtype=bool)
    _GF = _REF["genes"].astype(str)
    _ST = json.loads((ROOT / "experiments/467-strength-gated-programme/work/strength467.json").read_text())
    _vals = np.array([_ST["per_target"][str(i)][_ST["chosen_measure"]] for i in range(len(_ST["per_target"]))], float)
    STRONG = set(np.flatnonzero(_vals >= float(np.quantile(_vals, 0.33))).tolist())
    _d = ROOT / "data/donors/h1_2025_full"
    _perts = [ln.strip() for ln in open(_d / "perts.txt")]
    _genes = [ln.strip() for ln in open(_d / "genes.txt")]
    _rows = [i for i, p in enumerate(_perts) if p not in _H1T]
    _X = np.asarray(np.load(_d / "delta.npy", mmap_mode="r")[_rows], dtype=np.float64)
    _nz, _dn = (_X != 0).sum(0), (_X < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        _pd = np.where(_nz > 0, _dn / np.maximum(_nz, 1), np.nan)
    _pos = {g: j for j, g in enumerate(_genes)}
    _col = np.array([_pos.get(str(g), -1) for g in _GF])
    _Pa = np.where(_col >= 0, _pd[np.maximum(_col, 0)], np.nan)
    _ess = {ln.strip() for ln in open(ROOT / "data/donors/k562_essential/perts.txt") if ln.strip()}
    _ne = np.array([_perts[i] not in _ess for i in _rows])
    _sd = _X.std(0)
    _sd[_sd == 0] = 1.0
    _zb = (_X[_ne] / _sd).mean(0)
    del _X
    _z = np.where(_col >= 0, _zb[np.maximum(_col, 0)], np.nan)
    _ncov = int((np.isfinite(_Pa[_KEEPH]) & (np.abs(_Pa[_KEEPH] - 0.5) >= 0.05)).sum())
    _zq = float(np.sort(np.abs(np.nan_to_num(_z[_KEEPH])))[::-1][_ncov - 1])
    _Pz = np.where(np.isfinite(_z), np.clip(0.5 - 0.05 * _z / _zq, 0.001, 0.999), np.nan)
    _cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    _names = list(_cp["bank_names"].astype(str))
    _PC = _cp["P_banks"].astype(float)[[_names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    _cpos = {g: j for j, g in enumerate(_cp["genes"].astype(str))}
    _ccol = np.array([_cpos.get(str(g), -1) for g in _GF])
    _Pc = np.full(len(_GF), np.nan)
    _Pc[_ccol >= 0] = _PC[_ccol[_ccol >= 0]]

    def _hs(P, tau):
        sg = np.where(np.isfinite(P) & (np.abs(P - 0.5) >= tau), np.where(P > 0.5, -1.0, 1.0), 0.0)
        return np.ascontiguousarray(sg[_KEEPH])
    HK["strong"], HK["weak"] = _hs(_Pz, 0.05), _hs(_Pc, 0.05)
    REC["hook_signs"] = {k: int((v != 0).sum()) for k, v in HK.items()}


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    CUR[0] = i
    if BLOCK_KEY is not None:
        assert BLOCK_KEY.shape[1] == int(np.asarray(emit_c.keep).sum()), "key is not on the gated axis"
    if HOOK:
        M.MDE_ADD_SIGN = HK["strong"] if i in STRONG else HK["weak"]
    return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)


M.carrier_block = carrier_block

if __name__ == "__main__":
    (WORK / "run528.json").write_text(json.dumps(REC, indent=1))
    print(f"[528] {json.dumps(REC)}", flush=True)
    M.main()
