#!/usr/bin/env python3
"""537 — real production-emitter H1 runs for the DE combinations the 10-01 … 10-06 board reads open.

435's harness (run496.py with RUN426_EMITTER=prod: production rc_emit_k 0f56209d, CallSet(400, 2200, 1.9), CLEAN4 q80 block
prior, coarse MDE; `--zs dn_k250_zs` adds the incumbent's zs pass on the same blocks) = 446bj2z's DE construction on H1,
with the variables of 528 (call-propensity block, call-set size) and 496/517/527 (the uncalled-gene hook) combined:

MODE = <key>[-k<keep>][-kept][-n<cpm>][+c4h|+mzg][-r<r>][-s][-bs<g>]
  <key>    : cpm (the shipped CPM block, block_key unset: the verbatim code path) | lcap (528's key: log10 control CPM +
             22 x excess dropout, work/keys528.npz of experiments/528-callprop-block)
  -k<keep> : block_keep (top-CPM block genes kept before the key fills the rest; default 0 = full swap)
  -kept    : block_steer "kept" (528's shipped form: only the call set's CPM part keeps CLEAN4 steering)
  -n<cpm>  : CallSet(400, <cpm>, 1.9) (523's VCC_RC_CALLSET_CPM); default 2200
  +c4h[τ]  : 496's CLEAN4 hook, tau τ (default 0.05 = 496a's VCC_DE_MDE_SIGN=clean4:0.05; 517's c4t sign source), every target
  +mzg     : 527's hook (446bj2zmz): meanz_ne sign at tau .05 on 467's strong targets, CLEAN4 .05 on the weak third
  -r<r>    : attenuate the hook's source: flip its sign on a fixed random fraction r of the gated axis, u = default_rng(51717)
             (mzgrid517's emulation; coverage unchanged); +mzg flips z, +c4h flips P -> 1 - P
  -s       : the zero-information floor for an unseen line: permute the hook's signs among its covered genes on the gated
             axis (coverage set, count and up/down balance kept, gene identity destroyed), permutation default_rng(51718 +
             seed), i.e. a different one per emission seed (+c4h: added 2026-09-30 15:50Z; +mzg: 18:55Z, the strong source's
             z among |z| >= its tau-.05 cut and the weak third's CLEAN4 P, each permuted with the same generator)
  -bs<g>   : a bulk-only step (§2's re-draw test, added 16:10Z): the bulk's per-target scale s_i (410b's key-scaled bulk,
             `--bulk ps`) raised to the power g and re-normalised to the same non-panel energy; the DE design is untouched,
             so only the ballast (its K per target, hence the shared emission stream) can move pass 2
The emitter is the production module (src/vcc2026/rc_emit_k.py, sha asserted 0f56209d…), i.e. byte-identical to 528's
k_emit528v3.py; block_key / block_keep / block_steer are injected per target exactly as run528.py does.
Known answers (drive537 checks them before any new number is read):
  cpm        seed 0 == 528's h1_cpm stats (== 435's production run), pre-zs and zs
  lcap-kept  seed 0 == 528's h1_lcap-kept stats, pre-zs and zs
  cpm+c4h    seed 0 pre-zs == 517's h1_c4-c4t0.05 stats (== 496's real hook run)
  cpm+mzg    seed 0 == 528's h1_cpm+mzg stats
Outputs: experiments/537-arm-pipeline/work/h1_<MODE>/ only.
usage: uv run python experiments/537-arm-pipeline/run537.py <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 \
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
os.environ["RUN426_EMITTER"] = "prod"      # run496 imports the production rc_emit_k as K (tag `_prodem`)
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
MODE = sys.argv.pop(1)
m_ = re.fullmatch(r"(cpm|lcap)(?:-k([0-9]+))?(-kept)?(?:-n([0-9]+))?(\+c4h|\+mzg)?([0-9]\.[0-9]+)?(?:-r([0-9.]+))?(-s)?(?:-bs([0-9.]+))?", MODE)
assert m_, MODE
KEYNAME = m_.group(1)
KEEP = int(m_.group(2) or 0)
STEER = "kept" if m_.group(3) else "key"
N_CPM = int(m_.group(4) or 2200)
HOOK = (m_.group(5) or "").lstrip("+") or None
HOOK_TAU = float(m_.group(6) or 0.05)          # +c4h<tau>: the CLEAN4 hook's confidence threshold (default 0.05 = 496a's)
R_FLIP = float(m_.group(7) or 0.0)
assert not (m_.group(6) and HOOK != "c4h"), "a tau suffix is defined for +c4h only"
assert not (STEER == "kept" and KEYNAME == "cpm"), "-kept needs a key"
assert not (R_FLIP > 0 and HOOK not in ("mzg", "c4h")), "-r attenuates a hook source only"
SHUF = bool(m_.group(8))
assert not (SHUF and HOOK not in ("c4h", "mzg")), "-s shuffles a hook source only"
BULK_G = float(m_.group(9)) if m_.group(9) else None
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; its WORK is redirected below)
from vcc2026 import rc_emit_k as KPROD  # noqa: E402

EMITTER_SHA = hashlib.sha256(Path(KPROD.__file__).read_bytes()).hexdigest()
assert EMITTER_SHA.startswith("0f56209d57e8385a"), f"emitter is not the committed production 0f56209d: {EMITTER_SHA}"
assert M.K is KPROD, "run496 did not import the production emitter"

WORK = HERE / "work" / f"h1_{MODE}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
REF = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
REF_TARGETS = list(REF["targets"].astype(str))
KEEPH = np.asarray(REF["keep"], dtype=bool)
GF = REF["genes"].astype(str)
BLOCK_KEY = None
if KEYNAME != "cpm":
    with np.load(ROOT / "experiments/528-callprop-block/work/keys528.npz", allow_pickle=False) as z:
        assert list(z["targets"].astype(str)) == REF_TARGETS, "key targets are not in the reference order"
        BLOCK_KEY = np.ascontiguousarray(z[KEYNAME], dtype=np.float64)
CUR = [None]
REC = {"mode": MODE, "key": KEYNAME, "block_keep": KEEP, "block_steer": STEER, "n_cpm": N_CPM, "hook": HOOK,
       "hook_tau": HOOK_TAU if HOOK == "c4h" else None, "r_flip": R_FLIP, "sign_shuffle": SHUF, "bulk_step_g": BULK_G,
       "key_sha256": None if BLOCK_KEY is None else hashlib.sha256(BLOCK_KEY.tobytes()).hexdigest(),
       "emitter_sha256": EMITTER_SHA}


class KWrap:
    """run496 calls `K.rc_hybrid_block_k` and pops `K.LAST_PRIOR_REC`: the production module, with the block key of the
    current target injected (nothing injected in `cpm` mode, so the verbatim code path runs)."""
    LAST_PRIOR_REC = KPROD.LAST_PRIOR_REC

    @staticmethod
    def rc_hybrid_block_k(*a, **kw):
        if BLOCK_KEY is not None:
            kw["block_key"] = BLOCK_KEY[CUR[0]]
            kw["block_keep"] = KEEP
            if STEER != "key":
                kw["block_steer"] = STEER
        return KPROD.rc_hybrid_block_k(*a, **kw)


M.K = KWrap
_carrier = M.carrier_block
if N_CPM != 2200:
    M.H.CALLSET = (400, N_CPM, 1.9)

HK: dict[str, np.ndarray] = {}
STRONG: set[int] = set()


def _clean4_P_full():
    """517's / 478's CLEAN4 P (mean of the four banks' P_down), NaN where uncovered (run517.clean4_P_full, copied)."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in GF])
    out = np.full(len(GF), np.nan)
    out[col >= 0] = P[col[col >= 0]]
    return out


def _hs(P, tau):
    sg = np.where(np.isfinite(P) & (np.abs(P - 0.5) >= tau), np.where(P > 0.5, -1.0, 1.0), 0.0)
    return np.ascontiguousarray(sg[KEEPH])


def _shuffle_cov(P, tau, rng, tag):
    """-s: permute P among the covered genes of the gated axis (signs move, coverage and balance stay); records to REC."""
    Pk = P[KEEPH].copy()
    cov = np.flatnonzero(np.isfinite(Pk) & (np.abs(Pk - 0.5) >= tau))
    before = _hs(P, tau)
    Pk[cov] = Pk[rng.permutation(cov)]
    out = P.copy()
    out[KEEPH] = Pk
    after = _hs(out, tau)
    assert np.array_equal(before != 0, after != 0) and before.sum() == after.sum(), "shuffle moved coverage or balance"
    REC[f"shuffle_{tag}_n_covered"] = int(cov.size)
    REC[f"shuffle_{tag}_sign_agreement"] = float((before[cov] == after[cov]).mean())
    return out


if HOOK == "c4h":
    _Pc4 = _clean4_P_full()
    if R_FLIP > 0:
        # the same emulation as mzgrid517's: flip the source (P -> 1 - P) on a fixed random fraction r of the gated axis
        _u = np.random.default_rng(51717).random(int(KEEPH.sum()))
        _Pk = _Pc4[KEEPH].copy()
        _Pk = np.where(_u < R_FLIP, 1.0 - _Pk, _Pk)
        _Pc4 = _Pc4.copy()
        _Pc4[KEEPH] = _Pk
        REC["r_flip_n_flipped_gated"] = int((_u < R_FLIP).sum())
    if SHUF:
        # the zero-information floor: permute P among the covered genes of the gated axis (signs move, coverage stays)
        _seed = int(sys.argv[2])
        _Pk = _Pc4[KEEPH].copy()
        _cov = np.flatnonzero(np.isfinite(_Pk) & (np.abs(_Pk - 0.5) >= HOOK_TAU))
        _before = _hs(_Pc4, HOOK_TAU)
        _Pk[_cov] = _Pk[np.random.default_rng(51718 + _seed).permutation(_cov)]
        _Pc4 = _Pc4.copy()
        _Pc4[KEEPH] = _Pk
        _after = _hs(_Pc4, HOOK_TAU)
        assert np.array_equal(_before != 0, _after != 0) and _before.sum() == _after.sum(), "shuffle moved coverage or balance"
        REC["shuffle_seed"] = 51718 + _seed
        REC["shuffle_n_covered"] = int(_cov.size)
        REC["shuffle_sign_agreement"] = float((_before[_cov] == _after[_cov]).mean())
    HK["strong"] = HK["weak"] = _hs(_Pc4, HOOK_TAU)
elif HOOK == "mzg":
    # run528's +mzg block (= run517's `c4-mzg0.05`), copied; the optional -r flip is mzgrid517's emulation
    _H1T = set(REF_TARGETS)
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
    _col = np.array([_pos.get(str(g), -1) for g in GF])
    _Pa = np.where(_col >= 0, _pd[np.maximum(_col, 0)], np.nan)
    _ess = {ln.strip() for ln in open(ROOT / "data/donors/k562_essential/perts.txt") if ln.strip()}
    _ne = np.array([_perts[i] not in _ess for i in _rows])
    _sd = _X.std(0)
    _sd[_sd == 0] = 1.0
    _zb = (_X[_ne] / _sd).mean(0)
    del _X
    _z = np.where(_col >= 0, _zb[np.maximum(_col, 0)], np.nan)
    _ncov = int((np.isfinite(_Pa[KEEPH]) & (np.abs(_Pa[KEEPH] - 0.5) >= 0.05)).sum())
    _zq = float(np.sort(np.abs(np.nan_to_num(_z[KEEPH])))[::-1][_ncov - 1])
    if R_FLIP > 0:
        _u = np.random.default_rng(51717).random(int(KEEPH.sum()))
        _zk = _z[KEEPH].copy()
        _zk = np.where(_u < R_FLIP, -_zk, _zk)
        _z = _z.copy()
        _z[KEEPH] = _zk
        REC["r_flip_n_flipped_gated"] = int((_u < R_FLIP).sum())
    _Pz = np.where(np.isfinite(_z), np.clip(0.5 - 0.05 * _z / _zq, 0.001, 0.999), np.nan)
    _Pw = _clean4_P_full()
    if SHUF:
        _rng = np.random.default_rng(51718 + int(sys.argv[2]))
        REC["shuffle_seed"] = 51718 + int(sys.argv[2])
        _Pz = _shuffle_cov(_Pz, 0.05, _rng, "strong")
        _Pw = _shuffle_cov(_Pw, 0.05, _rng, "weak")
    HK["strong"], HK["weak"] = _hs(_Pz, 0.05), _hs(_Pw, 0.05)
if HK:
    REC["hook_signs"] = {k: int((v != 0).sum()) for k, v in HK.items()}
    REC["hook_signs_sha256"] = {k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in HK.items()}
    REC["n_strong"] = len(STRONG)


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    CUR[0] = i
    if BLOCK_KEY is not None:
        assert BLOCK_KEY.shape[1] == int(np.asarray(emit_c.keep).sum()), "key is not on the gated axis"
    if HOOK:
        M.MDE_ADD_SIGN = HK["strong"] if i in STRONG else HK["weak"]
    return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)


M.carrier_block = carrier_block

if BULK_G is not None:
    assert "--bulk" in sys.argv and sys.argv[sys.argv.index("--bulk") + 1] == "ps", "-bs steps the ps bulk"
    _ps_bulk = M.ps_bulk

    def ps_bulk_step(eff, genes, targets):
        """410b's key-scaled bulk with its per-target scale raised to BULK_G, same non-panel energy (a bulk-only step)."""
        out = _ps_bulk(eff, genes, targets)
        panel = set(map(str, targets))
        nonp = np.array([str(g) not in panel for g in genes])
        n0 = np.sqrt(np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp]))
        n1 = np.sqrt(np.einsum("ij,ij->i", out[:, nonp], out[:, nonp]))
        s = np.where(n0 > 0, n1 / np.maximum(n0, 1e-300), 1.0) ** BULK_G
        c = float(np.sqrt((n0 ** 2).sum() / ((s * n0) ** 2).sum()))
        REC["bulk_step_scale"] = dict(g=BULK_G, s_ps_min=float((n1 / np.maximum(n0, 1e-300)).min()),
                                      s_ps_max=float((n1 / np.maximum(n0, 1e-300)).max()),
                                      s_new_min=float((s * c).min()), s_new_max=float((s * c).max()), c=c)
        (WORK / "run537.json").write_text(json.dumps(REC, indent=1))
        return eff * (s * c)[:, None]

    M.ps_bulk = ps_bulk_step

if __name__ == "__main__":
    (WORK / "run537.json").write_text(json.dumps(REC, indent=1))
    print(f"[537] {json.dumps(REC)}", flush=True)
    M.main()
