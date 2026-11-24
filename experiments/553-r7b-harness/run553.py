#!/usr/bin/env python3
"""553 — the surface-parametric real-emitter DE harness: 537's run537.py with its H1 inputs replaced by a SURFACE's.

What stays identical (imported, never copied): run496.main — the per-target loop (profile / spread / rank axis / ballast K /
carrier block / zs pass / block_stats), 410b's ballast machinery, 426's MDE add, the production emitter
src/vcc2026/rc_emit_k.py (sha asserted 0f56209d...), the CallSet semantics (run400's CALLSET), evaluate323.block_stats (the
prediction's own Wilcoxon against the controls) and evaluate323.score (H1 only).
What the surface supplies (lib553.surface_spec):
  * the control reference (r, c, emit_c): 323's reference.npz on h1; controlref553.npz (lib553.controlref on
    context_R.h5ad, the same functions) on r7b. Any number of controls (the picks draw 400 of them, as production).
  * the emitter's library classes: the context's control-library q25/q50/q75 (production's `de_sign_classes`); on h1 this
    is run400's CLASSES3 constant exactly (asserted).
  * the carriers: the surface's cap15 effects (323's construction; effects553.py on r7b), sha-checked like native_cap15.
  * the bulk: 410b's key-scaled bulk (run496.ps_bulk's code) on the surface's own whitening operator and columns.
  * the call-propensity key (-kept / lcap): 528's keys528.npz `lcap` on h1 (as 537); on r7b the October form,
    528's callblock528.py on the surface's controls (log10 control CPM + 22 x excess dropout, median-rescaled to H1).
  * the CLEAN4 hook's P: 371's consensus prior on the surface's gene axis (537's _clean4_P_full with GF = the axis).
Scoring: on h1 the original evaluate323.score runs inline (metrics json comparable with 528 / 537); on a SEALED surface
the score is a stub (raw = null) — nothing reads the truth. The access (access553.py) scores the stats files later.

MODE grammar = run537's: <key>[-k<keep>][-kept][-n<cpm>][+c4h[tau]|+mzg][-r<r>][-s][-bs<g>]  (+mzg: h1 only)
usage: uv run python experiments/553-r7b-harness/run553.py <surface> <MODE> base <seed> --bulk ps --mde-c 0.1 \
         --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs [--bulk-cache live|use] [--class-rule production|flex]
         [--exclude-targets A,B] [--tag <suffix>]
  surfaces: h1 (the KA surface), h1cb (h1 through R7b's code path), h1r7 (h1 cut to R7b's shape), r7b (as built),
            r7bd<floor> (R7b's controls above a human-library floor: the candidate amendment, prediction side)
  --bulk-cache use : 410b's key from work/<surface>/bulk553.npz (written and verified by a live run; KA-B: bit-equal)
  --class-rule flex: see below (a deviation; R7b as built only)      --exclude-targets: not emitted (panel kept whole)
Outputs: experiments/553-r7b-harness/work/<surface>_<MODE><tag>/ only (tag: -flex, -x<targets> or --tag).
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
SURFACE = sys.argv.pop(1)
MODE = sys.argv.pop(1)
# 553 options (removed before run496's argparse): --bulk-cache live|use (default live: 410b's key computed from the operator;
# live also writes / verifies work/<surface>/bulk553.npz; `use` reads it and never loads the operator), --tag <suffix>
# --class-rule production|flex (default production): the emitter's library classes. production = the context's control-
# library q25 / q50 / q75 (make_submission_pm_cd4.py's `de_sign_classes`). flex = the context's median x H1's relative
# levels (16078.75 / 20190.5, 1, 25330.5 / 20190.5): a DEVIATION for surfaces whose quartiles are far apart (R7b: q75/q25
# 7.24 vs <= 1.95 on every Flex context), where the production placement fails (RESULTS.md); on h1 it IS production
# (KA-flex). flex outputs go to work/<surface>_<MODE>-flex unless --tag says otherwise.
BULK_CACHE = "live"
TAG_SUFFIX = None
CLASS_RULE = "production"
EXCLUDE = ()            # --exclude-targets A,B: not emitted (the panel -- operator, carriers, bulk key -- stays whole)
for _opt in ("--bulk-cache", "--tag", "--class-rule", "--exclude-targets"):
    _hit = [_j for _j, _x in enumerate(sys.argv) if _x == _opt or _x.startswith(_opt + "=")]
    if _hit:
        # `--opt value` or `--opt=value` (the form argparse needs for a value that begins with a dash, e.g. --tag=-flex)
        _i = _hit[0]
        if sys.argv[_i] == _opt:
            _v = sys.argv[_i + 1]
            del sys.argv[_i:_i + 2]
        else:
            _v = sys.argv[_i].split("=", 1)[1]
            del sys.argv[_i]
        if _opt == "--bulk-cache":
            assert _v in ("live", "use"), _v
            BULK_CACHE = _v
        elif _opt == "--class-rule":
            assert _v in ("production", "flex"), _v
            CLASS_RULE = _v
        elif _opt == "--exclude-targets":
            EXCLUDE = tuple(sorted(x for x in _v.split(",") if x))
        else:
            TAG_SUFFIX = _v
if TAG_SUFFIX is None:
    TAG_SUFFIX = ("-flex" if CLASS_RULE == "flex" else "") + (("-x" + "+".join(EXCLUDE)) if EXCLUDE else "")
m_ = re.fullmatch(r"(cpm|lcap)(?:-k([0-9]+))?(-kept)?(?:-n([0-9]+))?(\+c4h|\+mzg)?([0-9]\.[0-9]+)?(?:-r([0-9.]+))?(-s)?(?:-bs([0-9.]+))?", MODE)
assert m_, MODE
KEYNAME = m_.group(1)
KEEP = int(m_.group(2) or 0)
STEER = "kept" if m_.group(3) else "key"
N_CPM = int(m_.group(4) or 2200)
HOOK = (m_.group(5) or "").lstrip("+") or None
HOOK_TAU = float(m_.group(6) or 0.05)
R_FLIP = float(m_.group(7) or 0.0)
assert not (m_.group(6) and HOOK != "c4h"), "a tau suffix is defined for +c4h only"
assert not (STEER == "kept" and KEYNAME == "cpm"), "-kept needs a key"
assert not (R_FLIP > 0 and HOOK not in ("mzg", "c4h")), "-r attenuates a hook source only"
SHUF = bool(m_.group(8))
assert not (SHUF and HOOK not in ("c4h", "mzg")), "-s shuffles a hook source only"
BULK_G = float(m_.group(9)) if m_.group(9) else None
SEED = int(sys.argv[2])
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; its WORK and inputs are redirected below)
from vcc2026 import rc_emit  # noqa: E402
from vcc2026 import rc_emit_k as KPROD  # noqa: E402
import lib553 as S  # noqa: E402

EMITTER_SHA = hashlib.sha256(Path(KPROD.__file__).read_bytes()).hexdigest()
assert EMITTER_SHA == S.EMITTER_SHA, f"emitter is not the committed production 0f56209d: {EMITTER_SHA}"
assert M.K is KPROD, "run496 did not import the production emitter"
SPEC = S.surface_spec(SURFACE)
if HOOK == "mzg" and SURFACE != "h1":
    raise SystemExit("+mzg is generalised for h1 only: its strong/weak split is 467's H1 strength file and its source "
                     "excludes H1's targets; see RESULTS.md (gaps)")

WORK = HERE / "work" / f"{SURFACE}_{MODE}{TAG_SUFFIX}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
_write_json = M.R.write_json


def _finite(v):
    """non-finite floats -> None, recursively (JSON-only: run496 writes diagnostics with allow_nan=False, and a surface where
    some target carries no MDE record makes a diagnostic median NaN; never touches the stats arrays)"""
    if isinstance(v, float):
        return v if np.isfinite(v) else None
    if isinstance(v, dict):
        return {k: _finite(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_finite(x) for x in v]
    return v


def write_json_finite(path, value):
    return _write_json(path, _finite(value))


M.R.write_json = write_json_finite

# ---------------------------------------------------------------------------------------------------------------- inputs
# the small control arrays now; the sorted control blocks (cs / emit_cs: 3.2 GB on h1) only when run496.main asks for the
# reference, i.e. AFTER ps_bulk has released the whitening operator (run496's own order, so the peak RSS is run537's)
CR = S.load_controlref(SPEC, need_cs=False)
GF = CR["genes"].astype(str)
KEEPH = np.asarray(CR["keep"], dtype=bool)
TARGETS_ALL = CR["targets"].astype(str)
assert set(EXCLUDE) <= set(TARGETS_ALL), f"excluded targets not on the surface: {set(EXCLUDE) - set(TARGETS_ALL)}"
KEEP_T = ~np.isin(TARGETS_ALL, list(EXCLUDE))
TARGETS = TARGETS_ALL[KEEP_T]
N_CTRL = int(CR["n_controls"])
R_DICT = {"genes": CR["genes"], "targets": TARGETS, "libraries": CR["libraries"], "b_ctrl": CR["b_ctrl"]}
assert SURFACE != "h1" or not EXCLUDE, "h1 is never subset (its inline score reads the whole reference)"
if SURFACE == "h1":
    # h1 only: the truth half of 323's reference for evaluate323.score's inline known-answer score
    with np.load(SPEC["controlref"], allow_pickle=False) as z:
        R_DICT.update(perts=z["perts"], sums=z["sums"], counts=z["counts"], jk=z["jk"])
CLASSES = S.library_classes(R_DICT["libraries"])
if SURFACE == "h1":
    assert CLASSES == tuple(M.H.CLASSES3), (CLASSES, M.H.CLASSES3)
if CLASS_RULE == "flex":
    CLASSES = S.flex_classes(R_DICT["libraries"])
M.H.CLASSES3 = CLASSES


def surface_reference():
    """evaluate323.load_reference's return (r, scorer control block, emitter control block) for the surface."""
    with np.load(SPEC["controlref"], allow_pickle=False) as z:
        cs, ecs = z["cs"], z["emit_cs"]
    c = rc_emit.RcControlBlock(keep=CR["keep"], cs=cs, tie=CR["tie"], ref_lin=CR["ref_lin"], n_controls=N_CTRL, min_cpm=5)
    ec = rc_emit.RcControlBlock(keep=c.keep, cs=ecs, tie=CR["emit_tie"], ref_lin=CR["emit_ref_lin"], n_controls=N_CTRL,
                                min_cpm=5)
    return R_DICT, c, ec


M.E.load_reference = surface_reference


def surface_cap15(b):
    """run400.native_cap15's check on the surface's effects: the cap15 matrix must hash to its compose record."""
    comp = json.loads(Path(SPEC["compose"]).read_text())
    with np.load(SPEC["effects"], allow_pickle=False) as z:
        eff = z[SPEC["effects_key"]]
        if b.array_sha256(eff) != comp["arms"][SPEC["effects_key"]]["final_sha256"]:
            raise ValueError("cap15 final effect hash drift vs its compose record")
        if list(z["targets"].astype(str)) != list(TARGETS_ALL):
            raise ValueError("the effects' targets are not the surface's target list")
        return eff[KEEP_T].copy(), z["genes"].copy(), z["targets"][KEEP_T].copy()


M.H.native_cap15 = surface_cap15


class SurfaceOp:
    """run496.H1Op on the surface's operator: evals / evecs of the whitener, the operator columns on the official axis."""
    def __init__(self):
        with np.load(SPEC["operator"]) as z:
            self.evals = np.asarray(z["evals"], float)
            self.evecs = np.asarray(z["evecs"], float)
        with np.load(SPEC["op_columns"], allow_pickle=False) as z:
            self.cols = np.asarray(z["columns"], dtype=np.int64)
        assert self.evecs.shape[0] == len(self.cols), (self.evecs.shape, len(self.cols))

    def _check_axis(self, wcols):
        assert np.array_equal(np.asarray(wcols), self.cols)


BULK_FILE = HERE / "work" / ("h1" if SURFACE == "h1cb" or SURFACE.startswith("h1r7") else "r7b" if SURFACE.startswith("r7bd")
                              else SURFACE) / "bulk553.npz"


def ps_bulk(eff, genes, targets):
    """run496.ps_bulk verbatim with the surface's operator: 410's key weights, clip(0.3, 3), energy-preserving over the
    non-panel genes. `live` computes the key from the operator (and writes / verifies BULK_FILE); `use` applies BULK_FILE's
    (w, c) after checking its provenance (targets, the cap15 matrix's sha, the operator's sha) -- the same arithmetic."""
    import key410  # noqa: PLC0415  (410's file, read-only)
    eff_kept = eff
    if EXCLUDE:
        # the key and the energy constant are panel constructs: computed on the whole panel, the kept rows returned
        with np.load(SPEC["effects"], allow_pickle=False) as z:
            eff = np.asarray(z[SPEC["effects_key"]]).copy()
        targets = TARGETS_ALL
        assert np.array_equal(eff[KEEP_T], eff_kept)
    eff_sha = S.arr_sha(eff)
    op_sha = REC["inputs_sha256"]["operator"]
    if BULK_CACHE == "use":
        with np.load(BULK_FILE, allow_pickle=False) as z:
            assert list(z["targets"].astype(str)) == list(map(str, targets)), "bulk cache: targets differ"
            assert str(z["eff_sha256"]) == eff_sha, "bulk cache: built for another carrier matrix"
            assert str(z["operator_sha256"]) == op_sha, "bulk cache: built with another operator"
            w, c = np.asarray(z["w"], dtype=np.float64), float(z["c"])
            key, carried, nd, med = z["key"], z["carried"], z["n_donors"], float(z["median"])
        REC["bulk_cache"] = dict(mode="use", file=str(BULK_FILE.relative_to(ROOT)), sha256=S.sha256(BULK_FILE))
    else:
        op = SurfaceOp()
        key, carried, nd = key410.shared_energy_key(list(map(str, targets)), M.DONORS, op, op.cols, cols=None)
        w, med = key410.weights(key, carried, 0.3, 3.0)
        del op
        panel = set(map(str, targets))
        nonp = np.array([str(g) not in panel for g in genes])
        e2 = np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp])
        c = float(np.sqrt(e2.sum() / (w * w * e2).sum()))
        if BULK_FILE.exists():
            with np.load(BULK_FILE, allow_pickle=False) as z:
                same = (np.array_equal(z["w"], w) and float(z["c"]) == c and str(z["eff_sha256"]) == eff_sha
                        and str(z["operator_sha256"]) == op_sha)
            assert same, "the live bulk key differs from the cached one"
            REC["bulk_cache"] = dict(mode="live, cache verified equal", file=str(BULK_FILE.relative_to(ROOT)))
        else:
            BULK_FILE.parent.mkdir(parents=True, exist_ok=True)
            tmp = BULK_FILE.with_suffix(".tmp.npz")
            np.savez(tmp, targets=np.asarray(list(map(str, targets))), w=w, c=np.float64(c), key=key, carried=carried,
                     n_donors=nd, median=np.float64(med), eff_sha256=np.array(eff_sha), operator_sha256=np.array(op_sha))
            tmp.replace(BULK_FILE)
            REC["bulk_cache"] = dict(mode="live, cache written", file=str(BULK_FILE.relative_to(ROOT)))
    REC["bulk_key"] = dict(n_carried=int(np.sum(carried)), n_donors_hist={int(k): int(v) for k, v in zip(*np.unique(nd, return_counts=True))},
                           median=float(med), w_min=float(w.min()), w_max=float(w.max()), energy_c=c,
                           w_sha256=S.arr_sha(w))
    if EXCLUDE:
        return eff_kept * (w[KEEP_T] * c)[:, None]
    return eff * (w * c)[:, None]


M.ps_bulk = ps_bulk

if SURFACE != "h1":
    # every surface but h1 is scored only by access553.py (a sealed one only after the owner's go)
    def sealed_score(r, c, sums, jks, ps, qs, lfcs):
        """A sealed surface is never scored here: the stats files go to access553.py after the owner's go."""
        return {m: None for m in S.MEMBERS}, [], None
    M.E.score = sealed_score

# ---------------------------------------------------------------------------------------------------------------- the key
BLOCK_KEY = None
KEY_SRC = None
if KEYNAME != "cpm":
    if SPEC["key"] == "keys528":
        KEY_SRC = ROOT / "experiments/528-callprop-block/work/keys528.npz"
        with np.load(KEY_SRC, allow_pickle=False) as z:
            assert list(z["targets"].astype(str)) == list(TARGETS), "key targets are not in the reference order"
            BLOCK_KEY = np.ascontiguousarray(z[KEYNAME], dtype=np.float64)
    else:
        # the October form: callblock528.py's artifact from the surface's own controls (key on the official axis)
        KEY_SRC = HERE / "work" / SURFACE / f"callblock_{SURFACE}.npz"
        rec = json.loads(Path(str(KEY_SRC) + ".json").read_text())
        assert S.sha256(KEY_SRC) == rec["sha256"], "callblock artifact sha drift"
        with np.load(KEY_SRC, allow_pickle=False) as z:
            assert list(z["genes"].astype(str)) == list(GF), "callblock artifact axis differs from the surface's"
            (kname,) = [k for k in z.files if k.startswith("key_")]
            kfull = np.asarray(z[kname], dtype=np.float64)
        assert KEYNAME == "lcap", "the October key is the lcap form"
        kg = np.ascontiguousarray(kfull[KEEPH])
        assert np.all(np.isfinite(kg)), "key not finite on the gated axis"
        BLOCK_KEY = np.ascontiguousarray(np.tile(kg, (len(TARGETS), 1)))
CUR = [None]
REC = {"surface": SURFACE, "mode": MODE, "seed": SEED, "key": KEYNAME, "block_keep": KEEP, "block_steer": STEER,
       "n_cpm": N_CPM, "hook": HOOK, "hook_tau": HOOK_TAU if HOOK == "c4h" else None, "r_flip": R_FLIP,
       "sign_shuffle": SHUF, "bulk_step_g": BULK_G, "classes": CLASSES, "n_controls": N_CTRL, "n_targets": len(TARGETS),
       "n_gated": int(KEEPH.sum()),
       "key_source": None if KEY_SRC is None else str(KEY_SRC.relative_to(ROOT)),
       "key_sha256": None if BLOCK_KEY is None else hashlib.sha256(BLOCK_KEY.tobytes()).hexdigest(),
       "emitter_sha256": EMITTER_SHA,
       "bulk_cache_mode": BULK_CACHE, "tag_suffix": TAG_SUFFIX, "class_rule": CLASS_RULE, "excluded_targets": list(EXCLUDE),
       "inputs_sha256": {k: S.sha256(Path(SPEC[k])) for k in ("controlref", "effects", "operator", "op_columns")}}


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

# ---------------------------------------------------------------------------------------------------------------- hooks
HK: dict[str, np.ndarray] = {}
STRONG: set[int] = set()


def _clean4_P_full():
    """517's / 478's CLEAN4 P (mean of the four banks' P_down) on the surface's axis, NaN where uncovered (run537's)."""
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
        _u = np.random.default_rng(51717).random(int(KEEPH.sum()))
        _Pk = _Pc4[KEEPH].copy()
        _Pk = np.where(_u < R_FLIP, 1.0 - _Pk, _Pk)
        _Pc4 = _Pc4.copy()
        _Pc4[KEEPH] = _Pk
        REC["r_flip_n_flipped_gated"] = int((_u < R_FLIP).sum())
    if SHUF:
        _Pk = _Pc4[KEEPH].copy()
        _cov = np.flatnonzero(np.isfinite(_Pk) & (np.abs(_Pk - 0.5) >= HOOK_TAU))
        _before = _hs(_Pc4, HOOK_TAU)
        _Pk[_cov] = _Pk[np.random.default_rng(51718 + SEED).permutation(_cov)]
        _Pc4 = _Pc4.copy()
        _Pc4[KEEPH] = _Pk
        _after = _hs(_Pc4, HOOK_TAU)
        assert np.array_equal(_before != 0, _after != 0) and _before.sum() == _after.sum(), "shuffle moved coverage or balance"
        REC["shuffle_seed"] = 51718 + SEED
        REC["shuffle_n_covered"] = int(_cov.size)
        REC["shuffle_sign_agreement"] = float((_before[_cov] == _after[_cov]).mean())
    HK["strong"] = HK["weak"] = _hs(_Pc4, HOOK_TAU)
elif HOOK == "mzg":
    # run537's +mzg block verbatim (h1 only; see the guard above)
    _H1T = set(TARGETS)
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
        _rng = np.random.default_rng(51718 + SEED)
        REC["shuffle_seed"] = 51718 + SEED
        _Pz = _shuffle_cov(_Pz, 0.05, _rng, "strong")
        _Pw = _shuffle_cov(_Pw, 0.05, _rng, "weak")
    HK["strong"], HK["weak"] = _hs(_Pz, 0.05), _hs(_Pw, 0.05)
if HK:
    REC["hook_signs"] = {k: int((v != 0).sum()) for k, v in HK.items()}
    REC["hook_signs_down"] = {k: int((v < 0).sum()) for k, v in HK.items()}
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

_ps_bulk = M.ps_bulk
if BULK_G is not None:
    assert "--bulk" in sys.argv and sys.argv[sys.argv.index("--bulk") + 1] == "ps", "-bs steps the ps bulk"

    def ps_bulk_step(eff, genes, targets):
        """run537's bulk-only step, on the surface's key-scaled bulk."""
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
        (WORK / "run553.json").write_text(json.dumps(REC, indent=1))
        return eff * (s * c)[:, None]

    M.ps_bulk = ps_bulk_step

if __name__ == "__main__":
    (WORK / "run553.json").write_text(json.dumps(REC, indent=1))
    print(f"[553] {json.dumps(REC)}", flush=True)
    M.main()
    (WORK / "run553.json").write_text(json.dumps(REC, indent=1))
