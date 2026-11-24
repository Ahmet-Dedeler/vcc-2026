#!/usr/bin/env python3
"""599 — H1 real-emitter runs of the per-target bulk sharpening tier rule (tercile vs tercilec) on a coverage-thinned key.

run612's pattern: import run496 (435's construction = 446bj2z's DE design: production emitter 0f56209d, CallSet(400, 2200,
1.9), CLEAN4 q80 block prior, coarse MDE; `--zs dn_k250_zs` adds the zs pass), redirect its WORK, and swap three things:
  M.ps_bulk          -> the mode's bulk (common599.mode_bulk: the cached 410 key, the mode's carried set and clip),
                        energy-preserving; it checks that run496 passes the cap15 matrix it was built from;
  M.E.load_reference -> the same loader, keeping a handle on the reference it returns (no second 3 GB load);
  M.ballast_block    -> the same block, with profile_bulk replaced by the builder's own
                        profile_from_effect(eff_bk[i], b_ctrl, min_cpm=10, sharpen_q=q_i, protect=panel): production's
                        ballast render (make_submission_pm_cd4.py, the profile_bulk call). Before replacing it, it asserts
                        that the profile run496 passed equals the unsharpened profile recomputed here, bit for bit.
MODE in {ka, s175, t, tc} (common599.MODES). ka swaps nothing that matters (its bulk == ps_bulk, q None): KA vs 435.
The q / carried arrays are checked against prep599.npz (the registered design).
usage: uv run python experiments/599-tercilec-h1-de/run599.py <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 \
         --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
os.environ["RUN426_EMITTER"] = "prod"
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
MODE = sys.argv.pop(1)
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import common599 as C  # noqa: E402
import run496 as M  # noqa: E402  (read-only; WORK redirected below)
from vcc2026 import rc_emit_k as KPROD  # noqa: E402

assert MODE in C.MODES, MODE
EMITTER_SHA = hashlib.sha256(Path(KPROD.__file__).read_bytes()).hexdigest()
assert EMITTER_SHA.startswith("0f56209d57e8385a"), f"emitter is not the committed production 0f56209d: {EMITTER_SHA}"
assert M.K is KPROD, "run496 did not import the production emitter"
WORK = HERE / "work" / f"h1_{MODE}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = WORK

B = M.H.builder389()
EFF, GENES, TARGETS = M.H.native_cap15(B)
KZ = C.fit_key(TARGETS)
EB, W_, CARRIED, CC = C.mode_bulk(MODE, EFF, GENES, TARGETS, KZ)
QT, AUDIT = C.mode_q(MODE, EB, EFF, GENES, TARGETS, CARRIED)
with np.load(C.WORK / "prep599.npz", allow_pickle=False) as P:
    assert np.array_equal(P[f"carried_{MODE}"], CARRIED), "carried set differs from the registered prep"
    assert np.array_equal(P[f"w_{MODE}"], W_), "weights differ from the registered prep"
    if QT is not None:
        assert np.array_equal(P[f"q_{MODE}"], QT), "per-target q differs from the registered prep"
PANEL = set(map(str, TARGETS))
STATE: dict = {}
REC = {"mode": MODE, "spec": C.MODES[MODE], "emitter_sha256": EMITTER_SHA, "key_sha256": C.sha256(C.KEYFILE),
       "builder598_sha256": C.BUILDER598_SHA, "builder_live_sha256": C.sha256(ROOT / "scripts/make_submission_pm_cd4.py"),
       "n_carried": int(CARRIED.sum()), "c": CC, "audit": {k: v for k, v in AUDIT.items() if k != "moved_vs_tercile"},
       "q": None if QT is None else QT.tolist(), "per_target": {}}

_ps_bulk = M.ps_bulk


def ps_bulk(eff, genes, targets):
    assert np.array_equal(eff, EFF) and list(map(str, genes)) == list(map(str, GENES)) \
        and list(map(str, targets)) == list(map(str, TARGETS)), "run496 passed another cap15 matrix"
    return EB.copy()


_load_reference = M.E.load_reference


def load_reference():
    out = _load_reference()
    STATE["r"] = out[0]
    return out


_ballast_block = M.ballast_block


def ballast_block(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, k_force=None):
    r = STATE["r"]
    if "eff_bk" not in STATE:
        STATE["eff_bk"] = M.R.map_axis(EB, GENES, r["genes"])
        STATE["panel"] = np.array([str(g) in PANEL for g in r["genes"]])
    pb0 = B.profile_from_effect(STATE["eff_bk"][i], r["b_ctrl"], min_cpm=10.)
    assert np.array_equal(pb0, profile_bulk), f"target {i}: run496's bulk profile != the recomputed unsharpened one"
    if QT is not None:
        profile_bulk = B.profile_from_effect(STATE["eff_bk"][i], r["b_ctrl"], min_cpm=10., sharpen_q=float(QT[i]),
                                             protect=STATE["panel"])
    block, dg = _ballast_block(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, k_force=k_force)
    REC["per_target"][str(i)] = {"q": None if QT is None else float(QT[i]), "K": int(dg["ballast_K"]),
                                 "carried": bool(CARRIED[i])}
    return block, dg


M.ps_bulk = ps_bulk
M.E.load_reference = load_reference
M.ballast_block = ballast_block

if __name__ == "__main__":
    print(f"[599] {json.dumps({k: v for k, v in REC.items() if k not in ('per_target', 'q')})}", flush=True)
    M.main()
    seed = [a for a in sys.argv[1:] if a.isdigit()][0]
    (WORK / f"run599_s{seed}.json").write_text(json.dumps(REC, indent=1))
