#!/usr/bin/env python3
"""625f — ORACLE-direction calibration of 625e's exchange rate: the same real-emitter H1 runs as flip625.py, but the flip
set is the genes of ONE CLEAN4 octile where H1's own in-sample programme DISAGREES with CLEAN4. That flip IMPROVES gate
signs, whereas 625d's whole-octile flips degrade them, so the pair asks whether the exchange rate (overall per unit of
source gate-sign change) is symmetric. Leaky by design: a calibration of delivery, never a lever.
Programme: P_down_H1(g) = fraction of the targets with g in their f 1.0 reference DE pool (521's pools, own gene excluded)
where H1's lfc < 0, over genes in >= 3 pools. Flip set: octile k (modules625_v8 `octile`) AND CLEAN4 covered AND
(P_down_H1 > 1/2) != (P_C4 > 1/2).
(flip625.py's original docstring follows.)
625d — real production-emitter H1 runs that FLIP one CLEAN4 octile in the BLOCK prior and the hook together, the
block-prior share of the module readout's value. Spec agreed with [d4cd5f] (2026-10-02): a few octiles x 2-3 seeds,
`--dump 85`, with the dumps in this directory so 630's UP pass can be stacked on the same blocks (priced with and without UP).

Base = 612's `c4kg` (435's construction: production rc_emit_k 0f56209d, CallSet(400, 2200, 1.9), CLEAN4 q80 block prior,
coarse MDE, `--zs dn_k250_zs`; the hook passes CLEAN4's sign (tau 0.05) on each target's WEAK-key half, 0 on the strong half).
OCT = none | 0..7:
  none : exactly 612's c4kg (KNOWN ANSWER: stats must equal experiments/612-key-gated-hook/work/h1_c4kg bit for bit)
  k    : genes in CLEAN4 octile k (625's modules625_v8.npz, by name) get the opposite sign in BOTH the block prior
         (run400.prior_signs_full) and the hook source.
Outputs: experiments/625-module-readout/oct<OCT>/work/h1/ (stats, metrics) and .../oct<OCT>/work/blocks/ (dumps).
usage: uv run python experiments/625-module-readout/flip625.py <OCT> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 \
         --mde-min-count 0.3 --mde-coarse --zs dn_k250_zs [--dump 85]
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODE = sys.argv.pop(1)                            # o<k>: oracle-disagreeing genes of octile k
assert MODE[:1] == "o" and MODE[1:] in {str(k) for k in range(8)}, MODE
OCT = MODE[1:]
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
ODIR = HERE / f"oracle{OCT}"                    # per-mode tree: stats in oracle<k>/work/h1/, dumps in oracle<k>/work/blocks/
WORK = ODIR / "work" / "h1"                       # (630's stack driver reads --src oct<k>/work and looks in --src/h1)
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = ODIR                                     # run496 dumps to M.HERE / "work/blocks" (21dd4c's trap: NOT WORK)
REF = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
KEEPH = np.asarray(REF["keep"], dtype=bool)
GF = REF["genes"].astype(str)
V8 = np.load(HERE / "work/modules625_v8.npz", allow_pickle=True)
OCTMAP = dict(zip(V8["genes"].astype(str), V8["octile"].astype(int)))
OCT_FULL = np.array([OCTMAP.get(str(g), -1) for g in GF])
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as _L  # noqa: E402
import table521 as _TB  # noqa: E402
_S = _L.load_state(0)
_pool = _TB.regime_pool(_S, "1.0")[0] & ~_S["own"]
_rl = np.asarray(_S["rl"], float)
_npool = _pool.sum(0)
_pdn = np.where(_npool >= 3, ((_rl < 0) & _pool).sum(0) / np.maximum(_npool, 1), np.nan)
PDN_H1 = dict(zip(np.asarray(_S["genes"]).astype(str), _pdn))
del _S, _pool, _rl


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
_ph1 = np.array([PDN_H1.get(str(g), np.nan) for g in GF])
FLIP_FULL = (OCT_FULL == int(OCT)) & np.isfinite(P_FULL) & np.isfinite(_ph1) & ((_ph1 > 0.5) != (P_FULL > 0.5))
FLIP_NAMES = set(GF[FLIP_FULL].tolist())
HOOK_FULL = np.where(np.isfinite(P_FULL) & (np.abs(P_FULL - 0.5) >= 0.05), np.where(P_FULL > 0.5, -1.0, 1.0), 0.0)
HOOK_FULL = np.where(FLIP_FULL, -HOOK_FULL, HOOK_FULL)                 # octile k: the hook's sign flipped too
HOOK = np.ascontiguousarray(HOOK_FULL[KEEPH])

_prior = M.H.prior_signs_full


def prior_signs_full(genes_full):
    out, cov = _prior(genes_full)
    g = np.asarray([str(x) for x in genes_full])
    flip = np.array([x in FLIP_NAMES for x in g])
    out = np.where(flip, -out, out)
    return out, cov


M.H.prior_signs_full = prior_signs_full
REC = {"mode": MODE, "oct": OCT, "n_oracle_flip": int(FLIP_FULL.sum()), "n_octile": int((OCT_FULL == int(OCT)).sum()), "emitter_sha256": EMITTER_SHA, "n_flip_full": int(FLIP_FULL.sum()),
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
    (WORK / "oracle625.json").write_text(json.dumps(REC, indent=1))
    print(f"[625f] {json.dumps(REC)}", flush=True)
    M.main()
