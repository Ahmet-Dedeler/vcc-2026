#!/usr/bin/env python3
"""612 — real production-emitter H1 runs of the KEY-GATED uncalled-gene sign hook (council 600, [21dd4c]).

537's harness pattern (run537.py, read-only): 435's construction (production rc_emit_k 0f56209d, CallSet(400, 2200, 1.9),
CLEAN4 q80 block prior, coarse MDE; `--zs dn_k250_zs` adds the incumbent's zs pass) with the emitter's `mde_add_sign`
injected per target. 496's plain CLEAN4 hook passes the CLEAN4 sign (tau 0.05) on every complement gene; the KEY-GATED hook
passes it only where the target's design key is weak, and 0 (= keep the delivered sign, i.e. the key's) elsewhere:

MODE = c4kg[-q<q>] | c4kgall | idkg[-q<q>] | idkgd (= idkg, own dir; for --dump)
  c4kg    : weak half = |key| below the target's q-quantile (default 0.5) of its NONZERO |key| on the gated axis; the weak
            genes take CLEAN4's sign where |P_down - 0.5| >= 0.05, others 0
  c4kgall : the gate forced open (every gene "weak") = 537's `cpm+c4h` exactly -> KNOWN ANSWER (stats bit-equal)
  idkg    : 622's recommended form ([7772dc]): the weak half (and key-less genes) take sign(S), S = z(P_C4 - 0.5) +
            0.5 z(I) with Flex-only identity references (A/B/C), identity alone where CLEAN4 is undefined (no tau: S's own
            coverage, ~10.8k gated genes vs the hook's 2,550); work/idsign612.npz, built by idsign612.py
Outputs: experiments/612-key-gated-hook/work/h1_<MODE>/ only.
usage: uv run python experiments/612-key-gated-hook/run612.py <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 \
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
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
os.environ["RUN426_EMITTER"] = "prod"
os.environ.pop("RUN496_EMITTER", None)
sys.dont_write_bytecode = True
MODE = sys.argv.pop(1)
m_ = re.fullmatch(r"(c4kg|c4kgall|idkg|idkgd|idkg3|idall3|c4kgd|idkg3d|idall3d)(?:-q([0-9.]+))?", MODE)
assert m_, MODE
# *d modes = the same construction into its own dir, for --dump 85 ([d4cd5f]'s 630 stacks); KA: stats == the scored mode's
KIND = {"idkgd": "idkg", "idkg3": "idkg", "idall3": "idall", "idkg3d": "idkg", "c4kgd": "c4kg",
        "idall3d": "idall"}.get(m_.group(1), m_.group(1))
# idkg3 = the identity source in 604 v3's rule (no identity-only extension; idsign612v3.py), [d4cd5f] 10-02 ~13:55Z
# idall3 = the same v3 source on EVERY uncalled gene, no key gate (612/PREDICTION-idall3.md, registered before this mode)
IDFILE = "idsign612v3.npz" if m_.group(1) in ("idkg3", "idall3", "idkg3d", "idall3d") else "idsign612.npz"
Q = float(m_.group(2) or 0.5)
assert 0.0 < Q <= 1.0
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; WORK redirected below)
from vcc2026 import rc_emit_k as KPROD  # noqa: E402

EMITTER_SHA = hashlib.sha256(Path(KPROD.__file__).read_bytes()).hexdigest()
assert EMITTER_SHA.startswith("0f56209d57e8385a"), f"emitter is not the committed production 0f56209d: {EMITTER_SHA}"
assert M.K is KPROD, "run496 did not import the production emitter"
WORK = HERE / "work" / f"h1_{MODE}"
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.HERE = WORK   # run496's --dump writes M.HERE / "work/blocks" (its only post-import use of HERE): keep dumps in this run's dir
REF = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
KEEPH = np.asarray(REF["keep"], dtype=bool)
GF = REF["genes"].astype(str)


def _clean4_P_full():
    """run537._clean4_P_full, copied: the mean of the four CLEAN4 banks' P_down, NaN where uncovered."""
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


HOOK_C4 = _hs(_clean4_P_full(), 0.05)                      # == run537's HK["weak"] for +c4h at tau 0.05
if KIND in ("idkg", "idall"):
    with np.load(HERE / "work" / IDFILE, allow_pickle=False) as z:
        assert list(z["genes"].astype(str)) == list(GF), f"{IDFILE} is not on the reference gene axis"
        WEAK_SRC = np.ascontiguousarray(np.asarray(z["sign"], dtype=np.float64)[KEEPH])
        ID_SHA = hashlib.sha256((HERE / "work" / IDFILE).read_bytes()).hexdigest()
    # 622's delivery: S's sign on its own coverage (no tau restriction), key-gated like c4kg
else:
    WEAK_SRC, ID_SHA = HOOK_C4, None
REC = {"mode": MODE, "kind": KIND, "q": Q, "emitter_sha256": EMITTER_SHA, "idsign_sha256": ID_SHA,
       "hook_c4_n": int((HOOK_C4 != 0).sum()), "hook_c4_sha256": hashlib.sha256(HOOK_C4.tobytes()).hexdigest(),
       "weak_src_n": int((WEAK_SRC != 0).sum()), "weak_src_agree_c4": float((WEAK_SRC[HOOK_C4 != 0] == HOOK_C4[HOOK_C4 != 0]).mean()),
       "per_target": {}}
_carrier = M.carrier_block


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    # run496 passes the design key on the FULL axis (np.where(gate, eff_de[i], 0)); the emitter reads key[control.keep]
    # (rc_emit_k: key_gated = key[keep]), and mde_add_sign lives on that gated axis -> gate on the same slice.
    k = np.asarray(key, dtype=np.float64)[np.asarray(emit_c.keep, dtype=bool)]
    assert k.shape == HOOK_C4.shape, f"gated key {k.shape} vs hook {HOOK_C4.shape}"
    if KIND == "c4kgall":
        sign = HOOK_C4
    elif KIND == "idall":
        sign = WEAK_SRC                                   # ungated: the v3 source wherever it signs
    else:
        a = np.abs(k)
        nz = a > 0
        cut = float(np.quantile(a[nz], Q)) if nz.any() else 0.0
        weak = (~nz) | (a < cut) if Q < 1.0 else np.ones_like(nz)
        sign = np.ascontiguousarray(np.where(weak, WEAK_SRC, 0.0))
        REC["per_target"][str(i)] = {"n_weak": int(weak.sum()), "n_sign": int((sign != 0).sum()), "cut": cut}
    M.MDE_ADD_SIGN = np.ascontiguousarray(sign)
    return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)


M.carrier_block = carrier_block

if __name__ == "__main__":
    (WORK / "run612.json").write_text(json.dumps({k: v for k, v in REC.items() if k != "per_target"}, indent=1))
    print(f"[612] {json.dumps({k: v for k, v in REC.items() if k != 'per_target'})}", flush=True)
    M.main()
    (WORK / "run612.json").write_text(json.dumps(REC, indent=1))
