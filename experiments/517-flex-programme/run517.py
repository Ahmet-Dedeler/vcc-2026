#!/usr/bin/env python3
"""517 — REAL emitter runs on H1 for the programme-sourced hook (confirms price517's analytic numbers, and measures what the
analytic path cannot: leaked calls / jac and the |z| realisation's reach). 496's harness (`run496.py`, RUN496_EMITTER=hook:
the production emitter edb5b8dc + patch_emitter496.py, whose `mde_add_sign` is read per block) with a per-target block
prior (467's) and a per-target complement sign vector. Nothing of 496's is written: outputs go to work/h1_<mode>/.

MODE = <block>-<hook>[-r<r>][-q<q>]
  block : c4 (CLEAN4 everywhere = 435) | g33 (Arc on 467's strong 57, CLEAN4 on the weak 28 = 467's gated33) | arc (462l)
          | mz (sign of meanz_ne on every target) | mzg (meanz_ne on strong, CLEAN4 on weak)
  hook  : none | c4t<tau> (496's CLEAN4 hook) | arct<tau> (Arc, every target) | arcg<tau> (Arc on strong, CLEAN4 .05 on weak)
          | mzt<tau> / mzg<tau> (517's meanz_ne estimator: mean z-scored delta over Arc's non-essential rows, as a pseudo-P
            whose tau-.05 coverage equals Arc's P_down's on the gated axis; every target / gated)
  -r<r> : Arc with the fixed random fraction r of genes' signs flipped (grid517's mask), in the block and the hook alike
  -q<q> : protect_q (default 0.8), e.g. q0.7
usage: uv run python experiments/517-flex-programme/run517.py <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 \
         --mde-min-count 0.3 --mde-coarse
Known answers (drive517.sh checks them before anything is read): c4-none == 435's production run (stats bit-equal, seed 0);
g33-none == 467's gated33 (seed 0); c4-c4t0.05 == 496's hook run (seed 0).
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")      # before any import that sets it (run400/run496)
os.environ["RUN496_EMITTER"] = "hook"
MODE = sys.argv.pop(1)
m_ = re.fullmatch(r"(c4|g33|arc|mz|mzg)-(none|c4t[0-9.]+|arct[0-9.]+|arcg[0-9.]+|mzt[0-9.]+|mzg[0-9.]+)(?:-r([0-9.]+))?(?:-q([0-9.]+))?", MODE)
assert m_, MODE
BLOCK, HOOK = m_.group(1), m_.group(2)
R_FLIP = float(m_.group(3) or 0.0)
PQ = float(m_.group(4) or 0.8)
sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
import numpy as np  # noqa: E402
import run496 as M  # noqa: E402  (read-only; its WORK is redirected below)

WORK = HERE / "work" / (f"h1_{MODE}" + ("_zs" if "--zs" in sys.argv else ""))   # a zs run re-emits its pre-pass blocks
WORK.mkdir(parents=True, exist_ok=True)
M.WORK = WORK
M.H.WORK = WORK
M.PROTECT_Q = PQ
_clean4 = M.H.prior_signs_full
REF = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
H1_TARGETS = set(REF["targets"].astype(str))
KEEP = np.asarray(REF["keep"], dtype=bool)
ST = json.loads((ROOT / "experiments/467-strength-gated-programme/work/strength467.json").read_text())
_vals = np.array([ST["per_target"][str(i)][ST["chosen_measure"]] for i in range(len(ST["per_target"]))], float)
STRONG = set(np.flatnonzero(_vals >= float(np.quantile(_vals, 0.33))).tolist())     # 467's q33 (= le1 on H1): 57
PRI: dict[str, np.ndarray] = {}
HK: dict[str, np.ndarray] = {}


def arc_pdown_full(genes_full):
    """Arc's P_down on the full H1 axis, h1-test's rows out (467 / 462l's rule), with grid517's flip mask applied."""
    d = ROOT / "data/donors/h1_2025_full"
    delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
    genes = [l.strip() for l in open(d / "genes.txt")]
    rows = [i for i, p in enumerate(perts) if p not in H1_TARGETS]
    X = np.asarray(delta[rows], dtype=np.float64); nz = (X != 0).sum(0); dn = (X < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: j for j, g in enumerate(genes)}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    P = np.where(col >= 0, pd[np.maximum(col, 0)], np.nan)
    if R_FLIP > 0:
        # grid517 / hookgrid517's mask: one uniform per GATED gene (default_rng(51717)), flip where u < r
        u = np.random.default_rng(51717).random(int(KEEP.sum()))
        Pk = P[KEEP]
        Pk = np.where(u < R_FLIP, 1.0 - Pk, Pk)
        P = P.copy(); P[KEEP] = Pk
    return P


def meanz_ne_P_full(genes_full, P_arc):
    """517's estimator (variant517 / artifact517): mean over Arc's non-essential, non-h1-test rows of delta / sd (sd over all
    non-h1-test rows), as a pseudo-P with Arc P_down's tau-.05 coverage on the GATED axis (variant517's matching)."""
    ess = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
    d = ROOT / "data/donors/h1_2025_full"
    delta = np.load(d / "delta.npy"); perts = [l.strip() for l in open(d / "perts.txt")]
    genes = [l.strip() for l in open(d / "genes.txt")]
    rows = [i for i, p in enumerate(perts) if p not in H1_TARGETS]
    X = np.asarray(delta[rows], dtype=np.float64)
    ne = np.array([perts[i] not in ess for i in rows])
    sd = X.std(0); sd[sd == 0] = 1.0
    zb = (X[ne] / sd).mean(0)
    pos = {g: j for j, g in enumerate(genes)}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    z = np.where(col >= 0, zb[np.maximum(col, 0)], np.nan)
    n_cov = int((np.isfinite(P_arc[KEEP]) & (np.abs(P_arc[KEEP] - 0.5) >= 0.05)).sum())
    z_q = float(np.sort(np.abs(np.nan_to_num(z[KEEP])))[::-1][n_cov - 1])
    return np.where(np.isfinite(z), np.clip(0.5 - 0.05 * z / z_q, 0.001, 0.999), np.nan)


def clean4_P_full(genes_full):
    """478 / price517's CLEAN4 P (mean of the four banks' P_down), NaN where uncovered."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.full(len(genes_full), np.nan)
    out[col >= 0] = P[col[col >= 0]]
    return out


def prior_signs_full(genes_full):
    base, cov = _clean4(genes_full)
    Pa = arc_pdown_full(genes_full)
    PRI["c4"] = base
    PRI["arc"] = np.where(np.isfinite(Pa), np.where(Pa > 0.5, -1.0, 1.0), -1.0)
    Pz = meanz_ne_P_full(genes_full, arc_pdown_full(genes_full) if R_FLIP > 0 else Pa) if (BLOCK.startswith("mz") or HOOK.startswith("mz")) else None
    if Pz is not None:
        PRI["mz"] = np.where(np.isfinite(Pz), np.where(Pz > 0.5, -1.0, 1.0), -1.0)
    Pc = clean4_P_full(genes_full)

    def hs(P, tau):
        s = np.where(np.isfinite(P) & (np.abs(P - 0.5) >= tau), np.where(P > 0.5, -1.0, 1.0), 0.0)
        return np.ascontiguousarray(s[KEEP])
    if HOOK.startswith("c4t"):
        HK["strong"] = HK["weak"] = hs(Pc, float(HOOK[3:]))
    elif HOOK.startswith("arct"):
        HK["strong"] = HK["weak"] = hs(Pa, float(HOOK[4:]))
    elif HOOK.startswith("arcg"):
        HK["strong"], HK["weak"] = hs(Pa, float(HOOK[4:])), hs(Pc, 0.05)
    elif HOOK.startswith("mzt"):
        HK["strong"] = HK["weak"] = hs(Pz, float(HOOK[3:]))
    elif HOOK.startswith("mzg"):
        HK["strong"], HK["weak"] = hs(Pz, float(HOOK[3:])), hs(Pc, 0.05)
    print(f"[517] mode {MODE}: block {BLOCK}, hook {HOOK}, r {R_FLIP}, protect_q {PQ}; strong targets {len(STRONG)}; "
          f"hook signs strong {int((HK.get('strong', np.zeros(1)) != 0).sum())} weak {int((HK.get('weak', np.zeros(1)) != 0).sum())}",
          flush=True)
    return base, cov


_carrier = M.carrier_block


def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    strong = i in STRONG
    M.H.PRIOR_FULL = {"c4": PRI["c4"], "arc": PRI["arc"], "g33": PRI["arc"] if strong else PRI["c4"],
                      "mz": PRI.get("mz"), "mzg": PRI.get("mz") if strong else PRI["c4"]}[BLOCK]
    M.MDE_ADD_SIGN = None if HOOK == "none" else (HK["strong"] if strong else HK["weak"])
    return _carrier(profile, libs, emit_rng, spread, key, i, emit_c)


M.H.prior_signs_full = prior_signs_full
M.carrier_block = carrier_block

if __name__ == "__main__":
    M.main()
