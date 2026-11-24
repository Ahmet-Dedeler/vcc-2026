"""193 Phase 3 — `mse` for the operator maps, through 096's `claim` route (rule 13's instrument).

Rule 12 exists because of a19: a whitener ridge changes the effect vector's SHAPE, `mse` reads
the profile directly, and a19's two-surface screen never computed it and lost `-0.01352`.  This
is that computation, for exactly the family a19 belonged to.

162 §3 measured that a cross-line validator cannot price `mse` at all (LOO RMS 0.112 against a
predict-zero null of 0.044, permutation p = 0.96), so PLOLO's `nmse_row` is not admissible and
096's decomposition is used instead:

    mse - 1 = ( s^2 A - 2 s X + Kc - min(s^2 B, claim) ) / sumD

with `X` held at the `rX = 1` null, which is what rule 13 credits a shape-only arm.

THE OBJECT AND THE OPERATOR (rule 11)
-------------------------------------
096's cached object is **v11's delivered effect** on the real challenge panel, and v11 whitened
with `data/donors/whitener_union.npz` (sha `5608e807`) -- NOT with a15's `union_pruned`.  The
operator ratio `V diag(f_alt/f_prod) V^T` is only meaningful against the operator that did the
whitening, so this script uses `union`, asserts its sha, and reports the arms as the map's
effect measured at v11's operator.  The two operators share the union axis and differ by
091's essential prune, so this is the map's effect on the same family of spectra, not on a
different object; it is nonetheless an extrapolation to a15's operator and is labelled as one.

`claim188.py`'s four known-answer `claim` gates are run verbatim and abort the script.
"""
from __future__ import annotations

import gc
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{REPO}/src")
sys.path.insert(0, f"{REPO}/experiments/096-mse-floor")
sys.path.insert(0, f"{REPO}/experiments/188-magnitude-calibration")
sys.path.insert(0, HERE)
import mech_alpha_chal as M                              # noqa: E402
from claim188 import AB, claim_of, KEEP_SPREAD, KEEP_TOP, MIN_CPM, PHI, S15, SCORE, TAIL  # noqa: E402,E501
import maps193                                           # noqa: E402
from screen193 import sha256_file                        # noqa: E402
from pathlib import Path                                 # noqa: E402

STATS = f"{REPO}/experiments/096-mse-floor/work/v11-stats.npz"
OP = Path(f"{REPO}/data/donors/whitener_union.npz")
OP_SHA = "5608e807230dea04"          # v11's operator; full digest asserted below
N_BASIS = 33787
GAMMA = 0.65
MAPS = ("prod(lam=0)", "ridge1e-06", "ridge1e-05", "ridge0.0001", "ridge0.007",
        "floorrel1e-06", "cliprank12000", "clipMP", "rie-LP-eta0.5", "tail^0.8",
        "tail^1.2", "head^1.2")


def main() -> None:
    dig = sha256_file(OP)
    assert dig.startswith(OP_SHA), f"operator sha {dig[:16]} != {OP_SHA}"
    D = M.load(STATS)
    cpn = dict(np.load(f"{REPO}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(f"{REPO}/experiments/096-mse-floor/work/board_solve.npz")
    sD, A0, B0 = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
    Kc, X0 = float(bs["v11_Kc"]), float(bs["X"])
    gate = {}
    for c in M.CTX:
        b = D[c]["bctrl"]
        expr = np.expm1(b)
        cpm = 1e6 * expr / expr.sum()
        gate[c] = (cpm > MIN_CPM) & (b > 0)
    raw = {c: D[c]["dh"] for c in M.CTX}

    print("--- claim188's known-answer gates, verbatim ---")
    bad = []
    for name, tail, kt, ks, pk, pt, want in (
            ("flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860),
            ("flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216),
            ("flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478),
            ("a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581)):
        got, _ = claim_of(D, raw, tail, kt, ks, pk, pt, cpn)
        ok = abs(got - want) / want < 0.03
        print(f"  {name:<24} claim {got:.5f}  expect {want:.5f} "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(name)
    if bad:
        sys.exit(f"ABORT: {bad}")

    wcols = np.load(f"{REPO}/experiments/162-retrodicting-validator/work/fold-kolf2.npz",
                    allow_pickle=True)["wcols"]
    with np.load(OP) as z:
        evals = np.asarray(z["evals"], float)
        V = np.asarray(z["evecs"], float)
    assert len(evals) == len(wcols), "operator axis vs wcols"
    fmaps = maps193.build_maps(evals, N_BASIS, gamma=GAMMA)
    fmaps.update(maps193.build_maps_b(evals, N_BASIS, gamma=GAMMA))
    f_prod = fmaps["prod(lam=0)"]
    print(f"\noperator {OP.name} sha {dig[:8]} p={len(evals)} n={N_BASIS}")

    # pre-project every context once; the ratio is all that varies
    proj, norm0 = {}, {}
    for c in M.CTX:
        blk = raw[c][:, wcols]
        norm0[c] = float(np.median(np.linalg.norm(blk, axis=1))) or 1.0
        proj[c] = blk @ V
    print(f"{'map':16s}{'A':>10}{'rA':>9}{'B':>10}{'claim':>9}{'raw mse':>10}"
          f"{'d raw':>10}{'score':>9}{'d overall':>11}{'branch':>13}")
    rec, base = {}, None
    for name in MAPS:
        ratio = fmaps[name] / f_prod
        eff = {}
        for c in M.CTX:
            t = (proj[c] * ratio) @ V.T
            t *= norm0[c] / (float(np.median(np.linalg.norm(t, axis=1))) or 1.0)
            e = raw[c].copy()
            e[:, wcols] = t
            eff[c] = np.where(gate[c], e, 0.0)
            del t, e
        cl, ch = claim_of(D, eff, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
        Aq, Bq = AB(D, eff, 1.0)
        mse = 1 + (S15**2 * Aq - 2 * S15 * X0 + Kc - min(S15**2 * Bq, cl)) / sD
        sc = SCORE(mse)
        if base is None:
            base = (Aq, cl, mse, sc)
        rec[name] = dict(A=Aq, rA=Aq / base[0], B=Bq, claim=cl, chat=ch, raw_mse=mse,
                         score=sc, d_raw=mse - base[2], d_overall=(sc - base[3]) / 6)
        print(f"{name:16s}{Aq:>10.5f}{Aq/base[0]:>9.5f}{Bq:>10.5f}{cl:>9.5f}{mse:>10.5f}"
              f"{mse-base[2]:>+10.5f}{sc:>9.5f}{(sc-base[3])/6:>+11.6f}"
              f"{('s^2 B' if S15**2*Bq <= cl else 'CLAIM CLIFF'):>13}", flush=True)
        del eff
        gc.collect()
    with open(os.path.join(WORK, "mse193.json"), "w") as fh:
        json.dump(dict(operator=str(OP.relative_to(Path(REPO))), op_sha256=dig, arms=rec),
                  fh, indent=1, default=float)
    print(f"\nwrote {WORK}/mse193.json")


if __name__ == "__main__":
    sys.exit(main())
