"""241 — `mse` for every ablated component, through 096's `claim` route, under each base.

`mse200.py` with two changes and no others:
  1. the renderer's `mass_fix` flag is a property of the BASE, and three bases are run —
     `MF` (200's, `minnorm`), `NOMF` (the shipped renderer) and `SHIPPED` (NOMF at the
     shipped CD4 dose f = 0.34);
  2. every ratio is taken against ITS OWN base, not against 200's `SHIP`, because a ratio
     against a base that is not the one being ablated is not a leave-one-out.

Rule 13: `mse` is priceable only through 096's `claim` route.  EMISSION arms leave the
per-gene panel totals bit-identical and are priced ABSOLUTELY; GEOMETRY arms move `A`, `B`,
`claim` and `X`, so `X` is held at the rX = 1 null (169 §3) and A/B/claim are transferred as
RATIOS — a PARTIAL price, labelled as one, exactly as 200 labelled it.

KNOWN-ANSWER GATES (abort): 164 §1's four `claim` values, the A/B round trip, and a15's raw
`mse` 0.88228 / `score_mse` 0.11112.

    BASES=MF,NOMF,SHIPPED uv run python mse241.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E200 = os.path.join(ROOT, "experiments/200-production-ablation")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/162-retrodicting-validator")
sys.path.insert(0, E200)
sys.path.insert(0, HERE)

from cell_eval2.moments import jackknife_correction          # noqa: E402
import mech_alpha_chal as M                                  # noqa: E402
from mse200 import (AB, S15, SCORE, claim_of, gate_of, geometry_effects,  # noqa: E402
                    mass_fix_row, KEEP_SPREAD, KEEP_TOP, PHI, TAIL)
from pds241 import MF, NOMF, SHIPPED, loo_arms               # noqa: E402

STATS = f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz"
os.makedirs(WORK, exist_ok=True)

# which LOO arms to price per base.  `BASE-dup` is a pds control only.
FULL = ["BASE", "noWhiten", "noPM", "noCD4", "alpha_v11", "weights_4_1_1", "noGate",
        "noPrune", "noMassFix", "addMassFix"]
# on the shipped CD4 dose only the rows whose `mse` is NOT censored need their own render;
# every censored row's leg is exactly -base_score/6 whatever its own rA/rB.
SMALL = ["BASE", "noPM", "noCD4", "noPrune", "addMassFix"]
PLAN = {"MF": ["BASE", "noWhiten", "noPM", "noCD4", "noPrune", "noMassFix", "alpha_v11",
               "weights_4_1_1", "noGate"],
        "NOMF": FULL, "SHIPPED": FULL}
BASE_CFG = {"MF": MF, "NOMF": NOMF, "SHIPPED": SHIPPED}

# 200's published geometry column (work/mse200.json), for the Phase-0 gate
KA_GEOM = {"noWhiten": -0.018520, "noPM": +0.000256, "noCD4": -0.000024,
           "alpha_v11": -0.018520, "noPrune": +0.000018, "noMassFix": -0.000073}


def main() -> None:
    t0 = time.time()
    bases = os.environ.get("BASES", "MF,NOMF,SHIPPED").split(",")
    D = M.load(STATS)
    cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
    sD, A0, B0 = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
    Kc, X0, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
    print(f"board solve: sumD {sD:.6f}  X {X0:.6f}  Kc {Kc:.6f}  cap K {K:.5f}")
    print(f"v11 A {A0:.6f}  B {B0:.6f} | a15 amplitude s = {S15:.5f}, "
          f"s^2 B = {S15**2 * B0:.5f}")
    raw = {c: D[c]["dh"] for c in M.CTX}

    print("\n--- known-answer gates (164 §1, 169 §3) ---")
    GATES = [("flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860, 0.03),
             ("flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216, 0.03),
             ("flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478, 0.03),
             ("a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581, 0.03)]
    for name, tail, kt, ks, pk, pt, want, tol in GATES:
        got, ch, cp = claim_of(D, raw, tail, kt, ks, pk, pt, cpn)
        ok = abs(got - want) / want < tol
        print(f"  {name:<24} claim {got:.5f}  expect {want:.5f}  "
              f"({100 * (got - want) / want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            sys.exit(f"ABORT: claim gate {name}")
    A_chk, B_chk = AB(D, raw, 1.0)
    print(f"  A/B round trip at s=1: A {A_chk:.6f} ({100*(A_chk/A0-1):+.3f}%)  "
          f"B {B_chk:.6f} ({100*(B_chk/B0-1):+.3f}%)")
    assert abs(A_chk / A0 - 1) < 5e-3 and abs(B_chk / B0 - 1) < 5e-3
    a15_claim, _, _ = claim_of(D, raw, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
    a15_raw = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, a15_claim)) / sD
    print(f"  a15 raw mse {a15_raw:.5f} (164: 0.88228)  score_mse {SCORE(a15_raw):.5f} "
          f"(164: 0.11112)")
    assert abs(a15_raw - 0.88228) < 3e-4 and abs(SCORE(a15_raw) - 0.11112) < 3e-4

    out_path = os.path.join(WORK, "mse241.json")
    out = json.loads(open(out_path).read()) if os.path.exists(out_path) else {}
    out["anchor"] = dict(sumD=sD, A0=A0, B0=B0, Kc=Kc, X0=X0, K=K, S15=S15,
                         a15_claim=a15_claim, a15_raw=a15_raw, a15_score=SCORE(a15_raw))

    # ------------------------------------------------------------- EMISSION-law ablations
    # These are base-independent in `claim` (the emission law is read off the same rendered
    # profile), but their SCORE depends on the base's own raw `mse`, so both are recorded.
    if "emission" not in out:
        EM = [("SHIP  tier8000 phi4", TAIL, KEEP_TOP, KEEP_SPREAD, PHI),
              ("-tier (flat 0.25) phi4", None, 0, 0.25, PHI),
              ("-phi  tier8000 phi0", TAIL, KEEP_TOP, KEEP_SPREAD, 0.0)]
        erec = {}
        print("\n--- EMISSION-law ablations (profile bit-identical: A, B, X pinned) ---")
        for label, tail, kt, ks, phi in EM:
            cl, ch, cp = claim_of(D, raw, tail, kt, ks, phi, phi, cpn)
            erec[label] = dict(claim=cl, chat=ch, cp=cp)
            print(f"  {label:<30} claim {cl:>9.5f}  C_hat {ch:>7.2f}  [{time.time()-t0:.0f}s]",
                  flush=True)
        out["emission"] = erec
        json.dump(out, open(out_path, "w"), indent=1, default=float)

    # ------------------------------------------------------------- GEOMETRY ablations
    import chalparts200
    from plolo import PROD_OPS, Operator
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    tgt_index = {t: i for i, t in enumerate(p["targets"])}
    order = {c: np.asarray([tgt_index[t] for t in D[c]["tgt"]], dtype=np.int64)
             for c in M.CTX}

    for tag in bases:
        arms = loo_arms(BASE_CFG[tag], tag)
        names = [f"{tag}:{a}" for a in PLAN[tag] if f"{tag}:{a}" in arms]
        rec = out.setdefault(f"geometry:{tag}", {})
        todo = [n for n in names if n not in rec]
        if not todo:
            print(f"\n[{tag}] geometry already done")
            continue
        by_op: dict[str, list[str]] = {}
        for a in todo:
            by_op.setdefault(arms[a]["op"], []).append(a)
        print(f"\n--- GEOMETRY, base {tag} (rX = 1 null; A/B/claim as RATIOS vs {tag}:BASE) ---")
        for opkey, names_op in by_op.items():
            op = Operator(PROD_OPS / f"whitener_{opkey}.npz")
            for a in names_op:
                cfg = arms[a]
                E = geometry_effects(cfg, p, c4, op)
                eff = {}
                for c in M.CTX:
                    b = D[c]["bctrl"]
                    g = gate_of(b, cfg["min_cpm"])
                    h = np.where(g[None, :], E[order[c]], 0.0)
                    if cfg["mass"] is not None:
                        h = np.vstack([mass_fix_row(h[i], b, g, cfg["mass"])
                                       for i in range(h.shape[0])])
                    eff[c] = h
                Aq, Bq = AB(D, eff, 1.0)
                cl, ch, cp = claim_of(D, eff, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn,
                                      min_cpm=cfg["min_cpm"])
                rec[a] = dict(A=Aq, B=Bq, claim=cl, chat=ch)
                print(f"    built {a:<22} A {Aq:.5f} B {Bq:.5f} claim {cl:.5f} "
                      f"[{time.time()-t0:.0f}s]", flush=True)
                json.dump(out, open(out_path, "w"), indent=1, default=float)
            del op
    json.dump(out, open(out_path, "w"), indent=1, default=float)
    print(f"\nwrote {out_path} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    sys.exit(main())
