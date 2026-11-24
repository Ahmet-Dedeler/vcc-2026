"""232 Phase 0 — the mandatory known-answer rows.  Nothing downstream runs if one fails.

KA-A  228's `nmae` decomposition and its efficiency statistic `e = dM/dU = 0.072` for 224's
      per-gene reliability gain, recomputed from 228's own stored `adel228.json`
      (three cancer folds, `@adel`, `wald89`).  This is the instrument §8.1d is written on
      and it is REUSED here, not reinvented.

KA-B  Our own PLOLO baseline against the shipped `h1union` stack at **Delta = 0.0e+00**:
      `cd4`-only and `base222` (CD4 @ 0.34 + H1 union @ 0.28) must reproduce 223's
      `focus223` numbers on all four folds through THIS code path, to 1e-12.

KA-C  204's finding that EQUAL donor weighting maximises the transfer cosine.  Two rungs:
      (i) 204's absolute equal-weight `pds` per fold at 204's own geometry
          (alpha 0.982, gamma 0.65, lam 0, CD4 @ 0.10, no H1 union, a23 render);
      (ii) 204's weighted arms re-scored on the same geometry -- `equal` must be the
          maximum of the macro cosine over the alternatives, with `share_4_1_1`
          (2.0, 0.5, 0.5) losing ~6.5% relative cosine.

    uv run python phase0_232.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common232 import (ALL_FOLDS, CANCER, KA_204_PDS, KA_STACK, WORK,  # noqa: E402
                       FoldBundle, load_sources, weighted_baseline)

ROOT = HERE.parents[1]
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")

# 204 §4.1's arms, the scalar ones nearest and furthest from the incumbent
W204 = {
    "equal": (1.000, 1.000, 1.000),
    "k115": (1.095, 0.952, 0.952),
    "h115": (0.952, 0.952, 1.095),
    "k085": (0.895, 1.053, 1.053),
    "share_1_2_2": (0.600, 1.200, 1.200),
    "share_4_1_1": (2.000, 0.500, 0.500),
    "a_self": (1.246, 1.032, 0.722),
    "minvar": (0.618, 1.208, 1.174),
}
# 204 RESULTS.md §4.1, `Dcos` column (fold-mean macro cosine vs `equal`)
KA_204_DCOS = {"k115": -0.000045, "h115": -0.000008, "k085": -0.000002,
               "share_1_2_2": -0.000367, "share_4_1_1": -0.001429,
               "a_self": -0.000327, "minvar": -0.000328}


def ka_a() -> dict:
    """228's `e = dM/dU`, from 228's stored numbers.  `nmae_U`/`nmae_M` are already the
    per-target means of `U_p/D_p` and `M_p/D_p`; the identity below verifies that."""
    d = json.loads((ROOT / "experiments/228-theta-gain-plolo/work/adel228.json").read_text())
    out = {}
    for arm in ("BASE", "gain_reg@0.25"):
        rows = [r for r in d if r["arm"] == arm and r["fold"] in CANCER]
        assert len(rows) == 3, f"{arm}: {len(rows)} cancer rows"
        U = float(np.mean([r["nmae_U"]["wald89"] for r in rows]))
        M = float(np.mean([r["nmae_M"]["wald89"] for r in rows]))
        nm = float(np.mean([r["nmae_pools"]["wald89"] for r in rows]))
        ident = 1.0 + U - 2.0 * M
        assert abs(ident - nm) < 5e-6, f"188 identity broken for {arm}: {ident} vs {nm}"
        out[arm] = dict(U_over_D=U, M_over_D=M, nmae=nm, identity=ident)
    dU = out["gain_reg@0.25"]["U_over_D"] - out["BASE"]["U_over_D"]
    dM = out["gain_reg@0.25"]["M_over_D"] - out["BASE"]["M_over_D"]
    e = dM / dU
    dn = out["gain_reg@0.25"]["nmae"] - out["BASE"]["nmae"]
    ok = abs(e - 0.072) < 0.002 and abs(dn - 0.0373) < 0.0005
    out.update(dU=dU, dM=dM, e=e, d_nmae=dn, expected_e=0.072, expected_d_nmae=0.0373,
               PASS=bool(ok))
    print(f"[KA-A] 228 wald89, 3 cancer folds, @adel: U/D {out['BASE']['U_over_D']:.5f} -> "
          f"{out['gain_reg@0.25']['U_over_D']:.5f}  dU {dU:+.5f}   "
          f"M/D {out['BASE']['M_over_D']:.5f} -> {out['gain_reg@0.25']['M_over_D']:.5f} "
          f"dM {dM:+.5f}   e = {e:.4f} (228: 0.072)   d_nmae {dn:+.4f} (228: +0.0373)  "
          f"{'PASS' if ok else 'FAIL'}", flush=True)
    assert ok, "KA-A failed: 228's decomposition does not reproduce"
    return out


def ka_b(bundle: FoldBundle) -> dict:
    held = bundle.held
    out = {}
    for stack in ("cd4", "base222"):
        e = bundle.score(bundle.final(stack=stack))
        exp_cos, exp_pds = KA_STACK[(held, stack)]
        dc, dp = abs(e["cos"] - exp_cos), abs(e["pds"] - exp_pds)
        ok = dc < 1e-12 and dp < 1e-12
        print(f"[KA-B/{held}] {stack:8s} cos {e['cos']:.15f} (D {dc:.1e})  "
              f"pds {e['pds']:.15f} (D {dp:.1e})  {'PASS' if ok else 'FAIL'}", flush=True)
        assert ok, f"KA-B FAILED for {held}/{stack}"
        out[stack] = dict(cos=e["cos"], pds=e["pds"], cos_agg=e["cos_agg"],
                          sgn_w=e["sgn_w"], d_cos=dc, d_pds=dp, PASS=True)
    return out


def ka_c(bundle: FoldBundle, src: dict) -> dict:
    """204's geometry, 204's weight vectors, on the fold's own source set."""
    held = bundle.held
    names = [n.replace("_pm50", "") for n in src["names"]]
    out = {}
    base_pre = np.asarray(bundle.fold["pre"], dtype=np.float64)
    for arm, w3 in W204.items():
        wmap = dict(zip(DONORS, w3))
        w = [wmap[n] for n in names]
        pre = weighted_baseline(src["E"], src["avail"], w, fallback=base_pre)
        dev = float(np.max(np.abs(pre - base_pre)))
        if arm == "equal":
            assert dev < 1e-6, f"{held}: equal weighting is not 162's `pre` (dev {dev:.2e})"
        e = bundle.score(bundle.final_204(pre))
        out[arm] = dict(pds=e["pds"], a15_pds=e["a15_pds"], cos=e["cos"],
                        cos_agg=e["cos_agg"], sgn_w=e["sgn_w"], weights=w,
                        max_dev_from_162_pre=dev)
        print(f"[KA-C/{held}/{arm:12s}] pds {e['pds']:.6f}  cos {e['cos']:.6f}  "
              f"dev {dev:.2e}", flush=True)
    exp = KA_204_PDS[held]
    got = out["equal"]["pds"]
    ok = abs(got - exp) < 5e-5
    print(f"[KA-C/{held}] equal-weight pds {got:.6f} vs 204's {exp:.5f} "
          f"(D {abs(got-exp):.2e})  {'PASS' if ok else 'FAIL'}", flush=True)
    out["_equal_pds_expected"] = exp
    out["_equal_pds_PASS"] = bool(ok)
    return out


def main() -> None:
    t0 = time.time()
    audit: dict = {"KA_A": ka_a(), "KA_B": {}, "KA_C": {}}
    out_path = WORK / "phase0_232.json"
    for held in ALL_FOLDS:
        b = FoldBundle(held)
        src = load_sources(held)
        assert tuple(src["targets"]) == tuple(b.targets), f"{held}: target axis mismatch"
        audit["KA_B"][held] = ka_b(b)
        audit["KA_C"][held] = ka_c(b, src)
        audit.setdefault("folds", {})[held] = dict(
            n_targets=len(b.targets), sources=src["names"],
            prep_ka_rel=src["ka_rel"], n_scored=int(b.sc.size),
            targets_with_no_covering_donor=int((src["votes"] == 0).sum()),
            target_coverage_hist={str(v): int((src["votes"] == v).sum())
                                  for v in np.unique(src["votes"])})
        out_path.write_text(json.dumps(audit, indent=1))
        del b, src

    # ---- KA-C rung (ii): is `equal` the maximum of the macro cosine?
    cos = {arm: float(np.mean([audit["KA_C"][f][arm]["cos"] for f in ALL_FOLDS]))
           for arm in W204}
    dcos = {a: cos[a] - cos["equal"] for a in W204 if a != "equal"}
    rel = {a: 100.0 * dcos[a] / cos["equal"] for a in dcos}
    n_below = sum(1 for v in dcos.values() if v < 0)
    print("\n[KA-C] fold-mean macro cosine, 204 geometry, all four folds")
    print(f"  equal          {cos['equal']:.6f}   (204 baseline)")
    for a in sorted(dcos, key=lambda x: dcos[x]):
        print(f"  {a:14s} {cos[a]:.6f}   Dcos {dcos[a]:+.6f} ({rel[a]:+.2f}%)   "
              f"204 reported Dcos {KA_204_DCOS[a]:+.6f}")
    ok = (n_below == len(dcos)) and (rel["share_4_1_1"] < -4.0)
    print(f"[KA-C] equal is the maximum over {len(dcos)} alternatives: "
          f"{n_below}/{len(dcos)} below;  share_4_1_1 rel {rel['share_4_1_1']:+.2f}% "
          f"(204: -6.53%)  {'PASS' if ok else 'FAIL'}", flush=True)
    audit["KA_C_summary"] = dict(cos=cos, dcos=dcos, rel_pct=rel, n_below=n_below,
                                 n_alternatives=len(dcos),
                                 reported_204=KA_204_DCOS, PASS=bool(ok))
    audit["seconds"] = time.time() - t0
    out_path.write_text(json.dumps(audit, indent=1))
    assert ok, "KA-C failed: equal weighting is not the cosine maximum here"
    print(f"\n-> {out_path}  ({time.time() - t0:.0f}s)  ALL KNOWN-ANSWER ROWS PASS")


if __name__ == "__main__":
    main()
