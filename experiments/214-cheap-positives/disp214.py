"""214 — how much bigger is PerturbMatch on h1-test than on the challenge?

Arm A's four DE members are measured on h1-test, and h1-test carries only TWO production
donors (`k562_gwps` + `orion_hek293t`; `h1_2025` is oracle leakage and `orion_hct116` is not
in 095's h1-test operator).  PerturbMatch enters the K562 slot only, so on a two-donor
consensus it owns 1/2 of the effect and on the shipped three-donor consensus it owns 1/3.
Every h1-test reading of PerturbMatch — 200's and this experiment's — therefore OVERSTATES
the shipped arm, and by how much is measurable rather than arguable.

Three measurements:
  1. the relative displacement `||E(pm=0.5) - E(pm=0)|| / ||E(pm=0.5)||` on both surfaces;
  2. the per-row median cosine between the two arms on both surfaces;
  3. PerturbMatch's target coverage on both surfaces.

    uv run python disp214.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
for p in ("experiments/162-retrodicting-validator", "experiments/200-production-ablation",
          "experiments/196-emitted-depth", "experiments/096-mse-floor", "src"):
    sys.path.insert(0, str(ROOT / p))
sys.path.insert(0, str(HERE))

import chalparts200                                                      # noqa: E402
import h1parts200                                                        # noqa: E402
import mse200                                                            # noqa: E402
from plolo import PROD_OPS, Operator                                     # noqa: E402
from pds214 import arm                                                   # noqa: E402


def rowstats(P, Q):
    d = P - Q
    rel = float(np.linalg.norm(d) / np.linalg.norm(P))
    per = np.linalg.norm(d, axis=1) / np.maximum(np.linalg.norm(P, axis=1), 1e-300)
    cs = ((P * Q).sum(1) / np.maximum(np.linalg.norm(P, axis=1) * np.linalg.norm(Q, axis=1),
                                      1e-300))
    return dict(rel_frobenius=rel, rel_row_median=float(np.median(per)),
                rel_row_min=float(per.min()), rel_row_max=float(per.max()),
                cos_median=float(np.median(cs)), cos_min=float(cs.min()),
                rows_moved=int((per > 1e-12).sum()), rows=int(P.shape[0]))


def main() -> None:
    out = {}
    # ---------------------------------------------------------------- challenge, 3 donors
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    op = Operator(PROD_OPS / "whitener_union_pruned.npz")
    Ec = {k: mse200.geometry_effects(arm(cd4_f=0.34, pm=v, mass=None), p, c4, op)
          for k, v in (("pm50", 0.50), ("pm00", 0.00))}
    del op
    out["challenge"] = rowstats(Ec["pm50"], Ec["pm00"])
    out["challenge"]["donors"] = list(p["src"])
    out["challenge"]["pm_available"] = int(np.asarray(p["pm_avail"]).sum())
    del Ec, p, c4

    # ---------------------------------------------------------------- h1-test, 2 donors
    ph = h1parts200.parts()
    c4h = h1parts200.cd4(ph["targets"])
    oph = h1parts200.H1Operator()
    art = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                  allow_pickle=True)
    norm = float(np.median(np.linalg.norm(art["H"].astype(np.float64), axis=1)))
    Eh = {k: h1parts200.effects_for(
        dict(pm=v, gamma=0.65, lam=0.0, cd4_f=0.34, cd4_frac=None, weights=None),
        ph, c4h, oph, norm) for k, v in (("pm50", 0.50), ("pm00", 0.00))}
    del oph
    out["h1test"] = rowstats(Eh["pm50"], Eh["pm00"])
    out["h1test"]["donors"] = list(ph["src"])
    out["h1test"]["pm_available"] = int(np.asarray(ph["pm_avail"]).sum())

    print("=== how far does dropping PerturbMatch move the effect matrix? ===")
    print(f"  {'surface':<12}{'donors':>8}{'targets':>9}{'PM-covered':>12}"
          f"{'||dE||/||E||':>14}{'row median':>12}{'row cos median':>16}")
    for k in ("challenge", "h1test"):
        r = out[k]
        print(f"  {k:<12}{len(r['donors']):>8}{r['rows']:>9}{r['pm_available']:>12}"
              f"{r['rel_frobenius']:>14.5f}{r['rel_row_median']:>12.5f}"
              f"{r['cos_median']:>16.5f}")
    ratio = out["h1test"]["rel_frobenius"] / out["challenge"]["rel_frobenius"]
    out["h1_over_challenge"] = ratio
    print(f"\n  h1-test / challenge displacement ratio = {ratio:.3f}")
    print("  (1/2 vs 1/3 of the consensus is 1.500; the measured ratio is the honest "
          "deflator for\n   every PerturbMatch DE number read on h1-test, this "
          "experiment's and 200's alike.)")
    json.dump(out, open(WORK / "disp214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'disp214.json'}")


if __name__ == "__main__":
    main()
