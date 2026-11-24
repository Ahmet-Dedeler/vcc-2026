"""190 — the three closing measurements.

  K  per-covered-target breakdown: the effect BEFORE dilution, and the dilution factor,
     so the coverage cap is stated in overall units and not as a per-covered-target effect.
  L  signed Δ`fid` for the regulon keys alone (comp190 §G's max was driven by the
     expression key, which is a different comparator).
  M  FALSIFIER 4b: is the five-member SHAPE of a CollecTRI direction correction specific
     to CollecTRI, or does ANY signed vector of the same L2 produce it?  If a shuffled
     network reproduces tydymy's cost side, then the cost side carries no information
     about the biology and the whole retrodiction rests on the benefit side, which does
     not reproduce at all.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(ROOT, "experiments", "172-sign-confidence-key"))
sys.path.insert(0, HERE)

from comp190 import reach_per_target  # noqa: E402
from key190 import agreement, key_mult  # noqa: E402
from run190 import LAM, build_all, to_overall  # noqa: E402


def main() -> None:
    import harness172

    T = harness172.H1Test()
    harness172.check(T, verbose=False)
    print("known-answer rows: OK (silent)")
    S = T.S
    H = T.H.astype(np.float64)
    absH = np.abs(H)
    pool = T.pool
    sT = np.sign(S.real_lfc)
    base = T.reach(absH)
    nets = build_all(T, H)
    out = {}

    # ---------------------------------------------------------------- K
    print("\n" + "=" * 100)
    print("K. PER-COVERED-TARGET BREAKDOWN — the number before dilution, and the dilution")
    print("=" * 100)
    pt0 = reach_per_target(S, absH)
    print(f"   {'network':22} {'cov tgts':>9} {'Δ on covered':>14} {'Δ on all 85':>13} "
          f"{'dilution':>10} {'panel dilution':>15}")
    kk = {}
    for name, N in nets.items():
        R = agreement(N, H)
        vals = [T.reach(key_mult(H, R, l)) for l in LAM]
        lam = LAM[int(np.argmax(vals))]
        pt = reach_per_target(S, key_mult(H, R, lam))
        cov = np.flatnonzero((N != 0).any(1))
        m = np.isfinite(pt) & np.isfinite(pt0)
        dcov = float(np.nanmean(pt[cov]) - np.nanmean(pt0[cov])) if cov.size else 0.0
        dall = float(np.mean(pt[m]) - np.mean(pt0[m]))
        kk[name] = dict(lam=lam, n_cov=int(cov.size), d_cov=dcov, d_all=dall)
        print(f"   {name:22} {cov.size:9d} {dcov:+14.5f} {dall:+13.5f} "
              f"{cov.size/85:10.3f} {'41/300 = 0.137' if 'ct-' in name else '300/300':>15}")
    out["per_covered_target"] = kk
    print("\n   CollecTRI covers 17/85 h1-test targets and 41/300 challenge panel targets,")
    print("   so any per-covered-target effect is multiplied by 0.137 on the board, not 0.20.")

    # ---------------------------------------------------------------- L
    print("\n" + "=" * 100)
    print("L. SIGNED Δ`fid` FOR THE REGULON KEYS ALONE")
    print("=" * 100)
    C = json.load(open(os.path.join(WORK, "comp190.json")))
    fid = C["fid"]
    ref = fid["|H|"]
    best_d = {}
    for kname, tab in fid.items():
        if kname in ("|H|", "expression"):
            continue
        ds = [tab[m][nr] - ref[m][nr] for m in tab for nr in tab[m] if nr != "precision"]
        best_d[kname] = {"max": max(ds), "min": min(ds)}
        print(f"   {kname:34} Δfid in [{min(ds):+.5f}, {max(ds):+.5f}]  "
              f"-> overall [{to_overall('fid', min(ds)):+.5f}, "
              f"{to_overall('fid', max(ds)):+.5f}]")
    print(f"   {'expression (the 189 comparator)':34} Δfid in "
          f"[{min(fid['expression'][m][nr] - ref[m][nr] for m in fid['expression'] for nr in fid['expression'][m] if nr != 'precision'):+.5f}, "
          f"{max(fid['expression'][m][nr] - ref[m][nr] for m in fid['expression'] for nr in fid['expression'][m] if nr != 'precision'):+.5f}]")
    out["fid_signed"] = best_d

    # ---------------------------------------------------------------- M
    print("\n" + "=" * 100)
    print("M. FALSIFIER 4b — is tydymy's COST side specific to CollecTRI?")
    print("   A shuffled network of identical topology and identical L2, same c = 0.25.")
    print("   10 seeds of each shuffle, against the real network.")
    print("=" * 100)
    P0 = S.pred_lfc.astype(np.float64)
    n0 = np.linalg.norm(P0, axis=1, keepdims=True)
    c = 0.25

    def five(N):
        u = N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)
        P = P0 + c * n0 * u
        sp = S.pred_lfc
        S.pred_lfc = P
        try:
            return dict(pds=S.pds(), nmae=S.nmae(), fid=S.fid(),
                        reach=S.reach(key=absH),
                        sgn=float(np.mean(np.sign(P)[pool] == sT[pool])))
        finally:
            S.pred_lfc = sp

    b = {"pds": S.pds(), "nmae": S.nmae(), "fid": S.fid(), "reach": base,
         "sgn": float(np.mean(np.sign(P0)[pool] == sT[pool]))}
    print(f"   {'correction':32} {'Δpds':>9} {'Δnmae':>9} {'Δfid':>9} {'Δreach':>9} {'Δsgnacc':>9}")
    mm = {}
    for name in ("ct-direct", "ct-prop-k25", "ct-prop-kNone"):
        r = five(nets[name])
        mm[name] = {"real": {k: r[k] - b[k] for k in b}}
        print(f"   {name + ' (real)':32} " +
              "".join(f"{r[k]-b[k]:+9.5f}" for k in ("pds", "nmae", "fid", "reach", "sgn")))
        for sh in ("targets", "signs", "tfs"):
            acc = []
            for s in range(10):
                rng = np.random.default_rng(2000 + s)
                Ns = build_all(T, H, rng=rng, shuffle=sh)[name]
                if not (Ns != 0).any():
                    continue
                acc.append(five(Ns))
            if not acc:
                continue
            mu = {k: float(np.mean([a[k] - b[k] for a in acc])) for k in b}
            sd = {k: float(np.std([a[k] - b[k] for a in acc], ddof=1)) for k in b}
            mm[name][sh] = {"mean": mu, "sd": sd}
            print(f"   {'  shuffle=' + sh:32} " +
                  "".join(f"{mu[k]:+9.5f}" for k in ("pds", "nmae", "fid", "reach", "sgn")))
    out["falsifier4b"] = mm
    out["blend_base"] = b

    with open(os.path.join(WORK, "final190.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nwrote {WORK}/final190.json")


if __name__ == "__main__":
    main()
