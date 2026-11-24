"""250 Phase 3 — the projector's cross-panel identity, on the PRODUCTION-shaped objects.

For each panel the CD4 residual as the builder holds it on the day (all 300 panel rows, the
panel's own re-extracted banks, complete-mean3, counts kept): fit V_k (gene-space top-k
right-singular directions of the centred carried rows) and report
  * principal angles between V_k(abc) and V_k(panel), and between the disjoint panels;
  * the fraction of each panel's centred energy captured by its OWN top-k basis vs by the
    OTHER panel's basis (the "top-5 energy fraction" the D/E/F census would read);
  * the spectrum shape (top-15 eigenvalue fractions) per panel;
  * fold-level bases (saved by sweep250) against the production bases.

KNOWN-ANSWER ROWS (abort on failure):
  KA-P1  250's A/B/C production-shaped residual == 249's staged `cd4-off.npz` bit-for-bit
         (mean_effect, gene_votes, target_condition_count, target order).
  KA-P2  `lib250.deflate(svd, 5)` on it reproduces 249's shipped deflation audit:
         n_carried 291, energy_removed_frac 0.19866046506100588,
         top_k_energy_frac_of_centred 0.1982359547102955, basis_u_k_sha256 f06e0262...

    /usr/bin/time -l uv run python basis250.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gate250 import wait_for_memory  # noqa: E402
import lib250 as L  # noqa: E402

PANELS = ("abc", "real", "pess", "efp", "efr")
KS = (3, 5, 10, 15)
SHIPPED = dict(n_carried=291, energy_removed_frac=0.19866046506100588,
               top_k_energy_frac_of_centred=0.1982359547102955,
               basis_u_k_sha256="f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c")


def main() -> None:
    t0 = time.time()
    print(f"[basis] gate clear: {wait_for_memory(tag=' basis')}", flush=True)
    res, bases, audits, centred = {}, {}, {}, {}
    for panel in PANELS:
        L.S.bind_panel(panel, L.resid218)
        targets = L.B.read_panel(panel)
        r = L.S.cd4_prepared(panel, targets)
        res[panel] = r
        rows, R, Rc, sup = L.centred_rows(r)
        centred[panel] = (rows, Rc)
        for k in KS:
            V, U, aud = L.gene_basis(r, k)
            bases[(panel, k)] = V
            audits[(panel, k)] = aud
        n5 = audits[(panel, 5)]
        print(f"[basis/{panel}] carried {n5['n_carried']}/300, top-5 energy "
              f"{n5['top_k_energy_frac_of_centred']:.4f}, top-3 "
              f"{audits[(panel, 3)]['top_k_energy_frac_of_centred']:.4f}, top-10 "
              f"{audits[(panel, 10)]['top_k_energy_frac_of_centred']:.4f}, spectrum "
              f"{[round(x, 4) for x in n5['spectrum_frac'][:8]]}", flush=True)

    # ---------------------------------------------------------------- known answers
    prod, prod_targets = L.production_cd4_abc()
    mine = res["abc"]
    ka = {}
    ka["KA-P1"] = dict(
        targets_equal=bool(tuple(L.B.read_panel("abc")) == prod_targets),
        mean_effect_equal=bool(np.array_equal(np.asarray(mine.mean_effect), np.asarray(prod.mean_effect))),
        votes_equal=bool(np.array_equal(np.asarray(mine.gene_votes), np.asarray(prod.gene_votes))),
        counts_equal=bool(np.array_equal(np.asarray(mine.target_condition_count),
                                         np.asarray(prod.target_condition_count))),
        max_abs_diff=float(np.abs(np.asarray(mine.mean_effect) - np.asarray(prod.mean_effect)).max()))
    d5, a5 = L.deflate(mine, "svd", 5)
    ka["KA-P2"] = dict(measured={k: a5[k] for k in ("n_carried", "energy_removed_frac",
                                                     "top_k_energy_frac_of_centred",
                                                     "basis_u_k_sha256")},
                       shipped=SHIPPED,
                       d_energy=abs(a5["energy_removed_frac"] - SHIPPED["energy_removed_frac"]),
                       d_topk=abs(a5["top_k_energy_frac_of_centred"]
                                  - SHIPPED["top_k_energy_frac_of_centred"]),
                       sha_equal=a5["basis_u_k_sha256"] == SHIPPED["basis_u_k_sha256"])
    ok1 = all(v for k, v in ka["KA-P1"].items() if k != "max_abs_diff")
    ok2 = (a5["n_carried"] == 291 and ka["KA-P2"]["d_energy"] < 1e-12
           and ka["KA-P2"]["d_topk"] < 1e-12)
    print(f"[basis] KA-P1 production residual bitwise: {ka['KA-P1']}  -> {'PASS' if ok1 else 'FAIL'}")
    print(f"[basis] KA-P2 shipped deflation audit: {json.dumps(ka['KA-P2'], indent=None)}  -> "
          f"{'PASS' if ok2 else 'FAIL'} (U_k sha {'EQUAL' if ka['KA-P2']['sha_equal'] else 'DIFFERS'})")
    if not (ok1 and ok2):
        raise SystemExit("[basis] KNOWN-ANSWER FAILURE")

    # ---------------------------------------------------------------- angles + cross energy
    pairs = [(a, b) for i, a in enumerate(PANELS) for b in PANELS[i + 1:]]
    angles, cross = {}, {}
    for k in KS:
        for a, b in pairs:
            angles[f"{a}|{b}|{k}"] = L.principal_angles_deg(bases[(a, k)], bases[(b, k)])
        for p in PANELS:
            rows, Rc = centred[p]
            ec = float((Rc * Rc).sum())
            own = float(((Rc @ bases[(p, k)]) ** 2).sum() / ec)
            for q in PANELS:
                proj = Rc @ bases[(q, k)]
                cross[f"{p}|by_{q}|{k}"] = dict(captured=float((proj * proj).sum() / ec), own=own,
                                                ratio=float((proj * proj).sum() / ec / own))
    # random-basis floor for the cross-capture statistic: k random orthonormal directions in
    # the panel's own support capture k / n_support of the energy in expectation
    floor = {}
    rng = np.random.default_rng(L.RAND_SEED)
    for p in PANELS:
        rows, Rc = centred[p]
        sup = np.flatnonzero(np.abs(Rc).sum(axis=0) > 0)
        ec = float((Rc * Rc).sum())
        vals = []
        for _ in range(3):
            Q, _r = np.linalg.qr(rng.standard_normal((sup.size, 5)))
            proj = Rc[:, sup] @ Q
            vals.append(float((proj * proj).sum() / ec))
        floor[p] = dict(n_support=int(sup.size), random5_captured=float(np.mean(vals)))

    # ---------------------------------------------------------------- fold-level bases
    fold_angles = {}
    for f in sorted(L.WORK.glob("basis250-*-*.npz")):
        _, panel, held = f.stem.split("-", 2)
        z = np.load(f)
        V5 = z["V5"]
        fold_angles[f"{panel}|{held}"] = dict(
            vs_own_production=L.principal_angles_deg(V5, bases[(panel, 5)]),
            vs_abc_production=L.principal_angles_deg(V5, bases[("abc", 5)]),
            audit=json.loads(str(z["audit5"])))

    for k in (5,):
        print(f"\nPRINCIPAL ANGLES (deg), V_{k}, production-shaped residuals")
        for a, b in pairs:
            print(f"  {a:5}|{b:5}  " + "  ".join(f"{x:6.1f}" for x in angles[f'{a}|{b}|{k}']))
        print(f"\nENERGY CAPTURED by a {k}-dim basis (fraction of the panel's centred energy)")
        print(f"  {'panel':6}" + "".join(f"{'by_' + q:>10}" for q in PANELS) + f"{'rand5':>10}")
        for p in PANELS:
            print(f"  {p:6}" + "".join(f"{cross[f'{p}|by_{q}|{k}']['captured']:10.4f}" for q in PANELS)
                  + f"{floor[p]['random5_captured']:10.4f}")
    out = dict(ka=ka, audits={f"{p}|{k}": v for (p, k), v in audits.items()},
               angles=angles, cross=cross, random_floor=floor, fold_angles=fold_angles,
               seconds=round(time.time() - t0, 1), peak_rss_gib=L.peak_gib())
    (L.WORK / "basis250.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(L.WORK / "basis250-production.npz",
                        **{f"V{k}_{p}": bases[(p, k)] for (p, k) in bases})
    print(f"\n[basis] done in {time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB -> work/basis250.json")


if __name__ == "__main__":
    main()
