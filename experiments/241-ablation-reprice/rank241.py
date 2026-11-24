"""241 — the ten-row leave-one-out table, BEFORE (200's renderer) and AFTER (the shipped one).

Sign convention throughout, identical to 200's: **positive = what removing the component
would cost us**, in overall board units.

Members and their instruments (unchanged from 200 except where noted):
  pds     PLOLO, 4 donor folds, challenge frames A/B/C, scaled by the board `pds` span and
          rule 13's k = 1.011.  The 3-cancer mean and `kolf2` are carried SEPARATELY
          (falsifier 3) and the 4-fold mean is the headline only because 200 published it.
  mse     096's `claim` route.  200 read the SCORE off the a15-anchored board solve, whose
          whole `mse` score is 0.11112 scaled; the stack we ship reads **0.17957** on the
          board, so a censored row is worth **0.029928** overall, not 0.018520.  Both are
          printed; the shipped one is the headline, deflated by R3's measured 0.730 as the
          sensitivity.
  DE 4    h1-test, emitted, matched seeds, `debench.Surface`, from `six241.py`.

Rows h1-test CANNOT measure are left blank and named, never carried at zero silently: the
CD4 residual's DE four (coverage 4/85 -> solved fraction exactly 0), the pruned basis (h1-test
has one operator) and `VCC_ALPHA` (the rebuild renormalises amplitude away by construction).
For those the PLOLO cosine delta and 218's rate card give a RANGE, printed beside the row and
labelled as an unidentified conversion (falsifier 4).
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

# the board's own row for the stack we ship (`h1union`, entry iUziqE1gewfcPKzYUDaB)
BOARD = dict(overall=0.21183096391231482, raw_mse=0.8172630319018631,
             score_mse=0.17956879145192475, score_pds=0.7453739835254144,
             score_nmae=0.11074241388261237, score_fid=-0.016476121473095937,
             score_reach=0.22293924023183273, score_jac=0.02883747585520066,
             raw_reach=0.27654544083789157)
SCORE = lambda raw: min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)   # noqa: E731
R3_DEFLATOR = 0.730            # 200 R3: the claim route's cliff is 1.37x the board's
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089

# (label, pds arm key, mse key, DE row label in six241, class)
COMPONENTS = [
    ("union whitening operator (γ 0.65)", "noWhiten", "noWhiten",
     "union whitening operator (γ 0.65)", "shape"),
    ("equal donor weights (1,1,1)", "weights_4_1_1", "weights_4_1_1",
     "equal donor weights (1,1,1)", "shape"),
    ("--overdispersion 4.0", None, "EMISSION:-phi  tier8000 phi0",
     "--overdispersion 4.0", "shape"),
    ("tiered spread policy (8000/1.0)", None, "EMISSION:-tier (flat 0.25) phi4",
     "tiered spread policy (8000/1.0)", "shape"),
    ("VCC_ALPHA 0.982 (vs v11's 1.0796)", "alpha_v11", "alpha_v11", None, "amp"),
    ("CPM gate (--min-cpm 10)", "noGate", "noGate", "CPM gate (--min-cpm 10)", "shape"),
    ("post-whitening CD4 residual", "noCD4", "noCD4", "post-whitening CD4 residual", "amp"),
    ("a23 mass fix (minnorm)", "MASSFIX", "MASSFIX", "a23 mass fix (minnorm)", "shape"),
    ("pruned basis (union_pruned)", "noPrune", "noPrune", None, "shape"),
    ("PerturbMatch K562 replacement (0.5)", "noPM", "noPM",
     "PerturbMatch K562 replacement (0.5)", "shape"),
]
# rows whose DE four are a STRUCTURAL zero on h1-test, not a measurement
STRUCTURAL_ZERO = {"post-whitening CD4 residual":
                   "CD4 coverage 4/85 -> solve_cd4_fraction returns exactly 0; "
                   "de241.reb_noCD4_nomf is BIT-IDENTICAL to de241.reb_ship_nomf"}


def mse_leg(mse, tag, key, anchor, base_raw):
    """Overall-unit `mse` leg for one LOO row on one base.

    Emission rows are priced ABSOLUTELY on `claim` (rule 13's validated class).  Geometry rows
    are the rX = 1 partial price, as ratios against their OWN base.  Both produce a delta in
    RAW mse, which is then applied to `base_raw` and scored with the board's own map.
    """
    S15, A0, B0, X0, Kc, sD = (anchor["S15"], anchor["A0"], anchor["B0"], anchor["X0"],
                               anchor["Kc"], anchor["sumD"])
    a15_claim = anchor["a15_claim"]
    if key.startswith("EMISSION:"):
        lab = key.split(":", 1)[1]
        em = mse["emission"]
        if lab not in em or "SHIP  tier8000 phi4" not in em:
            return None
        cl0 = em["SHIP  tier8000 phi4"]["claim"]
        cl1 = em[lab]["claim"]
        r0 = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, cl0)) / sD
        r1 = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, cl1)) / sD
        d_raw = r1 - r0
    else:
        g = mse.get(f"geometry:{tag}")
        if not g:
            return None
        bn, an = f"{tag}:BASE", f"{tag}:{key}"
        if bn not in g or an not in g:
            return None
        b, a = g[bn], g[an]
        rA, rB, rC = a["A"] / b["A"], a["B"] / b["B"], a["claim"] / b["claim"]
        r0 = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, a15_claim)) / sD
        r1 = 1 + (S15**2 * A0 * rA - 2 * S15 * X0 + Kc
                  - min(S15**2 * B0 * rB, a15_claim * rC)) / sD
        d_raw = r1 - r0
    sc = (SCORE(base_raw + d_raw) - SCORE(base_raw)) / 6.0
    return dict(d_raw=d_raw, leg=-sc, base_raw=base_raw,
                censored=SCORE(base_raw + d_raw) == 0.0)


def main() -> None:
    pds = json.loads(open(os.path.join(WORK, "report_pds241.json")).read())
    mse = json.loads(open(os.path.join(WORK, "mse241.json")).read())
    six = json.loads(open(os.path.join(WORK, "six241.json")).read())["table"]
    cos = json.loads(open(os.path.join(WORK, "cos241.json")).read())
    anchor = mse["anchor"]

    BASES = [("BEFORE", "MF", "noMassFix", anchor["a15_raw"]),
             ("AFTER",  "NOMF", "addMassFix", anchor["a15_raw"]),
             ("SHIPPED", "SHIPPED", "addMassFix", BOARD["raw_mse"])]

    out, table = {}, []
    for label, pkey, mkey, dekey, cls in COMPONENTS:
        row = dict(component=label, cls=cls)
        for when, tag, massname, base_raw in BASES:
            r = {}
            pk = massname if pkey == "MASSFIX" else pkey
            mk = massname if mkey == "MASSFIX" else mkey
            # `addMassFix` is an ADD-IN, not a leave-one-out: the mass fix is OFF in the
            # shipped stack, so "what having it is worth" is +(arm - base), not -(arm - base).
            sign = +1.0 if pk == "addMassFix" else -1.0
            r["direction"] = "add-in" if sign > 0 else "leave-one-out"
            # ---- pds
            if pk is None:
                r["pds"] = 0.0
                r["pds3"] = 0.0
                r["pds_kolf2"] = 0.0
                r["pds_note"] = "0 by construction (emission-only; board-verified twice)"
            elif pk in pds.get(tag, {}):
                d = pds[tag][pk]
                sgn = sign
                r["pds"] = sgn * d["overall4"]
                r["pds3"] = sgn * d.get("overall3", float("nan"))
                r["pds_kolf2"] = sgn * d.get("kolf2_overall", float("nan"))
                r["pds_signs"] = d["signs"]
                r["pds_nfold"] = d["n"]
            # ---- mse
            m = mse_leg(mse, tag, mk, anchor, base_raw) if mk else None
            if m is not None:
                r["mse"] = m["leg"] * (-sign if sign > 0 else 1.0)
                r["mse_d_raw"] = m["d_raw"]
                r["mse_censored"] = m["censored"]
                r["mse_deflated"] = m["leg"] * R3_DEFLATOR if m["censored"] else m["leg"]
            # ---- DE four
            dl = dekey
            if dl and dl in six:
                w = "BEFORE" if when == "BEFORE" else "AFTER"
                if w in six[dl]:
                    d = six[dl][w]
                    for k in ("reach", "fid", "jac", "jac_np", "nmae", "e",
                              "d_reach_raw", "reach_t"):
                        r[k] = d.get(k)
            # ---- cosine
            ck = massname if pkey == "MASSFIX" else pkey
            if ck and ck in cos.get(tag, {}).get("loo", {}):
                c = cos[tag]["loo"][ck]
                s = -sign
                r["dcos3"] = s * c["dcos3"] * (-1.0)      # positive = the component ADDS cosine
                r["rel3"] = s * c["rel3"] * (-1.0)
                r["de4_fill_lo"] = s * c["de4_fill_k1000"]
                r["de4_fill_hi"] = s * c["de4_fill_k1384"]
            # ---- sum
            if label in STRUCTURAL_ZERO:
                r["structural_zero"] = STRUCTURAL_ZERO[label]
                r["n_de_measured"] = 0
            legs = [r.get(k) for k in ("pds", "mse", "reach", "fid", "jac", "nmae")]
            r["n_priced"] = sum(1 for x in legs if x is not None)
            if label in STRUCTURAL_ZERO:
                r["n_priced"] = 2
            r["sum"] = sum(x for x in legs if x is not None)
            row[when] = r
        table.append(row)
        out[label] = row

    def table_print(when, title):
        print()
        print("=" * 152)
        print(title)
        print("=" * 152)
        print(f"  {'#':>2} {'component':<36}{'pds':>10}{'mse':>10}{'reach':>10}{'fid':>10}"
              f"{'jac':>10}{'nmae':>10}{'SUM':>11}{'n':>3}{'pds 3-can':>11}{'pds kolf2':>11}"
              f"{'Δcos %':>9}{'e':>8}")
        rows = sorted(table, key=lambda r: -r[when]["sum"])
        order = []
        for i, row in enumerate(rows, 1):
            r = row[when]
            cells = []
            for k in ("pds", "mse", "reach", "fid", "jac", "nmae"):
                v = r.get(k)
                cells.append(f"{v:+10.6f}" if v is not None else f"{'—':>10}")
            p3 = r.get("pds3")
            pk = r.get("pds_kolf2")
            rc = r.get("rel3")
            ee = r.get("e")
            print(f"  {i:>2} {row['component']:<36}" + "".join(cells)
                  + f"{r['sum']:>+11.6f}{r['n_priced']:>3}"
                  + (f"{p3:>+11.6f}" if p3 is not None and p3 == p3 else f"{'—':>11}")
                  + (f"{pk:>+11.6f}" if pk is not None and pk == pk else f"{'—':>11}")
                  + (f"{rc:>+9.2f}" if rc is not None and rc == rc else f"{'—':>9}")
                  + (f"{ee:>8.3f}" if ee is not None and ee == ee else f"{'—':>8}"))
            order.append(row["component"])
        return order

    o_before = table_print("BEFORE", "BEFORE — 200's renderer (VCC_MASS_FIX=minnorm), CD4 f=0.10")
    o_after = table_print("AFTER", "AFTER — THE SHIPPED RENDERER (mass fix unset), CD4 f=0.10")
    o_ship = table_print("SHIPPED",
                         "AFTER — shipped renderer AT THE SHIPPED DOSE f=0.34, mse on the "
                         "board's own score_mse 0.17957")

    print()
    print("=" * 100)
    print("RANK CHANGES")
    print("=" * 100)
    print(f"  {'component':<40}{'200 rank':>10}{'shipped-renderer rank':>24}"
          f"{'at f=0.34':>12}")
    for c in o_before:
        print(f"  {c:<40}{o_before.index(c)+1:>10}{o_after.index(c)+1:>24}"
              f"{o_ship.index(c)+1:>12}")

    json.dump(dict(table=out, order=dict(before=o_before, after=o_after, shipped=o_ship),
                   board=BOARD),
              open(os.path.join(WORK, "rank241.json"), "w"), indent=1, default=float)
    print("\nwrote work/rank241.json")


if __name__ == "__main__":
    sys.exit(main())
