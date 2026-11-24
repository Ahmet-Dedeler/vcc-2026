"""214 — the composed six-member price and rule 13's slot gate, marginal to `cd4f034`.

Rule 13 prices only two channels and compares against the realised RMS of the four it does
not price:

    predicted = PLOLO d_pds/6 x 1.011  +  096 d_mse/6      (amplitude/spread/phi arms)
                                       -  0.0146/6         (emission-LAW, profile-CHANGING)
    SHIP only if predicted > +0.0019 (shape-only)  or  > +0.0089 (amplitude-changing)

ARM CLASSES, measured not asserted:
  A  drop PerturbMatch -- a DONOR-BLEND change, i.e. profile class.  169 falsified the
     `reach` toll for that class at +0.0006, so no toll.  But `head214.py` measures the
     DELIVERED amplitude moving +3.10% on the median row norm and +4.28% on `A`, four times
     the 1.0% bound 210's mandate used, so the arm is AMPLITUDE-CHANGING and its floor is
     +0.0089.  Per 211's standing rule alpha is NOT re-solved: amplitude is closed, and
     re-solving it would make this a two-variable arm.
  B  Hamilton rounding -- emission-LAW and profile-CHANGING: `chal214.py` measures the
     per-gene panel-total digest differing on 900/900 challenge blocks, max |dT| = 1, every
     move UP.  It carries the -0.0146 toll and claims no exemption.

    uv run python price214.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
import recon196                                                          # noqa: E402
from ham196 import hamilton                                              # noqa: E402

TS = 50_000.0
TOLL = -0.0146 / 6          # rule 13's emission-law, profile-changing toll
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
D_COS_RATE = 0.3655         # 212's measured d(overall mse*)/d(cos) on PLOLO


def hamilton_amplitude():
    """Does arm B move the delivered amplitude?  Measured on cd4f034's own 900 blocks."""
    z = np.load(WORK / "totals-cd4f034.npz", allow_pickle=True)
    blocks = np.asarray(z["blocks"]).astype(str)
    out = []
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        sel = np.flatnonzero(np.array([b.split("|")[0] for b in blocks]) == ctx)
        b = D["b"]
        af, ah, mf, mh = [], [], [], []
        step = max(1, sel.size // 24)
        for j, k in enumerate(sel[::step]):
            R = recon196.Recon(z["totals"][k], D["r"], D["gate"], n=400,
                               seed=1000 + int(j * step))
            for T, A, Mm in ((R.totals(1.0), af, mf), (hamilton(R.x), ah, mh)):
                row = np.log1p(TS * T / T.sum()) - b
                A.append(float((row ** 2).sum()))
                Mm.append(float(np.linalg.norm(row)))
        out.append(dict(ctx=ctx, A_floor=float(np.mean(af)), A_ham=float(np.mean(ah)),
                        med_floor=float(np.median(mf)), med_ham=float(np.median(mh))))
    rA = float(np.mean([o["A_ham"] for o in out]) / np.mean([o["A_floor"] for o in out]))
    rM = float(np.mean([o["med_ham"] for o in out])
               / np.mean([o["med_floor"] for o in out]))
    return rA, rM, out


def main() -> None:
    rep = json.load(open(WORK / "report214.json"))
    mse = json.load(open(WORK / "mse214.json"))
    six = json.load(open(WORK / "six214.json"))
    chal = json.load(open(WORK / "chal214.json"))
    head = json.load(open(WORK / "head214.json"))

    # ---------------------------------------------------------------- the two amplitudes
    rA_ham, rM_ham, ham_rows = hamilton_amplitude()
    aall = head["amplitude"]["all genes"]
    rA_A = float(np.mean([x["A_A"] for x in aall]) / np.mean([x["A_inc"] for x in aall]))
    rM_A = float(np.mean([x["med_A"] for x in aall])
                 / np.mean([x["med_inc"] for x in aall]))
    print("=== ARM CLASS, MEASURED (rule 13's table + 211's standing rule) ===")
    print(f"  {'arm':<26}{'digest identical':>18}{'delivered A ratio':>20}"
          f"{'median|row| ratio':>20}{'class':>26}")
    print(f"  {'A  drop PerturbMatch':<26}{'n/a (not emission)':>18}{rA_A:>20.5f}"
          f"{rM_A:>20.5f}{'profile, AMPLITUDE-CHANGING':>26}")
    ident = f"{chal['class_test']['identical']}/{chal['class_test']['blocks']}"
    print(f"  {'B  Hamilton rounding':<26}{ident:>18}{rA_ham:>20.5f}{rM_ham:>20.5f}"
          f"{'emission-law, profile-CHANGING':>26}")
    cls = dict(A=dict(toll=0.0, gate=GATE_AMP, rA=rA_A, rM=rM_A),
               B=dict(toll=TOLL, gate=GATE_SHAPE if abs(rA_ham - 1) < 0.01 else GATE_AMP,
                      rA=rA_ham, rM=rM_ham))
    cls["AB"] = dict(toll=TOLL, gate=GATE_AMP, rA=rA_A * rA_ham, rM=rM_A * rM_ham)

    # ---------------------------------------------------------------- the six members
    P = rep["comp"]
    d_pds = {"A": P["A  drop PerturbMatch"]["overall"],
             "B": P["B  Hamilton rounding"]["overall"],
             "AB": P["A∘B both"]["overall"]}
    dcos = {"A": P["A  drop PerturbMatch"]["dcos_pct"] / 100.0
            * P["incumbent cd4f034"]["cos"],
            "B": P["B  Hamilton rounding"]["dcos_pct"] / 100.0
            * P["incumbent cd4f034"]["cos"],
            "AB": P["A∘B both"]["dcos_pct"] / 100.0 * P["incumbent cd4f034"]["cos"]}
    gp = mse["geometry_pairs"]
    d_mse_amp = {"A": gp["ARM A: drop PerturbMatch, f 0.34"]["d_overall"],
                 "B": mse["armB_d_overall"]}
    d_mse_amp["AB"] = d_mse_amp["A"] + d_mse_amp["B"]
    d_mse_dir = {k: D_COS_RATE * v for k, v in dcos.items()}

    # ---- arm B on the ARTIFACT family: the surface with a bit-identical KA row -----------
    INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)
    af, ah = six["artifact"]["art_floor"], six["artifact"]["art_ham"]
    B_art = dict(reach=(ah["reach"] - af["reach"]) * 1.290 * INV["reach"] / 6,
                 fid=(ah["fid"] - af["fid"]) * INV["fid"] / 6,
                 jac=0.0,
                 nmae=-(ah["nmae"] - af["nmae"]) * INV["nmae"] / 6)
    B_art["sum"] = sum(B_art.values())
    print(f"\n  arm B on the ARTIFACT family (a15's own delivered h1-test profile and 169's "
          f"recovered key --\n  the surface whose `floor` row is BIT-IDENTICAL to "
          f"`rc_callset`, i.e. 208's own instrument):")
    print(f"    reach {B_art['reach']:+.6f}  fid {B_art['fid']:+.6f}  "
          f"jac {B_art['jac']:+.6f}  nmae {B_art['nmae']:+.6f}  SUM {B_art['sum']:+.6f}")

    DEFLATOR = json.load(open(WORK / "disp214.json"))["h1_over_challenge"]
    print(f"\n  h1-test overstates a PerturbMatch change by x{DEFLATOR:.3f} "
          f"(`disp214.py`: 2 donors vs 3, so the K562 slot carries 1/2 instead of 1/3).\n"
          f"  Every PerturbMatch DE number below is given raw and deflated by that factor.")

    print("\n=== THE SIX MEMBERS, marginal to `cd4f034` (overall board units) ===")
    pools = ("native", "wald89", "matched")
    KEY = {"A": "reb_pm00_floor", "B": "reb_pm50_ham", "AB": "reb_pm00_ham"}
    rows = {}
    for pool in pools:
        print(f"\n--- DE pool `{pool}` ---")
        print(f"  {'member':<34}{'A':>13}{'B':>13}{'A∘B':>13}")
        de = {k: six["board"][f"{pool}|{KEY[k]}"] for k in ("A", "B", "AB")}
        lines = [
            ("`pds`   PLOLO 4-fold x 1.011", {k: d_pds[k] for k in KEY}),
            ("`mse`   096 claim (A/B ratios)", {k: d_mse_amp[k] for k in KEY}),
            ("`mse`   cosine channel (212)", {k: d_mse_dir[k] for k in KEY}),
            ("`reach` h1-test x 1.290", {k: de[k]["reach"] for k in KEY}),
            ("`fid`   h1-test precision", {k: de[k]["fid"] for k in KEY}),
            ("`jac`   S*, NP pinned 2,599", {k: de[k]["jac"] for k in KEY}),
            ("`nmae`  h1-test raw / span", {k: de[k]["nmae"] for k in KEY}),
        ]
        for lab, v in lines:
            print(f"  {lab:<34}" + "".join(f"{v[k]:>+13.6f}" for k in ("A", "B", "AB")))
        tot = {k: sum(v[k] for _, v in lines) for k in KEY}
        print(f"  {'SUM (six members)':<34}" +
              "".join(f"{tot[k]:>+13.6f}" for k in ("A", "B", "AB")))
        print(f"  {'SUM, reach at 0.4277':<34}" +
              "".join(f"{tot[k] - de[k]['reach'] + de[k]['reach_alt']:>+13.6f}"
                      for k in ("A", "B", "AB")))
        # deflate only the PerturbMatch-attributable part of the DE channel
        de_sum = {k: de[k]["sum"] for k in KEY}
        defl = {"A": de_sum["A"] / DEFLATOR, "B": de_sum["B"],
                "AB": de_sum["B"] + (de_sum["AB"] - de_sum["B"]) / DEFLATOR}
        tot_d = {k: tot[k] - de_sum[k] + defl[k] for k in KEY}
        print(f"  {'SUM, DE deflated x1/%.3f' % DEFLATOR:<34}" +
              "".join(f"{tot_d[k]:>+13.6f}" for k in ("A", "B", "AB")))
        print(f"  {'SUM, arm B from the ARTIFACT row':<34}" +
              "".join(f"{(tot[k] - de_sum[k] + B_art['sum']) if k == 'B' else float('nan'):>+13.6f}"
                      for k in ("A", "B", "AB")))
        print(f"  composition on all six: A + B = {tot['A']+tot['B']:+.6f}   measured "
              f"A∘B = {tot['AB']:+.6f}   difference {tot['AB']-tot['A']-tot['B']:+.6f}")
        rows[pool] = dict(total=tot, de=de, deflated=tot_d,
                          alt={k: tot[k] - de[k]["reach"] + de[k]["reach_alt"] for k in KEY})

    # ---------------------------------------------------------------- rule 13's gate
    print("\n=== RULE 13's COMPOSED SLOT GATE ===")
    print(f"  {'arm':<22}{'priced pds':>12}{'priced mse':>12}{'toll':>11}"
          f"{'composite':>12}{'floor':>9}{'verdict':>10}{'six-member (wald89)':>21}")
    verdict = {}
    for k, lab in (("A", "A  drop PerturbMatch"), ("B", "B  Hamilton rounding"),
                   ("AB", "A∘B both")):
        priced = d_pds[k] + d_mse_amp[k]
        comp = priced + cls[k]["toll"]
        g = cls[k]["gate"]
        v = "SHIP" if comp > g else "VETO"
        verdict[k] = dict(priced_pds=d_pds[k], priced_mse=d_mse_amp[k], toll=cls[k]["toll"],
                          composite=comp, gate=g, verdict=v,
                          six_native=rows["native"]["total"][k],
                          six_wald89=rows["wald89"]["total"][k],
                          six_matched=rows["matched"]["total"][k],
                          six_wald89_deflated=rows["wald89"]["deflated"][k],
                          six_native_deflated=rows["native"]["deflated"][k])
        print(f"  {lab:<22}{d_pds[k]:>+12.6f}{d_mse_amp[k]:>+12.6f}{cls[k]['toll']:>+11.6f}"
              f"{comp:>+12.6f}{g:>9.4f}{v:>10}{rows['wald89']['total'][k]:>+21.6f}")
    print("\n  (rule 13's floor IS the realised RMS of the four members it does not price, "
          "so the\n   six-member column is the best estimate and the composite is the gate; "
          "they are not added.)")

    json.dump(dict(classes=cls, pds=d_pds, dcos=dcos, mse_amp=d_mse_amp,
                   mse_dir=d_mse_dir, pools=rows, verdict=verdict,
                   hamilton_amplitude=ham_rows),
              open(WORK / "price214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'price214.json'}")


if __name__ == "__main__":
    main()
