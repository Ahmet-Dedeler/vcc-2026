"""178 — read the four DE members off the h1-test tables and convert them to board units.

`six173.py`, unchanged in every convention, pointed at 178's tables. Rule 9 is in force and is
handled exactly as 173 handled it: `fid` is read as PRECISION (`mean_p k/n_pred`, no `n_real`),
`jac` is rebuilt on the board's own `n_real`/`n_pred` from the overlap RATIO, and `reach` is
reported at native BH with its re-seeded twin beside it and used for nothing on its own.

165's finding applies to every row in this table: **h1-test reproduces board `reach` for
emission-law changes but NOT for profile changes**, and every arm here is a profile change. So
`reach` is printed with its twin spread and the reading leans on `fid` and `nmae`, which 173
resolved to twin spreads of 5e-5 and 4e-4.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import ANCHOR, Surface  # noqa: E402

ALPHA = 0.05
ORDER = ["KA_a15_raw", "rec_base_mf", "rb_a15", "rb_a15_b", "rb_g045", "rb_g045_b",
         "rb_g030", "rb_flip", "rb_flip_b", "rb_w211", "rb_w211_b", "rb_g045_flip",
         "rb_g045_w211", "rb_all3_g045", "rb_all3_g045_b", "rb_all3_g030"]
TWIN = {"rb_a15": "rb_a15_b", "rb_g045": "rb_g045_b", "rb_flip": "rb_flip_b",
        "rb_w211": "rb_w211_b", "rb_all3_g045": "rb_all3_g045_b"}
PARTS = {"rb_g045_flip": ("rb_g045", "rb_flip"),
         "rb_g045_w211": ("rb_g045", "rb_w211"),
         "rb_all3_g045": ("rb_g045", "rb_flip", "rb_w211"),
         "rb_all3_g030": ("rb_g030", "rb_flip", "rb_w211")}
# 157 section 1, both at f = 1 on this surface, against the board's own pair.
BOARD_REACH = dict(a8=0.1890, a15=0.1811)
H1_REACH = dict(a8=0.2029, a15=0.1828)
NMAE_B, NMAE_R = 1.0013, 0.400
NR, NP = 240.0, 4424.0                 # the board's own n_real / n_pred (173's convention)


def precision(s):
    out = []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        if sel.size == 0:
            continue
        out.append(float((np.sign(s.pred_lfc[p, sel]) ==
                          np.sign(s.real_lfc[p, sel])).sum()) / sel.size)
    return float(np.mean(out))


def overlap_ratio(s):
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        out.append(int((s.real_sig[p] & s.pred_sig[p]).sum()) / nr)
    return float(np.mean(out))


def main() -> None:
    paths = {os.path.basename(p)[len("de178."):-4]: p
             for p in sorted(glob.glob(os.path.join(HERE, "de178.*.npz")))}
    arms = [a for a in ORDER if a in paths] + [a for a in paths if a not in ORDER]
    if not arms:
        sys.exit("no de178.*.npz -- run realize178.py first")

    rec = {}
    print(f"{'arm':<17}{'n_pred':>8}{'reach':>9}{'precision':>11}{'fid raw':>9}"
          f"{'jac raw':>9}{'recall':>8}{'nmae':>9}{'pds':>9}{'sign':>8}")
    for a in arms:
        s = Surface(paths[a])
        nr, npd = s.counts()
        row = dict(n_pred=float(np.median(npd)), n_real=float(np.median(nr)),
                   precision=precision(s), fid=s.fid(), jac=s.jac(), nmae=s.nmae(),
                   pds=s.pds(), sign=s.sign_accuracy(), recall=overlap_ratio(s),
                   reach=s.reach())
        rec[a] = row
        print(f"{a:<17}{row['n_pred']:>8.0f}{row['reach']:>9.5f}{row['precision']:>11.5f}"
              f"{row['fid']:>9.5f}{row['jac']:>9.5f}{row['recall']:>8.5f}"
              f"{row['nmae']:>9.5f}{row['pds']:>9.5f}{row['sign']:>8.5f}", flush=True)

    # ---------------------------------------------------------------- known-answer gates
    print("\n--- known-answer gates ---")
    ok = True
    if "KA_a15_raw" in rec:
        g = rec["KA_a15_raw"]
        for label, got, want, tol in (("median n_pred", g["n_pred"], 6034, 0.5),
                                      ("reach f=1", g["reach"], 0.18283, 5e-4)):
            good = abs(got - want) <= tol
            ok &= good
            print(f"  164/157 a15 {label:<16} got {got:.5f}  expect {want:.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    ref173 = os.path.join(REPO, "experiments/173-sign-repair")
    if os.path.exists(os.path.join(ref173, "de173.base_mf.npz")) and "rec_base_mf" in rec:
        s173 = Surface(os.path.join(ref173, "de173.base_mf.npz"))
        for m, f in (("precision", precision), ("nmae", lambda x: x.nmae()),
                     ("reach", lambda x: x.reach())):
            got, want = rec["rec_base_mf"][m], f(s173)
            good = abs(got - want) < (5e-4 if m != "reach" else 5e-3)
            ok &= good
            print(f"  173 base_mf {m:<16} got {got:.5f}  expect {want:.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    if not ok:
        print("  *** A GATE FAILED -- everything below is void ***")

    # ---------------------------------------------------------------- twin spreads
    print("\n--- re-seeded twin spreads (the resolution floor of every delta below) ---")
    print(f"{'pair':<20}{'d precision':>13}{'d nmae':>10}{'d reach':>10}{'d recall':>10}")
    for a, b in TWIN.items():
        if a in rec and b in rec:
            print(f"{a + ' / ' + b:<20}"
                  f"{rec[a]['precision'] - rec[b]['precision']:>+13.5f}"
                  f"{rec[a]['nmae'] - rec[b]['nmae']:>+10.5f}"
                  f"{rec[a]['reach'] - rec[b]['reach']:>+10.5f}"
                  f"{rec[a]['recall'] - rec[b]['recall']:>+10.5f}")

    # ---------------------------------------------------------------- board units
    k_reach = ((BOARD_REACH["a8"] - BOARD_REACH["a15"]) /
               (H1_REACH["a8"] - H1_REACH["a15"]))
    span_r = ANCHOR["reach"][1] - ANCHOR["reach"][0]
    span_f = ANCHOR["fid"][1] - ANCHOR["fid"][0]
    span_j = ANCHOR["jac"][1] - ANCHOR["jac"][0]
    span_n = NMAE_B - NMAE_R
    base = "rb_a15"
    b0 = rec[base]
    print(f"\nboard/h1-test `reach` transfer factor from the a8<->a15 pair: {k_reach:.4f}")
    print(f"DELTAS vs `{base}` (the REBUILT a15 config), in OVERALL board units")
    print(f"{'arm':<17}{'d reach raw':>12}{'d reach ovr':>12}{'d prec':>9}{'d fid ovr':>10}"
          f"{'d recall':>10}{'d jac ovr':>10}{'d nmae':>9}{'d nmae ovr':>11}{'SUM DE':>10}")
    ovr = {}
    for a in arms:
        r = rec[a]
        dr = r["reach"] - b0["reach"]
        dr_ovr = dr * k_reach / span_r / 6
        df_ovr = (r["precision"] - b0["precision"]) / span_f / 6
        j0 = b0["recall"] * NR / (NP + NR - b0["recall"] * NR)
        j1 = r["recall"] * NR / (NP + NR - r["recall"] * NR)
        dj_ovr = (j1 - j0) / span_j / 6
        dn = r["nmae"] - b0["nmae"]
        dn_ovr = -dn / span_n / 6
        ovr[a] = dict(reach=dr_ovr, fid=df_ovr, jac=dj_ovr, nmae=dn_ovr,
                      de_sum=dr_ovr + df_ovr + dj_ovr + dn_ovr,
                      de_sum_no_reach=df_ovr + dj_ovr + dn_ovr)
        print(f"{a:<17}{dr:>12.5f}{dr_ovr:>+12.6f}"
              f"{r['precision'] - b0['precision']:>+9.5f}{df_ovr:>+10.6f}"
              f"{r['recall'] - b0['recall']:>+10.5f}{dj_ovr:>+10.6f}"
              f"{dn:>+9.5f}{dn_ovr:>+11.6f}{ovr[a]['de_sum']:>+10.6f}")

    # ---------------------------------------------------------------- composition
    print(f"\n=== composition on the DE members (deltas vs {base}) ===")
    print(f"{'composite':<16}{'channel':<10}{'composite':>12}{'sum parts':>12}{'ratio':>8}"
          f"  reading")
    for comp, members in PARTS.items():
        if comp not in ovr or any(m not in ovr for m in members):
            continue
        for ch in ("fid", "nmae", "jac", "de_sum_no_reach", "de_sum"):
            c = ovr[comp][ch]
            s = sum(ovr[m][ch] for m in members)
            ratio = c / s if abs(s) > 1e-12 else float("nan")
            read = ("super-additive" if ratio > 1.05 else
                    "sub-additive" if ratio < 0.95 else "additive")
            print(f"{comp:<16}{ch:<10}{c:>+12.6f}{s:>+12.6f}{ratio:>8.3f}  {read}")

    with open(os.path.join(WORK, "six178.json"), "w") as fh:
        json.dump({"raw": rec, "overall": ovr}, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'six178.json')}")


if __name__ == "__main__":
    sys.exit(main())
