"""185 Part E — the six-member ledger in board units, three readings per arm.

  pds   identity (85/85 digests == a15's)                          -> 0.000000
  mse   claim at the cap 0.24653 >= knee 0.23268 on the challenge  -> 0.000000 (measured)
  nmae  pinned by fixed T: every rc arm reads 0.96670 vs a15's twin mean; converted anyway
  reach h1-test twin-mean delta x 0.4277 (164/171 shrink) x 1/span   [NOT board-validated for
        emission-law arms: this instrument inverts a21's sign, Part D]
  fid   (i) h1-test delta, (ii) the board's five-point n_pred series (n_pred unchanged -> 0)
  jac   h1-test delta x 1/span x 0.44 (181 §1's board/h1-test jac ratio)
Rule 13 composite for an emission-law arm: 1.011*dpds/6 + dmse/6 - 0.0146/6.
"""
import json, numpy as np
J = json.load(open("work/six185.json"))
M = J["members"]
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)
SHRINK, JACR = 0.4277, 0.44
base = {k: 0.5 * (M["a15"][k] + M["a15_b"][k]) for k in ("reach", "fid", "prec", "jac", "nmae")}
print(f"a15 twin mean: reach {base['reach']:.5f} fid {base['fid']:.5f} jac {base['jac']:.5f} nmae {base['nmae']:.5f}")
print(f"{'arm':<12}{'d reach raw':>12}{'reach ovl':>11}{'d fid raw':>11}{'fid ovl':>9}{'d jac raw':>11}{'jac ovl':>9}"
      f"{'d nmae raw':>11}{'nmae ovl':>10}{'SUM(h1)':>9}{'rule13':>9}")
for a in ("rc_hist_a", "rc_hybrid", "rc_hist_s", "rc_lfc", "rc_k2000", "rc_k4000", "rc_perm", "a21_like", "a22_like"):
    if a not in M:
        continue
    r = M[a]["reach"] if a != "rc_hist_a" else 0.5 * (M["rc_hist_a"]["reach"] + M["rc_hist_a_b"]["reach"])
    dr = r - base["reach"]; ro = dr * SHRINK * INV["reach"] / 6
    df = M[a]["fid"] - base["fid"]; fo = df * INV["fid"] / 6
    dj = M[a]["jac"] - base["jac"]; jo = dj * INV["jac"] * JACR / 6
    dn = -(M[a]["nmae"] - base["nmae"]); no = dn * INV["nmae"] / 6
    print(f"{a:<12}{dr:>+12.5f}{ro:>+11.5f}{df:>+11.5f}{fo:>+9.5f}{dj:>+11.5f}{jo:>+9.5f}"
          f"{dn:>+11.5f}{no:>+10.5f}{ro+fo+jo+no:>+9.5f}{-0.0146/6:>+9.5f}")
print("\nreading 1 (rule 13, emission-law class): predicted = 0 + 0 - 0.00243 = -0.00243 < +0.0019 -> VETO")
print("reading 2 (h1-test transferred with the a8/a15 shrink; NOT validated for this class): SUM(h1) column")
