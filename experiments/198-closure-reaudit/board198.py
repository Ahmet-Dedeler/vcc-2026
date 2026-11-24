#!/usr/bin/env python3
"""198 — board-wide adjudication of the AgingLab GeroAI_v9 -> v11 reach move.

Three questions, all answerable from the 613-entry snapshot on disk:
  Q1  exact six-member decomposition of the +0.0354 overall in board units
  Q2  does the 186 amplitude law explain it? how many residual sigma is it?
  Q3  does `fid` (sign accuracy on the call set) move with it? if not, the sign
      INFORMATION did not change and the move is ordering/DELIVERY, not key quality
"""
from __future__ import annotations
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(exist_ok=True)

rows = json.loads((ROOT / "experiments/186-intel-0903/board/leaderboard_2026-09-03.json"
                   ).read_text())["entries"]

RAW = dict(pds="pds_cosine", mse="expr_mse_unbiased_capped_norm",
           nmae="de_wilcoxon_lfc_nmae",
           fid="de_wilcoxon_direction_fidelity_yield_raw",
           reach="de_wilcoxon_direction_reach_raw", jac="de_wilcoxon_sig_jaccard")
# docs/BOARD-LEDGER.md: scaled = slope*(raw - b); overall = mean of 6
SLOPE = dict(pds=2.2138, mse=-1.051, nmae=-1.6532, fid=3.3645, reach=1.1292, jac=2.7263)
DOVERALL = {k: v / 6 for k, v in SLOPE.items()}

# GeroAI pair, from experiments/186-intel-0903/board/trajectories.md lines 17-18
V9 = dict(pds=0.8273, mse=1.822, nmae=0.9216, fid=0.5050, reach=0.1672, jac=0.0284)
V11 = dict(pds=0.8196, mse=0.959, nmae=0.8923, fid=0.5135, reach=0.2674, jac=0.0287)
OVERALL_V9, OVERALL_V11 = 0.1545, 0.1899

print("=== Q1  GeroAI_v9 -> v11, six-member decomposition in OVERALL units ===")
print(f"{'member':7s} {'v9':>9s} {'v11':>9s} {'d raw':>10s} {'dOverall/dRaw':>14s} "
      f"{'d overall':>11s}")
tot = 0.0
for k in ("pds", "mse", "nmae", "fid", "reach", "jac"):
    d = V11[k] - V9[k]
    # mse is clamped: raw above the ~0.989 baseline scales to 0
    if k == "mse":
        s9 = min(0.0, SLOPE[k] * (V9[k] - 0.989)) if V9[k] > 0.989 else SLOPE[k] * (V9[k] - 0.989)
        s9 = max(0.0, SLOPE[k] * (V9[k] - 0.989))
        s11 = max(0.0, SLOPE[k] * (V11[k] - 0.989))
        do = (s11 - s9) / 6
    else:
        do = d * DOVERALL[k]
    tot += do
    print(f"{k:7s} {V9[k]:9.4f} {V11[k]:9.4f} {d:+10.4f} {DOVERALL[k]:+14.4f} {do:+11.5f}")
print(f"{'SUM':7s} {OVERALL_V9:9.4f} {OVERALL_V11:9.4f} "
      f"{'':10s} {'':14s} {tot:+11.5f}")
print(f"reported d overall = {OVERALL_V11 - OVERALL_V9:+.5f}   "
      f"reconstruction error {tot - (OVERALL_V11 - OVERALL_V9):+.5f}")
print(f"\n`reach` alone = {(V11['reach'] - V9['reach']) * DOVERALL['reach']:+.5f} overall "
      f"= {100 * (V11['reach'] - V9['reach']) * DOVERALL['reach'] / (OVERALL_V11 - OVERALL_V9):.0f}% "
      f"of the whole move.")

# ---------------------------------------------------------------- Q2
print("\n=== Q2  does the 186 amplitude law explain it? ===")
ok = [r for r in rows if all(r.get(v) is not None for v in RAW.values())
      and r[RAW["mse"]] > 0]
X = {k: np.array([r[v] for r in ok]) for k, v in RAW.items()}
print(f"n entries with all six raw members: {len(ok)}")

hi = X["pds"] > 0.80
lm, rc, pd_ = np.log(X["mse"][hi]), X["reach"][hi], X["pds"][hi]
A = np.column_stack([np.ones(hi.sum()), lm, pd_])
beta, *_ = np.linalg.lstsq(A, rc, rcond=None)
resid = rc - A @ beta
print(f"pds>0.80 cohort n={hi.sum()}: reach ~ 1 + log(mse) + pds")
print(f"  d reach / d log mse = {beta[1]:+.5f}   (186 quotes +0.0127)")
print(f"  residual sd = {resid.std(ddof=3):.5f}")
dlm = np.log(V11["mse"] / V9["mse"])
pred = beta[1] * dlm + beta[2] * (V11["pds"] - V9["pds"])
print(f"  GeroAI d log mse = {dlm:+.4f}, d pds = {V11['pds'] - V9['pds']:+.4f}")
print(f"  law PREDICTS d reach = {pred:+.5f}; OBSERVED {V11['reach'] - V9['reach']:+.5f}")
print(f"  unexplained = {V11['reach'] - V9['reach'] - pred:+.5f} "
      f"= {(V11['reach'] - V9['reach'] - pred) / resid.std(ddof=3):.1f} residual sigma")

# ---------------------------------------------------------------- Q3
print("\n=== Q3  did the SIGN INFORMATION move, or only the ORDERING? ===")
print("`fid` = fraction of YOUR called DE genes with the correct sign, x yield cap.")
print("`reach` = depth your ranking of the REFERENCE's DE genes stays >=90% sign-pure.")
print("A change in the underlying sign INFORMATION moves both. A change in the")
print("ORDERING/DELIVERY moves `reach` and leaves `fid` and `jac` ~flat.\n")
print(f"  GeroAI d fid  = {V11['fid'] - V9['fid']:+.4f}  "
      f"({100 * (V11['fid'] - V9['fid']) / V9['fid']:+.2f}%)")
print(f"  GeroAI d jac  = {V11['jac'] - V9['jac']:+.4f}  "
      f"({100 * (V11['jac'] - V9['jac']) / V9['jac']:+.2f}%)")
print(f"  GeroAI d reach= {V11['reach'] - V9['reach']:+.4f}  "
      f"({100 * (V11['reach'] - V9['reach']) / V9['reach']:+.2f}%)")
print(f"  GeroAI d pds  = {V11['pds'] - V9['pds']:+.4f}  "
      f"({100 * (V11['pds'] - V9['pds']) / V9['pds']:+.2f}%)")

# board-wide: how strongly does reach track fid at all?
r_rf = np.corrcoef(X["reach"][hi], X["fid"][hi])[0, 1]
print(f"\n  board pds>0.80 cohort: corr(raw reach, raw fid) = {r_rf:+.4f}")
print(f"  a 60% reach jump with fid +1.7% is NOT a sign-information change.")

# where do 0.1672 and 0.2674 sit in the board's reach distribution?
q = np.percentile(X["reach"], [5, 25, 50, 75, 90, 95, 99])
print(f"\n  board raw reach percentiles (n={len(ok)}): "
      f"5% {q[0]:.4f} 25% {q[1]:.4f} 50% {q[2]:.4f} 75% {q[3]:.4f} "
      f"90% {q[4]:.4f} 95% {q[5]:.4f} 99% {q[6]:.4f}")
for nm, v in (("GeroAI_v9", V9["reach"]), ("GeroAI_v11", V11["reach"]),
              ("us a8 (high-water)", 0.1890), ("us abc-knee (current)", 0.177303)):
    print(f"    {nm:22s} {v:.4f}  -> board percentile "
          f"{100 * (X['reach'] < v).mean():5.1f}")

# ---------------------------------------------------------------- Q4
print("\n=== Q4  what 172 ACTUALLY measured, in board units ===")
H1T = dict(constant=0.18021, production=0.18283, M1=0.22448, M2=0.23163, H=0.23533,
           oracle=0.64553)
print("h1-test `reach`, 172 §4 (all with the SAME predictions and SAME signs):")
for k, v in H1T.items():
    print(f"  {k:12s} {v:.5f}   d vs production {v - H1T['production']:+.5f}")
print("\nBoard shrink factor for a h1-test `reach` delta: 0.390 (a8/a15 pair) to")
print("0.4277 (164 §6)  -- docs/BOARD-LEDGER.md §4.")
for k in ("M1", "M2", "H"):
    d = H1T[k] - H1T["production"]
    print(f"  172's `{k}` key: h1-test {d:+.5f} -> board {d * 0.390:+.5f} to "
          f"{d * 0.4277:+.5f} raw reach = {d * 0.390 * DOVERALL['reach']:+.5f} to "
          f"{d * 0.4277 * DOVERALL['reach']:+.5f} overall")
print(f"\nAgingLab realised {V11['reach'] - V9['reach']:+.4f} raw board reach.")
best = (H1T["H"] - H1T["production"]) * 0.4277
print(f"Our best measured ordering counterfactual is {best:+.4f} raw board reach.")
print(f"AgingLab's move is {(V11['reach'] - V9['reach']) / best:.1f}x LARGER than the "
      f"entire span\nfrom our delivered order to `|H|`, the best key 172 could find.")

json.dump(dict(decomposition_total=tot, reach_share=(V11["reach"] - V9["reach"])
               * DOVERALL["reach"] / (OVERALL_V11 - OVERALL_V9),
               amp_slope=float(beta[1]), amp_resid_sd=float(resid.std(ddof=3)),
               amp_pred=float(pred), amp_unexplained_sigma=float(
                   (V11["reach"] - V9["reach"] - pred) / resid.std(ddof=3)),
               corr_reach_fid=float(r_rf)),
          open(OUT / "board198.json", "w"), indent=1)
