"""434g: an OCTOBER-only DE lever. The essential screens (K562-essential, RPE1, Nadig Jurkat) cover 0/300 A/B/C
targets but ~30 D/E/F targets each (DEF-RUNBOOK 5.-2). On H1's targets they cover, how much does taking the MDE add's
SIGN (and optionally the arithmetic-mean lfc's sign) from the covering donors' mean buy over our own sign?
Analytic nmae on h1-test's gates, reported on the covered targets and scaled to a 300-target panel."""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "work/allcells434.npz", allow_pickle=True)
L, T, D = z["L"].astype(float), z["T"].astype(float), z["D"].astype(float)
donors = [str(x) for x in z["donors"]]
r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool); ref_lin = np.asarray(r["ref_lin"], float)
ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
cnt = ref_lin * np.median(np.asarray(r["libraries"], float)) / 1e6
ESS = [donors.index(d) for d in ("k562_essential", "rpe1", "nadig_jurkat")]
cov = np.isfinite(D[ESS]).any(2).any(0)
tg = [k for k in range(L.shape[0]) if np.isfinite(T[k]).sum() >= 10]
tc = [k for k in tg if cov[k]]
add = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)[None, :]
S = np.sign(L)
with np.errstate(all="ignore"):
    dm_all = np.nanmean(D, 0)                  # all six donors that cover (p, g)
    dm_ess = np.nanmean(D[ESS], 0)             # essential screens only
def nmae(P, ks):
    return float(np.mean([np.abs(P[k, np.isfinite(T[k])] - T[k, np.isfinite(T[k])]).sum() / np.abs(T[k, np.isfinite(T[k])]).sum() for k in ks]))
def acc(sig, ks):
    a = [np.mean(np.sign(sig[k, np.isfinite(T[k])]) == np.sign(T[k, np.isfinite(T[k])])) for k in ks]
    return float(np.mean(a))
print(f"gated targets {len(tg)}, covered by an essential screen {len(tc)}")
base = nmae(L, tc); mde = L + S * add
print(f"covered targets: base nmae {base:.5f}; our gate sign acc {acc(L, tc):.3f}; donor-mean(all) {acc(np.nan_to_num(dm_all), tc):.3f}; ess-mean {acc(np.nan_to_num(dm_ess), tc):.3f}")
rows = [("MDE, our sign", mde)]
for nm, dm in (("all-donor mean", dm_all), ("essential mean", dm_ess)):
    sd = np.sign(np.nan_to_num(dm)); sd = np.where(sd == 0, S, sd)
    rows.append((f"MDE add signed by {nm}", L + sd * add))
    for w in (1.0, 2.0, 3.0):
        # blend: our lfc + w * donor mean (log2 units), then MDE add on the blend's sign
        bl = L + w * np.nan_to_num(dm) / np.log(2)
        rows.append((f"blend ours + {w:g}x {nm} (log2), + MDE on its sign", bl + np.sign(bl) * add))
for nm, P in rows:
    d = nmae(P, tc) - base
    print(f"  {nm:58s} covered d {d:+.5f}  (x {len(tc)}/{len(tg)} over H1: {d * len(tc) / len(tg):+.5f}; per 300-target panel with 45 covered: {d * 45 / 300:+.5f})")
