"""436 analytic (no emission): the donor-signed MDE add on H1's essential-covered targets.

KA-434g: 434's `ess434.py` numbers on its own cell table (base_s0_bps delivered lfc, six donors, three essential
screens): MDE our sign -0.0075, add signed by the six-donor mean -0.0242, by the essential mean -0.0272 (covered).
Then 436's sign table (`signs436.donor_signs`: four essential screens incl. Nadig HepG2 for coverage; the mean over
ALL covering tables, or the essential ones only; optional >= 2-table agreement) on the same arithmetic, split by
covered / uncovered and by K562-essential-list membership of the target.
Writes work/ana436.json."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import signs436 as S6

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
z = np.load(ROOT / "experiments/434-posterior-median-magnitude/work/allcells434.npz", allow_pickle=True)
L, T, D = z["L"].astype(float), z["T"].astype(float), z["D"].astype(float)
donors = [str(x) for x in z["donors"]]; targets = z["targets"].astype(str); genes = z["genes"].astype(str)
r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool)
ref_lin = np.asarray(r["ref_lin"], float); ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
cnt = ref_lin * np.median(np.asarray(r["libraries"], float)) / 1e6
add = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)[None, :]
Sg = np.sign(L)
tg = [k for k in range(L.shape[0]) if np.isfinite(T[k]).sum() >= 10]


def nmae_k(P, k):
    m = np.isfinite(T[k])
    return float(np.abs(P[k, m] - T[k, m]).sum() / np.abs(T[k, m]).sum())


def nmae(P, ks):
    return float(np.mean([nmae_k(P, k) for k in ks])) if ks else float("nan")


def acc(sig, ks):
    return float(np.mean([np.mean(np.sign(sig[k, np.isfinite(T[k])]) == np.sign(T[k, np.isfinite(T[k])])) for k in ks]))


out = {}
# ---- KA-434g (ess434.py verbatim arithmetic)
ESS = [donors.index(d) for d in S6.ESS3]
cov3 = np.isfinite(D[ESS]).any(2).any(0)
tc3 = [k for k in tg if cov3[k]]
with np.errstate(all="ignore"):
    dm_all = np.nanmean(D, 0); dm_ess = np.nanmean(D[ESS], 0)
base3 = nmae(L, tc3)
ka = {"n_gated": len(tg), "n_cov3": len(tc3), "mde_our": nmae(L + Sg * add, tc3) - base3}
for nm, dm in (("all6", dm_all), ("ess3", dm_ess)):
    sd = np.sign(np.nan_to_num(dm)); sd = np.where(sd == 0, Sg, sd)
    ka[nm] = nmae(L + sd * add, tc3) - base3
print("KA-434g:", json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in ka.items()}))
assert len(tc3) == 20 and abs(ka["mde_our"] - (-0.0075)) < 5e-5 and abs(ka["all6"] - (-0.0242)) < 5e-5 \
    and abs(ka["ess3"] - (-0.0272)) < 5e-5, ka
out["ka434g"] = ka

# ---- 436's sign tables
cov = S6.coverage(targets)
out["h1_coverage"] = {nm: {"all85": int(v.sum()), "gated58": int(v[tg].sum()),
                           "targets": [str(t) for t, c in zip(targets, v) if c]} for nm, v in cov.items()}
print("H1 coverage (85 / 58 gated):", {nm: (int(v.sum()), int(v[tg].sum())) for nm, v in cov.items()})
ess_list = set(S6._lines(S6.ESS_DIRS["k562_essential"] / "perts.txt"))
tc = [k for k in tg if cov["union"][k]]; tu = [k for k in tg if not cov["union"][k]]
tc_e = [k for k in tc if targets[k] in ess_list]; tc_n = [k for k in tc if targets[k] not in ess_list]
out["covered_gated"] = [str(targets[k]) for k in tc]
out["covered_gated_not_on_k562ess_list"] = [str(targets[k]) for k in tc_n]
print(f"gated covered {len(tc)} (on the K562-essential list {len(tc_e)}, not {len(tc_n)}: {[targets[k] for k in tc_n]}); uncovered {len(tu)}")
rows = {}
mde_our = L + Sg * add
for label, mode, ma in (("all", "all", 1), ("ess", "ess", 1), ("allx", "allx", 1), ("all_agree2", "all", 2), ("ess_agree2", "ess", 2)):
    s, info = S6.donor_signs(targets, genes, mode=mode, min_agree=ma)
    sd = np.where(s != 0, s, Sg)
    P = L + sd * add
    nb = {}
    for nm_, ks in (("covered", tc), ("covered_ess_list", tc_e), ("covered_nonlist", tc_n), ("uncovered", tu)):
        if not ks:
            continue
        nb[nm_] = {"n": len(ks), "base": nmae(L, ks), "d_mde_our": nmae(mde_our, ks) - nmae(L, ks),
                   "d_mde_donor": nmae(P, ks) - nmae(L, ks),
                   "increment": nmae(P, ks) - nmae(mde_our, ks),
                   "acc_ours": acc(L, ks), "acc_donor_signed": acc(np.where(s != 0, s, 0) + 1e-12 * Sg, ks)}
    # on the covered targets: share of gated-axis genes carrying a donor sign, and how many disagree with ours
    frac = float(np.mean([(s[k] != 0).mean() for k in tc])); dis = float(np.mean([((s[k] != 0) & (s[k] != Sg[k])).sum() for k in tc]))
    gate_cov = float(np.mean([np.mean(s[k, np.isfinite(T[k])] != 0) for k in tc]))
    per_t = {str(targets[k]): {"d_mde_our": nmae_k(mde_our, k) - nmae_k(L, k), "d_mde_donor": nmae_k(P, k) - nmae_k(L, k)} for k in tc}
    rows[label] = {"info": {k: v for k, v in info.items() if k != "tables_per_covered_target"}, "blocks": nb,
                   "signed_frac_gated_axis": frac, "disagree_genes_per_covered_target": dis,
                   "gate_genes_signed_frac": gate_cov, "per_target": per_t}
    c = nb["covered"]
    print(f"{label:11s} covered {c['n']}: base {c['base']:.5f}  MDE ours {c['d_mde_our']:+.5f}  donor-signed {c['d_mde_donor']:+.5f}  "
          f"increment {c['increment']:+.5f}; signed {frac:.3f} of axis, {gate_cov:.3f} of gate cells; disagree {dis:.0f}/target; "
          f"acc ours {c['acc_ours']:.3f} donor {c['acc_donor_signed']:.3f}")
    for nm_ in ("covered_ess_list", "covered_nonlist", "uncovered"):
        if nm_ in nb:
            b = nb[nm_]
            print(f"             {nm_:17s} n {b['n']:2d}: MDE ours {b['d_mde_our']:+.5f} donor {b['d_mde_donor']:+.5f} incr {b['increment']:+.5f}")
out["rows"] = rows
(HERE / "work").mkdir(exist_ok=True)
(HERE / "work/ana436.json").write_text(json.dumps(out, indent=1, default=float))
