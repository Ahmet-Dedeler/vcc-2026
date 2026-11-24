"""382 item 5: K562 all-row mean vector read by sign on the H1 gate.

Gate rebuilt per 345 known-answer recipe: reference_de.parquet, p_adj<0.05 &
finite nonzero lfc & feature != target; 58 targets with >=10 gate genes.
Vector: mean over ALL 9787 rows of data/donors/k562_gwps/delta.npy, mapped by
gene name. Weighted sign accuracy per target = sum|ref|*1[sign match]/sum|ref|;
macro = mean over 58; pooled = global sums. Sign convention checked first via
own-target negativity (must be ~99% negative).
"""
import json
import numpy as np
import pandas as pd

R323 = "experiments/323-current-family-validation/work"
parq = pd.read_parquet(f"{R323}/reference_de.parquet")
ref = np.load(f"{R323}/reference.npz")
uni = ref["genes"][ref["keep"]]
U = {g: i for i, g in enumerate(uni)}

kg = [l.strip() for l in open("data/donors/k562_gwps/genes.txt")]
kp = [l.strip() for l in open("data/donors/k562_gwps/perts.txt")]
kd = np.load("data/donors/k562_gwps/delta.npy", mmap_mode="r")
assert kd.shape == (len(kp), len(kg)), (kd.shape, len(kp), len(kg))
K = {g: j for j, g in enumerate(kg)}

# sign-convention check: own-target delta < 0 fraction over covered targets
kp_arr = np.array(kp)
own_neg, own_cov = [], 0
for t in set(kp_arr):
    if t in K:
        rows = np.where(kp_arr == t)[0]
        own_cov += len(rows)
        own_neg.append(float((kd[rows][:, K[t]] < 0).mean()))
print(f"own-target: {own_cov} rows, mean neg frac={np.mean(own_neg):.5f}")

# gate per target (eligible set)
g = parq[(parq["p_adj"] < 0.05) & np.isfinite(parq["log2_fold_change"])
         & (parq["log2_fold_change"] != 0) & (parq["feature"] != parq["target"])]
cnt = g.groupby("target").size()
elig = sorted(cnt[cnt >= 10].index)
print("eligible targets:", len(elig))
assert len(elig) == 58, len(elig)

# K562 all-mean vector on the K562 axis, then mapped to universe
meanv = np.asarray(kd.mean(axis=0), dtype=np.float64)   # (8248,)
m_universe = np.full(len(uni), np.nan)
n_map = 0
for gu, i in U.items():
    if gu in K:
        m_universe[i] = meanv[K[gu]]
        n_map += 1
print(f"universe genes mapped to K562 bank: {n_map}/{len(uni)}")

wacc, uacc, pool_num, pool_den, pool_n = [], [], 0.0, 0.0, 0
cover_fracs = []
for t in elig:
    rows = g[g["target"] == t]
    feats = rows["feature"].to_numpy()
    rl = rows["log2_fold_change"].to_numpy(float)
    ui = np.array([U.get(f, -1) for f in feats])
    ok = (ui >= 0) & np.isfinite(m_universe[ui]) & (m_universe[ui] != 0)
    cover_fracs.append(float(ok.mean()))
    rl, mv = rl[ok], m_universe[ui[ok]]
    match = (np.sign(mv) == np.sign(rl))
    w = np.abs(rl)
    wacc.append(float((w * match).sum() / w.sum()))
    uacc.append(float(match.mean()))
    pool_num += float((w * match).sum())
    pool_den += float(w.sum())
    pool_n += int(ok.sum())

out = {
    "own_target_rows": own_cov, "own_target_neg_frac": float(np.mean(own_neg)),
    "n_eligible": len(elig),
    "universe_mapped": n_map,
    "gate_gene_coverage_in_k562_mean": float(np.mean(cover_fracs)),
    "macro_weighted_acc": float(np.mean(wacc)),
    "macro_unweighted_acc": float(np.mean(uacc)),
    "pooled_weighted_acc": pool_num / pool_den,
    "pooled_pairs": pool_n,
    "pi_macro": 0.527, "pi_pooled": 0.547, "report345": 0.275,
}
print(json.dumps(out, indent=1))
json.dump(out, open("experiments/382-night-audit/audit5.json", "w"), indent=1)
