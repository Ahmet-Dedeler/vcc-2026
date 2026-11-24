"""345 — common-program structure. Steps 0-3 + 4a + 5(data prep): saved arrays only.

Step 0: known answers (343: 0.588 delivered / 0.600 designed weighted; 323 NMAE 0.950836).
Step 1: SVD energy of R and D_sup; singular-vector cosines.
Step 2: LOO common vectors from R (mean / sign-vote / |lfc|-weighted mean) + all-down + gene-shuffle null.
Step 3: sign(D + a*s*c) grid with implied + exact NMAE.
Step 4a: design's own LOO common vector, alone + grid.
Writes work/stage345.npz (R, Dsup, Dfull, gates, delivered) for the follow-up script.
"""
from __future__ import annotations
import json
import numpy as np, polars as pl
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
W = REPO / "experiments/323-current-family-validation/work"
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(345)

# ---------- load ----------
s = np.load(W / "stats-cap15.npz")
r = np.load(W / "reference.npz")
e = np.load(W / "effects.npz", mmap_mode="r")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
G = len(genes)
targets = s["targets"].astype(str)
assert (targets == r["targets"].astype(str)).all() and (e["targets"].astype(str) == targets).all()
T = len(targets)
gi = {g: i for i, g in enumerate(genes)}
ref = pl.read_parquet(W / "reference_de.parquet")
lhat_del = np.nan_to_num(np.asarray(s["lfc"], dtype=np.float64))
basal = np.asarray(r["ref_lin"], dtype=float)

# design mapped by name onto universe (343 signs343.py construction)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}
col = np.array([eg.get(g, -1) for g in genes])
miss = col < 0
print("universe genes absent from effects axis:", int(miss.sum()),
      [genes[i] for i in np.where(miss)[0]])
Dfull = np.asarray(e["cap15"], dtype=np.float64)[:, np.where(miss, 0, col)]
Dfull[:, miss] = 0.0

# ---------- gate + known answers ----------
gate_mask = np.zeros((T, G), dtype=bool)
L = np.zeros((T, G))  # reference lfc on gate, else 0  (= R)
elig = []
per_t = {}
for ti, t in enumerate(targets):
    sub = ref.filter(pl.col("target") == t)
    g = sub["feature"].to_numpy().astype(str)
    q = sub["p_adj"].to_numpy().astype(float)
    l = sub["log2_fold_change"].to_numpy().astype(float)
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    n = int(gate.sum())
    if n == 0:
        continue
    idx = np.array([gi[x] for x in g[gate]])
    gate_mask[ti, idx] = True
    L[ti, idx] = l[gate]
    n = int(gate.sum())
    if n >= 10:
        elig.append(ti)
        lh = lhat_del[ti, idx]
        Dd = np.abs(l[gate]).sum()
        per_t[t] = dict(n_gate=n,
                        nmae=float(np.abs(lh - l[gate]).sum() / Dd),
                        m=float(np.abs(lh).sum() / Dd),
                        pw_del=float((np.abs(l[gate]) * (np.sign(lh) == np.sign(l[gate]))).sum() / Dd),
                        pw_des=float((np.abs(l[gate]) * (np.sign(Dfull[ti, idx]) == np.sign(l[gate]))).sum() / Dd))
elig = np.array(elig)
E = len(elig)
print(f"eligible targets: {E} (343: 58)")
pw_del = float(np.mean([per_t[targets[i]]["pw_del"] for i in elig]))
pw_des = float(np.mean([per_t[targets[i]]["pw_des"] for i in elig]))
nmae = float(np.mean([per_t[targets[i]]["nmae"] for i in elig]))
m_glob = float(np.mean([per_t[targets[i]]["m"] for i in elig]))
print(f"KNOWN-ANSWER delivered weighted sign acc: {pw_del:.6f} (record 0.588...)")
print(f"KNOWN-ANSWER designed  weighted sign acc: {pw_des:.6f} (record 0.600...)")
print(f"KNOWN-ANSWER mean NMAE: {nmae:.6f} (323 recorded 0.950836)")
print(f"mean magnitude ratio m: {m_glob:.6f} (343: 0.205)")

R = L  # 85 x G, gate values else 0
Dsup = np.where(gate_mask, Dfull, 0.0)

res = {"known_answer": {"pw_delivered": pw_del, "pw_designed": pw_des,
                        "nmae_recomputed": nmae, "m_mean": m_glob,
                        "n_eligible": E,
                        "nmae_323_recorded": json.load(open(W / "metrics-cap15.json"))["raw"]["nmae"]},
       "missing_design_genes": [genes[i] for i in np.where(miss)[0]]}

# ---------- step 1: structure ----------
def energy_frac(M, ks=(1, 2, 3, 5, 10)):
    U, Sv, Vt = np.linalg.svd(M, full_matrices=False)
    tot = float((Sv ** 2).sum())
    fr = {int(k): float((Sv[:k] ** 2).sum() / tot) for k in ks}
    return fr, U, Sv, Vt, tot

frR, UR, SR, VR, totR = energy_frac(R)
frD, UD, SD, VD, totD = energy_frac(Dsup)
cosL = float(UR[:, 0] @ UD[:, 0])
cosR = float(VR[0] @ VD[0])
# mean-direction (all-ones left vector) energy share
u0 = np.full(T, 1 / np.sqrt(T))
share_mean_R = float(np.linalg.norm(u0 @ R) ** 2 / totR)
share_mean_D = float(np.linalg.norm(u0 @ Dsup) ** 2 / totD)
# pool down-share (155 comparison): unweighted + |lfc|-weighted over pooled eligible gate
pm = np.zeros(G, dtype=bool)  # not needed; pool directly
lp, wp = [], []
for ti in elig:
    m_ = gate_mask[ti]
    lp.append(L[ti, m_])
lp = np.concatenate(lp)
down_unw = float((lp < 0).mean())
down_w = float(np.abs(lp[lp < 0]).sum() / np.abs(lp).sum())
res["step1"] = {"energy_R": frR, "energy_Dsup": frD,
                "F2_R": totR, "F2_Dsup": totD,
                "cos_left1_signed": cosL, "cos_left1_abs": abs(cosL),
                "cos_right1_signed": cosR, "cos_right1_abs": abs(cosR),
                "mean_dir_share_R": share_mean_R, "mean_dir_share_D": share_mean_D,
                "pool_down_unweighted": down_unw, "pool_down_weighted": down_w,
                "pooled_gate_genes": int(lp.size)}
print("energy R:", frR, " Dsup:", frD)
print(f"cos left1 {cosL:+.4f} (|.|={abs(cosL):.4f})  right1 {cosR:+.4f} (|.|={abs(cosR):.4f})")
print(f"mean-dir share R {share_mean_R:.4f} D {share_mean_D:.4f}; down unw {down_unw:.4f} w {down_w:.4f}")

# ---------- step 2: LOO common vectors from R ----------
def wacc(ti, pred_sign):
    m_ = gate_mask[ti]
    l = L[ti, m_]
    w = np.abs(l)
    D_ = w.sum()
    ok = (pred_sign[m_] == np.sign(l))
    return float((w * ok).sum() / D_), float(ok.mean()), float(D_)

Rsgn = np.sign(R)
Rabs = np.abs(R)
out2 = {}
variants = {}
# mean over ALL other 84 (state it)
for name in ("mean", "vote", "wmean"):
    pw_l, pu_l, cov_l = [], [], []
    C = np.zeros_like(R)
    for ti in range(T):
        o = [j for j in range(T) if j != ti]
        if name == "mean":
            c = R[o].mean(axis=0)
        elif name == "vote":
            c = Rsgn[o].sum(axis=0)
        else:
            num = (R[o] * Rabs[o]).sum(axis=0)
            den = Rabs[o].sum(axis=0)
            c = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        C[ti] = c
        if ti in elig:
            a, u, _ = wacc(ti, np.sign(c))
            pw_l.append(a); pu_l.append(u)
            m_ = gate_mask[ti]
            cov_l.append(float(np.abs(L[ti, m_][c[m_] != 0]).sum() / np.abs(L[ti, m_]).sum()))
    variants[name] = C
    out2[name] = {"weighted_mean": float(np.mean(pw_l)), "unweighted_mean": float(np.mean(pu_l)),
                  "per_target_weighted": pw_l, "coverage_mass": float(np.mean(cov_l))}
# constant all-down on this gate
ad_w, ad_u = [], []
for ti in elig:
    a, u, _ = wacc(ti, np.full(G, -1.0))
    ad_w.append(a); ad_u.append(u)
out2["all_down"] = {"weighted_mean": float(np.mean(ad_w)), "unweighted_mean": float(np.mean(ad_u))}
# gene-shuffle null for the mean variant: 200 seeds, permute gene labels of c_{-t}
sh = []
Cmean = variants["mean"]
for sd in range(200):
    rs = np.random.default_rng(10_000 + sd)
    acc = []
    for ti in elig:
        cp = Cmean[ti][rs.permutation(G)]
        a, _, _ = wacc(ti, np.sign(cp))
        acc.append(a)
    sh.append(float(np.mean(acc)))
sh = np.array(sh)
out2["gene_shuffle_null_meanvec"] = {"mean": float(sh.mean()), "sd": float(sh.std()),
                                     "real": out2["mean"]["weighted_mean"],
                                     "z": float((out2["mean"]["weighted_mean"] - sh.mean()) / sh.std())}
res["step2"] = out2
for k, v in out2.items():
    if "weighted_mean" in v:
        print(f"L2 {k}: w {v['weighted_mean']:.4f} u {v['unweighted_mean']:.4f}")
print("shuffle null: mean %.4f sd %.4f z %.2f" % (sh.mean(), sh.std(), out2["gene_shuffle_null_meanvec"]["z"]))

np.savez_compressed(OUT / "stage345.npz", R=R, Dsup=Dsup, Dfull=Dfull,
                    gate=gate_mask, elig=elig, Lhat_del=lhat_del,
                    Cmean=variants["mean"], Cvote=variants["vote"], Cwmean=variants["wmean"],
                    genes=genes, targets=targets, basal=basal)
json.dump(res, open(OUT / "results_part1.json", "w"), indent=1)
print("wrote stage345.npz + results_part1.json")
