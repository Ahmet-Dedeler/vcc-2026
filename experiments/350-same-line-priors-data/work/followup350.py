"""350 follow-up: depth check, H1-truth common vector, CPM-gated cosines.

No X re-reads except obs columns (codes + UMI_count). Everything else is
arithmetic on work/*.npz + reference.npz.
"""
import h5py
import json
import numpy as np

WORK = "experiments/350-same-line-priors-data/work"
out = {}

# 1. library-size (depth) controls vs perturbed cells per bank
for name, path in [("jurkat", "data/external/nadig2025/jurkat_raw_singlecell.h5ad"),
                   ("hepg2", "data/external/nadig2025/hepg2_raw_singlecell.h5ad")]:
    f = h5py.File(path, "r")
    cats = np.array([x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
                     for x in f["obs"]["__categories"]["gene"][:]])
    codes = np.asarray(f["obs"]["gene"][:]).astype(np.int64)
    umi = np.asarray(f["obs"]["UMI_count"][:], dtype=np.float64)
    f.close()
    ctrl = cats[codes] == "non-targeting"
    out[name + "_lib_ctrl_mean"] = float(umi[ctrl].mean())
    out[name + "_lib_ctrl_median"] = float(np.median(umi[ctrl]))
    out[name + "_lib_pert_mean"] = float(umi[~ctrl].mean())
    out[name + "_lib_pert_median"] = float(np.median(umi[~ctrl]))
    out[name + "_lib_ratio"] = float(umi[ctrl].mean() / umi[~ctrl].mean())
    print(name, "ctrl mean/median", round(umi[ctrl].mean()), round(np.median(umi[ctrl])),
          "pert mean/median", round(umi[~ctrl].mean()), round(np.median(umi[~ctrl])))

# 2. H1 TRUTH common vector in log2-CPM-ratio space from reference.npz sums
r = np.load("experiments/323-current-family-validation/work/reference.npz")
sums, counts = r["sums"].astype(np.float64), r["counts"].astype(np.float64)
perts = np.array([str(x) for x in r["perts"]])
genes18080 = np.array([str(x) for x in r["genes"]])
ci = list(perts).index("non-targeting") if "non-targeting" in perts else len(perts) - 1
print("control row:", perts[ci], "counts:", counts[ci])
cpm_ctrl = 1e6 * sums[ci] / sums[ci].sum()
pm = np.delete(np.arange(len(perts)), ci)
cpm_p = 1e6 * sums[pm] / sums[pm].sum(axis=1, keepdims=True)
lfc_truth = np.log2((cpm_p + 1e-9) / (cpm_ctrl + 1e-9))  # (85, 18080)
truth_common = lfc_truth.mean(axis=0)
print("truth common: norm", round(float(np.linalg.norm(truth_common)), 4),
      "frac down", round(float((truth_common < 0).mean()), 4),
      "n perts", len(pm))

axis = [l.strip() for l in open("data/challenge/gene_names.csv")][1:]
aidx = {g: i for i, g in enumerate(axis)}
rows = [(aidx[g], j) for j, g in enumerate(genes18080) if g in aidx]
truth_ax = np.zeros(len(axis))
truth_ax[[x[0] for x in rows]] = truth_common[[x[1] for x in rows]]
tm = np.zeros(len(axis), dtype=bool)
tm[[x[0] for x in rows]] = True
print("truth mapped to axis:", tm.sum(), "of", len(axis))

j = np.load(WORK + "/jurkat_common.npz")
h = np.load(WORK + "/hepg2_common.npz")
ca_j, ca_h = j["vector"], h["vector"]
map_j, map_h = j["mapped"], h["mapped"]


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def signagree(a, b):
    m = (a != 0) & (b != 0)
    return float((np.sign(a[m]) == np.sign(b[m])).mean()), int(m.sum())


e = np.load("experiments/323-current-family-validation/work/effects.npz")
h1_common = e["cap15"].mean(axis=0)
for nm, a, b, m in [("truth-jurkat", truth_ax, ca_j, tm & map_j),
                    ("truth-hepg2", truth_ax, ca_h, tm & map_h),
                    ("truth-design", truth_ax, h1_common, tm)]:
    c, (s, n) = cos(a[m], b[m]), signagree(a[m], b[m])
    print(f"cos {nm}: {c:.4f}  sign {s:.4f} (n={n})")
    out["cos_" + nm] = c
    out["sign_" + nm] = [s, n]

# 3. control-CPM-gated cosines: where do the extreme common tails live?
# recompute bank control CPM from saved lfc? No — need control sums. Use propensity
# files + report instead: fraction of top-50 |common| genes with bank control CPM<5.
# Bank control CPM is recoverable: rerun is expensive; approximate via context-A-free
# route is wrong. Instead gate on H1-truth control CPM and on mapped magnitude:
# report cosine after dropping the top-1% |common| genes on each side.
for nm, a, b, m in [("jurkat-hepg2", ca_j, ca_h, map_j & map_h),
                    ("truth-jurkat", truth_ax, ca_j, tm & map_j)]:
    aa, bb = a[m], b[m]
    cut = int(0.01 * len(aa))
    drop = np.argsort(np.abs(aa))[-cut:]
    keep = np.ones(len(aa), dtype=bool)
    keep[drop] = False
    drop2 = np.argsort(np.abs(bb))[-cut:]
    keep2 = np.ones(len(bb), dtype=bool)
    keep2[drop2] = False
    print(f"{nm}: full cos {cos(aa,bb):.4f} | drop-top1%|a| {cos(aa[keep],bb[keep]):.4f} "
          f"| drop-top1%|b| {cos(aa[keep2],bb[keep2]):.4f} | drop-both {cos(aa[keep&keep2],bb[keep&keep2]):.4f}")
    out[nm + "_drop_top1pct"] = [cos(aa[keep], bb[keep]), cos(aa[keep & keep2], bb[keep & keep2])]

# 4. how much of each common vector's ENERGY sits in the top-50 |genes|?
for nm, v, m in [("jurkat", ca_j, map_j), ("hepg2", ca_h, map_h), ("truth", truth_ax, tm)]:
    vv = v[m]
    e50 = (np.sort(vv ** 2)[-50:].sum() / (vv ** 2).sum())
    print(f"{nm}: top-50|.| energy share {float(e50):.4f}")
    out[nm + "_top50_energy"] = float(e50)

# 5. per-target ceiling: how well does one common vector align with individual
# H1 truth effects, vs our design's per-target alignment?
shr = tm & map_j
Tp = lfc_truth[:, [x[1] for x in rows]]  # (85, 18533-full) -> restrict next
full_rows = np.array([x[0] for x in rows])
col_idx = [x[1] for x in rows]
Tfull = np.zeros((lfc_truth.shape[0], len(axis)))
Tfull[:, full_rows] = lfc_truth[:, col_idx]
Tshr = Tfull[:, shr]
Dshr = e["cap15"][:, shr]
Jshr = ca_j[shr]
cj = np.array([float(Tshr[i] @ Jshr / (np.linalg.norm(Tshr[i]) * np.linalg.norm(Jshr)))
               for i in range(Tshr.shape[0])])
cd = np.array([float(Tshr[i] @ Dshr[i] / (np.linalg.norm(Tshr[i]) * np.linalg.norm(Dshr[i])))
               for i in range(Tshr.shape[0])])
print(f"per-target cos(truth, JurkatCommon): mean {cj.mean():.4f} sd {cj.std():.4f} "
      f"frac>0 {(cj > 0).mean():.3f}")
print(f"per-target cos(truth, design):      mean {cd.mean():.4f} sd {cd.std():.4f} "
      f"frac>0 {(cd > 0).mean():.3f}")
out["pertarget_cos_truth_jurkatcommon"] = [float(cj.mean()), float(cj.std())]
out["pertarget_cos_truth_design"] = [float(cd.mean()), float(cd.std())]

json.dump(out, open(WORK + "/followup.json", "w"), indent=2)
print("DONE")
