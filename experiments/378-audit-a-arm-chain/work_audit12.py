"""378 items 1+2: artifact integrity + equality. Own code, no audited-script reads."""
import hashlib, json
import numpy as np
import pandas as pd

A375 = "experiments/375-ctx-a-arm/work/ctx_a_blend_rows_b035.npz"
A372 = "experiments/372-jurkat-calibration-artifact/work/ctx_a_calibrated.npz"
EXP_SHA = "500c24ab1f7c0cc6f81e7a56f1ff361415ac845f76efc17aa2c7e9854627a8bd"

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

z5 = np.load(A375, allow_pickle=True)
print("375 keys:", list(z5.files))
for k in z5.files:
    a = z5[k]
    print(f"  {k}: shape={a.shape} dtype={a.dtype}", ("val=%r" % a.item()) if a.size == 1 else "")

t = np.asarray(z5["targets"]).astype(str)
g = np.asarray(z5["genes"]).astype(str)
m = np.asarray(z5["matrix"])
beta = float(np.asarray(z5["beta"]))
mode = str(np.asarray(z5["mode"]).item())
ctx = str(np.asarray(z5["context"]).item())
nr = str(np.asarray(z5["norm_rule"]).item())
print(f"beta={beta} mode={mode} ctx={ctx} norm_rule={nr}")
print(f"matrix finite={np.isfinite(m).all()} shape={m.shape} dtype={m.dtype}")

off = pd.read_csv("data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
off_set = set(off)
panel = pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].astype(str).tolist()
panel_set = set(panel)
print(f"official genes={len(off)} panel targets={len(panel)}")

print("genes n=", len(g), "unique=", len(set(g.tolist())), "all_on_axis=", set(g.tolist()) <= off_set)
print("genes ∩ panel =", sorted(set(g.tolist()) & panel_set), "count=", len(set(g.tolist()) & panel_set))
print("targets n=", len(t), "unique=", len(set(t.tolist())), "all_in_panel=", set(t.tolist()) <= panel_set)
print("targets missing_from_panel =", sorted(set(t.tolist()) - panel_set))
print("panel missing_from_targets =", len(panel_set - set(t.tolist())), sorted(panel_set - set(t.tolist()))[:10])
print("targets duplicates:", len(t) - len(set(t.tolist())))
print("genes duplicates:", len(g) - len(set(g.tolist())))
print("sha375 =", sha(A375))
print("sha match:", sha(A375) == EXP_SHA)

z2 = np.load(A372, allow_pickle=True)
print("372 keys:", list(z2.files))
pr = np.asarray(z2["prediction_raw"])
mx = np.asarray(z2["matrix"])
print(f"372 beta={float(np.asarray(z2['beta']))} mode={str(np.asarray(z2['mode']).item())}")
d_raw = np.max(np.abs(m - pr))
d_mx = np.max(np.abs(m - mx))
print(f"max|375.matrix-372.prediction_raw| = {d_raw:.3e}")
print(f"max|375.matrix-372.matrix|        = {d_mx:.3e}")
print(f"equal_raw_exact = {d_raw == 0.0}")
print(f"targets identical 375vs372: {np.array_equal(np.asarray(z5['targets']).astype(str), np.asarray(z2['targets']).astype(str))}")
print(f"genes identical 375vs372: {np.array_equal(np.asarray(z5['genes']).astype(str), np.asarray(z2['genes']).astype(str))}")
# also check panel-order claim: is targets in panel order?
order_idx = [panel.index(x) for x in t]
print("targets in panel order:", order_idx == sorted(order_idx))
