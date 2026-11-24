"""357 known-answer: independently recompute 343's headline decomposition numbers
from 323's saved files (read-only). Expected (343 RESULTS.md): mean NMAE 0.950836
over eligible targets; m=0.205; pw(delivered)=0.588; called_share=0.346;
zero_share=0.0000; cf all-signs-right=0.803; cf best-scalar=0.928."""
import json
import numpy as np
import polars as pl
from pathlib import Path

REPO = Path("~/Code/virtual-cell-challenge-2026")
W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/357-jurkat-shrna-truth-proxy/work"

s = np.load(W / "stats-cap15.npz", allow_pickle=False)
r = np.load(W / "reference.npz", allow_pickle=False)
genes_all = r["genes"].astype(str)
keep = r["keep"].astype(bool)
genes = genes_all[keep]
targets = s["targets"].astype(str)
assert (targets == r["targets"].astype(str)).all()
print("targets:", len(targets), "universe genes:", len(genes), "keep sum:", int(keep.sum()))
assert len(targets) == 85 and int(keep.sum()) == 10779, "axis mismatch vs record (85 / 10,779)"

ref = pl.read_parquet(W / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}
lhat_all = np.asarray(s["lfc"], dtype=np.float64)
qhat_all = np.asarray(s["q"], dtype=np.float64)

rows = []
for ti, t in enumerate(targets):
    sub = ref.filter(pl.col("target") == t)
    g = sub["feature"].to_numpy().astype(str)
    q = sub["p_adj"].to_numpy().astype(float)
    l = sub["log2_fold_change"].to_numpy().astype(float)
    gate = (q < 0.05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10:
        continue
    idx = np.array([gi[x] for x in g[gate]])
    l = l[gate]
    lh = np.where(np.isfinite(lhat_all[ti, idx]), lhat_all[ti, idx], 0.0)
    qh = qhat_all[ti, idx]
    D = np.abs(l).sum()
    rows.append(dict(
        target=t,
        nmae=float(np.abs(lh - l).sum() / D),
        m=float(np.abs(lh).sum() / D),
        pw=float((np.abs(l) * (np.sign(lh) == np.sign(l))).sum() / D),
        zero_share=float(np.abs(l[np.abs(lh) < 1e-6]).sum() / D),
        called_share=float(np.abs(l[qh < 0.05]).sum() / D),
        cf_sign=float(np.abs(np.abs(lh) - np.abs(l)).sum() / D),
    ))

T = pl.DataFrame(rows)
means = {k: float(T[k].mean()) for k in ("nmae", "m", "pw", "zero_share",
                                         "called_share", "cf_sign")}
means["n_targets"] = T.height
print(json.dumps(means, indent=1))

expected = {"nmae": 0.950836, "m": 0.205, "pw": 0.588, "zero_share": 0.0000,
            "called_share": 0.346, "cf_sign": 0.803}
print("\ncheck vs 343 record:")
ok = True
for k, v in expected.items():
    d = abs(means[k] - v)
    tol = 0.002 if k != "zero_share" else 1e-4  # record prints 0.0000 (4dp)
    flag = "OK " if d <= tol else "FAIL"
    if d > tol:
        ok = False
    print(f"  {flag} {k}: mine {means[k]:.6f} vs {v} (|d|={d:.6f})")
print("KNOWN-ANSWER:", "REPRODUCED" if ok else "MISMATCH")
(OUT / "known_answer_343.json").write_text(json.dumps(
    {"mine": means, "expected": expected, "reproduced": ok}, indent=1))
