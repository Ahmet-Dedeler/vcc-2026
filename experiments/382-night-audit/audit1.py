"""382 item 1: recompute raw fid + nmae for 363 arms from saved stats + 323 reference.

Spec: docs/reference/vcc2026-metrics-BRIEF.md sections 3/4/7.
fid: F_p = k/max(n_pred,n_real); own target excluded; adjudicable = reference lfc
     non-NaN and nonzero (inf carries direction); ours significant = own q<0.05;
     directionless own call (lfc==0) counts as miss; 0/0 dropped from the mean.
nmae: gate = reference q<0.05 & finite ref lfc, target removed; gate>=10 genes;
     denom zero/nonfinite -> omitted; nonfinite own lfc read as 0.
pds/mse are NOT recomputable from these files (no design effects, no counts).
"""
import json
import numpy as np
import pandas as pd

R323 = "experiments/323-current-family-validation/work"
W363 = "experiments/363-block-down-sixmember/work"
ARMS = ["cap15", "block_down", "block_random"]
TOL = 1e-9

ref = np.load(f"{R323}/reference.npz")
genes = ref["genes"]
keep = ref["keep"]
uni = genes[keep]                      # (10779,) universe symbols
rt, rq, rl = ref["targets"], ref["q"], ref["lfc"]
assert rq.shape == rl.shape == (len(rt), len(uni)), (rq.shape, rl.shape, len(uni))

# own-target universe index per target (None if absent from universe)
tidx = []
for t in rt:
    w = np.where(uni == t)[0]
    tidx.append(int(w[0]) if len(w) else None)

# parquet cross-check of the gate definition
parq = pd.read_parquet(f"{R323}/reference_de.parquet")
print("parquet columns:", list(parq.columns), "rows:", len(parq))

out = {"arms": {}, "notes": []}
out["notes"].append(f"ref lfc: nan={int(np.isnan(rl).sum())}, "
                    f"inf={int(np.isinf(rl).sum())}, zero={int((rl == 0).sum())}")

for arm in ARMS:
    st = np.load(f"{W363}/stats-{arm}.npz")
    assert list(st["targets"]) == list(rt), "target order mismatch"
    assert st["q"].shape == (85, 10779), st["q"].shape
    oq, ol = st["q"], st["lfc"]
    out["notes"].append(
        f"{arm}: own lfc nan={int(np.isnan(ol).sum())} inf={int(np.isinf(ol).sum())} "
        f"zero={int((ol == 0).sum())}; own q<0.05 frac={float((oq < 0.05).mean()):.4f}")

    Fp, Np, n_pred_l, n_real_l, gate_l = [], [], [], [], []
    for i in range(len(rt)):
        adjud = ~np.isnan(rl[i]) & (rl[i] != 0)          # inf counts (spec)
        refsig = (rq[i] < 0.05) & adjud
        oursig = oq[i] < 0.05
        excl = np.ones(len(uni), bool)
        if tidx[i] is not None:
            excl[tidx[i]] = False
        n_real = int((refsig & excl).sum())
        pred = oursig & adjud & excl
        n_pred = int(pred.sum())
        k = int(((ol[i] * rl[i] > 0) & pred).sum())      # zero own-lfc -> miss
        den = max(n_pred, n_real)
        Fp.append(k / den if den else np.nan)
        n_pred_l.append(n_pred)
        n_real_l.append(n_real)
        # nmae gate: finite ref lfc only
        gate = refsig & np.isfinite(rl[i]) & excl
        gate_l.append(int(gate.sum()))
        if gate.sum() >= 10:
            denom = float(np.abs(rl[i][gate]).sum())
            if denom != 0 and np.isfinite(denom):
                own = np.where(np.isfinite(ol[i][gate]), ol[i][gate], 0.0)
                Np.append(float(np.abs(own - rl[i][gate]).sum()) / denom)

    fid = float(np.nanmean(Fp))
    nmae = float(np.nanmean(Np))
    rec = json.load(open(f"{W363}/metrics-{arm}.json"))["raw"]
    dfid = abs(fid - rec["fid"])
    dnmae = abs(nmae - rec["nmae"])
    out["arms"][arm] = {
        "fid_recomputed": fid, "fid_recorded": rec["fid"], "fid_absdiff": dfid,
        "fid_match": bool(dfid <= TOL),
        "nmae_recomputed": nmae, "nmae_recorded": rec["nmae"], "nmae_absdiff": dnmae,
        "nmae_match": bool(dnmae <= TOL),
        "n_fid_defined": int(np.sum(~np.isnan(Fp))),
        "n_nmae_eligible": len(Np),
        "n_pred_mean": float(np.mean(n_pred_l)),
        "n_real_mean": float(np.mean(n_real_l)),
        "gate_mean": float(np.mean(gate_l)),
    }

# parquet gate cross-check: 58 targets with >=10 gate genes?
g = parq[(parq["p_adj"] < 0.05) & np.isfinite(parq["log2_fold_change"])
         & (parq["log2_fold_change"] != 0)
         & (parq["feature"] != parq["target"])]
cnt = g.groupby("target").size()
out["parquet_gate"] = {"n_targets_ge10": int((cnt >= 10).sum()),
                       "mean_gate": float(cnt[cnt >= 10].mean())}

print(json.dumps(out, indent=1))
json.dump(out, open("experiments/382-night-audit/audit1.json", "w"), indent=1)
