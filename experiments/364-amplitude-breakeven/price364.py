#!/usr/bin/env python3
"""364 pricing: overall(p,a) on H1-local anchors + breakeven lines. Reads
breakeven.json, h1-anchors.json, sweep317.json. Writes priced364.json + prints tables."""
import json
import numpy as np

HERE = "~/Code/virtual-cell-challenge-2026/experiments/364-amplitude-breakeven/work"
BASE = "~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work"
S317 = "~/Code/virtual-cell-challenge-2026/experiments/317-operating-point/work/sweep317.json"

tab = json.load(open(f"{HERE}/breakeven.json"))["table"]
anchors = json.load(open(f"{BASE}/h1-anchors.json"))
s317 = json.load(open(S317))

B, Rp = anchors["baseline"], anchors["replicate"]
span = {m: Rp[m] - B[m] for m in B}
print("H1-local spans (r-b):", {m: round(span[m], 6) for m in span})

# current delivered raw on THIS surface (323 cap15)
cur = {"nmae": 0.95083646, "fid": 0.49755452, "reach": 0.25434726,
       "pds": 0.88963585, "mse": 0.97494624, "jac": 0.07475983}
# 317 raw deltas for pds/mse at a (recorded, sweep317.json)
d317 = {a: {m: s317["arms"][f"amp{int(a*100) if a != 1.45 else 145}"]["delta_raw"][m]
            for m in ("pds", "mse")} for a in (1.2, 1.45, 2.0)}
d317[1.0] = {"pds": 0.0, "mse": 0.0}
print("317 deltas used:", json.dumps(d317, indent=1))

P = ["0.6", "0.62", "0.65", "0.7", "0.75"]
A = [1.0, 1.2, 1.45, 2.0]

def scaled(d_raw, m):
    return d_raw / span[m]

# overall(p,a): mean of six scaled deltas vs current
overall = {}
for p in P:
    overall[p] = {}
    for a in A:
        d = {}
        d["nmae"] = scaled(tab[p][f"nmae_a{a}_mean"] - cur["nmae"], "nmae")
        d["fid"] = scaled(tab[p][f"fid_a{a}_mean"] - cur["fid"], "fid")
        d["reach"] = scaled(tab[p][f"reach_a{a}_mean"] - cur["reach"], "reach")
        d["pds"] = scaled(d317[a]["pds"], "pds")
        d["mse"] = scaled(d317[a]["mse"], "mse")
        d["jac"] = 0.0
        overall[p][a] = {"members": d, "overall": float(sum(d.values()) / 6)}
        # also with mse counted FREE (ceiling, 317 §3.2 style)
        overall[p][a]["overall_msefree"] = float((sum(d.values()) - d["mse"]) / 6)

print("\n== overall(p,a) H1-local (mean of 6 scaled deltas vs delivered) ==")
print(f"{'p':>6}" + "".join(f"{'a='+str(a):>14}" for a in A))
for p in P:
    print(f"{p:>6}" + "".join(f"{overall[p][a]['overall']:+.6f}" for a in A))
print("\n== with mse FREE (ceiling) ==")
for p in P:
    print(f"{p:>6}" + "".join(f"{overall[p][a]['overall_msefree']:+.6f}" for a in A))

print("\n== member legs at a=1 (pure sign gain, no amplitude) ==")
for p in P:
    m = overall[p][1.0]["members"]
    print(f"p={p}: " + " ".join(f"{k}={v:+.6f}" for k, v in m.items()))

print("\n== member legs of x1.2 vs x1 at each p (what amplitude buys/pays) ==")
for p in P:
    d1, d2 = overall[p][1.0]["members"], overall[p][1.2]["members"]
    legs = {k: d2[k] - d1[k] for k in d1}
    print(f"p={p}: " + " ".join(f"{k}={v:+.6f}" for k, v in legs.items())
          + f" | net={(sum(legs.values())/6):+.6f}")

# breakeven a*(p): interpolate overall(a) zero crossing per p
print("\n== breakeven a*(p): overall=0 crossing ==")
for p in P:
    xs = [overall[p][a]["overall"] for a in A]
    # find first crossing from below/above
    cross = None
    for i in range(len(A) - 1):
        if xs[i] >= 0 and xs[i + 1] <= 0 or xs[i] <= 0 and xs[i + 1] >= 0:
            f = -xs[i] / (xs[i + 1] - xs[i])
            cross = A[i] + f * (A[i + 1] - A[i])
            break
    print(f"p={p}: overall(a)={[round(x,6) for x in xs]} -> a*={cross}")

# breakeven p*(a): interpolate overall(p) zero crossing per a
print("\n== breakeven p*(a) ==")
pf = [float(p) for p in P]
for a in A:
    xs = [overall[p][a]["overall"] for p in P]
    cross = None
    for i in range(len(P) - 1):
        if (xs[i] <= 0 <= xs[i + 1]) or (xs[i] >= 0 >= xs[i + 1]):
            f = -xs[i] / (xs[i + 1] - xs[i])
            cross = pf[i] + f * (pf[i + 1] - pf[i])
            break
    print(f"a={a}: overall(p)={[round(x,6) for x in xs]} -> p*={cross}")

# pay per +0.01 accuracy: slope of overall(p, a) in p
print("\n== pay per +0.01 pooled accuracy (OLS slope over the 5 p points) ==")
for a in A:
    xs = np.array(pf); ys = np.array([overall[p][a]["overall"] for p in P])
    slope = float(((xs - xs.mean()) * (ys - ys.mean())).sum() / ((xs - xs.mean()) ** 2).sum())
    print(f"a={a}: {slope*0.01:+.6f} overall per +0.01")

# raw nmae table for the record
print("\n== raw nmae(p,a) means ==")
for p in P:
    print(f"p={p} pach={tab[p]['p_ach_mean']:.5f} " +
          " ".join(f"a{a}={tab[p][f'nmae_a{a}_mean']:.5f}" for a in A))
print("\n== raw fid(p) / reach(p) means (a-invariant, verified 0.0e+00) ==")
for p in P:
    print(f"p={p} fid={tab[p]['fid_a1.0_mean']:.6f} reach={tab[p]['reach_a1.0_mean']:.6f}")

json.dump({"overall": {p: {str(a): overall[p][a] for a in A} for p in P},
           "spans": span, "d317": d317, "current": cur},
          open(f"{HERE}/priced364.json", "w"), indent=1)
