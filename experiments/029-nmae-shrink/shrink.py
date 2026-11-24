"""029 — buy `nmae` with per-gene confidence shrinkage, and price it in `pds`.

`nmae` is the largest reachable gap left (us −0.020 scaled, Vivai +0.120, and our raw 1.0123
is *worse than predicting nothing*). It has a structural property nothing else on the board
has:

    nmae_p = sum_{g in DE} |pred_g − true_g| / sum_{g in DE} |true_g|

A gene where we predict **zero** contributes exactly `|true_g|`, i.e. a ratio of 1. A gene
where we predict confidently in the **wrong direction** contributes more than `|true_g|`.
So on genes we cannot call, silence is strictly safer than commitment — and `nmae` is the only
metric that pays for restraint.

That makes per-gene shrinkage toward zero a real `nmae` lever:

    effect_g  ->  effect_g * (conf_g / median conf)^lambda

with `conf` the un-whitened cross-donor geometric agreement (exp 013's signal, and exp 023's
finding that it survives whitening while `|effect|` does not).

The catch is that this is *per-gene reweighting*, which lost for `pds` in experiments 004 and
006 and did not come back under whitening (exp 022). But those runs asked whether it improves
`pds`. This asks a different question: **how much `pds` does a given amount of `nmae` cost?**
Both are measured here from the same predictions — `pds` exactly, `nmae` by the experiment 018
proxy, which is a direct L1 computation and was validated in the right direction by the
count-level run (it predicted whitening would cost `nmae`, and it did).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
from proxy6 import pds_cosine, topk_mask, nmae_curve, scale       # noqa: E402
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
import run as proxy_run                                            # noqa: E402

LAMBDAS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
MODE = __import__("sys").argv[1] if len(__import__("sys").argv) > 1 else "mask"
GAMMA = 0.65
K_REAL = 300
ALPHAS = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0, 1.3, 1.7, 2.2])

folds, common = proxy_run.build_folds("1:1")
print(f"panel {len(common)}  gamma={GAMMA}  K={K_REAL}\n")
print(f"{'lambda':>7} {'pds':>8} {'nmae*':>8} {'a*':>5} {'s_pds':>7} {'s_nmae':>7} "
      f"{'sum':>8}")
out = {}
for lam in LAMBDAS:
    acc = {k: [] for k in ("pds", "nmae", "a")}
    for h, f in folds.items():
        P = f["W"].apply(f["P"], GAMMA)
        if lam > 0:
            c = proxy_run.geo_conf(f["parts"], f["others"])
            if MODE == "nomask":
                # geometric magnitude WITHOUT the sign-agreement gate, to separate "weight by
                # consensus magnitude" from "delete genes the donors disagree about"
                A, B = f["parts"][f["others"][0]], f["parts"][f["others"][1]]
                c = np.sqrt(np.abs(A) * np.abs(B))
            med = np.median(c[c > 0]) or 1.0
            # Floor at a small positive value: experiment 012 established that anything
            # which can zero a whole prediction is a `pds` hazard, because a zero vector has
            # no direction and ties against everything.
            w = np.clip(c / med, 1e-3, None) ** lam
            P = P * w
        T = f["T"]
        R = topk_mask(T, K_REAL)
        acc["pds"].append(pds_cosine(P, T))
        Pn = P / (np.median(np.linalg.norm(P, axis=1)) + 1e-12)
        Tn = T / (np.median(np.linalg.norm(T, axis=1)) + 1e-12)
        cur = nmae_curve(Pn, Tn, R, ALPHAS)
        j = int(np.argmin(cur))
        acc["nmae"].append(cur[j]); acc["a"].append(ALPHAS[j])
    v = {k: float(np.mean(x)) for k, x in acc.items()}
    s_pds, s_nm = scale("pds", v["pds"]), scale("nmae", v["nmae"])
    out[f"{MODE} lambda={lam:g}"] = {**v, "s_pds": s_pds, "s_nmae": s_nm,
                              "sum": s_pds + s_nm}
    print(f"{lam:7.2f} {v['pds']:8.4f} {v['nmae']:8.4f} {v['a']:5.2f} {s_pds:7.3f} "
          f"{s_nm:7.3f} {s_pds + s_nm:8.4f}", flush=True)

best = max(out, key=lambda k: out[k]["sum"])
print(f"\nbest s_pds+s_nmae: {best} -> {out[best]['sum']:.4f}")
print("(the two metrics weigh equally in Overall, so their scaled sum is the right objective)")
Path(f"experiments/029-nmae-shrink/results_{MODE}.json").write_text(json.dumps(out, indent=2))
