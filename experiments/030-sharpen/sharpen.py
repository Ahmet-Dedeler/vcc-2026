"""030 — is experiment 029 cross-donor consensus, or just magnitude sharpening?

Experiment 029 multiplied each predicted effect by the cross-donor confidence raised to a
power, and macro PDS went 0.8017 -> 0.8405. But the confidence signal is the geometric mean of
the donors' own magnitudes, and the prediction is (roughly) their arithmetic mean — so
`P * c^lambda` is close to `sign(P) * |P|^(1+lambda)`. That is a **power transform of the
prediction's own magnitudes**, and it would explain the whole effect without any consensus
information being used at all.

The two hypotheses make different predictions, so they are separable:

  consensus   the gain needs *cross-donor* information; a self-power transform should do less
  sharpening  the gain is about the shape of the magnitude distribution; the self-power
              transform should reproduce it

This matters for what gets built. If it is sharpening, the mechanism is that `gene_norm` and
whitening between them leave the predicted magnitude distribution **flatter than the truth's**,
and the fix is a one-parameter exponent that could be set from the donors alone.

Also reported: `nmae` at a **fixed** alpha rather than at its own optimum. Experiment 029
reported `nmae*` at the minimising alpha, and at high lambda that minimiser goes to **zero** —
i.e. the proxy is saying "predict nothing". Our actual operating point is alpha ~0.92 of the
true effect magnitude (inverted from v3's board `mse`), and at that alpha the cost is real
rather than hypothetical.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
from proxy6 import pds_cosine, topk_mask, nmae_curve, scale   # noqa: E402
import run as proxy_run                                        # noqa: E402

GAMMA, K_REAL, ALPHA_OP = 0.65, 300, 0.92
ALPHAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.92, 1.2, 1.6, 2.2])
J_OP = int(np.argmin(np.abs(ALPHAS - ALPHA_OP)))

folds, common = proxy_run.build_folds("1:1")
print(f"panel {len(common)}  gamma={GAMMA}  nmae reported at alpha={ALPHA_OP} "
      f"(our operating point) and at its own optimum\n")


def variants(P, parts, others, q):
    A, B = parts[others[0]], parts[others[1]]
    geo = np.sqrt(np.abs(A) * np.abs(B))
    out = {}
    out[f"self q={q:g}"] = np.sign(P) * np.abs(P) ** q
    out[f"geo  q={q:g}"] = P * (geo / (np.median(geo[geo > 0]) or 1.0)) ** (q - 1)
    out[f"mask q={q:g}"] = (P * (geo / (np.median(geo[geo > 0]) or 1.0)) ** (q - 1)
                            * (np.sign(A) == np.sign(B)))
    return out


rows = {}
for q in (1.0, 1.5, 2.0, 2.5, 3.0):
    for h, f in folds.items():
        P0 = f["W"].apply(f["P"], GAMMA)
        T = f["T"]
        R = topk_mask(T, K_REAL)
        Tn = T / (np.median(np.linalg.norm(T, axis=1)) + 1e-12)
        for name, P in variants(P0, f["parts"], f["others"], q).items():
            if q == 1.0 and not name.startswith("self"):
                continue
            Pn = P / (np.median(np.linalg.norm(P, axis=1)) + 1e-12)
            cur = nmae_curve(Pn, Tn, R, ALPHAS)
            rows.setdefault(name, {}).setdefault("pds", []).append(pds_cosine(P, T))
            rows[name].setdefault("nmae_op", []).append(cur[J_OP])
            rows[name].setdefault("nmae_best", []).append(float(cur.min()))

print(f"{'variant':14} {'pds':>8} {'s_pds':>7} {'nmae@.92':>9} {'s_nmae':>7} "
      f"{'nmae_best':>10} {'sum@op':>8}")
out = {}
for name, v in rows.items():
    pds_v = float(np.mean(v["pds"]))
    nm_op = float(np.mean(v["nmae_op"]))
    nm_b = float(np.mean(v["nmae_best"]))
    s_p, s_n = scale("pds", pds_v), scale("nmae", nm_op)
    out[name] = dict(pds=pds_v, nmae_op=nm_op, nmae_best=nm_b, s_pds=s_p, s_nmae=s_n,
                     sum=s_p + s_n)
    print(f"{name:14} {pds_v:8.4f} {s_p:7.3f} {nm_op:9.4f} {s_n:7.3f} {nm_b:10.4f} "
          f"{s_p + s_n:8.4f}")
best = max(out, key=lambda k: out[k]["sum"])
print(f"\nbest s_pds+s_nmae at our operating alpha: {best} -> {out[best]['sum']:.4f}")
Path("experiments/030-sharpen/results.json").write_text(json.dumps(out, indent=2))
