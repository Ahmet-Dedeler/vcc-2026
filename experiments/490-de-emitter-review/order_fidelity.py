#!/usr/bin/env python3
"""490 F-order: does the REALISED p-value order inside the delivered call set follow the documented intent (descending
|H| = |key|, rc_emit.py:941 / rc_emit_k.py:886), and what is reach worth under each order?

Variants (calls, signs and lfc held at the delivered ones; only the ORDER the scorer reads is changed):
  realised   the delivered p / p_adj (known answer: reproduces the harness)
  intent     within the delivered call set, p and p_adj strictly increasing in the |key| rank (the design's order);
             outside it, unchanged
  z_target   within the call set, the order of the designed |z| TARGETS (pass-1 multiset re-keyed by |key|): identical
             to `intent` by construction, so it is not computed separately
Also reports, per target, Spearman(|key| rank, realised p rank) inside the call set.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import score_de  # noqa: E402


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse"
    z = np.load(HERE / "work" / f"recon_{tag}.npz", allow_pickle=False)
    d = {k: z[k] for k in z.files}
    T, G = d["targets"].astype(str), d["genes_gated"].astype(str)
    p, q, lfc = d["p_del"], d["q_del"], d["lfc_del"]
    key = np.abs(d["key"])
    called = q < 0.05
    p2, q2 = p.copy(), q.copy()
    rhos = []
    for i in range(len(T)):
        c = np.flatnonzero(called[i])
        # the design's order: descending |key|, ties by the stable argsort the emitter uses (gene index)
        order = c[np.argsort(-key[i, c], kind="stable")]
        vals = 1e-300 * 10.0 ** (np.arange(order.size) * (250.0 / max(order.size, 1)))   # strictly increasing, < 1e-50
        p2[i, order] = vals
        q2[i, order] = vals
        rk_key = np.empty(c.size)
        rk_key[np.argsort(-key[i, c], kind="stable")] = np.arange(c.size)
        rhos.append(spearmanr(rk_key, np.argsort(np.argsort(p[i, c], kind="stable"), kind="stable")).statistic)
    assert np.array_equal(q2 < 0.05, called)
    raw_r, per_r = score_de.score(T, G, p, q, lfc)
    raw_i, per_i = score_de.score(T, G, p2, q2, lfc)
    print("realised:", {k: round(v, 6) for k, v in raw_r.items()})
    print("intent  :", {k: round(v, 6) for k, v in raw_i.items()})
    print("intent - realised:", {k: round(raw_i[k] - raw_r[k], 6) for k in raw_r})
    print("Spearman(|key| rank, realised p rank) inside the call set: quantiles",
          np.round(np.quantile(rhos, [0, .1, .5, .9, 1]), 3))
    rr = np.array([per_r["reach"][t] for t in T])
    ri = np.array([per_i["reach"][t] for t in T])
    dlt = ri - rr
    print(f"per-target reach delta: mean {dlt.mean():+.5f}, targets up {int((dlt > 1e-12).sum())}, "
          f"down {int((dlt < -1e-12).sum())}")
    rng = np.random.default_rng(490)
    boots = [dlt[rng.integers(0, len(dlt), len(dlt))].mean() for _ in range(4000)]
    print("target bootstrap 95% CI of mean reach delta:", np.round(np.quantile(boots, [.025, .975]), 5))
    (HERE / "work" / f"order_fidelity_{tag}.json").write_text(json.dumps(dict(
        tag=tag, raw_realised=raw_r, raw_intent=raw_i, spearman_q=list(np.quantile(rhos, [0, .1, .5, .9, 1])),
        reach_delta_mean=float(dlt.mean()), reach_delta_ci=list(np.quantile(boots, [.025, .975]))), indent=1))


if __name__ == "__main__":
    main()
