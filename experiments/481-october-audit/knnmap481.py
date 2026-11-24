"""481 F3, distributional: predict the builder's strength DISTRIBUTION (not the conditional mean) on a panel without a
build, by sampling each target's strength from its k nearest reference targets in whitened-proxy space (whiten481's
base-only proxy, log1p). References: A/B/C (453 dump) and efr (465 dump). Validation: A/B/C-only references must
reproduce efr's classes and vice versa. Reports shares of strength 0, 1, <= 1 and the q33 cut, over 2,000 draws.
usage: uv run python knnmap481.py -> work/knnmap481.json (~2.5 GB)"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import whiten481 as W
HERE, ROOT = W.HERE, W.ROOT
gate = W.gate_D()
abc_panel = set(pd.read_csv(HERE / "work/panels/abc.csv")["target_gene"].astype(str))
banks = [W.L.bank(n) for n in W.L.MAIN]
full = W.L.E064.operator_columns(banks, W.L.OFFICIAL, abc_panel)
evals, evecs = W.load_op(W.OPS["abc_op"])
prox, panels = {}, {}
for name in ("abc", "efr", "h1_2025"):
    panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
    pre = W.L.Donors438(panel, W.L.MAIN).pre(W.L.MAIN, normrel=True)
    pc = np.asarray([W.L.OFF_POS[g] for g in panel if g in W.L.OFF_POS]); pre[:, pc] = 0.0
    X = W.whiten(pre, full, evals, evecs); X[:, pc] = 0.0
    prox[name] = np.log1p((np.abs(np.where(gate[None, :], X, 0.0)) > W.TAU).sum(1).astype(float)); panels[name] = panel
del evecs
known = {"abc": W.dump_strength(ROOT / "experiments/453-abc-dress-rehearsal/work/run/dump_h1x2_cd4x1.npz", panels["abc"], gate),
         "efr": W.dump_strength(ROOT / "experiments/465-both-switches-efr/work/run/dump_h1x2_cd4x1.npz", panels["efr"], gate)}
rng = np.random.default_rng(481)
K, D = 15, 2000
def predict(x, refs):
    rx = np.concatenate([prox[r] for r in refs]); ry = np.concatenate([known[r] for r in refs])
    nn = np.argsort(np.abs(x[:, None] - rx[None, :]) + 1e-9 * rng.random((x.size, rx.size)), axis=1)[:, :K]
    stats = []
    for _ in range(D):
        s = ry[nn[np.arange(x.size), rng.integers(0, K, x.size)]]
        stats.append(((s == 0).mean(), (s == 1).mean(), (s <= 1).mean(), np.quantile(s, 0.33)))
    a = np.asarray(stats)
    q = lambda v: [round(float(np.quantile(v, p)), 3) for p in (0.025, 0.5, 0.975)]
    return {"eq0": q(a[:, 0]), "eq1": q(a[:, 1]), "le1": q(a[:, 2]), "p_q33_cut_le1": round(float((a[:, 3] <= 1).mean()), 3)}
out = {}
for n in ("abc", "efr"):
    s = known[n]
    out[f"{n}_builder"] = {"eq0": round(float((s == 0).mean()), 3), "eq1": round(float((s == 1).mean()), 3),
                           "le1": round(float((s <= 1).mean()), 3), "q33": float(np.quantile(s, 0.33))}
out["validate_abc_refs_predict_efr"] = predict(prox["efr"], ["abc"])
out["validate_efr_refs_predict_abc"] = predict(prox["abc"], ["efr"])
ne = np.array([t not in W.ESS for t in panels["h1_2025"]])
out["h1_2025_nonessential"] = predict(prox["h1_2025"][ne], ["abc", "efr"])
out["h1_2025_all"] = predict(prox["h1_2025"], ["abc", "efr"])
for k, v in out.items():
    print(k, v, flush=True)
(HERE / "work/knnmap481.json").write_text(json.dumps(out, indent=1))
