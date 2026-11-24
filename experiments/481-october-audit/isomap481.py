"""481 F3: predict the builder's design strength on Arc's 2025 H1 panel from the whitened base-only proxy (whiten481),
through a monotone (isotonic) map of builder strength on proxy strength fitted on A/B/C (453 dump) + efr (465 dump),
with a target bootstrap for the weak shares. Validation: fit on one panel, predict the other's weak share.
usage: uv run python isomap481.py -> work/isomap481.json  (reruns whiten481's composition; ~2.5 GB)"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
import whiten481 as W
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
gate = W.gate_D()
abc_panel = set(pd.read_csv(HERE / "work/panels/abc.csv")["target_gene"].astype(str))
banks = [W.L.bank(n) for n in W.L.MAIN]
full = W.L.E064.operator_columns(banks, W.L.OFFICIAL, abc_panel)
evals, evecs = W.load_op(W.OPS["abc_op"])
prox, panels = {}, {}
for name in ("abc", "efr", "h1_2025"):
    panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
    pre = W.L.Donors438(panel, W.L.MAIN).pre(W.L.MAIN, normrel=True)
    pcols = np.asarray([W.L.OFF_POS[g] for g in panel if g in W.L.OFF_POS])
    pre[:, pcols] = 0.0
    X = W.whiten(pre, full, evals, evecs)
    X[:, pcols] = 0.0
    prox[name] = (np.abs(np.where(gate[None, :], X, 0.0)) > W.TAU).sum(1).astype(float)
    panels[name] = panel
del evecs
known = {"abc": W.dump_strength(ROOT / "experiments/453-abc-dress-rehearsal/work/run/dump_h1x2_cd4x1.npz", panels["abc"], gate),
         "efr": W.dump_strength(ROOT / "experiments/465-both-switches-efr/work/run/dump_h1x2_cd4x1.npz", panels["efr"], gate)}
def fit(names):
    x = np.concatenate([np.log1p(prox[n]) for n in names]); y = np.concatenate([known[n] for n in names])
    return IsotonicRegression(out_of_bounds="clip").fit(x, y)
def weak_stats(s):
    return {"frac_le1": round(float((s <= 1.5).mean()), 3), "q33": float(np.quantile(s, 0.33))}
out = {}
# cross-validation: fit on one known panel, predict the other's weak share
for tr, te in (("abc", "efr"), ("efr", "abc")):
    m = fit([tr]); pred = np.round(m.predict(np.log1p(prox[te])))
    out[f"fit_{tr}_predict_{te}"] = {"pred": weak_stats(pred), "actual": weak_stats(known[te])}
    print(f"fit {tr} -> {te}: pred {out[f'fit_{tr}_predict_{te}']['pred']} actual {out[f'fit_{tr}_predict_{te}']['actual']}")
m = fit(["abc", "efr"])
ne = np.array([t not in W.ESS for t in panels["h1_2025"]])
rng = np.random.default_rng(481)
for tag, mask in (("h1_2025_all", np.ones_like(ne)), ("h1_2025_nonessential", ne)):
    pred = np.round(m.predict(np.log1p(prox["h1_2025"][mask])))
    boots = [weak_stats(pred[rng.integers(0, pred.size, pred.size)])["frac_le1"] for _ in range(4000)]
    out[tag] = {**weak_stats(pred), "n": int(mask.sum()), "frac_le1_boot95": np.quantile(boots, [0.025, 0.975]).round(3).tolist()}
    print(tag, out[tag])
for n in ("abc", "efr"):
    out[f"{n}_builder"] = weak_stats(known[n])
(HERE / "work/isomap481.json").write_text(json.dumps(out, indent=1))
