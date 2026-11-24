"""481 F3, calibrated: the builder's design strength predicted for panels we have no build for.

For each panel: 438's pre-whitening composition (lib438.Donors438: the three production donors' bank rows, normrel cap
-1.5, the PM50 correction), the panel's columns zeroed, whitened with a 3-donor union operator at gamma 0.65 (float32
arithmetic, the builder's on-axis median preserve), times alpha 0.982. Strength = count of context-D render-gated
coordinates with |value| > tau (the prior gate's measure, tau 0.2). No CD4 / H1 residuals, no global preserve: a
base-only proxy, calibrated against the builder's own strength (453 dump on A/B/C, 465 dump on efr) before it is read
on the 2025 H1 panel (Arc's other panel; essential targets reported separately).

Operators: the A/B/C production operator (e2f65ecf, axis = union - A/B/C panel) for every panel, plus efr's own day-of
re-fit (465) on efr as the check that the operator's panel does not matter for this proxy.
usage: uv run python whiten481.py -> work/whiten481.json   (peak ~2 GB: float32 evecs 1.29 GB)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import lib438 as L  # noqa: E402

GAMMA, ALPHA, TAU = 0.65, 0.982, 0.2
OPS = {"abc_op": ROOT / "data/donors/whitener_union_pruned.npz",
       "efr_op": ROOT / "experiments/465-both-switches-efr/work/run/whitener_union_pruned_def.npz"}
ESS = {ln.strip() for ln in (ROOT / "data/donors/k562_essential/perts.txt").read_text().split("\n") if ln.strip()}


def gate_D() -> np.ndarray:
    z = np.load(HERE / "work/ctrl_def_D.npz", allow_pickle=False)
    tot = z["totals"]
    assert list(z["genes"].astype(str)) == list(L.OFFICIAL)
    return (1e6 * tot / tot.sum() > 10.0) & (np.log1p(5e4 * tot / tot.sum()) > 0)


def load_op(path: Path):
    with np.load(path) as z:
        return np.asarray(z["evals"], dtype=np.float64), np.asarray(z["evecs"], dtype=np.float32)


def whiten(pre: np.ndarray, wcols: np.ndarray, evals, evecs) -> np.ndarray:
    x = np.asarray(pre[:, wcols], dtype=np.float32)
    before = float(np.median(np.linalg.norm(x.astype(np.float64), axis=1)))
    w = (evals ** (-GAMMA / 2.0)).astype(np.float32)
    t = ((x @ evecs) * w) @ evecs.T
    after = float(np.median(np.linalg.norm(t.astype(np.float64), axis=1)))
    out = np.zeros_like(pre)
    out[:, wcols] = t.astype(np.float64) * (before / after)
    return out * ALPHA


def dump_strength(dump: Path, panel, gate) -> np.ndarray:
    z = np.load(dump, allow_pickle=False)
    assert [str(p) for p in z["perts"]] == list(panel)
    car = np.asarray(z["carriers"], dtype=np.float64)
    return (np.abs(np.where(gate[None, :], car, 0.0)) > TAU).sum(1).astype(float)


def main() -> None:
    t0 = time.time()
    gate = gate_D()
    abc_panel = set(pd.read_csv(HERE / "work/panels/abc.csv")["target_gene"].astype(str))
    banks = [L.bank(n) for n in L.MAIN]
    wcols_abc = L.E064.operator_columns(banks, L.OFFICIAL, abc_panel)
    pres = {}
    for name in ("abc", "efr", "h1_2025", "efp", "real199", "pess199"):
        panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
        dn = L.Donors438(panel, L.MAIN)
        pre = dn.pre(L.MAIN, normrel=True)
        pcols = np.asarray([L.OFF_POS[g] for g in panel if g in L.OFF_POS])
        pre[:, pcols] = 0.0
        pres[name] = (panel, pre)
        print(f"[whiten481] composed {name} ({time.time() - t0:.0f}s)", flush=True)
    out = {"panels": {}, "calibration": {}}
    strengths = {}
    for opname, path in OPS.items():
        evals, evecs = load_op(path)
        for name, (panel, pre) in pres.items():
            if opname == "efr_op" and name != "efr":
                continue
            pset = set(panel)
            wcols = (wcols_abc if opname == "abc_op"
                     else L.E064.operator_columns(banks, L.OFFICIAL, pset))
            wcols = np.asarray([c for c in wcols if L.OFFICIAL[c] not in pset], dtype=np.int64)
            if opname == "efr_op":
                assert wcols.size == evals.size, (wcols.size, evals.size)
                W = whiten(pre, wcols, evals, evecs)
            else:
                # the A/B/C operator on another panel: its axis includes that panel's genes; zeroed inputs, dropped outputs
                full = L.E064.operator_columns(banks, L.OFFICIAL, abc_panel)
                W = whiten(pre, full, evals, evecs)
                W[:, [L.OFF_POS[g] for g in panel if g in L.OFF_POS]] = 0.0
            s = (np.abs(np.where(gate[None, :], W, 0.0)) > TAU).sum(1).astype(float)
            strengths[(opname, name)] = s
            ne = np.array([t not in ESS for t in panel])
            out["panels"][f"{name}|{opname}"] = {
                "quartiles": np.quantile(s, [0, .25, .5, .75, 1]).tolist(),
                "frac_le1": round(float((s <= 1).mean()), 3), "frac_eq0": round(float((s == 0).mean()), 3),
                "q33": float(np.quantile(s, 0.33)),
                "nonessential_n": int(ne.sum()),
                "nonessential_frac_le1": round(float((s[ne] <= 1).mean()), 3) if ne.any() else None,
                "nonessential_q33": float(np.quantile(s[ne], 0.33)) if ne.any() else None,
                "median_row_norm_on_axis_before_whitening": round(float(np.median(np.linalg.norm(pre, axis=1))), 4)}
            print(f"  {name:8s} {opname}: {out['panels'][f'{name}|{opname}']}", flush=True)
        del evecs
    known = {"abc": dump_strength(ROOT / "experiments/453-abc-dress-rehearsal/work/run/dump_h1x2_cd4x1.npz",
                                  pres["abc"][0], gate),
             "efr": dump_strength(ROOT / "experiments/465-both-switches-efr/work/run/dump_h1x2_cd4x1.npz",
                                  pres["efr"][0], gate)}
    for name, sb in known.items():
        for opname in ("abc_op", "efr_op"):
            if (opname, name) not in strengths:
                continue
            s = strengths[(opname, name)]
            out["calibration"][f"{name}|{opname}"] = {
                "spearman_vs_builder": round(float(spearmanr(s, sb).statistic), 3),
                "builder_frac_le1": round(float((sb <= 1).mean()), 3), "proxy_frac_le1": round(float((s <= 1).mean()), 3),
                "builder_q33": float(np.quantile(sb, 0.33)), "proxy_q33": float(np.quantile(s, 0.33)),
                "agree_le1_class": round(float(((s <= 1) == (sb <= 1)).mean()), 3)}
            print(f"  calib {name} {opname}: {out['calibration'][f'{name}|{opname}']}", flush=True)
    (HERE / "work/whiten481.json").write_text(json.dumps(out, indent=1))
    print(f"[whiten481] done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
