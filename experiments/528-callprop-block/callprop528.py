#!/usr/bin/env python3
"""528 — KA-a and the capture screen: does a call-propensity block capture more of the truth's DE genes than the shipped
top-CPM block, on the EXACT production call set?

KA-a: 517's callprop517 numbers (2-fold by target on h1-test, block 2,200 among non-panel genes ranked by pooled CPM vs the
fitted propensity; capture of held-out gate genes 0.3386 / 0.3672, per-target 0.3896 / 0.4462) reproduced exactly from
528's cached features (feat528_H1.npz, itself equal to 517's excess_dropout to 0.0).

Then the production construction (lib521's reconstruction of 435's CallSet(400, 2200, 1.9): head = top 400 by |key| per
target, block = top 2,200 by the emitter's ref_lin among the rest of the 5-CPM universe) against blocks ranked by
propensity keys with the SAME head and size, per target:
  capture   = |block ∩ gate| / |gate|, pooled over held-out targets and per-target mean (gate = q < .05, finite non-zero
              lfc, not the target; 517's definition)
  jac_emul  = the jac member if the prescribed set C = head ∪ block were delivered exactly (lib521.jac_rows on its pool)
  prec      = |block ∩ pool| / |block|
Keys (all fitted on the OTHER half of the targets; 517's halves, default_rng(517)): 517's logistic on [lc, x, lc x, lc^2]
with lc = log10 pooled CPM (p517); the same on log10 ref_lin (p517r); one-parameter lin keys log10 ref_lin + beta x
(beta from the fit); partial swaps (bottom 10 / 25 / 50 % of the CPM block replaced by the top propensity genes not in it).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/callprop528.py
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402  (read-only import)

N_CPM = 2200


def logit_fit(F, y, n):
    """517's weighted-logistic fit of per-gene call rates y (n targets each)."""
    w = np.full(y.size, n)
    Xr = np.repeat(F, 2, 0)
    yr = np.tile([1, 0], y.size)
    wr = np.column_stack([y * w, (1 - y) * w]).ravel()
    return LogisticRegression(max_iter=2000).fit(Xr, yr, sample_weight=wr)


def feats(lc, x):
    return np.column_stack([lc, x, lc * x, lc ** 2])


def ka_517(H, z):
    genes, keep, targets = z["genes"].astype(str), np.asarray(z["keep"], bool), list(z["targets"].astype(str))
    q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
    kg = genes[keep]
    x = H["x"][keep]
    lc = np.log10(np.maximum(H["pooled_cpm"][keep], 1e-3))
    own = np.array([[g == t for g in kg] for t in targets])
    gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own
    tset = set(targets)
    ok = np.array([g not in tset for g in kg])
    rng = np.random.default_rng(517)
    perm = rng.permutation(len(targets))
    halves = (perm[: len(perm) // 2], perm[len(perm) // 2:])
    F = feats(lc, x)
    res = {}
    for fit, ev in (halves, halves[::-1]):
        y = gate[fit].mean(0)
        m = logit_fit(F, y, len(fit))
        prop = m.predict_proba(F)[:, 1]
        for name, score in (("cpm", lc), ("propensity", prop)):
            blk = np.zeros(kg.size, bool)
            blk[np.argsort(-np.where(ok, score, -np.inf))[:2200]] = True
            cap = gate[ev][:, blk].sum() / gate[ev].sum()
            per_t = np.mean([gate[t][blk].sum() / max(gate[t].sum(), 1) for t in ev if gate[t].sum() > 0])
            res.setdefault(name, []).append((float(cap), float(per_t)))
    got = {k: [float(np.mean([a for a, _ in v])), float(np.mean([b for _, b in v]))] for k, v in res.items()}
    want = json.load(open(ROOT / "experiments/517-flex-programme/work/callprop517.json"))
    want = {k: [float(np.mean([r[0] for r in want[k]])), float(np.mean([r[1] for r in want[k]]))] for k in ("cpm", "propensity")}
    ok_ = all(abs(got[k][i] - want[k][i]) < 1e-12 for k in want for i in range(2))
    return dict(got=got, want=want, ok=ok_), halves


def block_by(score, head_i, n=N_CPM):
    """top-n by `score` among non-head genes, stable (the emitter's tie rule on the CPM key)."""
    b = np.zeros(score.size, bool)
    b[np.argsort(-np.where(head_i, -np.inf, score), kind="stable")[:n]] = True
    return b


def main():
    H = dict(np.load(HERE / "work/feat528_H1.npz", allow_pickle=True))
    z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    ka, halves = ka_517(H, z)
    print("KA-a (517's callprop):", json.dumps(ka), flush=True)
    assert ka["ok"], "KA-a failed"
    S = L.load_state(0)
    keep = np.asarray(z["keep"], bool)
    targets = S["targets"]
    T, G = S["key"].shape
    ref_lin = S["ref_lin"]
    assert np.array_equal(np.argsort(-ref_lin, kind="stable"), np.argsort(-H["mean_cpm"][keep], kind="stable"))
    x = H["x"][keep]
    lcp = np.log10(np.maximum(H["pooled_cpm"][keep], 1e-3))
    lcr = np.log10(np.maximum(ref_lin, 1e-3))
    gate = (S["rq"] < 0.05) & np.isfinite(S["rl"]) & (S["rl"] != 0) & ~S["own"]
    pool = L.pools(S, 1.0)
    head = S["head"]
    # the production block, reconstructed (KA: equals lib521's)
    blk0 = np.array([block_by(ref_lin, head[i]) for i in range(T)])
    assert np.array_equal(blk0, S["block"]), "production block reconstruction differs from lib521"
    tset = set(targets)
    notpanel = np.array([g not in tset for g in S["genes"]])
    out = {"KA_a": ka}
    keys = {}
    for h, (fit, ev) in enumerate((halves, halves[::-1])):
        y = gate[fit].mean(0)
        m1 = logit_fit(feats(lcp, x), y, len(fit))
        m2 = logit_fit(feats(lcr, x), y, len(fit))
        m3 = logit_fit(np.column_stack([lcr, x]), y, len(fit))            # linear-logit, 2 features
        keys.setdefault("p517", {}).update({int(t): m1.predict_proba(feats(lcp, x))[:, 1] for t in ev})
        keys.setdefault("p517r", {}).update({int(t): m2.predict_proba(feats(lcr, x))[:, 1] for t in ev})
        beta = float(m3.coef_[0][1] / m3.coef_[0][0])                    # log10 ref_lin units per unit x
        keys.setdefault("lin2", {}).update({int(t): lcr + beta * x for t in ev})
        out[f"half{h}_lin2_beta"] = beta
        out[f"half{h}_p517r_coef"] = m2.coef_[0].tolist() + [float(m2.intercept_[0])]
    for b in (5.0, 10.0, 20.0, 40.0):
        keys[f"lin_b{b:g}"] = {t: lcr + b * x for t in range(T)}
    rows = {}

    def summarize(name, blocks):
        cap_num = sum(int((blocks[i] & gate[i]).sum()) for i in range(T))
        cap_den = int(gate.sum())
        per = [float((blocks[i] & gate[i]).sum() / gate[i].sum()) for i in range(T) if gate[i].sum() > 0]
        C = head | blocks
        jr = L.jac_rows(C, pool)
        prec = float(np.mean([(blocks[i] & pool[i]).sum() / blocks[i].sum() for i in range(T)]))
        swapped = float(np.mean([(blocks[i] & ~blk0[i]).sum() for i in range(T)]))
        lo_cpm = float(np.median(np.concatenate([ref_lin[blocks[i] & ~blk0[i]] for i in range(T)]))) if swapped else float("nan")
        rows[name] = dict(capture=cap_num / cap_den, capture_per_target=float(np.mean(per)), jac_emul=float(np.mean(jr)),
                          block_precision=prec, swapped_per_target=swapped, swapped_in_median_ref_lin=lo_cpm,
                          jac_rows=jr)
        print(f"{name:14s} capture {cap_num / cap_den:.4f} per-target {np.mean(per):.4f}  jac_emul {np.mean(jr):.5f}  "
              f"precision {prec:.4f}  swapped/target {swapped:7.1f}  swapped-in median ref_lin {lo_cpm:8.2f}", flush=True)
    summarize("cpm (shipped)", blk0)
    for k, kv in keys.items():
        summarize(k, np.array([block_by(kv[i], head[i]) for i in range(T)]))
    # partial swaps: keep the top (1 - f) of the CPM block by ref_lin, fill the rest with the top propensity genes (p517r)
    for f in (0.1, 0.25, 0.5):
        bl = []
        for i in range(T):
            nk = int(round(N_CPM * (1 - f)))
            kept = block_by(ref_lin, head[i], nk)
            bl.append(kept | block_by(keys["p517r"][i], head[i] | kept, N_CPM - nk))
        summarize(f"partial{int(f * 100)}_p517r", np.array(bl))
    # the bootstrap of the jac_emul difference (targets), for the keys
    rng = np.random.default_rng(528)
    base = rows["cpm (shipped)"]["jac_rows"]
    for k in rows:
        d = rows[k]["jac_rows"] - base
        bs = np.array([d[rng.integers(0, T, T)].mean() for _ in range(4000)])
        rows[k]["d_jac_emul"] = float(d.mean())
        rows[k]["d_jac_emul_ci"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
        rows[k].pop("jac_rows")
        print(f"  {k:14s} d jac_emul {d.mean():+.5f} [{np.quantile(bs, .025):+.5f}, {np.quantile(bs, .975):+.5f}]", flush=True)
    out["rows"] = rows
    (HERE / "work/callprop528.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
