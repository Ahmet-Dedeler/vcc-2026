"""262 -- PILOT (not a screen, not a price): one-variable read of the two hypotheses the diagnosis motivates.

    uv run python pilot262.py

(a) Agreement-keyed CD4 row weight.  Per target, the observable donor-donor agreement a_p (cosine between the
    fold's two un-whitened donor rows for that target, from the fold caches' truths), rank-transformed within the
    fold to r_p in [0, 1] (0 = most disagreeing); the deflated CD4 row is multiplied by
        w_p = exp(beta * (0.5 - r_p))         (median weight 1; rows without the covariate: w = 1)
    and the object is RE-SOLVED by `plolo.solve_fraction` at f = 0.34 (the builder's rule), merged with the union
    at its shipped fraction -- 261's `cd4svd5_parts` with the one extra multiplication.  beta = 0 is `cd4svd5`
    (KA row at 0.0e+00).  Controls: the same weight multiset on DERANGED rows (218's seed 20726), and the reverse
    key (-beta).  One fold per process would be the rule for a screen; this pilot stages the three cancer folds
    serially in one process (~5 GB peak, the operator released per fold).
(b) Headroom for a per-target donor choice on the anti-aligned group, from the arrays already on disk: for each
    cancer fold, each donor's agreement with the held truth (cross-line truth cosine), the mean of the two (the
    base's agreement proxy), the best single donor (oracle), and the donor an independent line (CD4) would pick.
Writes work/pilot262.json and prints.
"""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
import diag262 as G  # noqa: E402  (pins plolo.WORK to ours; L, X, plolo)

L, plolo = G.L, G.plolo
spec = importlib.util.spec_from_file_location("an", HERE / "analyse262.py")
an = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(an)
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BETAS = (1.0, 2.0, 4.0)
SEED = L.SHUF_SEED
out: dict = {"a": {}, "b": {}}


def derange(n, rng):
    perm = np.arange(n)
    while True:
        rng.shuffle(perm)
        if not np.any(perm == np.arange(n)):
            return perm


def weighted(d, w):
    eff = np.asarray(d.mean_effect, dtype=np.float64) * w[:, None]
    return G.E062.CD4Prepared(mode=d.mode + "|rowweight", mean_effect=eff, gene_votes=np.asarray(d.gene_votes).copy(),
                              target_condition_count=np.asarray(d.target_condition_count).copy(),
                              condition_details=d.condition_details)


print("=== (a) agreement-keyed CD4 row weight, resolved sizing; three cancer folds")
ref = {f: json.loads((WORK / f"diag262-{f}.json").read_text())["arms"]["cd4svd5"] for f in CANCER}
for held in CANCER:
    t0 = time.time()
    st = L.stage(held)
    targets = list(st["targets"])
    d, _ = L.deflate(st["cd4"], "svd", 5)
    a = an.TAB[held][1]["agree_donor_donor"]
    ok = np.isfinite(a)
    r = np.full(len(targets), 0.5)
    r[ok] = (np.argsort(np.argsort(a[ok])) + 0.5) / ok.sum()
    rows = {}

    def run(name, dw):
        fr, pu = L.solve_fraction(st["whitened"], dw, L.CD4_F)
        final = L.compose(st, [(dw, float(fr)), (st["union"], st["frac_x"])])
        e = L.score(st, final)
        rows[name] = dict(pds=float(e["pds"]), cos_agg=float(e["cos_agg"]), macro=float(e["a23_cos"]), fraction=float(fr),
                          pds_per_target=[float(x) for x in e["_pds_per_target_a23"]])
        print(f"  [{held}/{name}] pds {e['pds']:.6f} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} frac {fr:.5f}", flush=True)
        return e

    e0 = run("beta0", weighted(d, np.ones(len(targets))))
    ka = max(abs(e0[k] - ref[held][k]) for k in ("pds", "cos_agg", "a23_cos"))
    print(f"  KA beta=0 vs cd4svd5: max|d| {ka:.1e} -> {'PASS' if ka < 1e-9 else 'FAIL'}", flush=True)
    assert ka < 1e-9
    rng = np.random.default_rng(SEED)
    perm = derange(len(targets), rng)
    for beta in BETAS:
        w = np.exp(beta * (0.5 - r))
        run(f"beta{beta:g}", weighted(d, w))
        run(f"beta{beta:g}_shuf", weighted(d, w[perm]))
        run(f"beta{beta:g}_rev", weighted(d, np.exp(-beta * (0.5 - r))))
    base = rows["beta0"]
    for k, v in rows.items():
        v["d_pds"] = v["pds"] - base["pds"]
        v["rel_cos_agg"] = 100.0 * (v["cos_agg"] / base["cos_agg"] - 1.0)
        v["rel_macro"] = 100.0 * (v["macro"] / base["macro"] - 1.0)
    out["a"][held] = dict(ka=ka, rows=rows, seconds=round(time.time() - t0, 1))
    del st

print("\n  summary (three cancer folds, mean of per-fold values):")
print(f"  {'cell':14s} {'d pds':>10s} {'folds+':>7s} {'rel cos_agg %':>14s} {'rel macro %':>12s}")
summ = {}
for name in ["beta0"] + [f"beta{b:g}{s}" for b in BETAS for s in ("", "_shuf", "_rev")]:
    dp = [out["a"][f]["rows"][name]["d_pds"] for f in CANCER]
    rc = [out["a"][f]["rows"][name]["rel_cos_agg"] for f in CANCER]
    rm = [out["a"][f]["rows"][name]["rel_macro"] for f in CANCER]
    summ[name] = dict(d_pds=float(np.mean(dp)), folds_pos=int(sum(x > 0 for x in dp)), rel_cos_agg=float(np.mean(rc)),
                      rel_macro=float(np.mean(rm)), per_fold_d_pds=dp)
    print(f"  {name:14s} {np.mean(dp):+.6f} {sum(x > 0 for x in dp):>4d}/3 {np.mean(rc):+13.3f} {np.mean(rm):+11.3f}   per fold {['%+.5f' % x for x in dp]}")
out["a"]["summary"] = summ

print("\n=== (b) headroom for a per-target donor choice on the anti-aligned group (arrays on disk; no new folds)")
LINE = an.LINE
hb = {}
for held in CANCER:
    t, cov = an.TAB[held]
    donors = an.M[held]["meta"]["sources"]
    A = np.array([[an.pair_cos[(held, s)].get(x, np.nan) for s in donors] for x in t])      # (n, 2) donor-vs-held truth cos
    # CD4 as arbiter: CD4's agreement with donor s's row for target x == cd4_cos_truth measured on fold s
    C = np.array([[an.TAB[s][1]["cov|cd4_cos_truth"][an.tpos[s][x]] if x in an.tpos[s] else np.nan for s in donors] for x in t])
    ok = np.isfinite(A).all(1) & np.isfinite(C).all(1)
    mean_ag = A.mean(1); best = A.max(1)
    pick = np.argmax(C, axis=1)
    arb = A[np.arange(len(t)), pick]
    rng = np.random.default_rng(SEED)
    rnd = A[np.arange(len(t)), rng.integers(0, 2, len(t))]
    e = cov["cov|base_cos_truth"] <= 0
    band = (cov["cov|base_cos_truth"] > 0) & (cov["cov|base_cos_truth"] <= 0.05)
    rho_proxy = spearmanr(mean_ag[ok], cov["cov|base_cos_truth"][ok]).correlation
    res = {}
    for gname, gm in (("all", ok), ("group e (base·truth <= 0)", ok & e), ("marginal band", ok & band)):
        res[gname] = dict(n=int(gm.sum()), mean_of_two=float(mean_ag[gm].mean()), best_single=float(best[gm].mean()),
                          cd4_arbiter=float(arb[gm].mean()), random_pick=float(rnd[gm].mean()),
                          frac_best_above_0_05=float((best[gm] > 0.05).mean()), frac_mean_above_0_05=float((mean_ag[gm] > 0.05).mean()),
                          frac_arbiter_above_0_05=float((arb[gm] > 0.05).mean()),
                          arbiter_picks_best=float((pick[gm] == np.argmax(A, axis=1)[gm]).mean()))
    hb[held] = dict(donors=donors, rho_mean_agreement_vs_base_cos_truth=float(rho_proxy), groups=res)
    print(f"  [{held}] donors {donors}; Spearman(mean donor agreement, base·truth) = {rho_proxy:+.3f} (the proxy's validity)")
    for gname, rr in res.items():
        print(f"     {gname:28s} n {rr['n']:3d}: mean-of-two {rr['mean_of_two']:+.4f}  best-single {rr['best_single']:+.4f}  "
              f"CD4-arbiter {rr['cd4_arbiter']:+.4f}  random-pick {rr['random_pick']:+.4f}  | >0.05: mean {100*rr['frac_mean_above_0_05']:.0f} % "
              f"best {100*rr['frac_best_above_0_05']:.0f} % arbiter {100*rr['frac_arbiter_above_0_05']:.0f} % | arbiter==best {100*rr['arbiter_picks_best']:.0f} %")
out["b"] = hb
(WORK / "pilot262.json").write_text(json.dumps(out, indent=1, default=float))
