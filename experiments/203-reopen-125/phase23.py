"""203 Phase 2 + Phase 3 — the direction-only variant, `nmae`, and the mechanism.

PHASE 2 (the variant 125 itself proposed).  "122's magnitudes with 125's posterior
ordering" is, in production terms, a pure ORDERING arm: it changes only the key the DE
call is ranked by, so `pds`, `mse`, `nmae` and the cosine are bit-identical to the
incumbent by construction (§8.2 rule 2, board-verified twice).  It can therefore only pay
through `reach`/`fid`/`jac`.  Measured here: |truth|-weighted sign precision of the top-K
coordinates under four ranking keys, on the PLOLO folds' real truth, at K = 20% (125's own
depth) and K = 89 (188's measured reference BH depth).

PHASE 3 `nmae`.  188's identity, `nmae = 1 + (U - 2M)/D`, in log2 fold-change space, on a
SIGNIFICANCE-selected pool (`|lfc|*sqrt(cpm)`, 188's measured proxy for the reference's own
Wilcoxon BH set) at the reference's depth of 89 genes/row.  Two other pool constructions
are reported as sensitivity because 188 showed the construction can invert the sign.

MECHANISM.  Is the hierarchical layer adding information, or is the best setting simply the
one that behaves most like plain averaging?  Measured as the per-target cosine between each
posterior and the equal-mean estimator, against that setting's `pds` delta.

KNOWN-ANSWER ROWS
  KA1  the `a15` arm's `nmae*` and pool composition must reproduce across folds with the
       pool's median control CPM far above the gate median -- 188 measured the reference
       BH pool at ~10x the gate median; a pool that does not show that is not the
       reference's pool.
  KA2  `shrink6.78` must leave every RANKING key's precision EXACTLY unchanged (a positive
       rescale cannot reorder |H|), while its `nmae` moves -- the two channels separate.
  KA3  a random key must land at the base rate.
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

warnings.filterwarnings("ignore", message="All-NaN slice encountered")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "src"))

import est203  # noqa: E402
import plolo  # noqa: E402
from mapfit import to_lfc  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, solve_fraction  # noqa: E402

OUT = HERE / "work"
A15 = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10)
SHRINK = 6.7801
POOL_N = 89                 # 188: h1-test's median reference n_real
NMAE_SPAN = 0.60            # 188/169
RNG = np.random.default_rng(203)


def topn(key, n):
    out = np.zeros(key.shape, bool)
    for i in range(key.shape[0]):
        k = min(n, key.shape[1])
        out[i, np.argpartition(-key[i], k - 1)[:k]] = True
    return out


def umd(P, R, pool):
    U, M, D = [], [], []
    for i in range(P.shape[0]):
        s = np.flatnonzero(pool[i])
        h, r = np.nan_to_num(P[i, s]), R[i, s]
        d = float(np.abs(r).sum())
        if not np.isfinite(d) or d == 0:
            continue
        m = np.where(np.sign(h) == np.sign(r), np.minimum(np.abs(h), np.abs(r)), 0.0)
        U.append(np.abs(h).sum())
        M.append(m.sum())
        D.append(d)
    return np.array(U), np.array(M), np.array(D)


def key_precision(key, pred, truth, frac=None, n=None):
    """|truth|-weighted sign precision of the top-K coordinates chosen by `key`."""
    vals = []
    for i in range(key.shape[0]):
        k = n if n is not None else max(1, int(np.ceil(frac * key.shape[1])))
        k = min(k, key.shape[1])
        ch = np.argpartition(-key[i], k - 1)[:k]
        w = np.abs(truth[i, ch])
        ok = (np.sign(pred[i, ch]) != 0) & (np.sign(pred[i, ch]) == np.sign(truth[i, ch]))
        s = float(w.sum())
        if s > 0:
            vals.append(float(np.where(ok, w, 0.0).sum() / s))
    return float(np.mean(vals))


def main():
    acc, prec, mech = {}, {}, {}
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        z = np.load(OUT / f"sources-{held}.npz", allow_pickle=True)
        panel = est203.prepare(np.asarray(z["effects"], dtype=np.float64),
                               available=z["available"], support=z["support"])
        targets = [str(t) for t in z["targets"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        eq_eff, eq_conf, _ = est203.equal_mean_agreement(panel)
        bp_eff, bp_conf, _ = est203.nested_posterior(panel, targets,
                                                     est203.SETTINGS["eb_broad_prior"])
        cs_eff, cs_conf, _ = est203.nested_posterior(panel, targets,
                                                     est203.SETTINGS["eb_conservative"])
        arms = {"a15": pre, "eq122": eq_eff, f"shrink{SHRINK:.2f}": pre / SHRINK,
                "post125:eb_broad_prior": bp_eff, "post125:eb_conservative": cs_eff}

        # ---- mechanism: how close is each posterior to the plain equal mean?
        for nm, setting in est203.SETTINGS.items():
            e, _, _ = est203.nested_posterior(panel, targets, setting)
            c = []
            for i in range(e.shape[0]):
                a, b = e[i], eq_eff[i]
                na, nb = np.linalg.norm(a), np.linalg.norm(b)
                if na > 0 and nb > 0:
                    c.append(float(a @ b / (na * nb)))
            mech.setdefault(nm, []).append(float(np.median(c)))

        # ---- deliver in the challenge frame (context A; the three are within 1e-4)
        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
        pools = {
            "wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
            "lfc89": topn(np.abs(LR), POOL_N),
            "h89": topn(np.abs(truth_h[:, sc]), POOL_N),
        }
        for pk, m in pools.items():
            med = float(np.median(np.concatenate([cpm[m[i]] for i in range(m.shape[0])])))
            print(f"  [{held}] pool {pk:<8} median ctrl CPM in pool {med:9.2f} "
                  f"(gate median {np.median(cpm):.2f}, ratio {med/np.median(cpm):.1f}x)",
                  flush=True)

        delivered = {}
        for name, eff in arms.items():
            p = np.asarray(eff, dtype=np.float64) * A15["alpha"]
            wh = op.apply(p, fold["wcols"], A15["gamma"], A15["lam"])
            frac, _ = solve_fraction(wh, cd4, A15["cd4_f"])
            fin, _ = plolo.E062.add_post_whiten_residual(wh, cd4, frac, fold["wcols"])
            fin[:, fold["panel_positions"]] = 0.0
            H = deliver(fin, bA, gA, None)
            delivered[name] = H
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, m in pools.items():
                U, M, D = umd(LP, LR, m)
                acc.setdefault((name, pk), []).append((U, M, D))

        # ---- Phase 2: ranking keys, on the SAME delivered a15 profile
        Ha = delivered["a15"][:, sc]
        Ts = truth_h[:, sc]
        keys = {"|H| (production, 185/191)": np.abs(Ha),
                "122 confidence": eq_conf[:, sc],
                "125 posterior confidence": bp_conf[:, sc],
                "125 conservative confidence": cs_conf[:, sc],
                "random": RNG.random(Ha.shape)}
        for kn, kv in keys.items():
            prec.setdefault((kn, "top20%"), []).append(key_precision(kv, Ha, Ts, frac=0.20))
            prec.setdefault((kn, f"top{POOL_N}"), []).append(
                key_precision(kv, Ha, Ts, n=POOL_N))
        # KA2: a positive rescale cannot reorder |H|
        Hs = delivered[f"shrink{SHRINK:.2f}"][:, sc]
        prec.setdefault(("|H| after 6.78x shrink", "top20%"), []).append(
            key_precision(np.abs(Hs), Hs, Ts, frac=0.20))
        del op

    # ---------------------------------------------------------------- nmae table
    print(f"\n{'='*104}\nPHASE 3 -- `nmae` through 188's identity, lfc space, "
          f"significance-selected pool\n{'='*104}")
    print(f"{'pool':<10}{'arm':<30}{'nmae*':>10}{'d nmae*':>11}{'dU/D':>10}{'dM/D':>10}"
          f"{'e':>8}{'d overall':>12}")
    out = {}
    names = ["a15", "eq122", f"shrink{SHRINK:.2f}",
             "post125:eb_broad_prior", "post125:eb_conservative"]
    for pk in ("wald89", "lfc89", "h89"):
        b = acc[("a15", pk)]
        U0 = np.concatenate([x[0] for x in b])
        M0 = np.concatenate([x[1] for x in b])
        D0 = np.concatenate([x[2] for x in b])
        nm0 = float(np.mean(1 + (U0 - 2 * M0) / D0))
        for name in names:
            v = acc[(name, pk)]
            U = np.concatenate([x[0] for x in v])
            M = np.concatenate([x[1] for x in v])
            nm = float(np.mean(1 + (U - 2 * M) / D0))
            du = float(np.mean((U - U0) / D0))
            dm = float(np.mean((M - M0) / D0))
            e = dm / du if abs(du) > 1e-12 else float("nan")
            ov = -(nm - nm0) / NMAE_SPAN / 6.0
            out[f"{pk}|{name}"] = dict(nmae=nm, dnmae=nm - nm0, dU=du, dM=dm, e=e,
                                       overall=ov)
            print(f"{pk:<10}{name:<30}{nm:>10.5f}{nm-nm0:>+11.5f}{du:>+10.5f}"
                  f"{dm:>+10.5f}{e:>8.3f}{ov:>+12.6f}")

    # ---------------------------------------------------------------- phase 2
    print(f"\n{'='*104}\nPHASE 2 -- ordering keys on the IDENTICAL delivered a15 profile "
          f"(|truth|-weighted sign precision)\n{'='*104}")
    print(f"{'key':<34}{'top-20%':>12}{'vs |H|':>10}{'top-89':>12}{'vs |H|':>10}")
    ref20 = float(np.mean(prec[("|H| (production, 185/191)", "top20%")]))
    ref89 = float(np.mean(prec[("|H| (production, 185/191)", f"top{POOL_N}")]))
    keyout = {}
    for kn in ("|H| (production, 185/191)", "122 confidence", "125 posterior confidence",
               "125 conservative confidence", "random"):
        a = float(np.mean(prec[(kn, "top20%")]))
        b = float(np.mean(prec[(kn, f"top{POOL_N}")]))
        keyout[kn] = dict(top20=a, top89=b, d20=a - ref20, d89=b - ref89)
        print(f"{kn:<34}{a:>12.5f}{a-ref20:>+10.5f}{b:>12.5f}{b-ref89:>+10.5f}")
    ka2 = float(np.mean(prec[("|H| after 6.78x shrink", "top20%")]))
    print(f"\nKA2  |H| key after a 6.78x shrink: {ka2:.8f} vs unshrunk {ref20:.8f} "
          f"(diff {ka2-ref20:+.2e}; a positive rescale cannot reorder |H|)")

    # ---------------------------------------------------------------- mechanism
    print(f"\n{'='*104}\nMECHANISM -- is the hierarchical layer adding information, or is the "
          f"best setting the one\nthat behaves most like plain averaging?\n{'='*104}")
    tab = json.loads((OUT / "table203.json").read_text())
    print(f"{'setting':<24}{'median cos(posterior, equal mean)':>36}{'d pds (4-fold)':>17}")
    pts = []
    for nm in est203.SETTINGS:
        c = float(np.mean(mech[nm]))
        dp = tab[f"post125:{nm}"]["d_pds"]
        pts.append((c, dp))
        print(f"{nm:<24}{c:>36.5f}{dp:>+17.6f}")
    c = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts])
    r = float(np.corrcoef(c, y)[0, 1])
    print(f"\nPearson r(similarity-to-equal-mean, d pds) over the four frozen settings = {r:+.4f}")

    (OUT / "phase23.json").write_text(json.dumps(
        dict(nmae=out, keys=keyout, ka2_shrink_key=ka2,
             mechanism={k: float(np.mean(v)) for k, v in mech.items()},
             mechanism_r=r), indent=1))
    print(f"\nwrote {OUT/'phase23.json'}")


if __name__ == "__main__":
    main()
