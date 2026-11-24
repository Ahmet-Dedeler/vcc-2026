#!/usr/bin/env python3
"""537 — P9: what the CLEAN4 hook is worth on the block as a function of its sign accuracy, and what that implies for an
unseen line (D/E/F), where CLEAN4's accuracy is unknown.

x = the realized sign accuracy of the hook's source on H1's covered gate genes (h1-test reference: q < .05, finite non-zero
lfc, not the target), in two populations: the block's swapped-in genes (where the synergy acts) and the rest (neither
block: the hook-alone territory). Blocks as in foldsign537.py (CPM block vs lc + 22 x, 2,200 genes, gene level), so x is
defined exactly as the public folds' accuracies in work/foldsign537.json.
Sources (exactly run537.py's: the mean of CLEAN4's four banks' P_down, sign where |P - .5| >= .05, on the gated axis):
  real  (`+c4h`, 3 seeds)            — one x
  r 0.17 (`+c4h-r0.17`)              — P -> 1 - P on the fixed flip draw default_rng(51717) < 0.17: one x
  shuffle (`+c4h-s`, seed s)         — P permuted among covered genes with default_rng(51718 + s): one x per seed
y = the H1 increment on the block (lcap-kept+c4h* − lcap-kept) and alone (cpm+c4h* − base), per seed and regime, from
work/price537.json (price537.py's rows).
Fit y = a + b (x − 0.5) by least squares over the (variant, seed) points, per regime; predict at the folds' accuracies
(K562 with its own banks held out, HCT116, HEK293T; Orion folds carry a technical blob: relative direction only) and at
H1's. Writes work/dvalue537.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E528 = ROOT / "experiments/528-callprop-block"
TAU = 0.05
REG = ("1.0", "0.5", "0.35", "cmp")
MEM = ("fid", "jac", "reach", "nmae")


def block(score, excl, n=2200):
    b = np.zeros(score.size, bool)
    b[np.argsort(-np.where(excl, -np.inf, score), kind="stable")[:n]] = True
    return b


def main():
    ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    keep = ref["keep"].astype(bool)
    GF = ref["genes"].astype(str)
    gk = GF[keep]
    targets = ref["targets"].astype(str)
    q, lfc = np.asarray(ref["q"], float), np.asarray(ref["lfc"], float)
    gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~np.equal.outer(targets, gk)
    # the blocks (foldsign537's H1 construction)
    H = dict(np.load(E528 / "work/feat528_H1.npz", allow_pickle=True))
    pos = {g: i for i, g in enumerate(H["genes"].astype(str))}
    colh = np.array([pos[g] for g in gk])
    lc = np.log10(np.maximum(H["mean_cpm"][colh], 1e-3))
    excl = np.isin(gk, targets)
    b0, b1 = block(lc, excl), block(lc + 22.0 * H["x"][colh], excl)
    pops = {"swapped_in": b1 & ~b0, "rest": ~b0 & ~b1 & ~excl,
            # the hook's whole territory: on the block (`kept`: all but the steered CPM part) and alone (all but the CPM block)
            "territory_on_block": ~(b0 & b1) & ~excl, "territory_alone": ~b0 & ~excl}
    # run537's CLEAN4 P on the full axis (mean of the four banks; NaN where any bank is missing)
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cpos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([cpos.get(str(g), -1) for g in GF])
    Pf = np.full(len(GF), np.nan)
    Pf[col >= 0] = P[col[col >= 0]]
    Pk = Pf[keep]

    def signs(Pg):
        return np.where(np.isfinite(Pg) & (np.abs(Pg - 0.5) >= TAU), np.where(Pg > 0.5, -1.0, 1.0), 0.0)

    variants = {"real": signs(Pk)}
    u = np.random.default_rng(51717).random(int(keep.sum()))
    variants["r0.17"] = signs(np.where(u < 0.17, 1.0 - Pk, Pk))
    cov = np.flatnonzero(np.isfinite(Pk) & (np.abs(Pk - 0.5) >= TAU))
    for s in (0, 1, 2):
        Ps = Pk.copy()
        Ps[cov] = Pk[np.random.default_rng(51718 + s).permutation(cov)]
        variants[f"shuf_s{s}"] = signs(Ps)
    tr = np.sign(lfc)
    xacc = {}
    for v, sg in variants.items():
        xacc[v] = {}
        for pn, M in pops.items():
            G = gate & M[None, :] & (sg != 0)[None, :]
            xacc[v][pn] = dict(acc=float((tr[G] == np.broadcast_to(sg, G.shape)[G]).mean()), n=int(G.sum()))
    print("realized sign accuracy on H1's covered gate genes (" + " | ".join(pops) + "):")
    for v, d in xacc.items():
        print(f"   {v:8s} " + " | ".join(f"{d[pn]['acc']:.3f} (n {d[pn]['n']})" for pn in pops))
    fs = json.loads((HERE / "work/foldsign537.json").read_text())
    def fold_acc(cls):
        """per-class accuracies, plus the two territories pooled by covered pairs (foldsign537's classes)."""
        out = {pn: cls[pn]["clean4_acc"] for pn in ("swapped_in", "rest")}
        for pn, parts in (("territory_on_block", ("swapped_in", "swapped_out", "rest")), ("territory_alone", ("swapped_in", "rest"))):
            n = sum(cls[c]["covered"] for c in parts)
            out[pn] = sum(cls[c]["clean4_acc"] * cls[c]["covered"] for c in parts) / n
        return out
    folds = {f: fold_acc(fs[f]["classes"]) for f in ("k562_gwps", "orion_hct116", "orion_hek293t")}
    folds["H1 (foldsign537)"] = fold_acc(fs["H1"]["classes"])
    pr = json.loads((HERE / "work/price537.json").read_text())["comparisons"]
    series = {
        "on the block": {"real": "P1 CLEAN4 hook on the cp base", "r0.17": "P8 CLEAN4 hook r 0.17 on the cp base",
                         "shuf": "P9 sign-shuffled CLEAN4 hook on the cp base"},
        "alone": {"real": "P1 CLEAN4 hook on 446bj2z (= 496a's knob on the zs base)",
                  "r0.17": "P8 CLEAN4 hook attenuated r 0.17 on 446bj2z", "shuf": "P9 sign-shuffled CLEAN4 hook on 446bj2z (zero-information floor)"},
    }
    out = {"x": xacc, "folds": folds, "fits": {}}
    for sname, labs in series.items():
        if any(labs[k] not in pr for k in labs):
            print(f"{sname}: prices missing -- run price537.py first")
            continue
        out["fits"][sname] = {}
        for pn in pops:
            out["fits"][sname][pn] = {}
            for reg in REG:
                X, Y, YM = [], [], []
                for k, lab in labs.items():
                    rec = pr[lab]
                    for s, ps in zip(rec["seeds"], rec[reg]["per_seed"]):
                        x = xacc[f"shuf_s{s}" if k == "shuf" else k][pn]["acc"]
                        X.append(x)
                        Y.append(ps["value"])
                        YM.append([ps["d"][m] for m in MEM])
                X, Y, YM = np.array(X), np.array(Y), np.array(YM)
                A = np.vstack([np.ones_like(X), X - 0.5]).T
                coef, res, *_ = np.linalg.lstsq(A, Y, rcond=None)
                resid = Y - A @ coef
                dof = max(len(Y) - 2, 1)
                s2 = float(resid @ resid / dof)
                cv = s2 * np.linalg.inv(A.T @ A)
                cm = np.linalg.lstsq(A, YM, rcond=None)[0]
                pred = {}
                for f, acc in folds.items():
                    xx = np.array([1.0, acc[pn] - 0.5])
                    pred[f] = dict(x=acc[pn], y=float(xx @ coef), se=float(np.sqrt(xx @ cv @ xx)),
                                   members=dict(zip(MEM, (xx @ cm).tolist())))
                out["fits"][sname][pn][reg] = dict(a=float(coef[0]), b=float(coef[1]), se_a=float(np.sqrt(cv[0, 0])),
                                                   se_b=float(np.sqrt(cv[1, 1])), n=len(Y), resid_sd=float(np.sqrt(s2)),
                                                   members_at_0_5=dict(zip(MEM, cm[0].tolist())), pred=pred)
            print(f"hook {sname}, x = accuracy on {pn}:")
            for reg in REG:
                r = out["fits"][sname][pn][reg]
                print(f"   f {reg:4s} y(0.5) {r['a']:+.5f} ± {r['se_a']:.5f}, slope {r['b']:+.4f} ± {r['se_b']:.4f} per unit "
                      f"accuracy (n {r['n']}, resid sd {r['resid_sd']:.5f}); at 0.5 members " + " ".join(
                          f"{m} {r['members_at_0_5'][m]:+.5f}" for m in MEM))
                print("        predicted: " + "; ".join(f"{f.split(' ')[0]} (x {p['x']:.3f}) {p['y']:+.5f} ± {p['se']:.5f}"
                                                     for f, p in r["pred"].items()))
    (HERE / "work/dvalue537.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
