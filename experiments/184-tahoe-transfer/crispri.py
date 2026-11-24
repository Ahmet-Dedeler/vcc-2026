"""184 Phase 2d -- the same measurement on OUR CRISPRi donors, so the learning curve's low-k end
is measured rather than extrapolated.

This is NOT phase 3.  No operator is transferred from Tahoe, no h1-test surface is touched, no
pds is computed and nothing is built.  It runs the identical phase-2 arms and nulls on the donor
tables in `data/donors/` (read-only), at the only k values our own data supports (2..5), so the
Tahoe curve and the CRISPRi curve can be read on one axis.

Comparator: 161 measured, on its `core5` surface (the same five lines), cos(orth, TRUE orth) of
0.0606 for its target-side family and 0.0592 for its consensus x context family.  A diagonal
basal gain landing near that says 161's feature set was not the limitation; landing far above it
says it was.

Writes work/crispri.json only.  Reads data/donors/** and
experiments/166-donor-acquisition/donors/** read-only.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common import HERE, SEED, WORK, unit
from phase2 import cos_rows, fit_gain, gain_features, orth, sign_acc_idx, top_index

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DONORS = {
    "k562_gwps": os.path.join(REPO, "data/donors/k562_gwps"),
    "rpe1": os.path.join(REPO, "data/donors/rpe1"),
    "nadig_jurkat": os.path.join(REPO, "data/donors/nadig_jurkat"),
    "orion_hct116": os.path.join(REPO, "data/donors/orion_hct116"),
    "orion_hek293t": os.path.join(REPO, "data/donors/orion_hek293t"),
    "nadig_hepg2": os.path.join(REPO, "experiments/166-donor-acquisition/donors/nadig_hepg2"),
}
CORE5 = ["k562_gwps", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"]
CORE6 = CORE5 + ["nadig_hepg2"]


def load_donor(path):
    genes = np.asarray([l.strip() for l in open(os.path.join(path, "genes.txt")) if l.strip()])
    perts = np.asarray([l.strip() for l in open(os.path.join(path, "perts.txt")) if l.strip()])
    delta = np.load(os.path.join(path, "delta.npy"), mmap_mode="r")
    bctrl = np.load(os.path.join(path, "b_ctrl.npy"))
    assert delta.shape == (len(perts), len(genes))
    return genes, perts, delta, bctrl


def build(names, cap=3000):
    meta = {n: load_donor(DONORS[n]) for n in names}
    genes = sorted(set.intersection(*[set(m[0].tolist()) for m in meta.values()]))
    perts = sorted(set.intersection(*[set(m[1].tolist()) for m in meta.values()]))
    genes = [g for g in genes if g not in set(perts)]      # drop each target's own coordinate
    rng = np.random.default_rng(SEED)
    if len(perts) > cap:
        perts = [perts[i] for i in sorted(rng.choice(len(perts), cap, replace=False))]
    R = np.zeros((len(perts), len(names), len(genes)), dtype=np.float32)
    B = np.zeros((len(names), len(genes)), dtype=np.float32)
    for j, n in enumerate(names):
        g, p, d, b = meta[n]
        gi = {x: i for i, x in enumerate(g)}
        pi = {x: i for i, x in enumerate(p)}
        ci = np.asarray([gi[x] for x in genes])
        ri = np.asarray([pi[x] for x in perts])
        R[:, j, :] = np.asarray(d[ri])[:, ci]
        B[j] = b[ci]
    R -= R.mean(0, keepdims=True)                          # the across-target mean is the zero
    print(f"{len(names)} lines x {len(perts)} shared targets x {len(genes)} genes")
    return R, B, np.asarray(genes), np.asarray(perts)


def run(R, B, names, ks, draws=6):
    T, L, G = R.shape
    rng = np.random.default_rng(SEED)
    rows = []
    for held in range(L):
        truth = R[:, held, :]
        idx = top_index(truth, 200)
        others = [l for l in range(L) if l != held]
        for k in ks:
            if k > len(others):
                continue
            nd = 1 if k == len(others) else draws
            for _ in range(nd):
                S = np.asarray(others if k == len(others)
                               else rng.choice(others, size=k, replace=False))
                M = R[:, S, :].mean(1)
                am = np.abs(np.take_along_axis(M, idx, 1))
                half = idx.shape[1] // 2
                weak = np.take_along_axis(idx, np.argpartition(am, half - 1, 1)[:, :half], 1)
                mu, sd = B[S].mean(0), B[S].std(0) + 1e-6
                Z = (B - mu) / sd
                gh = np.asarray([fit_gain(R[:, l, :], R[:, np.delete(S, i), :].mean(1))
                                 for i, l in enumerate(S)]) if k > 1 else None
                gbar = gh.mean(0)
                X = np.concatenate([gain_features(gbar, Z[l]) for l in S])
                y = np.concatenate(list(gh))
                cf = np.linalg.lstsq(X[:, :2], y, rcond=None)[0]
                cx = np.linalg.lstsq(X, y, rcond=None)[0]
                gb = gain_features(gbar, Z[held]) @ cx
                zn = Z[rng.choice([l for l in range(L) if l != held])]
                preds = {
                    "M": M,
                    "gain-mean": M * (gain_features(gbar, Z[held])[:, :2] @ cf),
                    "gain-basal": M * gb,
                    "gain-basal-clip": M * np.maximum(gb, 0.0),
                    "gain-basal/perm-line": M * (gain_features(gbar, zn) @ cx),
                    "gain-ORACLE": M * fit_gain(truth, M),
                }
                perm = rng.permutation(T)
                preds["M/perm-drug"] = M[perm]
                for arm, P in preds.items():
                    rows.append({
                        "held": names[held], "k": int(k), "arm": arm,
                        "sign": float(np.nanmean(sign_acc_idx(P, truth, idx))),
                        "sign_weak": float(np.nanmean(sign_acc_idx(P, truth, weak))),
                        "cos": float(np.nanmean(cos_rows(P, truth))),
                        "cos_orth": float(np.nanmean(cos_rows(orth(P, M), orth(truth, M)))),
                    })
        print(f"  held {names[held]:16s} done")
    return rows


def summarise(rows, ks):
    import collections
    agg = collections.defaultdict(list)
    for r in rows:
        agg[(r["k"], r["arm"])].append(r)
    arms = ["M", "gain-mean", "gain-basal", "gain-basal-clip", "gain-basal/perm-line",
            "gain-ORACLE", "M/perm-drug"]
    for m in ["sign", "sign_weak", "cos", "cos_orth"]:
        print(f"  --- {m}")
        print(f"  {'arm':24s}" + "".join(f"k={k:<9}" for k in ks))
        for a in arms:
            print(f"  {a:24s}" + "".join(
                f"{np.mean([x[m] for x in agg[(k, a)]]):+.4f}   " if agg[(k, a)] else "  --      "
                for k in ks))


def main():
    out = {}
    for tag, names, ks in [("core5", CORE5, [2, 3, 4]), ("core6", CORE6, [2, 3, 5])]:
        print(f"===== {tag}")
        R, B, genes, perts = build(names)
        rows = run(R, B, names, ks)
        summarise(rows, ks)
        out[tag] = {"lines": names, "n_targets": int(R.shape[0]), "n_genes": int(R.shape[2]),
                    "ks": ks, "rows": rows}
    with open(os.path.join(WORK, "crispri.json"), "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
