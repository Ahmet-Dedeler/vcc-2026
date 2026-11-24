"""487 decision 1d: the in-mix lift with the PRODUCTION composition (not the unwhitened mix3 of lift487.py).

465's stack dump (`465/work/run/dump_h1x2_cd4x1.npz`, the October path on the efr synthetic panel: whitened 3-donor
composition + H1 x2 + CD4 x1) holds our bulk rows and DE carriers for 300 October-like targets. For the targets VIPerturb
K562 also screened, compare cos(ours, truth) with cos(unit(ours) + w * unit(K562-GWPS row), truth): the same-line mode for
a K562 context, whose bank is already one of our three donors. Same for XuCao HEK293 against Orion HEK293T (few targets),
and, as the new-information control, the RPE1 / Jurkat / HepG2 banks are NOT used (no independent truth).
Statistic as in lift487.py: centred cosine (each matrix centred over the evaluated targets), panel genes and each target's
own gene excluded; ratio of means with a target bootstrap. Rows are read one at a time from the mmap (RSS < 1 GB).
usage: uv run python experiments/487-october-decision-review/lift487b.py -> work/lift487b.json
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
D = ROOT / "data/donors"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
ESS = set((D / "k562_essential/perts.txt").read_text().split())
WS = (0.5, 1, 2, 4)
u = lambda A: A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
cos = lambda A, B: np.sum(u(A) * u(B), 1)


def bank(path):
    return (np.load(path / "delta.npy", mmap_mode="r"), [l.strip() for l in open(path / "perts.txt")],
            [l.strip() for l in open(path / "genes.txt")])


def rows_on(b, ts, genes):
    delta, perts, g = b
    pi = {p: i for i, p in enumerate(perts)}
    gi = {x: i for i, x in enumerate(g)}
    cols = np.array([gi.get(x, -1) for x in genes])
    X = np.zeros((len(ts), len(genes)))
    for r, t in enumerate(ts):
        row = np.asarray(delta[pi[t]], dtype=np.float64)
        X[r] = np.where(cols >= 0, row[np.maximum(cols, 0)], 0.0)
    X[~np.isfinite(X)] = 0.0
    return X


def run(truth_name, truth, same, dump, key, rng, nb=2000):
    perts = list(dump["perts"].astype(str))
    genes_all = list(dump["genes"].astype(str))
    ts = [t for t in perts if t in set(truth[1]) and t in set(same[1])]
    panel = set(perts)
    genes = [x for x in genes_all if x not in panel and x in set(truth[2]) and x in set(same[2])]
    gi = {x: i for i, x in enumerate(genes_all)}
    cols = np.array([gi[x] for x in genes])
    pi = {p: i for i, p in enumerate(perts)}
    ours = np.asarray(dump[key])[[pi[t] for t in ts]][:, cols]
    T = rows_on(truth, ts, genes)
    S = rows_on(same, ts, genes)
    Tc, Oc, Sc = T - T.mean(0), ours - ours.mean(0), S - S.mean(0)
    ess = np.array([t in ESS for t in ts])
    e = np.linalg.norm(Tc, axis=1)
    masks = {"all": np.ones(len(ts), bool), "non-essential": ~ess,
             "non-essential, truth energy top half": (~ess) & (e >= np.quantile(e[~ess], 0.5))}
    res = {"targets": len(ts), "genes": len(genes), "essential": int(ess.sum())}
    c0 = cos(Oc, Tc)
    for name, m in masks.items():
        idx = np.flatnonzero(m)
        r = {"n": int(m.sum()), "ours": round(float(c0[idx].mean()), 4), "same": round(float(cos(Sc, Tc)[idx].mean()), 4)}
        for w in WS:
            cw = cos(u(Oc) + w * u(Sc), Tc)
            bs = [cw[j].mean() / c0[j].mean() for j in (idx[rng.integers(0, len(idx), len(idx))] for _ in range(nb))]
            r[f"w{w:g}"] = {"ratio": round(float(cw[idx].mean() / c0[idx].mean()), 3),
                            "ci95": [round(float(np.quantile(bs, .025)), 3), round(float(np.quantile(bs, .975)), 3)]}
        res[name] = r
        print(f"  {truth_name:10s} {key:9s} {name:38s} n {r['n']:4d} ours {r['ours']:.4f} same {r['same']:.4f} "
              + " ".join(f"w{w:g} x{r[f'w{w:g}']['ratio']:.2f} [{r[f'w{w:g}']['ci95'][0]:.2f},{r[f'w{w:g}']['ci95'][1]:.2f}]" for w in (0.5, 2)),
              flush=True)
    return res


def main():
    rng = np.random.default_rng(4872)
    dump = np.load(ROOT / "experiments/465-both-switches-efr/work/run/dump_h1x2_cd4x1.npz", allow_pickle=True)
    out = {}
    vip, gwps = bank(D / "viperturb_k562"), bank(D / "k562_gwps")
    xu, hek = bank(T225 / "xucao_hek293"), bank(D / "orion_hek293t")
    for key in ("bulk", "carriers"):
        out[f"K562_{key}"] = run("VIPerturb", vip, gwps, dump, key, rng)
        out[f"HEK293_{key}"] = run("XuCao", xu, hek, dump, key, rng) if len(set(xu[1]) & set(dump["perts"].astype(str))) >= 20 else None
    (HERE / "work/lift487b.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
