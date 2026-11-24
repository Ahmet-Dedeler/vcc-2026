"""487 decision 1f: the in-mix same-line lift with a PRODUCTION-LIKE whitened composition (the settling check).

lift487.py used an unwhitened unit-mean of the three donors as "ours". The builder's composition is whitened (union
operator, gamma 0.65), and whitening may strip line-specific signal that a raw same-line row then restores. Here "ours"
is 481's base-only proxy of the production composition (whiten481.py's recipe): 438's pre-whitening composition
(Donors438.pre: the three production donors, normrel cap -1.5, the PM50 correction), the targets' own columns zeroed,
whitened with the A/B/C production operator (e2f65ecf) at gamma 0.65 with the on-axis median preserve, x 0.982.
The same-line row is the RAW bank delta (what def_sameline stores). Blend as the builder: unit(ours) + w unit(row).
Statistic: cosine to the independent same-line truth, per target, on the operator axis minus the targets' genes, both
centred over the evaluated targets (lift487's statistic) and uncentred; ratio of means with a target bootstrap.
Cases: HEK293 (truth XuCao, row Orion HEK293T; the Orion lab is also HCT116's) and K562 (truth VIPerturb, row K562-GWPS;
800 non-essential targets, seed 20260925).
KA: the composition for the A/B/C panel reproduces 481's proxy median pre-whitening row norm (whiten481.json abc|abc_op
4.2401; the builder's own meta reads 4.1637 for the full composition with its residuals).
Calibration on the same construction: H1 with Arc's H1 screen as the truth and Pan's rows (new information), clean and
with 450f's sign flips p 0.2 / 0.3 (emitter prices +0.0934 / +0.0268 / +0.0007 at w 2 on h1-test).
Memory: operator float32 1.29 GB; rows for <= 800 targets. Peak ~3 GB. No writes outside this directory
(bytecode off; 438's bank statistics are read from its cache, which exists for all four banks).
usage: PYTHONDONTWRITEBYTECODE=1 uv run python experiments/487-october-decision-review/whitemix487.py -> work/whitemix487.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import csv  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import lib438 as L  # noqa: E402

GAMMA, ALPHA = 0.65, 0.982
OP = ROOT / "data/donors/whitener_union_pruned.npz"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
ESS = set((ROOT / "data/donors/k562_essential/perts.txt").read_text().split())
WS = (0.5, 1, 2, 4)


def whiten(pre, wcols, evals, evecs):
    x = np.asarray(pre[:, wcols], dtype=np.float32)
    before = float(np.median(np.linalg.norm(x.astype(np.float64), axis=1)))
    w = (evals ** (-GAMMA / 2.0)).astype(np.float32)
    t = ((x @ evecs) * w) @ evecs.T
    after = float(np.median(np.linalg.norm(t.astype(np.float64), axis=1)))
    out = np.zeros_like(pre)
    out[:, wcols] = t.astype(np.float64) * (before / after)
    return out * ALPHA


def raw_rows(path, targets):
    perts = [l.strip() for l in open(path / "perts.txt")]
    genes = [l.strip() for l in open(path / "genes.txt")]
    pi = {p: i for i, p in enumerate(perts)}
    gi = {g: i for i, g in enumerate(genes)}
    cols = np.array([gi.get(g, -1) for g in L.OFFICIAL])
    delta = np.load(path / "delta.npy", mmap_mode="r")
    X = np.zeros((len(targets), len(L.OFFICIAL)))
    M = np.zeros((len(targets), len(L.OFFICIAL)), bool)
    for r, t in enumerate(targets):
        v = np.asarray(delta[pi[t]], dtype=np.float64)
        X[r] = np.where(cols >= 0, v[np.maximum(cols, 0)], 0.0)
        M[r] = cols >= 0
    X[~np.isfinite(X)] = 0.0
    return X, M


u = lambda A: A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
cos = lambda A, B: np.sum(u(A) * u(B), 1)


def evaluate(ours, same, truth, keep, rng, label, nb=2000):
    res = {}
    for centred in (True, False):
        O, S, T = ours[:, keep], same[:, keep], truth[:, keep]
        if centred:
            O, S, T = O - O.mean(0), S - S.mean(0), T - T.mean(0)
        c0 = cos(O, T)
        r = {"ours": round(float(c0.mean()), 4), "same": round(float(cos(S, T).mean()), 4)}
        n = len(c0)
        for w in WS:
            cw = cos(u(O) + w * u(S), T)
            bs = []
            for _ in range(nb):
                j = rng.integers(0, n, n)
                bs.append(cw[j].mean() / c0[j].mean())
            r[f"w{w:g}"] = {"ratio": round(float(cw.mean() / c0.mean()), 3),
                            "ci95": [round(float(np.quantile(bs, .025)), 3), round(float(np.quantile(bs, .975)), 3)]}
        res["centred" if centred else "uncentred"] = r
        print(f"  {label:44s} {'centred  ' if centred else 'uncentred'} n {n:4d} ours {r['ours']:.4f} same {r['same']:.4f} "
              + " ".join(f"w{w:g} x{r[f'w{w:g}']['ratio']:.2f} [{r[f'w{w:g}']['ci95'][0]:.2f},{r[f'w{w:g}']['ci95'][1]:.2f}]" for w in WS),
              flush=True)
    return res


def main():
    rng = np.random.default_rng(48761)
    abc_panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
    banks = [L.bank(n) for n in L.MAIN]
    wcols = np.asarray(L.E064.operator_columns(banks, L.OFFICIAL, abc_panel), dtype=np.int64)
    with np.load(OP) as z:
        evals, evecs = np.asarray(z["evals"], dtype=np.float64), np.asarray(z["evecs"], dtype=np.float32)
    assert evals.size == wcols.size, (evals.size, wcols.size)
    # KA: the A/B/C panel's pre-whitening median row norm (481 F6: 4.1637)
    abc = sorted(abc_panel)
    pre_abc = L.Donors438(abc, L.MAIN).pre(L.MAIN, normrel=True)
    pre_abc[:, [L.OFF_POS[g] for g in abc if g in L.OFF_POS]] = 0.0
    ka = float(np.median(np.linalg.norm(pre_abc, axis=1)))
    print(f"KA: A/B/C pre-whitening median row norm {ka:.4f} (481 whiten481.json abc|abc_op: 4.2401)", flush=True)
    assert abs(ka - 4.2401) < 5e-4, ka
    del pre_abc
    out = {"ka_abc_median_norm": ka}
    # ---- calibration: H1, truth = Arc's H1 screen (all cells), row = Pan 2025 hESC (new information), flips
    sys.path.insert(0, str(HERE))
    import lift487 as LF
    donor_sets = [set(l.strip() for l in open(ROOT / f"data/donors/{n}/perts.txt")) for n in L.MAIN]
    h1p = ROOT / "data/donors/h1_2025_full"
    ts = sorted(set(l.strip() for l in open(h1p / "perts.txt")) & LF.pan_targets() & donor_sets[0] & donor_sets[1] & donor_sets[2])
    pre = L.Donors438(ts, L.MAIN).pre(L.MAIN, normrel=True)
    pcols = [L.OFF_POS[g] for g in ts if g in L.OFF_POS]
    pre[:, pcols] = 0.0
    ours = whiten(pre, wcols, evals, evecs)
    del pre
    truth, tm = raw_rows(h1p, ts)
    P, pts, pg = LF.pan_bank(ts)
    assert pts == ts
    pgi = {g: i for i, g in enumerate(pg)}
    pcol = np.array([pgi.get(g, -1) for g in L.OFFICIAL])
    pan = np.where(pcol[None, :] >= 0, P[:, np.maximum(pcol, 0)], 0.0)
    pan[~np.isfinite(pan)] = 0.0
    keep = np.zeros(len(L.OFFICIAL), bool)
    keep[wcols] = True
    keep &= tm[0] & (pcol >= 0)
    keep[pcols] = False
    ess = np.array([t in ESS for t in ts])
    print(f"\nH1 calibration: {len(ts)} targets ({int(ess.sum())} essential), {int(keep.sum())} genes", flush=True)
    cal = {"targets": len(ts), "genes": int(keep.sum())}
    for pflip in (0.0, 0.2, 0.3):
        row = pan if pflip == 0 else LF.flip(pan, pflip, 0)
        cal[f"flip{pflip:g}"] = {"all": evaluate(ours, row, truth, keep, rng, f"H1 Pan flip {pflip:g}, all"),
                                 "non-essential": evaluate(ours[~ess], row[~ess], truth[~ess], keep, rng, f"H1 Pan flip {pflip:g}, non-essential")}
    out["H1 calibration: truth Arc H1, row Pan (new information)"] = cal
    donor_perts = [set(l.strip() for l in open(ROOT / f"data/donors/{n}/perts.txt")) for n in L.MAIN]
    cases = {"HEK293: truth XuCao, row Orion HEK293T": (T225 / "xucao_hek293", ROOT / "data/donors/orion_hek293t", None, L.MAIN),
             "K562: truth VIPerturb, row K562-GWPS": (ROOT / "data/donors/viperturb_k562", ROOT / "data/donors/k562_gwps", 800, L.MAIN),
             # the same truths with the bank's own line LEFT OUT of ours: the lift a new-information cancer bank gets
             "HEK293 as new information (ours = K562 + HCT116)": (T225 / "xucao_hek293", ROOT / "data/donors/orion_hek293t", None,
                                                                   ("k562_gwps", "orion_hct116")),
             "K562 as new information (ours = HCT116 + HEK293T)": (ROOT / "data/donors/viperturb_k562", ROOT / "data/donors/k562_gwps", 800,
                                                                    ("orion_hct116", "orion_hek293t"))}
    for name, (tpath, spath, cap, srcs) in cases.items():
        tperts = set(l.strip() for l in open(tpath / "perts.txt"))
        ts = sorted(tperts & donor_perts[0] & donor_perts[1] & donor_perts[2])
        if cap is not None:
            ts = [t for t in ts if t not in ESS]
            ts = sorted(np.random.default_rng(20260925).choice(ts, size=min(cap, len(ts)), replace=False).tolist())
        pre = L.Donors438(ts, list(srcs)).pre(list(srcs), normrel=True)
        pcols = [L.OFF_POS[g] for g in ts if g in L.OFF_POS]
        pre[:, pcols] = 0.0
        ours = whiten(pre, wcols, evals, evecs)
        del pre
        truth, tm = raw_rows(tpath, ts)
        same, sm = raw_rows(spath, ts)
        keep = np.zeros(len(L.OFFICIAL), bool)
        keep[wcols] = True
        keep &= tm[0] & sm[0]
        keep[pcols] = False
        ess = np.array([t in ESS for t in ts])
        print(f"\n{name}: {len(ts)} targets ({int(ess.sum())} essential), {int(keep.sum())} genes", flush=True)
        res = {"targets": len(ts), "genes": int(keep.sum())}
        res["all"] = evaluate(ours, same, truth, keep, rng, "all targets")
        if (~ess).sum() >= 20 and ess.sum() >= 20:
            idx = np.flatnonzero(~ess)
            res["non-essential"] = evaluate(ours[idx], same[idx], truth[idx], keep, rng, "non-essential")
        # strength strata by the truth row's centred energy (the rows matter most for mse on strong targets)
        tc = truth[:, keep] - truth[:, keep].mean(0)
        e = np.linalg.norm(tc, axis=1)
        for lab, m in (("truth energy top half", e >= np.median(e)), ("truth energy bottom half", e < np.median(e))):
            idx = np.flatnonzero(m)
            res[lab] = evaluate(ours[idx], same[idx], truth[idx], keep, rng, lab)
        out[name] = res
    (HERE / "work/whitemix487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
