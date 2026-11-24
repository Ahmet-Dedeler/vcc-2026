"""159 SS5 -- the ceiling on a context-INDEPENDENT donor-transfer model, by construction.

A context-independent model emits one effect matrix and applies it to every context.  The best
such matrix, given infinitely many donor lines, is the ACROSS-LINE MEAN effect for each target:
it is the conditional expectation of a line's effect given only the target's identity.  So the
ceiling is measurable without any model at all --

    predictor for line A  :=  mean over the k OTHER lines of their measured effect for target t
    ceiling               :=  the k -> infinity limit of that predictor's score on line A

and the k = 1..5 curve says how far from the limit a finite donor pool is.  Two things fall out
that algebra alone would not give: the limit is reached almost immediately (the curve is flat
past k = 2, because averaging donors removes donor NOISE but not line-SPECIFIC biology), and
the limiting cosine is directly comparable to the board's implied cosine.

`sse*` is 156's endpoint: min over a global amplitude of ||g*P - T||^2 / ||T||^2 = 1 - cos_agg^2,
the uncorrected `expr_mse` ratio at its own best alpha.  It is a PROXY for the scored `mse`
(the jackknife and rho terms are emitter properties) and is flagged as one.

`pds` is `experiments/153-de-bench/debench.py:Surface.pds` / 156's `endpoints`, recomputed here
so the panel size can be held at 300 -- pds depends on the number of decoys and the challenge
has 300, while the surfaces here carry 1,000+.

SURFACE NOTE (AGENTS.md SS11.1 rule 9 and this experiment's `reconcile.py`): everything is run
on the DONOR tables, not on the thinned Orion harnesses, because the harnesses carry a
cell-pool offset that puts 87% of their effect energy in a target-independent blob and drives
their decoy cosine to +0.179 where the same two lines' donor tables read +0.001.
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

SEED = 20260902
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
PANEL = 300

LINES = {
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
    "k562_gwps": "data/donors/k562_gwps",
    "k562_essential": "data/donors/k562_essential",
    "rpe1": "data/donors/rpe1",
    "nadig_jurkat": "data/donors/nadig_jurkat",
}
# k562_essential is the SAME cell line as k562_gwps; keep it out of cross-line consensus
SAME_LINE_GROUP = {"k562_gwps": "K562", "k562_essential": "K562"}


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def load(rel):
    p = os.path.join(REPO, rel)
    return (np.nan_to_num(np.load(os.path.join(p, "delta.npy")).astype(np.float32)),
            [l.strip() for l in open(os.path.join(p, "perts.txt")) if l.strip()],
            [l.strip() for l in open(os.path.join(p, "genes.txt")) if l.strip()])


def pds(P, T):
    A, B = unit(P.astype(np.float64)), unit(T.astype(np.float64))
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S)
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


def endpoints(P, T):
    A, B = unit(P.astype(np.float64)), unit(T.astype(np.float64))
    own = float((A * B).sum(1).mean())
    a = float((P * P).sum())
    x = float((P * T).sum())
    dd = float((T * T).sum())
    sse = 1.0 - x * x / (a * dd) if a > 0 else 1.0
    G = A @ A.T
    rr = float(G[~np.eye(len(G), dtype=bool)].mean())
    return {"pds": pds(P, T), "cos": own, "sse": sse, "rr": rr}


def main() -> None:
    rng = np.random.default_rng(SEED)
    raw = {k: load(v) for k, v in LINES.items()}
    perts = set.intersection(*[set(v[1]) for v in raw.values()])
    genes = set.intersection(*[set(v[2]) for v in raw.values()])
    ts = sorted(perts)
    gs = sorted(genes - perts)          # pds_cosine drops the panel target coordinates
    print(f"{len(LINES)}-line intersection: {len(ts)} targets x {len(gs)} genes")

    E = {}
    for k, (d, pl, gl) in raw.items():
        ip = {p: i for i, p in enumerate(pl)}
        ig = {g: i for i, g in enumerate(gl)}
        E[k] = d[np.asarray([ip[t] for t in ts])][:, np.asarray([ig[g] for g in gs])
                                                  ].astype(np.float64)
    # every line is centred on its own across-target mean: that mean IS the metric's zero
    # (the "predict the context-mean response for all targets" baseline), so no model can be
    # credited for it and no ceiling should include it.
    C = {k: v - v.mean(0, keepdims=True) for k, v in E.items()}

    out = {"n_targets": len(ts), "n_genes": len(gs), "panel": PANEL, "seed": SEED,
           "lines": list(LINES)}

    # panel subsample: pds depends on the decoy count and the challenge has 300
    sel = np.sort(rng.choice(len(ts), min(PANEL, len(ts)), replace=False))

    print(f"\n{'=' * 108}\nA. CONSENSUS-OF-k DONORS -> the context-independent ceiling "
          f"(panel of {len(sel)} targets)")
    print(f"{'reference line':<16}{'k':>3}{'donors':<44}{'pds':>8}{'cos':>9}{'sse*':>9}")
    rowsA = []
    for ref in LINES:
        others = [k for k in LINES if SAME_LINE_GROUP.get(k, k) != SAME_LINE_GROUP.get(ref, ref)]
        T = C[ref][sel]
        for k in range(1, len(others) + 1):
            best = None
            combos = list(itertools.combinations(others, k))
            if len(combos) > 10:
                combos = [combos[i] for i in rng.choice(len(combos), 10, replace=False)]
            vals = []
            for combo in combos:
                # equal weight on unit-normalised rows: AGENTS.md SS4.2, equal donor weights win
                P = np.mean([unit(C[o])[sel] for o in combo], axis=0)
                vals.append((endpoints(P, T), combo))
            e = {m: float(np.mean([v[0][m] for v in vals])) for m in vals[0][0]}
            best = max(vals, key=lambda v: v[0]["pds"])
            rowsA.append(dict(ref=ref, k=k, n_combo=len(combos), mean=e,
                              best_pds=best[0]["pds"], best_combo=list(best[1])))
            lbl = "+".join(best[1]) if k <= 2 else f"mean of {len(combos)} {k}-subsets"
            print(f"{ref:<16}{k:>3} {lbl:<43}{e['pds']:8.4f}{e['cos']:9.5f}{e['sse']:9.6f}")
    out["consensus"] = rowsA

    print(f"\n{'=' * 108}\nB. THE SAME-LINE POSITIVE CONTROL "
          f"(K562 gwps vs K562 essential -- one cell line, two libraries)")
    T = C["k562_gwps"][sel]
    e = endpoints(unit(C["k562_essential"])[sel], T)
    print(f"  predict K562-gwps from K562-essential:      pds {e['pds']:.4f}  "
          f"cos {e['cos']:.5f}  sse* {e['sse']:.6f}")
    out["same_line_control"] = e
    ecross = endpoints(unit(C["orion_hct116"])[sel], T)
    print(f"  predict K562-gwps from HCT116 (cross-line): pds {ecross['pds']:.4f}  "
          f"cos {ecross['cos']:.5f}  sse* {ecross['sse']:.6f}")
    out["same_line_control_crossref"] = ecross

    print(f"\n{'=' * 108}\nC. SYNTHETIC pds-vs-TRANSFER-LEVEL CURVE\n"
          f"  P_t = tau * s_t/||s_t|| + sqrt(1-tau^2) * z_t, z_t a random combination of the\n"
          f"  OTHER targets' residuals orthogonalised against s_t -- so the decoy geometry is\n"
          f"  the real one and only the target-specific information is dialled.")
    curve = {}
    for ref in ("orion_hct116", "k562_gwps"):
        S = unit(C[ref][sel])
        n, g = S.shape
        W = rng.normal(size=(n, n))
        np.fill_diagonal(W, 0.0)
        Z = W @ S
        Z = Z - (Z * S).sum(1, keepdims=True) * S          # orthogonalise against own truth
        Z = unit(Z)
        rows = []
        print(f"\n  reference {ref}   ({n} targets, {g} genes)")
        print(f"    {'tau':>6}{'pds':>9}{'cos to truth':>14}{'sse*':>10}")
        for tau in (0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0):
            P = tau * S + np.sqrt(max(1 - tau ** 2, 0.0)) * Z
            e = endpoints(P, C[ref][sel])
            rows.append(dict(tau=tau, **e))
            print(f"    {tau:6.2f}{e['pds']:9.4f}{e['cos']:14.5f}{e['sse']:10.6f}")
        curve[ref] = rows
    out["synthetic_curve"] = curve

    with open(os.path.join(WORK, "ceiling.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'ceiling.json')}")


if __name__ == "__main__":
    main()
