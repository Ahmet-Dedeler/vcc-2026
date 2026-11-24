"""099 arm 5 -- add the candidate donor AFTER the whitening operator, as production adds CD4.

`arms.py` at `gamma = 0.65` collapses every added-donor arm (M 0.05-0.51). `diag1.py` shows why:
the operator is fitted on the BASE donors, so directions a new donor introduces are low-variance
in that covariance and are amplified by `lambda^(-gamma/2)`; after the median-norm renormalisation
the amplified component crowds out the signal, and the surviving energy inside the CPM>5 gate
falls from 6.52 to 2.22 (kolf2) and 1.72 (cd4_rest).

That is a property of the FROZEN operator, not necessarily of the donor: production would refit
the whitener if it adopted a fourth donor. Refitting is a new, unvalidated operator, so instead
this measures the route production ALREADY uses for an auxiliary source -- experiment 062's
post-whitening residual, at a fraction of the whitened effect's norm:

    E = whiten(mix(base donors));  E <- E + f * ||E||_med / ||V||_med * V

with `V` the candidate's own effect on the target axis, so `f` is a share of amplitude and the
operator never sees the new donor at all.
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
import pandas as pd
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT, ALPHA

CAND = ("kolf2", "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "perturbmatch", "feng_ipsc",
        "orion_hct116", "h1_2025")
# `cd4_mean3` reproduces experiment 062's production residual: each condition rescaled to the
# donors' reference median norm, then AVERAGED over the three, then added post-whitening.
# `cd4_mean3_gn` adds 062's `gene_norm` (divide by each gene's mean |effect|, floored at the
# 25th percentile, rescaled to preserve the median row norm), which is production's `CD4_MODE`.
# `cd4_mean3_gn_prod` adds the two policies the shipped builder enforces on top of that:
#   * every VCC panel coordinate is zeroed in the residual (`apply_frozen_cd4_residual`
#     asserts it), and
#   * a target missing any one of the three conditions receives NO residual
#     (`freeze_complete_mean3`).
# Without those two the arm is not buildable, so it is the one that should be quoted.
PSEUDO = ("cd4_mean3", "cd4_mean3_gn", "cd4_mean3_gn_prod")
FRACS = (0.0, 0.05, 0.10, 0.25, 0.50, 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--gamma", type=float, default=0.65)
    a = ap.parse_args()
    C = Ctx(a.harness)
    cands = [c for c in CAND if c not in C.banned]
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in C.banned)[:2]
    union = tuple(dict.fromkeys(base + tuple(cands)))
    m = C.set_donors(union)
    cov = C.coverage()
    print(f"# {a.harness}, POST-whitening add, gamma={a.gamma}, base={base}")
    print("# coverage: " + ", ".join(f"{n} {cov[n]}/{C.P}" for n in union))

    w0 = {n: (1.0 if n in base else 0.0) for n in union}
    # the whitened base mixture, BEFORE the final amplitude renormalisation
    rows = [m.effect(p, 1.0, 0.0, weights=dict(w0)) for p in C.perts]
    E0 = np.vstack(rows).astype(np.float64)
    if a.gamma > 0:
        before = float(np.median(np.linalg.norm(E0[:, C.wcols], axis=1)))
        E0[:, C.wcols] = C.W.apply(E0[:, C.wcols], a.gamma, 0.0, renorm=False)
        E0[:, C.wcols] *= before / (
            float(np.median(np.linalg.norm(E0[:, C.wcols], axis=1))) or 1.0)
    nE = float(np.median(np.linalg.norm(E0, axis=1)))

    def finish(E):
        E = E * (ALPHA * m.effect_scale / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0))
        return C.dh_exact(E)

    dh0 = finish(E0)
    X0, A0, W0 = C.XAW(dh0)
    k0 = C.kappa(dh0)
    print(f"\n# base: kappa={k0:.5f} pds={C.pds(dh0):.5f} nmae={C.nmae(dh0):.5f}")
    print(f"  {'candidate':<15}{'cov':>5}" + "".join(f"{'f='+str(f):>9}" for f in FRACS)
          + f"{'argmax':>8}")
    def raw_of(c):
        V = np.zeros_like(E0)
        for i, p in enumerate(C.perts):
            j = m.pert_index[c].get(p)
            if j is None or j >= m.s[c].shape[0]:
                continue
            V[i, m.tcols[c]] = m.s[c][j][m.dcols[c]]
        n = float(np.median(np.linalg.norm(V, axis=1))) or 1.0
        return V * (nE / n)

    cd4n = [c for c in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr") if c in cands]

    out = {}
    for c in cands + list(PSEUDO if len(cd4n) == 3 else ()):
        if c == "cd4_mean3":
            V = sum(raw_of(x) for x in cd4n) / 3.0
        elif c.startswith("cd4_mean3_gn"):
            V = sum(raw_of(x) for x in cd4n) / 3.0
            if c.endswith("_prod"):
                have = np.ones(len(C.perts), bool)
                for x in cd4n:
                    have &= np.array([(m.pert_index[x].get(p) is not None
                                       and m.pert_index[x][p] < m.s[x].shape[0])
                                      for p in C.perts])
                V[~have] = 0.0
                panel = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")
                            ["target_gene"].astype(str))
                pcols = np.array([i for i, g in enumerate(m.target_genes) if str(g) in panel])
                V[:, pcols] = 0.0
                print(f"    [{c}] {int(have.sum())}/{len(C.perts)} targets carry all three "
                      f"conditions; {pcols.size} panel coordinates zeroed")
            sc = np.abs(V).mean(0)
            live = sc[sc > 0]
            sc = np.maximum(sc, np.percentile(live, 25.0) if live.size else 1e-9)
            before = float(np.median(np.linalg.norm(V, axis=1)))
            V = V / sc
            V *= before / (float(np.median(np.linalg.norm(V, axis=1))) or 1.0)
        else:
            V = raw_of(c)
        nV = float(np.median(np.linalg.norm(V, axis=1))) or 1.0
        V *= nE / nV
        cov.setdefault(c, cov.get(cd4n[0], 0) if c.startswith("cd4_mean") else 0)
        ms, rows_ = [], {}
        for f in FRACS:
            dh = finish(E0 + f * V)
            X, A, W = C.XAW(dh)
            ms.append(float(np.sign(X) * np.sqrt(W / W0)))
            rows_[str(f)] = dict(X=X, A=A, W=W, M=ms[-1], kappa=C.kappa(dh),
                                 pds=C.pds(dh), nmae=C.nmae(dh))
        out[c] = rows_
        print(f"  {c:<15}{cov[c]:>5}" + "".join(f"{v:>9.3f}" for v in ms)
              + f"{FRACS[int(np.argmax(ms))]:>8}")
    p = ROOT / f"experiments/099-donors/work/post_{a.harness}_g{a.gamma}.json"
    json.dump(dict(harness=a.harness, gamma=a.gamma, base=list(base), kappa0=k0,
                   coverage=cov, arms=out), open(p, "w"), indent=1)
    print(f"-> {p}")


if __name__ == "__main__":
    main()
