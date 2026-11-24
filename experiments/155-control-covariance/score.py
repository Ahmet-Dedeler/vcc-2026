"""155 — score the control-covariance candidate against the incumbent, on one fixed frame.

Endpoints (all from 153's `debench.Surface`, no redefinition):
  sign acc      per-gene sign accuracy over the reference's DE pool
  head purity   sign accuracy in the top decile by |predicted effect| -- what `reach` eats
  reach         de_wilcoxon_direction_reach_raw, keyed on |predicted effect|
  pds           pds_cosine
  nmae          de_wilcoxon_lfc_nmae (lower is better)

NULLS.  Sign accuracy has a floor that is NOT 0.5: if the reference's DE pool is
sign-imbalanced, a constant all-down vector scores the imbalance for free and is worth
nothing.  So every candidate is reported next to
  * `null-alldown`  -- one constant vector, no target specificity at all;
  * `null-permuted` -- the candidate's own rows dealt to the wrong targets.
A candidate only carries information if it beats BOTH.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench import HDR, ROW, load_frame, score  # noqa: E402

s, inc, gdst, pok = load_frame()
GM = np.zeros(s.genes.size, bool)
GM[gdst] = True


def onto_frame(M: np.ndarray) -> np.ndarray:
    """A (132 x 10261) coefficient matrix restricted to the scored columns."""
    P = np.zeros_like(s.pred_lfc)
    P[:, gdst] = M[:, gdst]
    return P


def report(name, P, key=None):
    r = score(s, P, key=key)
    print(ROW.format(name=name, **r))
    return r


def main() -> None:
    print(f"frame: {pok.sum()} perturbations x {gdst.size} genes")
    down = np.array([(s.real_lfc[p][s.real_sig[p] & s.adj[p]] < 0).mean()
                     for p in range(len(s.perts))])
    pool = np.concatenate([s.real_lfc[p][s.real_sig[p] & s.adj[p]] for p in range(len(s.perts))])
    print(f"reference DE pool: {pool.size} gene-calls, {100*(pool<0).mean():.1f}% down "
          f"(per-perturbation mean {100*down.mean():.1f}%)")
    print()
    print(HDR)
    base = report("incumbent (104 baseline, g=0.65)", inc)

    # --- nulls -------------------------------------------------------------------------
    rng = np.random.default_rng(0)
    alldown = onto_frame(np.tile(-np.ones(s.genes.size), (len(s.perts), 1)))
    alldown[:, gdst] *= rng.random(gdst.size)          # break ties so the decile is defined
    report("null: all-down constant", alldown)

    # --- candidates --------------------------------------------------------------------
    cands = {}
    for k in (0, 20, 50):
        z = np.load(os.path.join(HERE, "work", f"coef_k{k}.npz"), allow_pickle=True)
        assert [str(g) for g in z["genes"]] == list(s.genes)
        assert [str(t) for t in z["targets"]] == list(s.perts)
        for tag, M in (("beta", z["beta"]), ("corr", z["corr"])):
            P = onto_frame(-M.astype(np.float64))       # knockdown lowers T -> effect = -beta
            cands[f"{tag} k={k}"] = P
            report(f"control-cov -{tag}, k={k} PCs", P)

    print()
    print("permutation null (candidate rows dealt to the wrong targets, 5 seeds):")
    for name in ("beta k=20", "corr k=20"):
        P = cands[name]
        accs, hps = [], []
        for seed in range(5):
            perm = np.random.default_rng(seed).permutation(len(s.perts))
            while (perm == np.arange(len(s.perts))).any():
                perm = np.random.default_rng(seed + 1000).permutation(len(s.perts))
            r = score(s, P[perm])
            accs.append(r["sign_acc"])
            hps.append(r["head_purity"])
        print(f"  {name:14} sign acc {np.mean(accs):.4f} +- {np.std(accs):.4f}   "
              f"head purity {np.mean(hps):.4f} +- {np.std(hps):.4f}")

    return base


if __name__ == "__main__":
    main()
