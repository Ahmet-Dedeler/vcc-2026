"""159 -- consistency check on the reliability estimator, before anything is built on it.

The cosine split-half r is an estimator of rho = ||d||^2/(||d||^2+sigma^2) only under a
ratio-of-expectations approximation.  A second, moment-based estimator uses no such step:

    e1 = d + eps1,  e2 = d + eps2   (disjoint cells, disjoint control halves)
    E||e1 - e2||^2      = 2 sigma_m^2
    E||(e1+e2)/2||^2    = ||d||^2 + sigma_m^2 / 2
 => ||d||^2 = E||mean||^2 - E||diff||^2 / 4
 => rho(m)  = ||d||^2 / (||d||^2 + sigma_m^2) with sigma_m^2 = E||diff||^2 / 2
    rho(2m) = ||d||^2 / E||mean||^2          (the pooled 2m-cell estimate)

If the two disagree the disattenuation factor is wrong and every number downstream moves.
Also prints the energy budget (shared programme vs target-specific vs noise), which is what
the ceiling in SS5 is actually built from.
"""
from __future__ import annotations

import json
import os

import numpy as np

from reliability import load, prof, unit

SEED = 20260902
WORK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work")
N_DRAW = 12


def main() -> None:
    out = {}
    for name in ("hct116", "hek293t", "h1test"):
        p = os.path.join(WORK, f"shards-{name}.npz")
        if not os.path.exists(p):
            continue
        d = load(name)
        kt, kc, mask = d["kt"], d["kc"], d["mask"]
        rng = np.random.default_rng(SEED + 3)
        c1 = prof(d["csum"][: kc // 2].sum(0))[mask]
        c2 = prof(d["csum"][kc // 2:].sum(0))[mask]
        cf = prof(d["csum"].sum(0))[mask]
        j = kt // 2
        dif2 = men2 = 0.0
        cosr = 0.0
        for _ in range(N_DRAW):
            perm = rng.permutation(kt)
            e1 = prof(d["tsum"][:, perm[:j]].sum(1))[:, mask] - c1[None, :]
            e2 = prof(d["tsum"][:, perm[j:2 * j]].sum(1))[:, mask] - c2[None, :]
            dif2 += float((np.linalg.norm(e1 - e2, axis=1) ** 2).mean())
            men2 += float((np.linalg.norm(0.5 * (e1 + e2), axis=1) ** 2).mean())
            cosr += float((unit(e1) * unit(e2)).sum(1).mean())
        dif2, men2, cosr = dif2 / N_DRAW, men2 / N_DRAW, cosr / N_DRAW
        n_half = j * float(d["tcell"].mean())

        sig2_half = dif2 / 2.0                 # noise energy of one m-cell effect
        d2 = men2 - dif2 / 4.0                 # true effect energy
        rho_m = d2 / (d2 + sig2_half)          # reliability at m cells  (compare to cosr)
        rho_2m = d2 / men2                     # reliability at 2m = the surface's own depth

        E = prof(d["tsum"].sum(1))[:, mask] - cf[None, :]
        mline = E.mean(0)
        obs2 = float((np.linalg.norm(E, axis=1) ** 2).mean())
        m2 = float(mline @ mline)
        s2_obs = float((np.linalg.norm(E - mline[None, :], axis=1) ** 2).mean())
        sig2_full = sig2_half / 2.0            # noise at 2m cells
        s2_true = s2_obs - sig2_full           # target-specific energy, noise removed

        rec = dict(
            n_cells_per_side=n_half, n_targets=int(E.shape[0]), n_genes=int(len(mask)),
            cos_splithalf_at_m=cosr, moment_rho_at_m=rho_m,
            spearman_brown_from_cos=2 * cosr / (1 + cosr), moment_rho_at_2m=rho_2m,
            E_obs_sq=obs2, E_shared_sq=m2, E_resid_obs_sq=s2_obs,
            E_noise_sq_at_2m=sig2_full, E_resid_true_sq=s2_true,
            shared_frac_of_true=m2 / (m2 + s2_true),
            resid_frac_of_true=s2_true / (m2 + s2_true),
            true_frac_of_obs=(m2 + s2_true) / obs2,
            resid_reliability_at_2m=s2_true / s2_obs,
            implied_within_line_cos=m2 / (m2 + s2_true),
        )
        out[name] = rec
        print(f"\n=== {name}  ({rec['n_targets']} targets, {rec['n_genes']} genes, "
              f"{n_half:.0f}+{n_half:.0f} cells)")
        print(f"  reliability at {n_half:.0f} cells : cosine {cosr:.4f}   "
              f"moment {rho_m:.4f}   (ratio {rho_m / cosr:.3f})")
        print(f"  reliability at {2 * n_half:.0f} cells: Spearman-Brown(cos) "
              f"{rec['spearman_brown_from_cos']:.4f}   moment {rho_2m:.4f}   "
              f"(ratio {rho_2m / rec['spearman_brown_from_cos']:.3f})")
        print(f"  energy per target at {2 * n_half:.0f} cells:  observed {obs2:9.2f}"
              f" = shared {m2:8.2f} + target-specific {s2_true:8.2f} + noise {sig2_full:8.2f}")
        print(f"  -> shared / target-specific split of the TRUE effect: "
              f"{rec['shared_frac_of_true']:.3f} / {rec['resid_frac_of_true']:.3f}"
              f"   (implies within-line pairwise cos {rec['implied_within_line_cos']:.3f})")
        print(f"  -> reliability of the RESIDUAL alone at {2 * n_half:.0f} cells: "
              f"{rec['resid_reliability_at_2m']:.4f}")

    with open(os.path.join(WORK, "diag.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'diag.json')}")


if __name__ == "__main__":
    main()
