"""159 §3 -- cross-line transfer, measured AS A FUNCTION OF DEPTH, then disattenuated.

THE ESTIMATOR AND ITS ASSUMPTIONS
---------------------------------
Observed effect = truth + noise, d_obs = d + eps, with eps (i) zero-mean, (ii) independent of
d, (iii) independent between the two lines and between disjoint cell sets.  Then

    E cos(d_A + eps_A, d_B + eps_B) ~= cos(d_A, d_B) * sqrt(rho_A * rho_B),
    rho = ||d||^2 / (||d||^2 + sigma^2)   -- exactly the split-half r of SS2.

Two consequences that decide how this experiment reads:

1. The attenuation factor sqrt(rho_A rho_B) multiplies the SAME-target and the
   DIFFERENT-target cosine identically, because it lives entirely in the denominators and the
   noise does not know which pair it is in.  So the RATIO same/diff -- and therefore the
   verdict "is a target's own effect more informative than a random target's" -- is
   attenuation-INVARIANT to first order.  The brief's expectation that noise suppresses the
   excess more than the baseline is true in the RESIDUAL space (where the target-specific
   vector's own SNR is what attenuates) and false in the FULL-vector space.  Both are measured
   here so the two readings can be separated.
2. The invariance in (1) is testable, not assumed: sub-sampling the same cells to 15/30/45/60/
   75/150 cells per target sweeps sqrt(rho_A rho_B) over a factor of ~4.  If same/diff is flat
   over that sweep, the model holds and the donor-table matrix (SS4), where no cell-level split
   is possible, can be read through the ratio with no disattenuation at all.

Where the assumption is weakest: rho varies across targets and is CORRELATED across lines (a
strong perturbation is strong in both), so E[sqrt(rho_A,t rho_B,t)] on the diagonal exceeds
E[sqrt(rho_A,t)]E[sqrt(rho_B,u)] off it.  That biases the observed excess UP, not down.  The
per-target disattenuation variant below removes it; the aggregate variant does not.
"""
from __future__ import annotations

import json
import os

import numpy as np

from reliability import KEYS, fit, load, prof, unit  # noqa: F401

SEED = 20260902
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
N_DRAW = 12
N_BOOT = 2000


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def effects(d, idx, shard_sel, ctrl_sel, mask):
    """Pseudobulk effect for targets `idx` from shards `shard_sel`, vs control `ctrl_sel`."""
    b = prof(d["tsum"][idx][:, shard_sel].sum(1))[:, mask]
    c = prof(d["csum"][ctrl_sel].sum(0))[mask]
    return b - c[None, :]


def stats(EA, EB):
    """same / diff / excess in the full space and in the per-line mean-removed space."""
    out = {}
    for space in ("full", "resid"):
        A = EA - EA.mean(0, keepdims=True) if space == "resid" else EA
        B = EB - EB.mean(0, keepdims=True) if space == "resid" else EB
        a, b = unit(A), unit(B)
        S = a @ b.T
        same = np.diag(S).copy()
        diff = offdiag_mean(S)
        out[space] = {"same": float(same.mean()), "same_median": float(np.median(same)),
                      "same_sd": float(same.std(ddof=1)), "diff": diff,
                      "excess": float(same.mean() - diff),
                      "ratio": float(same.mean() / diff) if diff != 0 else np.nan,
                      "within_A": offdiag_mean(a @ a.T), "within_B": offdiag_mean(b @ b.T),
                      "per_target_same": same}
    return out


def splithalf(d, idx, j, mask, rng, n_draw, kc):
    """r(m = j*shard) per target, disjoint cell halves AND disjoint control halves."""
    kt = d["kt"]
    c1, c2 = list(range(kc // 2)), list(range(kc // 2, kc))
    acc = {"full": np.zeros(len(idx)), "resid": np.zeros(len(idx))}
    for _ in range(n_draw):
        perm = rng.permutation(kt)
        s1, s2 = perm[:j].tolist(), perm[j:2 * j].tolist()
        e1 = effects(d, idx, s1, c1, mask)
        e2 = effects(d, idx, s2, c2, mask)
        acc["full"] += (unit(e1) * unit(e2)).sum(1)
        r1, r2 = e1 - e1.mean(0, keepdims=True), e2 - e2.mean(0, keepdims=True)
        acc["resid"] += (unit(r1) * unit(r2)).sum(1)
    return {k: v / n_draw for k, v in acc.items()}


def main() -> None:
    rng = np.random.default_rng(SEED + 7)
    A, B = load("hct116"), load("hek293t")
    assert np.array_equal(A["genes"], B["genes"]), "gene axes differ"
    shared = np.asarray(sorted(set(A["targets"]) & set(B["targets"])))
    ia = np.asarray([int(np.flatnonzero(A["targets"] == s)[0]) for s in shared])
    ib = np.asarray([int(np.flatnonzero(B["targets"] == s)[0]) for s in shared])
    gpos = {g: i for i, g in enumerate(A["genes"])}
    drop = {gpos[s] for s in shared if s in gpos}
    mask = np.asarray([i for i in range(len(A["genes"])) if i not in drop])
    kt, kc = A["kt"], A["kc"]
    cps = float(A["tcell"].mean())
    print(f"shared targets {len(shared)} | scored genes {len(mask)} | "
          f"{kt} shards of {cps:.0f} cells | control {A['ccell'].sum():.0f} / "
          f"{B['ccell'].sum():.0f} cells")

    out = {"n_shared": int(len(shared)), "n_scored_genes": int(len(mask)), "seed": SEED,
           "cells_per_shard": cps, "n_draw": N_DRAW, "depth": {}}

    # ---- depth sweep: does the same/diff RATIO survive a 10x change in attenuation? -----
    print(f"\n{'m/target':>9} {'space':>6} {'same':>8} {'diff':>8} {'excess':>9} "
          f"{'ratio':>7} {'r_A':>7} {'r_B':>7} {'sqrt(rArB)':>10} {'exc/atten':>10}")
    jlist = list(range(1, kt + 1))
    for j in jlist:
        m = j * cps
        # cross-line cosines: independent shard draws per line, averaged
        accum = None
        for _ in range(N_DRAW if j < kt else 1):
            sa = rng.permutation(kt)[:j].tolist() if j < kt else list(range(kt))
            sb = rng.permutation(kt)[:j].tolist() if j < kt else list(range(kt))
            EA = effects(A, ia, sa, list(range(kc)), mask)
            EB = effects(B, ib, sb, list(range(kc)), mask)
            s = stats(EA, EB)
            if accum is None:
                accum = {sp: {k: (v.copy() if isinstance(v, np.ndarray) else v)
                              for k, v in s[sp].items()} for sp in s}
            else:
                for sp in s:
                    for k, v in s[sp].items():
                        accum[sp][k] = accum[sp][k] + v
        n = N_DRAW if j < kt else 1
        for sp in accum:
            for k in accum[sp]:
                accum[sp][k] = accum[sp][k] / n

        # reliability at the SAME depth (needs 2 disjoint sets, so only j <= kt/2 is direct)
        if j <= kt // 2:
            ra = splithalf(A, ia, j, mask, np.random.default_rng(SEED + 100 + j), N_DRAW, kc)
            rb = splithalf(B, ib, j, mask, np.random.default_rng(SEED + 200 + j), N_DRAW, kc)
        else:
            ra = rb = None

        rec = {"m": m}
        for sp in ("full", "resid"):
            v = accum[sp]
            rec[sp] = {k: (float(x) if not isinstance(x, np.ndarray) else x.tolist())
                       for k, x in v.items()}
            if ra is not None:
                rA, rB = float(ra[sp].mean()), float(rb[sp].mean())
                at = float(np.sqrt(max(rA, 1e-9) * max(rB, 1e-9)))
                rec[sp].update(r_A=rA, r_B=rB, atten=at,
                               excess_disatten=v["excess"] / at,
                               same_disatten=v["same"] / at, diff_disatten=v["diff"] / at)
                rec[sp]["per_target_rA"] = ra[sp].tolist()
                rec[sp]["per_target_rB"] = rb[sp].tolist()
                print(f"{m:9.0f} {sp:>6} {v['same']:8.4f} {v['diff']:8.4f} "
                      f"{v['excess']:+9.5f} {v['ratio']:7.4f} {rA:7.4f} {rB:7.4f} "
                      f"{at:10.4f} {v['excess'] / at:+10.5f}")
            else:
                print(f"{m:9.0f} {sp:>6} {v['same']:8.4f} {v['diff']:8.4f} "
                      f"{v['excess']:+9.5f} {v['ratio']:7.4f} {'':7} {'':7} {'':10} {'':10}")
        out["depth"][f"{m:.0f}"] = rec

    # ---- reliability at the surfaces' own depth (150) from the fitted law ---------------
    ms = [j * cps for j in range(1, kt // 2 + 1)]
    rho150 = {}
    print(f"\nreliability at the surfaces' own depth (150 cells), by fitted law:")
    for lbl, d, idx, off in (("A_hct116", A, ia, 100), ("B_hek293t", B, ib, 200)):
        rho150[lbl] = {}
        for sp in ("full", "resid"):
            rs = [float(splithalf(d, idx, j, mask, np.random.default_rng(SEED + off + j),
                                  N_DRAW, kc)[sp].mean()) for j in range(1, kt // 2 + 1)]
            row = {"r_curve": rs}
            for model in ("linear", "power"):
                par, r2, f = fit(ms, rs, model)
                row[model] = {**par, "r2": r2, "rho_150": float(f(150)),
                              "rho_400": float(f(400))}
            # Spearman-Brown from the directly measured r(75): 150 cells is exactly 2x 75
            row["spearman_brown_150"] = 2 * rs[-1] / (1 + rs[-1])
            rho150[lbl][sp] = row
            print(f"  {lbl:10} {sp:>6}  r(75)={rs[-1]:.4f}  SB(150)={row['spearman_brown_150']:.4f}"
                  f"  linear(150)={row['linear']['rho_150']:.4f}"
                  f"  power(150)={row['power']['rho_150']:.4f}"
                  f"  | power(400)={row['power']['rho_400']:.4f}"
                  f"  linear(400)={row['linear']['rho_400']:.4f}")
    out["rho150"] = rho150

    # ---- full-depth disattenuation + target-cluster bootstrap ---------------------------
    EA = effects(A, ia, list(range(kt)), list(range(kc)), mask)
    EB = effects(B, ib, list(range(kt)), list(range(kc)), mask)
    full = stats(EA, EB)
    rA_t = {sp: np.asarray(rho150["A_hct116"][sp]["r_curve"]) for sp in ("full", "resid")}

    # per-target r at 150 cells, Spearman-Brown from the per-target r(75)
    pt = {}
    for lbl, d, idx, off in (("A", A, ia, 100), ("B", B, ib, 200)):
        pt[lbl] = splithalf(d, idx, kt // 2, mask,
                            np.random.default_rng(SEED + off + kt // 2), N_DRAW, kc)
    sb = lambda r: 2 * r / (1 + r)                                          # noqa: E731

    print(f"\n{'=' * 96}\nFULL DEPTH (150 cells/target), disattenuated, "
          f"target-cluster bootstrap {N_BOOT} draws")
    boot = {}
    for sp in ("full", "resid"):
        same_t = np.asarray(full[sp]["per_target_same"])
        a, b = unit(EA - EA.mean(0, keepdims=True) if sp == "resid" else EA), None
        b = unit(EB - EB.mean(0, keepdims=True) if sp == "resid" else EB)
        S = a @ b.T
        rA = np.clip(sb(np.asarray(pt["A"][sp])), 1e-4, 0.999)
        rB = np.clip(sb(np.asarray(pt["B"][sp])), 1e-4, 0.999)
        agg_at = float(np.sqrt(rA.mean() * rB.mean()))
        est = {
            "same_obs": float(same_t.mean()), "diff_obs": offdiag_mean(S),
            "excess_obs": float(same_t.mean() - offdiag_mean(S)),
            "ratio_obs": float(same_t.mean() / offdiag_mean(S)),
            "rA_mean": float(rA.mean()), "rB_mean": float(rB.mean()), "atten_agg": agg_at,
            "same_dis": float(same_t.mean()) / agg_at,
            "diff_dis": offdiag_mean(S) / agg_at,
            "excess_dis": float(same_t.mean() - offdiag_mean(S)) / agg_at,
        }
        # per-target disattenuation (removes the rho-correlation bias described in the header)
        pt_at = np.sqrt(rA * rB)
        est["same_dis_pertarget"] = float(np.mean(same_t / pt_at))
        n = len(shared)
        bs = {k: [] for k in ("same_dis", "diff_dis", "excess_dis", "ratio_obs", "excess_obs")}
        for _ in range(N_BOOT):
            k = rng.integers(0, n, n)
            sm = float(np.mean(np.diag(S)[k]))
            Sk = S[np.ix_(k, k)]
            df = float(Sk[~np.eye(len(k), dtype=bool)].mean())
            at = float(np.sqrt(rA[k].mean() * rB[k].mean()))
            bs["same_dis"].append(sm / at)
            bs["diff_dis"].append(df / at)
            bs["excess_dis"].append((sm - df) / at)
            bs["excess_obs"].append(sm - df)
            bs["ratio_obs"].append(sm / df)
        for k, v in bs.items():
            v = np.asarray(v)
            est[f"{k}_ci"] = [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
        boot[sp] = est
        print(f"\n  [{sp}]  attenuation sqrt(rA*rB) = {agg_at:.4f}  "
              f"(rA {rA.mean():.4f}, rB {rB.mean():.4f})")
        print(f"    same    obs {est['same_obs']:+.4f} -> disatten {est['same_dis']:+.4f} "
              f"CI [{est['same_dis_ci'][0]:+.4f}, {est['same_dis_ci'][1]:+.4f}]")
        print(f"    diff    obs {est['diff_obs']:+.4f} -> disatten {est['diff_dis']:+.4f} "
              f"CI [{est['diff_dis_ci'][0]:+.4f}, {est['diff_dis_ci'][1]:+.4f}]")
        print(f"    EXCESS  obs {est['excess_obs']:+.5f} "
              f"CI [{est['excess_obs_ci'][0]:+.5f}, {est['excess_obs_ci'][1]:+.5f}] "
              f"-> disatten {est['excess_dis']:+.5f} "
              f"CI [{est['excess_dis_ci'][0]:+.5f}, {est['excess_dis_ci'][1]:+.5f}]")
        print(f"    ratio same/diff (attenuation-invariant) {est['ratio_obs']:.4f} "
              f"CI [{est['ratio_obs_ci'][0]:.4f}, {est['ratio_obs_ci'][1]:.4f}]")
    out["full_depth"] = {sp: {k: v for k, v in boot[sp].items()} for sp in boot}
    out["full_depth_raw"] = {sp: {k: v for k, v in full[sp].items()
                                  if not isinstance(v, np.ndarray)} for sp in full}

    with open(os.path.join(WORK, "crossline.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'crossline.json')}")


if __name__ == "__main__":
    main()
