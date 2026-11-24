"""160b — the payoff test, and the confound test.

`subspace.py` established that a line's own control basis captures 8-11x more of its
target-specific response residual than another line's basis, and 150-340x more than an
envelope-matched random subspace. Two things must follow before that means anything.

## 1. THE CONFOUND (run first, because it can kill the result)

"Line-specific" and "same-experiment-specific" are not the same claim. The control cells and the
perturbation cells of one harness come from the same experiment, so they share batch, depth and
pool structure. If the captured residual is technical rather than biological, the other line's
basis would fail for the trivial reason that it is a different experiment.

Controls applied here:
  * `resid_lib`  — regress log library size out of the control cells before the basis is built.
  * `resid_lib_mito` — also regress out the mitochondrial fraction.
  * `shuffle_targets` — RETRACTED, it is vacuous. `||R V^T||^2 / ||R||^2` is a sum of
    independent per-row terms, so it is invariant to any permutation of R's rows BY
    CONSTRUCTION. It returned exactly `own` to four decimals on both lines, which is a theorem
    and not evidence. Kept in the output as a reminder, per AGENTS.md 11.1 rule 10's corollary.
    The control that can actually fail is `discrim` below.
  * `discrim` — the real target-specificity test. Project each residual row into the subspace and
    measure the mean row-row cosine of the projections. A subspace that only carries the shared
    programme makes every projection point the same way (cosine -> 1) and cannot distinguish
    targets, however much energy it captures. Compared against the row-row cosine of the
    residual itself.
  * `split_half_cells` — basis from one disjoint half of the control cells, capture measured on
    the residual as usual; and basis from half A tested against basis from half B, which bounds
    how much of the basis is reproducible at all.

## 2. THE PAYOFF

Capture of the TRUTH in a subspace is necessary, not sufficient: our prediction must carry
signal INSIDE that subspace. Tested on the a23 delivered effect, pushed through the production
renderer with the mass fix (experiment 156 — scoring the pre-render effect measures a different
object), against exact-space truth:

  project   P -> P V^T V                     keep only the context subspace
  blend     P -> (1-w) P + w P V^T V         amplify it without discarding the rest
  amplify   P -> P + g P V^T V               same family, unnormalised

and always against the same operation using the OTHER line's basis, so a gain has to be
context-specific to count.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "156-delivered-geometry"))
from massfix import mass_fix, production_effect                      # noqa: E402
from screen import endpoints, render                                 # noqa: E402
from subspace import (KS, REPO, SEED, SURFACES, WORK, captured,      # noqa: E402
                      control_cells, random_basis, top_basis)

os.makedirs(WORK, exist_ok=True)


def resid_basis(X, k, covars, seed=SEED):
    """Control basis after projecting explicit technical covariates out of the cells."""
    Z = np.column_stack([np.ones(len(X))] + covars)
    Q, _ = np.linalg.qr(Z)
    Xr = X - Q @ (Q.T @ X)
    return top_basis(Xr, k, seed=seed)


def main() -> None:
    rng = np.random.default_rng(SEED)
    store, out = {}, {}

    for name, (h5, truth) in SURFACES.items():
        Xfull = control_cells(os.path.join(REPO, h5))
        z = np.load(os.path.join(REPO, truth), allow_pickle=True)
        genes = np.asarray([str(g) for g in z["genes"]])
        tg = np.asarray([str(t) for t in z["targets"]])
        gpos = {g: i for i, g in enumerate(genes)}
        drop = {gpos[t] for t in tg if t in gpos}
        mask = np.asarray([i for i in range(len(genes)) if i not in drop])
        E = z["effect"].astype(np.float64)[:, mask]
        X = Xfull[:, mask].astype(np.float64)
        lib = np.log(np.maximum(np.expm1(Xfull).sum(1), 1.0))
        mito = np.asarray([i for i, g in enumerate(genes) if g.startswith("MT-")])
        mfrac = (np.expm1(Xfull[:, mito]).sum(1)
                 / np.maximum(np.expm1(Xfull).sum(1), 1.0)) if len(mito) else np.zeros(len(Xfull))
        store[name] = dict(X=X, E=E, mask=mask, targets=tg,
                           lib=lib, mfrac=mfrac,
                           sd=X.std(0), nmito=len(mito))
        print(f"[{name}] {X.shape[0]} cells x {X.shape[1]} scored genes, "
              f"{E.shape[0]} targets, {len(mito)} MT- genes")

    K = 50
    print(f"\n{'=' * 100}\n1. CONFOUND TESTS — captured fraction of the residual at k={K}")
    print(f"  {'line':<9} {'own':>8} {'-lib':>8} {'-lib-mito':>10} {'halfA':>8} "
          f"{'halfB':>8} {'A~B':>7} {'shuffled':>9} {'other':>8} {'random':>8}")
    for name, d in store.items():
        other = [o for o in store if o != name][0]
        R = d["E"] - d["E"].mean(0, keepdims=True)
        V = top_basis(d["X"], K)
        Vl = resid_basis(d["X"], K, [d["lib"][:, None]])
        Vlm = resid_basis(d["X"], K, [d["lib"][:, None], d["mfrac"][:, None]])
        h = len(d["X"]) // 2
        idx = rng.permutation(len(d["X"]))
        Va = top_basis(d["X"][idx[:h]], K)
        Vb = top_basis(d["X"][idx[h:]], K)
        ab = float(np.linalg.norm(Va @ Vb.T, "fro") ** 2 / K)   # subspace overlap in [0,1]
        Rs = R[rng.permutation(R.shape[0])]                     # vacuous: see the docstring
        Vo = top_basis(store[other]["X"], K)
        d["V"], d["Vother"] = V, Vo

        def rr(M):
            U = M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-300)
            G = U @ U.T
            return float(G[~np.eye(len(G), dtype=bool)].mean())

        row = dict(own=captured(R, V), lib=captured(R, Vl), libmito=captured(R, Vlm),
                   halfA=captured(R, Va), halfB=captured(R, Vb), ab_overlap=ab,
                   shuffled=captured(Rs, V), other=captured(R, Vo),
                   random=float(np.mean([captured(R, random_basis(d["sd"], K, rng))
                                         for _ in range(3)])),
                   rr_resid=rr(R), rr_proj_own=rr((R @ V.T) @ V),
                   rr_proj_other=rr((R @ Vo.T) @ Vo))
        out[f"confound_{name}"] = row
        print(f"  {name:<9} {row['own']:8.4f} {row['lib']:8.4f} {row['libmito']:10.4f} "
              f"{row['halfA']:8.4f} {row['halfB']:8.4f} {ab:7.3f} {row['shuffled']:9.4f} "
              f"{row['other']:8.4f} {row['random']:8.4f}")
        print(f"  {'':<9} DISCRIM row-row cosine: residual {row['rr_resid']:+.4f} -> "
              f"projected(own) {row['rr_proj_own']:+.4f}, projected(other) "
              f"{row['rr_proj_other']:+.4f}   (-> 1.0 means non-discriminative)")
    print("  NOTE `shuffled` deals the residual rows to the WRONG targets. If it matches `own`,")
    print("       the captured subspace is target-INDEPENDENT and cannot carry target-specific")
    print("       information, however context-specific it is.")

    print(f"\n{'=' * 100}\n2. PAYOFF — does projecting our PREDICTION onto the context basis help?")
    for name in store:
        d = store[name]
        other = [o for o in store if o != name][0]
        D, T, b, b_ref, gate, mask_p, ptar = production_effect(name)
        assert np.array_equal(mask_p, d["mask"]), "gene masks differ between scripts"
        base_eff = mass_fix(D, b, gate)
        base = endpoints(render(base_eff, b)[:, d["mask"]], T[:, d["mask"]])
        print(f"\n  --- {name}: baseline (a23 delivered) pds {base['pds']:.5f} "
              f"cos {base['cos']:.5f} sse* {base['sse']:.6f}")
        print(f"  {'arm':<30} {'basis':<7} {'pds':>9} {'d pds':>9} {'cos':>9} "
              f"{'d cos %':>9} {'sse*':>10}")
        rows = []
        Dm = D[:, d["mask"]]
        for bname, V in (("own", d["V"]), ("other", d["Vother"])):
            for k in (5, 20, 50, 200):
                Vk = V[:k]
                proj = (Dm @ Vk.T) @ Vk
                for label, Pm in (
                        (f"project k{k}", proj),
                        (f"blend w0.5 k{k}", 0.5 * Dm + 0.5 * proj),
                        (f"amplify g1.0 k{k}", Dm + 1.0 * proj)):
                    P = D.copy()
                    P[:, d["mask"]] = Pm
                    Pg = np.where(gate[None, :], P, 0.0)
                    e = endpoints(render(mass_fix(Pg, b, gate), b)[:, d["mask"]],
                                  T[:, d["mask"]])
                    rows.append(dict(arm=label, basis=bname, **e,
                                     d_pds=e["pds"] - base["pds"]))
                    print(f"  {label:<30} {bname:<7} {e['pds']:9.5f} "
                          f"{e['pds'] - base['pds']:+9.5f} {e['cos']:9.5f} "
                          f"{100 * (e['cos'] / base['cos'] - 1):+8.1f}% {e['sse']:10.6f}")
        out[f"payoff_{name}"] = dict(base=base, rows=rows)

    with open(os.path.join(WORK, "transfer.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'transfer.json')}")


if __name__ == "__main__":
    sys.exit(main())
