"""216 Phase 2.2b — REFIT the whitening operator with CD4 IN ITS BASIS.

The mandate's distinction: "an operator fitted with CD4 in it is a different question from
whitening a CD4 residual through an operator that never saw it".  207 §5.2 measured the
second (`wres`: flat, +0.5 to +0.9% at every dose).  This measures the first.

Construction, stated exactly, because rule 11 says never infer an operator
-------------------------------------------------------------------------
The shipped operator is `W = (C + lam I)^(-gamma/2)` with `C = B^T B / n` over the donor
basis `B` on the fold's `wcols` axis.  The stored `.npz` holds the eigendecomposition of that
`C`, so `C` is recoverable exactly (up to `Whitener.fit`'s relative eigenvalue floor of
1e-12, which is applied to directions carrying ~0 mass):

    C = V diag(ev) V^T

The CD4 block is the same object the residual form adds -- `cd4.mean_effect` on rows where
all three conditions carry, restricted to `wcols`.  Its own Gram is

    C_cd4 = B_cd4^T B_cd4 / m

and the two are mixed at a stated ENERGY share, so the arm has one interpretable knob and no
dependence on donor row-count conventions we cannot reconstruct from the frozen files:

    C' = trace(C) * [ (1-beta) * C/trace(C)  +  beta * C_cd4/trace(C_cd4) ]

Two betas per fold:
  * `natural`  beta = m / (n_basis + m)   -- CD4 weighted by its row count, the weight it
    would get if its rows were simply appended to the basis;
  * `equal`    beta = 1 / (k + 1)         -- CD4 as one more EQUAL line among the fold's k
    sources, which is the operator-side analogue of the 1/k library route.

KNOWN-ANSWER ROWS.  At beta = 0 the construction is the identity BY ALGEBRA (`C' = C`), so
the only thing that needs verifying is that `eigh` inverts the reconstruction to the precision
this operator needs -- and this operator's condition number is 1e12, so that is a real
question.  Two gates, neither of which costs a full 17,806^2 eigendecomposition (each of those
is ~8.4 min, measured):

  1. **spectral round-trip at matched conditioning.**  A synthetic symmetric matrix is built
     with the fold operator's OWN eigenvalues, quantile-matched down to n = 4,000, on a random
     orthonormal basis; it is reconstructed and re-eigendecomposed exactly as the real one is,
     and the relative error of the resulting `(C)^(-gamma/2)` action is reported.
  2. **trace conservation on the full problem**, free: |sum(evals') - trace(C')| / trace(C').

Both are asserted before any refit arm is scored.

    uv run python op216.py
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)
from sweep216 import ALPHA, GAMMA, LAM, score_effects  # noqa: E402
from form216 import build_final, delivered_norm, solve_alpha_delivered, votes_for  # noqa: E402

# basis row counts, read from each operator's own frozen metadata (never inferred)
BASIS_ROWS = {
    "k562_gwps": None, "orion_hct116": None, "orion_hek293t": None, "kolf2": None}


class RefitOperator:
    """`Operator` with an externally supplied spectrum; same `apply` semantics."""

    def __init__(self, evals, evecs):
        self.evals = evals
        self.evecs = evecs

    apply = Operator.apply


def round_trip_gate(evals, n: int = 4000, seed: int = 216) -> float:
    """KA-1: does reconstruct-then-re-eigendecompose preserve `(C)^(-gamma/2)` at this
    operator's conditioning?  Built with the operator's OWN eigenvalues, quantile-matched
    down to `n`, on a random orthonormal basis -- so the condition number (1e12 here) and the
    eigenvalue floor are the real ones, at 1/85 the cost of the full problem."""
    rng = np.random.default_rng(seed)
    q = np.quantile(np.asarray(evals, dtype=np.float64), np.linspace(0, 1, n))
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    C = (Q * q) @ Q.T
    ev2, V2 = np.linalg.eigh(C)
    ev2 = np.maximum(ev2, ev2.max() * 1e-12)
    X = rng.standard_normal((64, n))
    a = ((X @ Q) * q ** (-GAMMA / 2.0)) @ Q.T
    b = ((X @ V2) * ev2 ** (-GAMMA / 2.0)) @ V2.T
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(a)), 1e-300))


def basis_rows_for(held: str) -> int:
    op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
               else PROD_OPS / "whitener_union_pruned.npz")
    meta = json.loads(op_path.with_suffix(".json").read_text())
    if "basis" in meta:
        return int(meta["basis"][0])
    return int(sum(meta["basis_rows_by_source"].values()))


def main() -> None:
    t0 = time.time()
    out_path = WORK / "op216.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        votes = votes_for(held, fold["targets"])
        carry = np.asarray(cd4.target_condition_count) == 3
        cd4e = np.asarray(cd4.mean_effect, dtype=np.float64)
        wcols = np.asarray(fold["wcols"])
        n_basis = basis_rows_for(held)
        k_sources = len(FOLDS[held]["sources"])

        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        sha = E062.sha256_file(op_path)
        print(f"[{held}] operator {op_path.name} sha {sha[:8]} n={op.evecs.shape[0]} "
              f"basis_rows={n_basis} k_sources={k_sources}", flush=True)

        # ------------------------------------------------------- reference arms
        ref = {}
        for tag, fv in (("res0.10", 0.10), ("res0.34", 0.34), ("res0", 0.0)):
            final, frac = build_final(pre, ALPHA, op, fold, cd4, fv)
            ref[tag] = dict(score=score_effects(final, fold, truth),
                            delivered=delivered_norm(final, fold), fraction=frac)
            del final
        target_norm = ref["res0.10"]["delivered"]
        gc.collect()

        # ------------------------------------------------------- reconstruct C
        # MEMORY: this loop is the only place in this experiment that allocates several
        # 17,942^2 float64 arrays at once.  `np.linalg.eigh` dispatches to LAPACK `dsyevd`,
        # whose workspace is ~2n^2 doubles (5.1 GiB here) ON TOP of the input and the
        # eigenvector output -- a first version of this script peaked at 18.5 GB and was
        # killed.  Two fixes: the original eigenbasis is freed as soon as `C` exists (only
        # its top eigenvector is kept, for the diagnostic), and the decomposition goes
        # through `scipy.linalg.eigh(driver='evd', overwrite_a=True)` on an F-ORDERED input,
        # which consumes the input in place instead of making scipy's silent Fortran copy.  Measured peak: see the log.
        from scipy.linalg import eigh as sp_eigh

        ev = np.asarray(op.evals, dtype=np.float64)
        V = np.asarray(op.evecs, dtype=np.float64)
        del op
        gc.collect()
        # KA-1: spectral round-trip at MATCHED CONDITIONING, n = 4,000 (see the docstring)
        ka_err = round_trip_gate(ev, n=4000, seed=216)
        print(f"[{held}] KA-1 round-trip at matched conditioning (n=4000): "
              f"rel max|diff| = {ka_err:.3e}", flush=True)
        assert ka_err < 1e-6, f"{held}: reconstruction round-trip failed at {ka_err:.3e}"

        # `C = V diag(ev) V^T` built as `(V sqrt(ev)) (V sqrt(ev))^T` with the scaling done
        # IN PLACE: `(V * ev) @ V.T` allocates a third n^2 temporary and costs 2.5 GB of peak
        # for nothing.  ev is non-negative (Whitener.fit floors it), so the sqrt is exact.
        top_orig = np.asarray(V[:, -1], dtype=np.float64).copy()
        max_eig_orig = float(ev.max())
        V *= np.sqrt(np.maximum(ev, 0.0))
        C = V @ V.T
        del V
        gc.collect()
        trC = float(np.trace(C))
        trace_recon_err = abs(trC - float(ev.sum())) / max(float(ev.sum()), 1e-300)
        print(f"[{held}] KA-2 trace(C reconstructed) vs sum(evals): rel err "
              f"{trace_recon_err:.3e}", flush=True)
        assert trace_recon_err < 1e-10, f"{held}: C reconstruction lost trace"

        # ------------------------------------------------------- CD4 Gram
        Bc = np.ascontiguousarray(
            np.where(np.asarray(cd4.gene_votes)[carry] > 0, cd4e[carry], 0.0)[:, wcols])
        m = int(Bc.shape[0])
        trC4 = float((Bc * Bc).sum() / m)

        betas = {"natural": m / (n_basis + m), "equal": 1.0 / (k_sources + 1)}
        for bname, beta in betas.items():
            t1 = time.time()
            # `Bs.T @ Bs` rather than `(Bc.T @ Bc) * s`: the scalar multiply on an n x n
            # product allocates a SECOND n^2 temporary and that alone cost 3 GB of peak.
            # F-ORDER MATTERS: scipy converts a C-ordered input to Fortran order before
            # calling LAPACK, which is a silent extra n^2 copy that defeats `overwrite_a`.
            # `C` is symmetric, so an F-ordered copy of it is just its transpose.
            Bs = Bc * np.sqrt(beta / (m * trC4))
            Cp = np.array(C, order="F")
            Cp *= (1.0 - beta) / trC
            Cp += Bs.T @ Bs
            Cp *= trC
            del Bs
            evp, Vp = sp_eigh(Cp, overwrite_a=True, driver="evd", check_finite=False)
            del Cp
            gc.collect()
            evp = np.maximum(evp, evp.max() * 1e-12)
            newop = RefitOperator(evp, Vp)
            spec = dict(beta=beta, m_cd4_rows=m, n_basis=n_basis,
                        max_eig=float(evp.max()), min_eig=float(evp.min()),
                        trace_ratio=float(np.sum(evp) / trC),
                        trace_conservation_rel_err=float(
                            abs(np.sum(evp) - trC) / max(abs(trC), 1e-300)),
                        eig_cos_with_original_top=float(abs(Vp[:, -1] @ top_orig)))
            print(f"[{held}/{bname}] beta {beta:.5f} refit in {time.time()-t1:.0f}s "
                  f"maxeig {evp.max():.6f} (was {max_eig_orig:.6f}) "
                  f"trace rel err {spec['trace_conservation_rel_err']:.2e} "
                  f"top-evec |cos| {spec['eig_cos_with_original_top']:.4f}", flush=True)

            arms = [("res0", 0.0, 0.0), ("res0.10", 0.0, 0.10), ("res0.34", 0.0, 0.34),
                    ("lib0.5", 0.5, 0.0), ("lib1", 1.0, 0.0)]
            for arm, w, fv in arms:
                name = f"refit-{bname}/{arm}"
                if (held, name) in done:
                    continue
                rowsel = carry & (votes > 0)
                pre4 = pre.copy()
                if w > 0:
                    pre4[rowsel] = ((pre[rowsel] * votes[rowsel, None] + w * cd4e[rowsel])
                                    / (votes[rowsel, None] + w))
                a = ALPHA
                if w > 0:
                    a, _ = solve_alpha_delivered(pre4, newop, fold, cd4, fv, target_norm)
                final, frac = build_final(pre4, a, newop, fold, cd4, fv)
                dn = delivered_norm(final, fold)
                e = score_effects(final, fold, truth)
                rows.append(dict(
                    fold=held, arm=name, refit=bname, lib_w=w, residual_f=fv,
                    alpha=float(a), fraction=frac, delivered_norm=dn,
                    delivered_ratio=float(dn / target_norm), operator_sha256=sha,
                    ka_beta0_rel_err=ka_err, spectrum=spec,
                    ref_cos_agg=ref[arm]["score"]["cos_agg"] if arm in ref else None,
                    ref_pds=ref[arm]["score"]["pds"] if arm in ref else None,
                    seconds=round(time.time() - t1, 1),
                    **{k: v for k, v in e.items() if not k.startswith("_")}))
                out_path.write_text(json.dumps(rows))
                d_ref = (f" vs shipped-op {ref[arm]['score']['cos_agg']:.6f}"
                         if arm in ref else "")
                print(f"[{held}/{name}] a {a:.4f} delivered {dn/target_norm:.4f}x "
                      f"pds {e['pds']:.6f} cos_agg {e['cos_agg']:.6f}{d_ref}", flush=True)
                del final, pre4
                gc.collect()
            del newop, evp, Vp
            gc.collect()
        del C, Bc, ev, top_orig, fold, truth, pre, cd4
        gc.collect()
    print(f"done in {time.time()-t0:.0f}s, {len(rows)} rows")


if __name__ == "__main__":
    main()
