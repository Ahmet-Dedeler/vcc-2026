"""What we can measure about our OWN submission on the REAL challenge panel.

The reference perturbation profiles are withheld, so `X` (alignment) and `R` (real amplitude)
are not computable. Everything else in

    MSE_c - 1 = ( A_c - 2 X_c + Kc_c - rho_c * claim_c ) / sumD_c

IS computable, because `data/challenge/context_{A,B,C}.h5ad` are the reference bundle's own
control cells:

    A_c      from the submission's group sums + the real control profile   EXACT
    Kc_c     from the real control's jackknife                             EXACT
    B_c      from the submission's group sums alone                        EXACT
    C_hat_p  from the submission's cells                                   EXACT
    claim_c  = sum_p w_p min(C_hat_p, C_p) -- needs the reference's C_p    ESTIMATED
               (proxy: jackknife of random 400-cell subsets of that context's control cells,
                which have the same cell count and the same 20k median depth)

    uv run python mech_challenge.py <stats.npz> [label]
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/096-mse-floor")
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402
from cell_eval2.moments import jackknife_correction  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

ROOT = "~/Code/virtual-cell-challenge-2026"
TS = 50_000.0
CTX = ("A", "B", "C")


def cp_proxy(n_draw=12, n_cells=400, seed=0):
    """`C_p` for a 400-cell, 20k-UMI group in each context, from that context's own controls.

    The reference's perturbation groups have exactly 400 cells at the same median depth, so a
    400-cell draw from the control pool earns the same order of correction. It is a PROXY: real
    perturbation groups are slightly more dispersed than controls, so this is a mild
    UNDER-estimate of C_p, which makes `min(C_hat, C_p)` an under-estimate too -- conservative
    in the direction that matters (it under-states the credit we are allowed).
    """
    import anndata as ad
    rng = np.random.default_rng(seed)
    out = {}
    for c in CTX:
        A = ad.read_h5ad(f"{ROOT}/data/challenge/context_{c}.h5ad")
        X = A.X.tocsr()
        vals = []
        for _ in range(n_draw):
            rows = np.sort(rng.choice(X.shape[0], n_cells, replace=False))
            sub = csr_matrix(X[rows])
            vals.append(jackknife_correction(sub, np.zeros(n_cells, np.intp), 1, TS, chunk=64)[0])
        out[c] = np.asarray(vals)
        print(f"  C_p proxy, context {c}: mean {out[c].mean():.5f}  sd {out[c].std():.5f} "
              f"over {n_draw} draws of {n_cells} control cells", flush=True)
    return out


def main(stats_path, label="submission", cp=None):
    z = np.load(stats_path, allow_pickle=True)
    perts = z["perts"].astype(str)
    sums, counts, jk, genes = z["sums"], z["counts"], z["jk"], z["genes"].astype(str)
    ctrl = np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    assert np.array_equal(genes, ctrl["genes"].astype(str)), "gene axis mismatch"
    G = genes.size
    gpos = {g: i for i, g in enumerate(genes)}
    pm_all = bulk_lognorm_means(sums.astype(np.float64), TS)

    rows = []
    for c in CTX:
        sel = np.flatnonzero(np.char.startswith(perts, f"{c}|"))
        tgt = np.asarray([p.split("|", 1)[1] for p in perts[sel]])
        bctrl = bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
        Cctrl = float(ctrl[f"jk_{c}"])
        cols = {k: gpos[t] for k, t in enumerate(tgt) if t in gpos}
        w = _row_weights(sel.size, G, cols)
        pm = pm_all[sel]
        dh = pm - bctrl[None, :]
        a = np.einsum("ij,ij->i", dh, dh)
        for k, col in cols.items():
            a[k] -= dh[k, col] ** 2
        A = float((a * w).sum())
        B = _across_pert_budget(pm, cols, w)
        Chat = jk[sel]
        claim_uncapped = float((Chat * w).sum())
        Kc = float((np.full(sel.size, Cctrl) * w).sum())
        r = dict(ctx=c, n=sel.size, A=A, B=B, Chat_sum=claim_uncapped, Kc=Kc,
                 Cctrl=Cctrl, resolved=len(cols), Chat_med=float(np.median(Chat)))
        if cp is not None:
            Cp = float(np.mean(cp[c]))
            capped = np.minimum(Chat, Cp)
            r["Cp_proxy"] = Cp
            r["claim"] = float((capped * w).sum())
            r["cap_binds"] = int((Chat > Cp).sum())
            r["K_proxy"] = float((np.full(sel.size, Cp) * w).sum())
            r["rho"] = min(1.0, B / r["claim"]) if r["claim"] > 0 else 1.0
        rows.append(r)

    print(f"\n===== {label}: measured on the real challenge panel =====")
    hdr = ("ctx", "n", "A", "B", "B/claim", "sum w*Chat", "claim", "Kc", "C_ctrl",
           "Cp_proxy", "cap binds")
    print(f"  {hdr[0]:>4}{hdr[1]:>5}{hdr[2]:>12}{hdr[3]:>12}{hdr[4]:>10}{hdr[5]:>12}"
          f"{hdr[6]:>12}{hdr[7]:>10}{hdr[8]:>9}{hdr[9]:>10}{hdr[10]:>11}")
    for r in rows:
        cl = r.get("claim", float("nan"))
        print(f"  {r['ctx']:>4}{r['n']:>5}{r['A']:>12.6f}{r['B']:>12.6f}"
              f"{(r['B']/cl if cl==cl else float('nan')):>10.3f}{r['Chat_sum']:>12.6f}"
              f"{cl:>12.6f}{r['Kc']:>10.6f}{r['Cctrl']:>9.5f}"
              f"{r.get('Cp_proxy', float('nan')):>10.5f}{r.get('cap_binds', -1):>11d}")
    tot = {k: float(np.mean([r[k] for r in rows])) for k in ("A", "B", "Kc", "Chat_sum")}
    print(f"  {'mean':>4}{'':>5}{tot['A']:>12.6f}{tot['B']:>12.6f}{'':>10}{tot['Chat_sum']:>12.6f}")
    return rows


if __name__ == "__main__":
    cp = cp_proxy()
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "submission", cp=cp)
    np.savez(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy.npz",
             **{k: v for k, v in cp.items()})
