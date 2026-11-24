"""423 exploration: rank-free capacity structure of the emitted base blocks (dumped by run423.py --dump)."""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.argv = [sys.argv[0]] + sys.argv[1:]
import run423 as RN  # sets paths
E = RN.E
r, cc, emit_c = E.load_reference()
keep = cc.keep; cs = cc.cs; m = cs.shape[1]
ref_lin = cc.ref_lin
files = sorted((HERE / "work/explore").glob("block_s0_*.npz"))
out = []
for f in files:
    z = np.load(f)
    import scipy.sparse as sp
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    n_car = int(z["n_car"]); tgt = str(z["target"])
    L = X.sum(1).astype(np.float64); scale = L / 1e6
    Xk = X[:, keep]
    V = np.log1p(Xk / scale[:, None])            # scorer's exact values (400, G)
    A = (Xk / scale[:, None]).mean(0)
    lfc = np.log2((A + 1e-9) / (ref_lin + 1e-9))
    G = Xk.shape[1]
    above_max = np.zeros(G, int); tied = np.zeros(G, int); zero = np.zeros(G, int)
    up_cap_cpm = np.zeros(G); dn_cap_cpm = np.zeros(G)
    up_inf = np.zeros(G, bool)
    ball_gap_lo = np.zeros(G); ball_gap_hi = np.zeros(G); ball_frac = np.zeros(G)
    for j in range(G):
        c = cs[j]; v = V[:n_car, j]
        lt = np.searchsorted(c, v, "left"); rt = np.searchsorted(c, v, "right")
        tie = rt > lt
        tied[j] = tie.sum(); zero[j] = (Xk[:n_car, j] == 0).sum()
        above = lt == m
        above_max[j] = above.sum()
        hi = np.where(rt < m, c[np.minimum(rt, m - 1)], np.inf)
        lo = np.where(lt > 0, c[np.maximum(lt - 1, 0)], -np.inf)
        s = scale[:n_car]
        # count headroom (float estimate; exact version in rp423)
        xhi = np.where(np.isfinite(hi), np.floor(np.expm1(hi) * s - 1e-9), np.inf)
        xlo = np.where(np.isfinite(lo), np.floor(np.expm1(lo) * s) + 1, 0)
        x = Xk[:n_car, j]
        upc = np.where(tie, 0, np.maximum(xhi - x, 0)); dnc = np.where(tie, 0, np.maximum(x - np.maximum(xlo, 0), 0))
        fin = np.isfinite(upc)
        up_cap_cpm[j] = (upc[fin] / s[fin]).sum() / 400; up_inf[j] = (~fin).any()
        dn_cap_cpm[j] = (dnc / s).sum() / 400
        # ballast
        vb = V[n_car:, j]; ltb = np.searchsorted(c, vb, "left"); rtb = np.searchsorted(c, vb, "right")
        hib = c[rtb[0]] if rtb[0] < m else np.inf; lob = c[ltb[0] - 1] if ltb[0] > 0 else -np.inf
        ball_gap_hi[j] = (np.expm1(hib) - np.expm1(vb[0])) if np.isfinite(hib) else np.inf
        ball_gap_lo[j] = (np.expm1(vb[0]) - np.expm1(lob)) if np.isfinite(lob) else np.inf
        ball_frac[j] = ltb[0] / m
    fl = np.isfinite(lfc) & (A > 0)
    for k in (1.4,):
        At = ref_lin * 2.0 ** (k * lfc)
        dA = np.where(fl, At - A, 0)
        need_up = dA > 0; need_dn = dA < 0
        feas_up = np.where(up_inf, 1.0, np.minimum(1, up_cap_cpm / np.maximum(dA, 1e-12)))
        feas_dn = np.minimum(1, dn_cap_cpm / np.maximum(-dA, 1e-12))
    rec = dict(target=tgt, n_car=n_car, G=G,
               above_max_any=float((above_max > 0).mean()), above_max_cells_med=float(np.median(above_max)),
               tied_frac=float(tied.sum() / (G * n_car)), zero_frac=float(zero.sum() / (G * n_car)),
               feas_up_med=float(np.median(feas_up[need_up])), feas_up_mean=float(feas_up[need_up].mean()),
               feas_dn_med=float(np.median(feas_dn[need_dn])), feas_dn_mean=float(feas_dn[need_dn].mean()),
               n_up=int(need_up.sum()), n_dn=int(need_dn.sum()), up_inf=float(up_inf.mean()),
               ball_gap_lo_med=float(np.median(ball_gap_lo)), ball_gap_hi_med=float(np.median(ball_gap_hi)),
               ball_pos_med=float(np.median(ball_frac)),
               lfc_abs_med=float(np.median(np.abs(lfc[fl]))))
    # top-|lfc| genes detail
    top = np.argsort(-np.where(fl, np.abs(lfc), -1))[:200]
    rec["top200_feas_up"] = float(np.mean(feas_up[top][need_up[top]])) if need_up[top].any() else None
    rec["top200_feas_dn"] = float(np.mean(feas_dn[top][need_dn[top]])) if need_dn[top].any() else None
    rec["top200_above_max_any"] = float((above_max[top] > 0).mean())
    rec["top200_up_frac"] = float((lfc[top] > 0).mean())
    print(json.dumps(rec), flush=True)
    out.append(rec)
    np.savez(HERE / f"work/explore/cap_{f.stem}.npz", lfc=lfc, A=A, up_cap=up_cap_cpm, dn_cap=dn_cap_cpm, up_inf=up_inf,
             above_max=above_max, tied=tied, zero=zero, ball_gap_lo=ball_gap_lo, ball_gap_hi=ball_gap_hi, ball_frac=ball_frac)
