"""196 Phase 1b — can the depth knob be COMPENSATED?

The naive depth sweep moves two things the mandate says must be held fixed or explicitly
re-optimised: DE power (`n_pred` explodes) and the jackknife credit (`C_hat` falls through the
`min(C_hat, C_p)` cap, which is what actually prices `mse`).  Both are functions of the emitted
cell-to-cell variance, and the quasi emitter has exactly one knob for that at fixed per-gene
totals -- the gamma-Poisson `overdispersion` phi (Var = mu + phi mu^2).

Raising phi raises `C_hat` and lowers Wilcoxon power, i.e. it moves BOTH quantities back the
way depth moved them.  So the question is whether ONE phi restores both at a given depth.  If
it does, the compensated depth arm keeps whatever `pds` the finer profile buys and pays nothing
for it; if it does not, depth is not a free parameter at any setting.

Measured on the challenge's own control cells, so `n_pred` and `C_hat` are both exact.

    MAXT=6 uv run python comp196.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import false_discovery_control, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/092-call-set-truncation"))
import recon196                                                # noqa: E402
import vcc2026.rc_emit as RC                                   # noqa: E402
from de_kernel import wilcoxon_z                               # noqa: E402
from vcc2026.donor_transfer import emit_counts                 # noqa: E402

TS = 50_000.0
n, DISPERSION = 400, 0.1
MAXT = int(os.environ.get("MAXT", "6"))
GRID = [(1.0, 4.0), (2.0, 4.0), (2.0, 8.0), (2.0, 16.0), (2.0, 32.0), (2.0, 64.0),
        (4.0, 4.0), (4.0, 16.0), (4.0, 64.0), (4.0, 256.0)]
if os.environ.get("GRID"):
    GRID = [tuple(float(y) for y in x.split(":")) for x in os.environ["GRID"].split(",")]


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction
    t0 = time.time()
    res, cp = {}, {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        Ad = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        M = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        ids = np.asarray(Ad.obs["ntc_id"].astype(str).values)
        cp[ctx] = float(np.mean([
            float(jackknife_correction(sp.csr_matrix(M[np.flatnonzero(ids == nm)]),
                                       np.zeros(int((ids == nm).sum()), np.intp), 1, TS,
                                       chunk=64)[0])
            for nm in dict.fromkeys(ids.tolist())]))
        cb = RC.control_block_from_counts(M, min_cpm=5.0)
        del Ad, M
        keep, cs, tie, ref_lin = cb.keep, cb.cs, cb.tie, cb.ref_lin
        G = D["genes"].size
        sel = np.linspace(0, D["sums"].shape[0] - 1, MAXT).astype(int)
        for i in sel:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
            sv = np.full(G, 1.0)
            sv[order[:min(8000, int(D["gate"].sum()))]] = 0.25
            for m, phi in GRID:
                rng = np.random.default_rng(90_000 + int(i))
                v, c, ip = emit_counts(R.profile(m), np.full(n, R.emit_L(m)), DISPERSION,
                                       rng, mode="quasi", spread=sv, overdispersion=phi)
                blk = sp.csr_matrix((v.astype(np.float64), c.astype(np.int64),
                                     np.asarray(ip)), shape=(n, G))
                ch = float(jackknife_correction(blk, np.zeros(n, np.intp), 1, TS, chunk=64)[0])
                B = np.asarray(blk.todense())
                rs = B.sum(1)
                rs[rs == 0] = 1.0
                nb = np.log1p(B[:, keep] * (1e6 / rs[:, None]))
                z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
                q = false_discovery_control(2 * norm.sf(np.abs(z)), method="bh")
                lfc = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
                call = q < 0.05
                si = np.flatnonzero(D["genes"][keep] == D["targets"][i])
                if si.size:
                    call[si[0]] = False
                res.setdefault((m, phi), []).append(dict(
                    ctx=ctx, chat=ch, cp=cp[ctx], n_pred=int(call.sum()),
                    up=float((lfc[call] > 0).mean()) if call.any() else np.nan,
                    nnz=int((B > 0).sum()), max_cell=int(B.sum(1).max())))
                del blk, B, nb
        print(f"  ctx {ctx} done [{time.time()-t0:.0f}s]", flush=True)
        del cb, cs, tie

    print(f"\nCOMPENSATION SCREEN — a15's operating point is (m=1, phi=4): "
          f"n_pred 4,436 (full panel), C_hat/C_p 1.094 (191)\n")
    print(f"{'depth':>7}{'phi':>7}{'n_pred':>9}{'vs a15':>9}{'Chat/C_p':>10}{'claim/cap':>11}"
          f"{'rows<cap':>10}{'up frac':>9}{'nnz/block':>12}{'max cell':>11}")
    base_np = float(np.median([r["n_pred"] for r in res[(1.0, 4.0)]]))
    out = []
    for m, phi in GRID:
        rs = res[(m, phi)]
        npd = float(np.median([r["n_pred"] for r in rs]))
        ratio = float(np.mean([r["chat"] / r["cp"] for r in rs]))
        claim = float(np.mean([min(r["chat"], PRED_TRACE_CAP_K * r["cp"]) / r["cp"]
                               for r in rs]))
        below = int(sum(r["chat"] < PRED_TRACE_CAP_K * r["cp"] for r in rs))
        up = float(np.mean([r["up"] for r in rs]))
        nnz = float(np.mean([r["nnz"] for r in rs]))
        mx = max(r["max_cell"] for r in rs)
        out.append(dict(depth=m, phi=phi, n_pred=npd, chat_ratio=ratio, claim_ratio=claim,
                        rows_below_cap=below, up=up, nnz=nnz, max_cell=mx))
        print(f"{m:>7g}{phi:>7g}{npd:>9,.0f}{npd/base_np:>9.3f}{ratio:>10.3f}"
              f"{claim:>11.4f}{below:>10d}{up:>9.4f}{nnz:>12,.0f}{mx:>11,.0f}")
    json.dump(dict(rows=out, C_p=cp, n_rows=len(res[(1.0, 4.0)])),
              open(WORK / "comp196.json", "w"), indent=1, default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
