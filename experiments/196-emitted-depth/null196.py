"""196 — the ZERO-EFFECT NULL at each emitted depth.

177 §0 diagnosed the quasi emitter's central DE defect: it places all 400 cells at ONE library
size, so the emitted CPM lives on a lattice the real controls do not, and the scorer's
Mann-Whitney is not null under the null -- ~19% of gated genes clear BH at ZERO predicted
effect, ~88% of them called UP.  That is the artefact a21/a22 spent two board slots on.

Emitted depth changes the lattice spacing (finer at higher L) and the emitted variance (larger
absolute, smaller relative), so it is a candidate fix for the same defect.  This measures it
directly: emit the CONTROL profile itself -- zero predicted effect on every gene -- at each
depth and count what the scorer calls.

The `m = 1` row is the known-answer row: it must reproduce 177's ~19% / ~88% UP.

    MAXT=6 uv run python null196.py
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

n, DISPERSION, PHI = 400, 0.1, 4.0
MAXT = int(os.environ.get("MAXT", "6"))
DEPTHS = tuple(float(x) for x in os.environ.get("DEPTHS", "1,2,4,8,16").split(","))


def main() -> None:
    t0 = time.time()
    res = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        Ad = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        M = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        cb = RC.control_block_from_counts(M, min_cpm=5.0)
        del Ad, M
        keep, cs, tie, ref_lin = cb.keep, cb.cs, cb.tie, cb.ref_lin
        G = D["genes"].size
        sel = np.linspace(0, D["sums"].shape[0] - 1, MAXT).astype(int)
        for i in sel:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            x0 = R.K * D["r"]                       # the CONTROL profile, same normaliser
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
            sv = np.full(G, 1.0)
            sv[order[:min(8000, int(D["gate"].sum()))]] = 0.25
            for m in DEPTHS:
                Lm = R.emit_L(m)
                tot = np.floor(m * x0)
                prof = (tot + 0.5) / (n * Lm)
                rng = np.random.default_rng(90_000 + int(i))
                v, c, ip = emit_counts(prof, np.full(n, Lm), DISPERSION, rng, mode="quasi",
                                       spread=sv, overdispersion=PHI)
                B = np.zeros((n, G), np.int64)
                for r_ in range(n):
                    s_ = slice(ip[r_], ip[r_ + 1])
                    B[r_, c[s_]] = v[s_]
                assert np.array_equal(B.sum(0), tot.astype(np.int64))
                rs = B.sum(1).astype(np.float64)
                rs[rs == 0] = 1.0
                nb = np.log1p(B[:, keep] * (1e6 / rs[:, None]))
                z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
                q = false_discovery_control(2 * norm.sf(np.abs(z)), method="bh")
                lfc = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
                call = q < 0.05
                res.setdefault(m, []).append(dict(
                    ctx=ctx, gated=int(keep.sum()), n_call=int(call.sum()),
                    frac=float(call.mean()),
                    up=float((lfc[call] > 0).mean()) if call.any() else np.nan,
                    median_abs_lfc=float(np.median(np.abs(lfc))),
                    median_lfc=float(np.median(lfc))))
                del B, nb
        print(f"  ctx {ctx} done [{time.time()-t0:.0f}s]", flush=True)
        del cb, cs, tie

    print("\nZERO-EFFECT NULL — emit the CONTROL profile, count what the scorer calls\n"
          "(177 §0 known-answer row at m=1: ~19% of gated genes called, ~88% of them UP)\n")
    print(f"{'depth':>7}{'gated':>8}{'called':>9}{'% of gate':>11}{'% UP':>9}"
          f"{'median |lfc|':>14}{'median lfc':>13}")
    out = []
    for m in DEPTHS:
        rs = res[m]
        row = dict(depth=m, gated=float(np.mean([r["gated"] for r in rs])),
                   called=float(np.median([r["n_call"] for r in rs])),
                   frac=float(np.mean([r["frac"] for r in rs])),
                   up=float(np.mean([r["up"] for r in rs])),
                   mabs=float(np.mean([r["median_abs_lfc"] for r in rs])),
                   med=float(np.mean([r["median_lfc"] for r in rs])))
        out.append(row)
        print(f"{m:>7g}{row['gated']:>8,.0f}{row['called']:>9,.0f}{row['frac']*100:>10.2f}%"
              f"{row['up']*100:>8.1f}%{row['mabs']:>14.5f}{row['med']:>+13.5f}")
    json.dump(out, open(WORK / "null196.json", "w"), indent=1, default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
