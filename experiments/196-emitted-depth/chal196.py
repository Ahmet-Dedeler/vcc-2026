"""196 — the CHALLENGE side of the depth sweep: format contract, DE power, and the
zero-direction call count (Phase 2).

Everything here is measured against the challenge's OWN control cells, which per BRIEF §3 are
the comparison group on both sides, so `n_pred`, the delivered lfc and the call set are exact
rather than proxies.  There is no ground truth on A/B/C, so `reach`/`fid`/`jac`/`nmae` are NOT
computed here -- they come from h1-test (`sweep196.py`).

What this file produces
  1. `nnz` and `max counts per cell` per emitted block at each depth, extrapolated to the full
     900-block package and priced against the 4.75e9 stored-entry cap and the 1e6 per-cell cap
     (AGENTS.md §2).  Falsifier 1.
  2. `n_pred` at each depth -- the DE-power drift the depth knob causes.
  3. Phase 2: among the genes a15 CALLS, how many carry a predicted lfc of exactly zero, how
     many sit outside the builder's own 10-CPM effect gate (so their lfc is pure emission
     artefact rather than prediction), and what the resulting guaranteed-miss bill is.

    MAXT=20 DEPTHS=1,2,4,8,16 uv run python chal196.py
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
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/092-call-set-truncation"))
import recon196                                            # noqa: E402
import vcc2026.rc_emit as RC                                # noqa: E402
from de_kernel import wilcoxon_z                            # noqa: E402
from vcc2026.donor_transfer import emit_counts              # noqa: E402

WORK = HERE / "work"
n, PHI, DISPERSION = 400, 4.0, 0.1
DEPTHS = tuple(float(x) for x in os.environ.get("DEPTHS", "1,2,4,8,16").split(","))
ARMS = tuple(os.environ.get("ARMS", "quasi").split(","))
MAXT = int(os.environ.get("MAXT", "20"))
NNZ_CAP = 4_750_000_000
CELL_CAP = 1_000_000
# The SHIPPED packages' true full-panel stored-entry counts (191 §0), used to anchor the
# extrapolation so the m = 1 row is an exact known-answer row rather than an estimate.
SHIPPED_NNZ = {"quasi": 2_377_814_190, "rc_hybrid": 1_647_149_337}


def to_block(vals, cols, ptr, G):
    B = np.zeros((n, G), np.int64)
    for r in range(n):
        sl = slice(ptr[r], ptr[r + 1])
        B[r, cols[sl]] = vals[sl]
    return B


def main() -> None:
    t0 = time.time()
    out = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        A = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        X = A.X.tocsr() if not sp.isspmatrix_csr(A.X) else A.X
        cb = RC.control_block_from_counts(X, min_cpm=5.0)
        del A, X
        keep, cs, tie, ref_lin = cb.keep, cb.cs, cb.tie, cb.ref_lin
        G_full = D["genes"].size
        # the two gates: the SCORER's DE universe vs the BUILDER's effect gate
        scorer = keep
        builder = D["gate"]
        both = scorer & builder
        dead = scorer & ~builder            # scored by DE, zero predicted effect by construction
        print(f"\nctx {ctx}: {cb.n_controls} controls, scorer gate {int(scorer.sum())}, "
              f"builder gate {int(builder.sum())}, scored-but-uneffected {int(dead.sum())} "
              f"[{time.time()-t0:.0f}s]", flush=True)

        sel = np.linspace(0, D["sums"].shape[0] - 1, MAXT).astype(int)
        rows = {f"{a}_m{m:g}": [] for a in ARMS for m in DEPTHS}
        shipped_nnz = D["nnz"][sel].astype(np.float64)
        cbrc = cb if "rc_hybrid" in ARMS else None
        for j, i in enumerate(sel):
            Rec = recon196.Recon(D["sums"][i], D["r"], builder, n=n, seed=1000 + int(i))
            assert np.array_equal(Rec.totals(1.0), D["sums"][i])
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G_full), -np.where(scorer, eff, -np.inf)))
            sv = np.full(G_full, 1.0)
            sv[order[:min(8000, int(scorer.sum()))]] = 0.25
            for m in DEPTHS:
                Lm, prof, want = Rec.emit_L(m), Rec.profile(m), Rec.totals(m)
                for arm in ARMS:
                  rng = np.random.default_rng(90_000 + int(i))
                  if arm == "quasi":
                    vals, cols, ptr = emit_counts(prof, np.full(n, Lm), DISPERSION, rng,
                                                  mode="quasi", spread=sv, overdispersion=PHI)
                    B = to_block(vals, cols, ptr, G_full)
                  else:
                    keyv = np.where(builder, np.log1p(50_000.0 * prof) - D["b"], 0.0)
                    B, _dg = RC.rc_hybrid_block(
                        prof, np.full(n, Lm), DISPERSION, rng,
                        pass2_rng=np.random.default_rng(90_000 + int(i)), spread=sv,
                        overdispersion=PHI, control=cbrc, key=keyv,
                        jitter_seed=1000 + int(i), emit_counts=emit_counts)
                  assert np.array_equal(B.sum(0), want)
                  rs = B.sum(1)
                  s = rs.astype(np.float64)
                  s[s == 0] = 1.0
                  nb = np.log1p(B[:, keep] * (1e6 / s[:, None]))
                  z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
                  p = 2 * norm.sf(np.abs(z))
                  q = false_discovery_control(p, method="bh")
                  lfc = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
                  call = q < 0.05
                  # exclude the perturbation's own target gene, BRIEF §4
                  self_i = np.flatnonzero(D["genes"][keep] == D["targets"][i])
                  if self_i.size:
                      call[self_i[0]] = False
                  dead_k = dead[keep]
                  rows[f"{arm}_m{m:g}"].append(dict(
                      nnz=int((B > 0).sum()), max_cell=int(rs.max()), min_cell=int(rs.min()),
                      cell_sd=float(rs.std()), emit_L=Lm, S=int(want.sum()),
                      shipped_nnz=float(shipped_nnz[j]), n_pred=int(call.sum()),
                      up=float((lfc[call] > 0).mean()) if call.any() else float("nan"),
                      zero_lfc=int((call & (lfc == 0.0)).sum()),
                      tiny_lfc=int((call & (np.abs(lfc) < 1e-3)).sum()),
                      dead_called=int((call & dead_k).sum()),
                      dead_frac=float((call & dead_k).sum() / max(int(call.sum()), 1)),
                      zero_total=int((call & (B.sum(0)[keep] == 0)).sum()),
                  ))
                  del B, nb
            if (j + 1) % 5 == 0:
                print(f"  {ctx} {j+1}/{sel.size} [{time.time()-t0:.0f}s]  " +
                      " ".join(f"{k}={v[-1]['n_pred']}" for k, v in rows.items()),
                      flush=True)
        out[ctx] = {k: {j: float(np.median([r[j] for r in v])) for j in v[0]}
                    for k, v in rows.items()}
        out[ctx]["_gates"] = dict(scorer=int(scorer.sum()), builder=int(builder.sum()),
                                  dead=int(dead.sum()), n_blocks=int(D["sums"].shape[0]))
        del cb, cs, tie
    json.dump(out, open(WORK / "chal196.json", "w"), indent=1)

    print("\n\nFORMAT CONTRACT — the m=1 row is ANCHORED on the shipped package's true "
          f"stored-entry count\n(cap {NNZ_CAP:,} entries, {CELL_CAP:,} counts/cell); deeper "
          "rows scale it by the MEASURED\nnnz ratio on the same sampled blocks, so the "
          "sampling error cancels out of the ratio.")
    fmt = {}
    for arm in ARMS:
        base = float(np.mean([out[c][f"{arm}_m1"]["nnz"] for c in "ABC"]))
        print(f"\n--- emitter `{arm}` (shipped package nnz {SHIPPED_NNZ[arm]:,} = "
              f"{100*SHIPPED_NNZ[arm]/NNZ_CAP:.2f}% of cap) ---")
        print(f"{'depth':>7}{'nnz/block':>13}{'nnz ratio':>11}{'package nnz':>16}"
              f"{'% of cap':>10}{'max counts/cell':>17}{'% cell cap':>12}{'n_pred':>9}")
        for m in DEPTHS:
            k = f"{arm}_m{m:g}"
            per = float(np.mean([out[c][k]["nnz"] for c in "ABC"]))
            ratio = per / base
            pkg = SHIPPED_NNZ[arm] * ratio
            mx = max(out[c][k]["max_cell"] for c in "ABC")
            npd = float(np.mean([out[c][k]["n_pred"] for c in "ABC"]))
            fmt[k] = dict(nnz_per_block=per, nnz_ratio=ratio, package_nnz=pkg,
                          pct_cap=100 * pkg / NNZ_CAP, max_cell=mx,
                          pct_cell_cap=100 * mx / CELL_CAP, n_pred=npd)
            print(f"{m:>7g}{per:>13,.0f}{ratio:>11.4f}{pkg:>16,.0f}"
                  f"{100*pkg/NNZ_CAP:>9.2f}%{mx:>17,.0f}{100*mx/CELL_CAP:>11.2f}%"
                  f"{npd:>9,.0f}")
    ship = float(np.mean([np.mean([out[c][f"quasi_m1"]["shipped_nnz"]] ) for c in "ABC"]))
    print(f"\n  sampling check: re-emitted nnz/block at m=1 {base:,.0f} against the shipped "
          f"blocks' own {ship:,.0f} ({100*(base/ship-1):+.3f}%)")
    json.dump(dict(per_context=out, format=fmt), open(WORK / "chal196.json", "w"),
              indent=1, default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
