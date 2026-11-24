"""230 Phase 1b — EMIT THE CELLS.  Does a21/a22's emission law move the profile at all?

magnitude230.py argues from `emit_counts`'s source that `mode="resample_mix"` recomputes the
SAME `totals = floor(n L p)` the quasi law does and repairs to it, so a21 and a22 -- rule 13's
only two calibration points -- cannot have moved the delivered per-gene panel totals.  An
argument is not a measurement.  This file emits real cells through the production `emit_counts`
on real challenge control cells and reads the delivered totals back.

Four laws, one block, one seed, one profile:

  quasi         a15 / the incumbent                                    <- known-answer row
  x=0.14        a21: `VCC_RESAMPLE_X=0.14`, top 14% of |effect| parametric, the rest
                resampled from real control cells  (LOG 2026-08-31)
  x=0.14 inv    the complementary mask, i.e. the head resampled -- a22's family, carried
                because the LOG's two descriptions of which half a22 resampled disagree
  phi=0         a8's channel: `overdispersion` 4.0 -> 0.0, the other emission-law knob
                with a board point, to show the identity is not special to resample_mix

Also measured, because it is the channel these arms DID move: the emitted per-cell library
distribution.  `quasi` puts every cell at exactly `L`; the resampling law inherits the real
control libraries.  That is the "cell channel" of LOG 2026-09-01.

    uv run python emit230.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
import recon196                                                        # noqa: E402
from vcc2026.donor_transfer import emit_counts                         # noqa: E402

TS = 50_000.0
N = 400
PHI = 4.0
SPREAD_KEEP, SPREAD_TAIL, KEEP_TOP = 0.25, 1.0, 8000
N_BLOCKS = int(6)


def sha16(a) -> str:
    return hashlib.sha256(np.asarray(a, dtype=np.int64).tobytes()).hexdigest()[:16]


def load_controls(context: str, max_cells: int = 4000, seed: int = 0):
    """The context's real control cells, as a dense float64 block of `max_cells` rows."""
    with h5py.File(ROOT / f"data/challenge/context_{context}.h5ad", "r") as fh:
        g = fh["X"]
        shape = tuple(int(v) for v in g.attrs["shape"])
        indptr = g["indptr"][:]
        rng = np.random.default_rng(seed)
        rows = np.sort(rng.choice(shape[0], min(max_cells, shape[0]), replace=False))
        blocks, ptr = [], [0]
        for r in rows:
            lo, hi = int(indptr[r]), int(indptr[r + 1])
            blocks.append((g["indices"][lo:hi], g["data"][lo:hi]))
            ptr.append(ptr[-1] + (hi - lo))
        idx = np.concatenate([b[0] for b in blocks]).astype(np.int32)
        dat = np.concatenate([b[1] for b in blocks]).astype(np.float64)
    M = sp.csr_matrix((dat, idx, np.asarray(ptr, dtype=np.int64)),
                      shape=(len(rows), shape[1]))
    return M


def spread_vector(gate_rank: np.ndarray, G: int) -> np.ndarray:
    """a15's tiered spread: the top `KEEP_TOP` gated genes at 0.25, the tail at 1.0."""
    s = np.full(G, SPREAD_TAIL)
    s[gate_rank[:KEEP_TOP]] = SPREAD_KEEP
    return s


def main() -> None:
    t0 = time.time()
    D = recon196.load_chal("A")
    b, r, gate = D["b"], D["r"], D["gate"]
    G = b.size
    C = load_controls("A")
    libs_all = np.asarray(C.sum(axis=1)).ravel()
    print(f"controls {C.shape}, library median {np.median(libs_all):.0f}, "
          f"CV {libs_all.std()/libs_all.mean():.3f}", flush=True)
    ctrl_profile = np.asarray(C.sum(axis=0)).ravel()
    ctrl_profile = ctrl_profile / ctrl_profile.sum()

    rows, cellstats = [], {}
    for bi in range(N_BLOCKS):
        R = recon196.Recon(D["sums"][bi], r, gate, n=N, seed=4200 + bi)
        # Hand `emit_counts` the profile whose own floor reproduces the artifact's totals.
        prof = R.profile(1.0)
        Lp = R.emit_L(1.0)
        pick_rng = np.random.default_rng(90_000 + bi)
        picks = pick_rng.choice(C.shape[0], N, replace=False)
        cells = np.asarray(C[picks].todense(), dtype=np.float64)
        libs = libs_all[picks]
        # `emit_counts` recomputes L = int(min(median(lib_sizes), max_counts)); feed it a
        # constant vector so L is exactly the reconstruction's, as production's own pinning
        # of `L` to the picked control libraries would give.
        lib_sizes = np.full(N, Lp, dtype=np.float64)

        # a15's ordering axis, and the two param masks a21/a22 use.
        eff_proxy = np.log((prof + 1e-12) / (r + 1e-12))          # delivered log fold change
        order = np.argsort(-np.abs(np.where(gate, eff_proxy, 0.0)), kind="stable")
        k14 = int(round(0.14 * G))
        m_a21 = np.zeros(G, dtype=bool)
        m_a21[order[:k14]] = True                     # top 14% parametric  (a21, LOG 08-31)
        m_a22 = ~m_a21                                # the complement      (a22's family)
        gate_rank = order[np.isin(order, np.flatnonzero(gate))]
        spread = spread_vector(gate_rank, G)

        laws = {
            "quasi":      dict(mode="quasi", spread=spread, overdispersion=PHI),
            "a8_phi0":    dict(mode="quasi", spread=spread, overdispersion=0.0),
            "a21_x0.14":  dict(mode="resample_mix", spread=spread, overdispersion=PHI,
                               control_counts=cells, control_profile=ctrl_profile,
                               param_mask=m_a21),
            "a22_inv":    dict(mode="resample_mix", spread=spread, overdispersion=PHI,
                               control_counts=cells, control_profile=ctrl_profile,
                               param_mask=m_a22),
        }
        ref_T = None
        for name, kw in laws.items():
            rng = np.random.default_rng(777)                    # SAME seed for every law
            vals, cols, ptr = emit_counts(prof, lib_sizes, 0.1, rng, mode=kw.pop("mode"),
                                          **kw)
            B = sp.csr_matrix((vals.astype(np.float64), cols, ptr), shape=(N, G))
            T = np.asarray(B.sum(axis=0)).ravel().astype(np.int64)
            cl = np.asarray(B.sum(axis=1)).ravel()
            if name == "quasi":
                ref_T = T
            rows.append(dict(
                block=str(D["targets"][bi]), law=name, digest=sha16(T),
                same_as_quasi=int(sha16(T) == sha16(ref_T)),
                n_gene_diff=int(np.count_nonzero(T - ref_T)),
                max_abs_dT=int(np.abs(T - ref_T).max()),
                relL1=float(np.abs(T - ref_T).sum() / max(ref_T.sum(), 1)),
                lib_cv=float(cl.std() / max(cl.mean(), 1e-9)),
                lib_min=float(cl.min()), lib_max=float(cl.max()),
                nnz=int(B.nnz),
                frac_resampled=float(1.0 - kw.get("param_mask", np.ones(G, bool)).mean()),
            ))
            cellstats.setdefault(name, []).append(float(cl.std() / max(cl.mean(), 1e-9)))
        # known-answer: quasi's own delivered totals must be floor(n L p) exactly
        assert np.array_equal(ref_T, R.totals(1.0)), f"quasi != floor(nLp) on block {bi}"
        print(f"  block {bi} {str(D['targets'][bi]):<10} " + "  ".join(
            "{}:{}".format(x["law"], "SAME" if x["same_as_quasi"]
                           else "DIFF({})".format(x["n_gene_diff"]))
            for x in rows[-4:]), flush=True)

    print(f"\n=== DELIVERED PER-GENE PANEL TOTALS, {N_BLOCKS} real challenge blocks ===")
    print(f"{'law':<12}{'digest == quasi':>17}{'genes differing':>17}{'max|dT|':>9}"
          f"{'relL1':>11}{'emitted lib CV':>16}{'frac resampled':>16}")
    summary = {}
    for law in ("quasi", "a8_phi0", "a21_x0.14", "a22_inv"):
        rs = [x for x in rows if x["law"] == law]
        s = dict(same=sum(x["same_as_quasi"] for x in rs), n=len(rs),
                 gene_diff=float(np.mean([x["n_gene_diff"] for x in rs])),
                 max_dT=int(max(x["max_abs_dT"] for x in rs)),
                 relL1=float(np.mean([x["relL1"] for x in rs])),
                 lib_cv=float(np.mean([x["lib_cv"] for x in rs])),
                 frac_resampled=float(np.mean([x["frac_resampled"] for x in rs])))
        summary[law] = s
        print(f"{law:<12}{s['same']:>10}/{s['n']:<6}{s['gene_diff']:>17.1f}"
              f"{s['max_dT']:>9}{s['relL1']:>11.3e}{s['lib_cv']:>16.4f}"
              f"{s['frac_resampled']:>16.4f}")

    ok = all(summary[l]["same"] == summary[l]["n"] for l in summary)
    json.dump(dict(rows=rows, summary=summary, all_identical=bool(ok),
                   secs=time.time() - t0), open(WORK / "emit230.json", "w"), indent=1)
    print(f"\n  every emission law delivered IDENTICAL per-gene panel totals: {ok}")
    print(f"  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
