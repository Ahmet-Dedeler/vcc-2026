"""179 gate 5 — the zero-effect null, driven through the BUILDER'S OWN emission path.

WHY IT IS RUN AGAIN HERE
------------------------
177 section 2 measured this property with a hand-written call into `emit_counts`. The gate the
package needs is different: it is that the property survives the path the build actually takes
-- `write_submission_from_effects` -> `profile_from_effect` -> `emission_spread_vector` ->
`emission_param_mask` -> `emit_counts`, with the control block drawn by the builder's own
`pick_rng` and passed through `assert_raw_control_counts`. A call-site mistake (the wrong
control profile, a rescaled block, an inverted mask) is invisible to 177's measurement and
fatal to the arm.

CONSTRUCTION
------------
The builder is handed an effect matrix of ZEROS at the a23 configuration, so `profile_from_
effect` renders exactly the context's own control profile and the predicted fold change is 1 on
every gene. Anything the emitted cells then call significant against the context's 18,400 real
controls is manufactured by the EMISSION LAW alone. Both laws are built from the same zero
effect, so the comparison is like-for-like.

KNOWN-ANSWER ROW (AGENTS.md 11.1 rule 10)
-----------------------------------------
The `quasi` arm is the shipped law and its answer is already known from two independent places:
157's closed form and 177 section 2 read **43.3% / 44.8% / 38.4%** BH-significant on contexts
A / B / C at `tail 1.0 / keep_top 8000 / keep_spread 0.25 / phi 4`, ~100% of the calls DOWN. If
this run does not land there, the instrument is wrong and the resample row means nothing.
The resample row's expectation from 177 is **0.0-0.2%** (1 / 18 / 19 genes).

Scoring is the scorer's own DE universe: mean per-cell CPM > 5 on the context's controls, exact
mid-rank Mann-Whitney of the 400 emitted cells against all 18,400 controls, BH at 0.05.

    uv run python null179.py
"""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import false_discovery_control, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))

# a23's environment, so the module-level constants are the incumbent's.
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

spec = importlib.util.spec_from_file_location(
    "builder179", ROOT / "scripts/make_submission_pm_cd4.py"
)
M = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = M
spec.loader.exec_module(M)

BUNDLE = ROOT / "data/challenge"
MIN_CPM_SCORER = 5.0
SHIPPED = dict(min_cpm=10.0, tail_spread=1.0, keep_top=8000, keep_spread=0.25,
               overdispersion=4.0)


class Ref:
    """Per-gene control CPM distribution; exact mid-rank Mann-Whitney of N emitted cells."""

    def __init__(self, csc, gated, m, n_emit):
        self.m = m
        self.nz = [np.sort(csc.data[csc.indptr[g]:csc.indptr[g + 1]]) for g in gated]
        self.n_zero = [m - a.size for a in self.nz]
        self.fac = np.sqrt(12.0 * n_emit * m / (m + n_emit + 1.0))

    def z(self, emitted):
        out = np.empty(emitted.shape[1])
        for k in range(emitted.shape[1]):
            e = emitted[:, k]
            a, nz0 = self.nz[k], self.n_zero[k]
            lo = np.searchsorted(a, e, "left")
            hi = np.searchsorted(a, e, "right")
            below = np.where(e == 0.0, 0, nz0 + lo)
            eq = np.where(e == 0.0, nz0, hi - lo)
            out[k] = ((below + 0.5 * eq) / self.m).mean()
        return (out - 0.5) * self.fac


def emitted_block(path: Path, n_genes: int):
    import h5py

    with h5py.File(path, "r") as handle:
        data = handle["X"]["data"][:]
        indices = handle["X"]["indices"][:]
        indptr = handle["X"]["indptr"][:]
        n_rows = int(handle["X"].attrs["shape"][0])
    block = np.zeros((n_rows, n_genes))
    for r in range(n_rows):
        sl = slice(int(indptr[r]), int(indptr[r + 1]))
        block[r, indices[sl]] = data[sl]
    return block


def main() -> None:
    t0 = time.time()
    genes = np.asarray(
        [line.split(",")[0] for line in
         (BUNDLE / "gene_names.csv").read_text().splitlines()[1:]]
    )
    assert genes.size == M.N_GENES, genes.size
    # one perturbation per context; the target's identity is irrelevant at zero effect, so a
    # single real panel symbol is used and the effect handed to the builder is exactly zero.
    target = "ADNP" if "ADNP" in set(genes) else str(genes[0])
    rows = []
    for ctx in ("A", "B", "C"):
        with tempfile.TemporaryDirectory(dir=WORK) as directory:
            root = Path(directory)
            bundle = root / "bundle"
            bundle.mkdir()
            os.symlink(BUNDLE / f"context_{ctx}.h5ad", bundle / f"context_{ctx}.h5ad")
            spec_ctx = M.ChallengeSpec(
                genes=genes, perts=np.asarray([target]), contexts=[ctx],
                cells_per_pert=400, manifest={},
            )
            effects = np.zeros((1, M.N_GENES))
            # Two of these are pure known-answer rows with NO keep block, so their expected
            # value is unambiguous (177 section 1's gate table, measured against 157's
            # closed form): flat 0.25/phi4 -> 19.3 / 21.8 / 14.8 %, flat 1.00/phi4 ->
            # 77.3 / 72.1 / 74.2 %. The shipped row carries a keep block whose membership
            # depends on which CPM gate defines the axis, so it is reported as the arm's
            # own baseline rather than pinned to a third-decimal expectation.
            arms = (
                ("KA_flat_0.25_phi4",
                 dict(min_cpm=10.0, tail_spread=0.25, keep_top=0, overdispersion=4.0)),
                ("KA_flat_1.00_phi4",
                 dict(min_cpm=10.0, tail_spread=1.0, keep_top=0, overdispersion=4.0)),
                ("quasi", dict(SHIPPED)),
                ("resample_mix",
                 dict(SHIPPED, emitter="resample_mix", resample_keep_top=8000)),
            )
            paths = {}
            for law, kwargs in arms:
                out = root / f"{law}.h5ad"
                sink = io.StringIO()
                with contextlib.redirect_stdout(sink):
                    M.write_submission_from_effects(
                        spec_ctx, effects, bundle, out, max_rss_gib=0.0, **kwargs,
                    )
                paths[law] = out
                print(f"  [{ctx}] emitted {law} [{time.time() - t0:.0f}s]", flush=True)

            import anndata as ad
            controls = ad.read_h5ad(BUNDLE / f"context_{ctx}.h5ad")
            X = controls.X.tocsr().astype(np.float64)
            m = X.shape[0]
            lib = np.asarray(X.sum(1)).ravel()
            lib[lib == 0] = 1.0
            csc = (sp.diags(1e6 / lib) @ X).tocsc()
            cpm = np.asarray(csc.mean(0)).ravel()
            gated = np.flatnonzero(cpm > MIN_CPM_SCORER)
            ref = Ref(csc, gated, m, 400)
            del csc, X, controls

            for law, path in paths.items():
                block = emitted_block(path, M.N_GENES)
                total = block.sum(1)
                total[total == 0] = 1.0
                emitted = block[:, gated] * (1e6 / total[:, None])
                z = ref.z(emitted)
                pv = 2 * norm.sf(np.abs(z))
                q = false_discovery_control(pv, method="bh")
                sig = q < 0.05
                rows.append(dict(
                    context=ctx, law=law, gated=int(gated.size),
                    n_sig=int(sig.sum()), rate=float(sig.mean()),
                    frac_up=float((z[sig] > 0).mean()) if sig.any() else float("nan"),
                    median_abs_z=float(np.median(np.abs(z))),
                    median_Fbar=float(np.median(z / ref.fac + 0.5)),
                ))
                print(f"  [{ctx}] {law:<13} gated {gated.size}  BH sig {int(sig.sum()):>5} "
                      f"({sig.mean():.2%})  {rows[-1]['frac_up']:.1%} UP  "
                      f"med|z| {rows[-1]['median_abs_z']:.2f}", flush=True)
                del block, emitted
            del ref

    print("\n=== gate 5: zero-effect null on the built emission path ===")
    print(f"{'ctx':>4}{'law':>14}{'gated':>8}{'BH sig':>9}{'rate':>9}{'UP':>8}{'med|z|':>9}")
    for r in rows:
        print(f"{r['context']:>4}{r['law']:>14}{r['gated']:>8}{r['n_sig']:>9}"
              f"{r['rate']:>9.2%}{r['frac_up']:>8.1%}{r['median_abs_z']:>9.2f}")
    expect = {
        "KA_flat_0.25_phi4": {"A": 0.193, "B": 0.218, "C": 0.148},
        "KA_flat_1.00_phi4": {"A": 0.773, "B": 0.721, "C": 0.742},
    }
    print("\nKNOWN-ANSWER rows (177 section 1's gate table, from 157's closed form)")
    ok = True
    for law, wants in expect.items():
        got = {r["context"]: r["rate"] for r in rows if r["law"] == law}
        for ctx, want in wants.items():
            have = got.get(ctx, float("nan"))
            good = abs(have - want) < 0.04
            ok &= bool(good)
            print(f"  {law:<18} ctx {ctx}: expected {want:.1%}  measured {have:.1%}  "
                  f"{'OK' if good else 'DRIFTED'}")
    print(f"known-answer gate: {'PASS' if ok else 'FAIL'}")
    (WORK / "null179.json").write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
