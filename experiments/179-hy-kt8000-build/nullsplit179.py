"""179 gate 5, diagnosed — WHERE the zero-effect calls live, inside vs outside `param_mask`.

`null179.py` read 28.1 / 30.9 / 17.2 % BH-significant at zero predicted effect for
`--emitter resample_mix --resample-keep-top 8000`, against the 0.0-0.2 % the 179 mandate
quoted. That 0.0-0.2 % is 177 section 2's row for **`rs_all`** -- EVERY gene resampled. It was
never measured for the hybrid, and the hybrid cannot inherit it: 8,000 genes are still emitted
by the parametric quasi law at `keep_spread 0.25 / phi 4`, whose own zero-effect rate is the
`KA_flat_0.25_phi4` row (32.6 / 34.0 / 19.5 %).

So the falsifiable statement is not "the arm is null" but "the arm silenced the TAIL and left
the head alone". This splits the null by `param_mask` and tests exactly that:

    inside  `param_mask`  -> should match the quasi law's own rate on the same genes
    outside `param_mask`  -> should collapse to ~0, which is 177's `rs_all` number

If the outside rate did NOT collapse, the resampling is not working and the arm is dead. If the
inside rate moved, the parametric half is not the shipped law and the build is not this arm.

It also asserts the call-site contract the whole thing rests on, on the REAL rendered object:
at zero effect `profile_from_effect` must return the pooled control profile `expm1(b_ctrl)/sum`
bit-close, so `fold = 1` and the donated cells are exchangeable with the reference controls
(177 section 7 (b)).

    uv run python nullsplit179.py
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
sys.path.insert(0, str(ROOT / "src"))

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

_spec = importlib.util.spec_from_file_location(
    "builder179s", ROOT / "scripts/make_submission_pm_cd4.py"
)
M = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = M
_spec.loader.exec_module(M)

from null179 import Ref, emitted_block  # noqa: E402

BUNDLE = ROOT / "data/challenge"
MIN_CPM_SCORER = 5.0
SHIPPED = dict(min_cpm=10.0, tail_spread=1.0, keep_top=8000, keep_spread=0.25,
               overdispersion=4.0)


def main() -> None:
    t0 = time.time()
    genes = np.asarray(
        [line.split(",")[0] for line in
         (BUNDLE / "gene_names.csv").read_text().splitlines()[1:]]
    )
    rows = []
    for ctx in ("A", "B", "C"):
        import anndata as ad

        controls = ad.read_h5ad(BUNDLE / f"context_{ctx}.h5ad")
        X = controls.X.tocsr().astype(np.float64)
        m = X.shape[0]
        colsum = np.asarray(X.sum(0)).ravel()
        pooled = colsum / colsum.sum()
        b_ctrl = np.log1p(M.TS * pooled)
        lib = np.asarray(X.sum(1)).ravel()
        lib[lib == 0] = 1.0
        csc = (sp.diags(1e6 / lib) @ X).tocsc()
        cpm = np.asarray(csc.mean(0)).ravel()
        gated = np.flatnonzero(cpm > MIN_CPM_SCORER)
        ref = Ref(csc, gated, m, 400)
        del csc, X, controls

        zero = np.zeros(M.N_GENES)
        # contract (b), on the real rendered object
        profile = M.profile_from_effect(zero, b_ctrl, min_cpm=10.0)
        drift = float(np.max(np.abs(profile - pooled)))
        rel = float(np.max(np.abs(profile - pooled) / np.maximum(pooled, 1e-12)))
        pmask = M.emission_param_mask(zero, b_ctrl, min_cpm=10.0, resample_keep_top=8000)
        print(f"[{ctx}] rendered profile vs pooled control profile: max abs {drift:.3e}, "
              f"max rel {rel:.3e}   param_mask {int(pmask.sum())} genes", flush=True)

        with tempfile.TemporaryDirectory(dir=WORK) as directory:
            root = Path(directory)
            bundle = root / "bundle"
            bundle.mkdir()
            os.symlink(BUNDLE / f"context_{ctx}.h5ad", bundle / f"context_{ctx}.h5ad")
            spec_ctx = M.ChallengeSpec(
                genes=genes, perts=np.asarray(["ADNP"]), contexts=[ctx],
                cells_per_pert=400, manifest={},
            )
            out = {}
            for law, extra in (("quasi", {}),
                               ("resample_mix",
                                {"emitter": "resample_mix", "resample_keep_top": 8000})):
                path = root / f"{law}.h5ad"
                with contextlib.redirect_stdout(io.StringIO()):
                    M.write_submission_from_effects(
                        spec_ctx, zero[None, :], bundle, path, max_rss_gib=0.0,
                        **SHIPPED, **extra,
                    )
                block = emitted_block(path, M.N_GENES)
                total = block.sum(1)
                total[total == 0] = 1.0
                z = ref.z(block[:, gated] * (1e6 / total[:, None]))
                q = false_discovery_control(2 * norm.sf(np.abs(z)), method="bh")
                out[law] = (q < 0.05, z)
                del block

        inside = pmask[gated]
        for law, (sig, z) in out.items():
            r = dict(
                context=ctx, law=law, gated=int(gated.size),
                n_in=int(inside.sum()), n_out=int((~inside).sum()),
                sig_in=int(sig[inside].sum()), sig_out=int(sig[~inside].sum()),
                rate_in=float(sig[inside].mean()), rate_out=float(sig[~inside].mean()),
                up_out=float((z[~inside][sig[~inside]] > 0).mean())
                if sig[~inside].any() else float("nan"),
                med_abs_z_in=float(np.median(np.abs(z[inside]))),
                med_abs_z_out=float(np.median(np.abs(z[~inside]))),
                profile_drift_rel=rel,
            )
            rows.append(r)
            print(f"[{ctx}] {law:<13} IN mask {r['sig_in']:>5}/{r['n_in']} "
                  f"({r['rate_in']:.2%}, med|z| {r['med_abs_z_in']:.2f})   "
                  f"OUT {r['sig_out']:>5}/{r['n_out']} ({r['rate_out']:.2%}, "
                  f"med|z| {r['med_abs_z_out']:.2f})", flush=True)
        del ref

    print("\n=== gate 5, split by `param_mask` (zero predicted effect) ===")
    print(f"{'ctx':>4}{'law':>14}{'in mask':>9}{'rate in':>9}{'out':>8}{'rate out':>10}"
          f"{'UP out':>8}")
    for r in rows:
        print(f"{r['context']:>4}{r['law']:>14}{r['sig_in']:>9}{r['rate_in']:>9.2%}"
              f"{r['sig_out']:>8}{r['rate_out']:>10.2%}{r['up_out']:>8.1%}")
    print("\nWhat this reads if the resampling were ABSENT: `rate out` would equal the quasi "
          "row's, i.e. 60-80% at tail_spread 1.0 / phi 4.")
    (WORK / "nullsplit179.json").write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
