"""217 — is the delivered amplitude drop REAL, or an artefact of integer emission?

`profile210.py` measured the shipped artifacts at delivered median row norm **0.9937x** and
`A = mean|row|^2` **0.9797x** the incumbent's, against 207 section 6.2's PLOLO-frame prediction
of 1.0009x and 0.9904x.  Two things could produce that gap and they have opposite consequences:

  (a) the emitter's INTEGER quantisation, which adds a fixed noise energy `n^2` to every row of
      both artifacts.  It biases the RATIO towards 1, so it cannot manufacture a drop -- but if
      it were large it would mean the whole measurement is noise-dominated.  Analytic bound:
      the exact-total allocation leaves |dt| <= 1 count per gene, `dE_g = (TS/T) dt/(1 + TS p_g)`
      with `TS/T = 5e4/8e6 = 6.25e-3`, so `sd(dE_g) <= 1.8e-3/(1 + TS p_g)` and `n^2 ~ 0.03`
      against `A = 14.7`, i.e. **0.2%**.  Checked empirically below by restricting to genes whose
      counts are large enough that quantisation is negligible.
  (b) a genuine amplitude move in the rendered profile, which is what 207's instrument would
      have missed: `amp207.py` renders with `deliver(..., "minnorm")` -- experiment 156's
      MASS-FIX renderer -- while production ships `VCC_MASS_FIX` UNSET and uses the implicit
      renormalisation.  On 156's own numbers that correction is 4.0-4.5% of the row norm in a
      nearly-shared direction, so the two renderers need not agree on a 1% amplitude claim.

So: the same two statistics over four nested gene sets, from the least to the most
quantisation-immune.  If the ratio survives on the high-count sets, (a) is excluded.

    uv run python amp210.py <arm.h5ad> <incumbent.h5ad>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
BUNDLE = ROOT / "data/challenge"
N_GENES = 18533
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import TS  # noqa: E402

sys.path.insert(0, str(HERE))
from profile217 import control_frames, read_obs_blocks  # noqa: E402


def totals_of(path: Path, blocks, indptr) -> np.ndarray:
    out = np.zeros((len(blocks), N_GENES), dtype=np.float64)
    with h5py.File(path, "r") as handle:
        data, indices = handle["X"]["data"], handle["X"]["indices"]
        for k, (_key, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            out[k] = np.bincount(indices[lo:hi].astype(np.int64),
                                 weights=data[lo:hi].astype(np.float64),
                                 minlength=N_GENES)
    return out


def main() -> None:
    arm_path, base_path = Path(sys.argv[1]), Path(sys.argv[2])
    blocks, indptr, _ = read_obs_blocks(arm_path)
    b_blocks, b_indptr, _ = read_obs_blocks(base_path)
    keys = [k for k, _, _ in blocks]
    assert keys == [k for k, _, _ in b_blocks]
    a_tot = totals_of(arm_path, blocks, indptr)
    b_tot = totals_of(base_path, b_blocks, b_indptr)
    frames = control_frames()

    # pooled control CPM per context, for the "top-N by expression" sets
    ctrl_cpm = {}
    for ctx in ("A", "B", "C"):
        controls = ad.read_h5ad(BUNDLE / f"context_{ctx}.h5ad")
        pooled = np.asarray(controls.X.sum(0), dtype=np.float64).ravel()
        ctrl_cpm[ctx] = 1e6 * pooled / pooled.sum()
        del controls

    a_eff = np.empty_like(a_tot)
    b_eff = np.empty_like(b_tot)
    for k, (ctx, _t) in enumerate(keys):
        b_ctrl = frames[ctx]["b_ctrl"]
        a_eff[k] = np.log1p(TS * a_tot[k] / a_tot[k].sum()) - b_ctrl
        b_eff[k] = np.log1p(TS * b_tot[k] / b_tot[k].sum()) - b_ctrl

    # median incumbent panel total per (context, gene) -- the quantisation-immunity axis
    med_tot = {ctx: np.median(b_tot[[i for i, k in enumerate(keys) if k[0] == ctx]], axis=0)
               for ctx in ("A", "B", "C")}

    def masks(name):
        out = {}
        for ctx in ("A", "B", "C"):
            if name == "all genes":
                out[ctx] = np.ones(N_GENES, dtype=bool)
            elif name == "builder 10-CPM gate":
                out[ctx] = frames[ctx]["builder_gate"]
            elif name == "top 2,000 by control CPM":
                m = np.zeros(N_GENES, dtype=bool)
                m[np.argsort(-ctrl_cpm[ctx])[:2000]] = True
                out[ctx] = m
            elif name == "panel total >= 1,000 counts":
                out[ctx] = med_tot[ctx] >= 1000
            elif name == "panel total >= 10,000 counts":
                out[ctx] = med_tot[ctx] >= 10000
        return out

    names = ["all genes", "builder 10-CPM gate", "top 2,000 by control CPM",
             "panel total >= 1,000 counts", "panel total >= 10,000 counts"]
    print(f"{'gene set':<32}{'n (A/B/C)':>20}{'median|row| inc':>18}{'arm':>12}"
          f"{'ratio':>10}{'A inc':>12}{'arm':>12}{'ratio':>10}{'q-noise/A':>12}")
    report = {}
    for name in names:
        mk = masks(name)
        m = np.stack([mk[k[0]] for k in keys])
        ar = np.linalg.norm(np.where(m, a_eff, 0.0), axis=1)
        br = np.linalg.norm(np.where(m, b_eff, 0.0), axis=1)
        A_a, A_b = float(np.mean(ar ** 2)), float(np.mean(br ** 2))
        # analytic quantisation energy on this set, from the incumbent's own totals:
        #   sd(dE_g) = (TS/T) * sd(dt) / (1 + TS p_g),  sd(dt) = 1/sqrt(12)
        q = []
        for k, (ctx, _t) in enumerate(keys):
            T = b_tot[k].sum()
            p = b_tot[k] / T
            s = (TS / T) * (1.0 / np.sqrt(12.0)) / (1.0 + TS * p)
            q.append(float(np.sum((s ** 2)[mk[ctx]])))
        q_energy = float(np.mean(q))
        print(f"{name:<32}{'/'.join(str(int(mk[c].sum())) for c in 'ABC'):>20}"
              f"{np.median(br):>18.6f}{np.median(ar):>12.6f}"
              f"{np.median(ar) / np.median(br):>10.5f}"
              f"{A_b:>12.6f}{A_a:>12.6f}{A_a / A_b:>10.5f}{q_energy / A_b:>12.2%}")
        report[name] = dict(
            n=[int(mk[c].sum()) for c in "ABC"],
            median_inc=float(np.median(br)), median_arm=float(np.median(ar)),
            median_ratio=float(np.median(ar) / np.median(br)),
            A_inc=A_b, A_arm=A_a, A_ratio=float(A_a / A_b),
            quantisation_energy=q_energy, quantisation_frac_of_A=float(q_energy / A_b),
        )
    # 207 amp207.py's own conversion of an A move into raw mse
    r = report["all genes"]["A_ratio"]
    implied = (r - 1.0) * 0.2331583 / 0.229141
    print(f"\n207 amp207's conversion, applied to the SHIPPED artifacts:")
    print(f"  A ratio {r:.5f}  ->  implied d raw mse from the A channel alone "
          f"{implied:+.6f}")
    print(f"  207 predicted A ratio 0.99036  ->  implied {-0.009809283014426127:+.6f}")
    print(f"  difference (extra, UNCREDITED, favourable): "
          f"{implied - -0.009809283014426127:+.6f} raw")
    report["implied_d_raw_mse_A_channel"] = implied
    report["implied_d_raw_mse_A_channel_207"] = -0.009809283014426127
    (HERE / "work" / "amp217.json").write_text(json.dumps(report, indent=1))
    print(f"\nwrote {HERE / 'work' / 'amp217.json'}")


if __name__ == "__main__":
    main()
