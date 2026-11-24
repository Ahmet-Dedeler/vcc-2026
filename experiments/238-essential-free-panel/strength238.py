"""238 — the confound check 233 §6 ran, repeated for the essential-free draws.

Removing the common-essential list from the eligible pool does not only remove the rows the
October sources carry.  Essential genes have LARGE effects, so an essential-free panel could
differ from 199's in effect STRENGTH and in model-free cross-line TRANSFERABILITY as well --
either of which would move the retention number for reasons that have nothing to do with the
October asset.  233 §6 measured both on its two panels and found strength "not it" and
transferability "this is it"; the same two statistics are measured here so the essential-free
rows can be read against the same yardstick.

  ||truth||                 the held-out line's own effect norm per target -- effect STRENGTH
  ceiling  cos(s_i, s_j)    the same perturbation's agreement between the two SOURCE lines,
                            neither of which is the truth or the prediction.  It uses no model,
                            so it is a property of the panel: how cross-line reproducible its
                            targets are.  233 measured Arc's panel at 1.3-2.0x a random draw.

Runs off the fold caches `six238.py` already built, so it costs one npz read per fold.

    uv run python strength238.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import bind238  # noqa: E402,F401
import panels238 as P  # noqa: E402
import fold233  # noqa: E402

ESSENTIAL = {ln.strip() for ln in
             (P.ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines()
             if ln.strip()}
PANELS = ("abc", "real", "pess", "efr", "efp")


def main() -> None:
    out: dict = {}
    print(f"{'panel':6}{'fold':16}{'n':>5}{'||truth|| mean':>15}{'median':>10}"
          f"{'ceiling cos':>13}{'ceil (ess)':>12}{'ceil (non-ess)':>16}")
    for panel in PANELS:
        if not P.BUNDLES[panel].exists():
            continue
        for held in P.CANCER_FOLDS:
            cache = P.WORK / "folds" / f"{panel}__{held}.npz"
            if not cache.exists():
                continue
            fold = fold233.build(panel, held)
            sc = np.asarray(fold["scored"])
            T = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
            norms = np.linalg.norm(T, axis=1)
            srcs = fold["sources"]
            A = np.asarray(fold[f"eff__{srcs[0]}"], dtype=np.float64)[:, sc]
            B = np.asarray(fold[f"eff__{srcs[1]}"], dtype=np.float64)[:, sc]
            av = np.asarray(fold[f"avail__{srcs[0]}"], bool) & \
                np.asarray(fold[f"avail__{srcs[1]}"], bool)
            na = np.linalg.norm(A, axis=1)
            nb = np.linalg.norm(B, axis=1)
            ok = av & (na > 0) & (nb > 0)
            ceil = np.full(len(na), np.nan)
            ceil[ok] = (A[ok] * B[ok]).sum(1) / (na[ok] * nb[ok])
            ess = np.array([t in ESSENTIAL for t in fold["targets"]], dtype=bool)
            row = dict(
                n=int(T.shape[0]), n_essential=int(ess.sum()),
                truth_norm_mean=float(norms.mean()), truth_norm_median=float(np.median(norms)),
                truth_norm_mean_essential=(float(norms[ess].mean()) if ess.any() else None),
                truth_norm_mean_nonessential=(float(norms[~ess].mean())
                                              if (~ess).any() else None),
                ceiling=float(np.nanmean(ceil)),
                ceiling_essential=(float(np.nanmean(ceil[ess])) if ess.any() else None),
                ceiling_nonessential=(float(np.nanmean(ceil[~ess])) if (~ess).any() else None),
                n_ceiling_defined=int(ok.sum()))
            out[f"{panel}|{held}"] = row
            ce = ("     —" if row["ceiling_essential"] is None
                  else f"{row['ceiling_essential']:12.5f}")
            cn = ("       —" if row["ceiling_nonessential"] is None
                  else f"{row['ceiling_nonessential']:16.5f}")
            print(f"{panel:6}{held:16}{row['n']:>5}{row['truth_norm_mean']:15.4f}"
                  f"{row['truth_norm_median']:10.4f}{row['ceiling']:13.5f}{ce}{cn}")
        print()

    # panel-level means over the folds where the statistic is defined
    print(f"{'panel':6}{'mean ||truth||':>16}{'mean ceiling':>15}"
          f"{'ceiling vs abc':>17}")
    abc_c = np.mean([v["ceiling"] for k, v in out.items() if k.startswith("abc|")]) \
        if any(k.startswith("abc|") for k in out) else float("nan")
    summ = {}
    for panel in PANELS:
        vs = [v for k, v in out.items() if k.startswith(f"{panel}|")]
        if not vs:
            continue
        t = float(np.mean([v["truth_norm_mean"] for v in vs]))
        c = float(np.mean([v["ceiling"] for v in vs]))
        summ[panel] = dict(truth_norm=t, ceiling=c, ceiling_ratio=c / abc_c)
        print(f"{panel:6}{t:16.4f}{c:15.5f}{100 * c / abc_c:16.1f}%")
    out["summary"] = summ
    (P.WORK / "strength238.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {P.WORK / 'strength238.json'}")


if __name__ == "__main__":
    main()
