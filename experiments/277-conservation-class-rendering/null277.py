"""277 -- the PROPER null for the row selector: V1's construction on 24 size-matched RANDOM row sets.

The three seeds the frozen definition asks for are enough to say "the control did not beat the arm"; they
are NOT enough to say whether the conservation map picks rows where un-deflating the CD4 residual costs
more than an arbitrary row set does.  This file answers that: the identical construction (un-deflated CD4
on the selected rows, deflated elsewhere, the CD4 fraction re-solved once) on 24 random row sets of the
same size per fold, then the observed CLASS30 / AGREE30 / TOPCONS30 rows read as percentiles of that null.

    PYTHONDONTWRITEBYTECODE=1 uv run python null277.py --fold <fold>

Writes `work/null277-<fold>.json`.  KA: the seed-20260906 draw reproduces `render277`'s `C_rand1` row at
max|d| == 0.0.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib277 as Z  # noqa: E402

L = Z.L
N_DRAWS = 24


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    st = L.stage(held)
    targets = [str(t) for t in st["targets"]]
    n = len(targets)
    W = np.asarray(st["whitened"], dtype=np.float64)
    cd4 = st["cd4"]
    d, _ = L.deflate(cd4, "svd", 5)
    fr0, _ = L.solve_fraction(W, d, Z.CD4_F)
    base_parts = [(d, fr0), (st["union"], st["frac_x"])]
    e0 = L.score(st, L.compose(st, base_parts))
    p0 = np.asarray(e0["_pds_per_target_a23"], dtype=np.float64)

    cl = Z.c_lolo(held, targets)
    zp = np.load(Z.E276 / "work/predict276.npz", allow_pickle=True)
    dd_map = {str(t): float(v) for t, v in zip(zp[f"{held}|targets"], zp[f"{held}|dd"])}
    dd = np.asarray([dd_map.get(t, np.nan) for t in targets], dtype=np.float64)
    C30 = Z.bottom_class(cl, targets, 0.30)
    n_class = int(C30.sum())
    named = {"CLASS30": C30, "AGREE30": Z.bottom_class(dd, targets, 0.30),
             "AGREEN": Z.bottom_class(dd, targets, n_class / max(int(np.isfinite(dd).sum()), 1)),
             "TOPCONS30": Z.bottom_class(cl, targets, 0.30, top=True)}
    print(f"[{held}] {n} targets, CLASS30 {n_class}; V0 pds {e0['pds']:.6f} cos_agg {e0['cos_agg']:.6f} "
          f"({time.time()-t0:.0f}s, peak {Z.peak_gib():.2f} GiB)", flush=True)

    def run(mask, tag):
        m = Z.mix_rows(cd4, d, mask, tag)
        fr, _ = L.solve_fraction(W, m, Z.CD4_F)
        e = L.score(st, L.compose(st, [(m, fr), (st["union"], st["frac_x"])]))
        p = np.asarray(e["_pds_per_target_a23"], dtype=np.float64)
        dif = p - p0
        return dict(tag=tag, n_rows=int(mask.sum()), d_pds=float(dif.mean()),
                    per_row=float(dif[mask].sum() / max(int(mask.sum()), 1)),
                    rel_agg=float(e["cos_agg"] / e0["cos_agg"] - 1.0),
                    d_pds_on=float(dif[mask].mean()), pds=float(e["pds"]))

    obs = {k: run(m, k) for k, m in named.items()}
    draws = []
    for i in range(N_DRAWS):
        seed = 20260906 + i
        r = run(Z.random_class(cl, n_class, seed), f"rand{seed}")
        r["seed"] = seed
        draws.append(r)
        print(f"  draw {i+1:2d}/{N_DRAWS} seed {seed}: d pds {r['d_pds']:+.5f}, per edited row "
              f"{r['per_row']:+.5f}", flush=True)

    ka = abs(draws[0]["d_pds"] - (json.loads((Z.WORK / f"render277-{held}.json").read_text())
                                  ["arms"]["C_rand1"]["pds"] - e0["pds"]))
    out = dict(fold=held, n_targets=n, n_class=n_class, v0=dict(pds=e0["pds"], cos_agg=e0["cos_agg"]),
               observed=obs, draws=draws,
               ka_rand1_vs_render277=float(ka),
               null=dict(d_pds=dict(mean=float(np.mean([x["d_pds"] for x in draws])),
                                    sd=float(np.std([x["d_pds"] for x in draws], ddof=1)),
                                    min=float(np.min([x["d_pds"] for x in draws])),
                                    max=float(np.max([x["d_pds"] for x in draws]))),
                         per_row=dict(mean=float(np.mean([x["per_row"] for x in draws])),
                                      sd=float(np.std([x["per_row"] for x in draws], ddof=1)))),
               seconds=round(time.time() - t0, 1), peak_rss_gib=Z.peak_gib())
    for k, o in obs.items():
        lo = float(np.mean([x["d_pds"] <= o["d_pds"] for x in draws]))
        lo_pr = float(np.mean([x["per_row"] <= o["per_row"] for x in draws]))
        o["null_frac_at_or_below_d_pds"] = lo
        o["null_frac_at_or_below_per_row"] = lo_pr
        print(f"  {k:10s} d pds {o['d_pds']:+.5f} (null mean {out['null']['d_pds']['mean']:+.5f} sd "
              f"{out['null']['d_pds']['sd']:.5f}; {100*lo:.0f} % of draws at or below) | per edited row "
              f"{o['per_row']:+.5f} ({100*lo_pr:.0f} % at or below)", flush=True)
    print(f"  KA seed-20260906 draw == render277's C_rand1: |d| {ka:.1e} "
          f"(tol 1e-12; the two paths sum the same per-target array in a different order) "
          f"-> {'PASS' if ka < 1e-12 else 'FAIL'}", flush=True)
    (Z.WORK / f"null277-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    if ka >= 1e-12:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
