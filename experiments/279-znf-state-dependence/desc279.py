"""279 item 5 — the reliability-corrected test and the noise-vs-state-dependence descriptor.

Reads work/held279.npz + work/rel279.npz.  Everything here is declared in PREREGISTRATION §3
(`C_rel`) and §7 (the descriptor, labelled POST-HOC in advance).

    PYTHONDONTWRITEBYTECODE=1 uv run python desc279.py

Writes work/desc279.json and appends the reliability columns to work/znf279.csv.
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/279-znf-state-dependence"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(HERE))
import analyse279 as A  # noqa: E402

SEED = 20260906


def main() -> None:
    t0 = time.time()
    zh = np.load(WORK / "held279.npz", allow_pickle=True)
    zr = np.load(WORK / "rel279.npz", allow_pickle=True)
    genes = [str(g) for g in zh["genes"]]
    c = np.asarray(zh["c_held"], np.float64)
    expr = np.asarray(zh["expr"], np.float64)
    mag = np.asarray(zh["mag"], np.float64)
    cls, groups = A.classify(genes)
    is_z, is_n = cls == "ZNF", cls == "NONZNF"
    de, dm = A.deciles(expr), A.deciles(mag)
    cell = np.where((de >= 0) & (dm >= 0), de * 10 + dm, -1)

    c_rel = np.asarray(zr["c_rel"], np.float64)
    c_rel3 = np.asarray(zr["c_rel3"], np.float64)
    c_full4 = np.asarray(zr["c_full4"], np.float64)
    rbar = np.asarray(zr["rbar"], np.float64)
    rkeys = [str(k) for k in zr["keys"]]
    out: dict = {}

    # ---- 1. is the ZNF response measured less reliably? ---------------------------------------
    out["reliability_by_source"] = {}
    for k in rkeys:
        r = np.asarray(zr[f"r_half|{k}"], np.float64)
        f = np.isfinite(r)
        x, y = r[is_z & f], r[is_n & f]
        out["reliability_by_source"][k] = dict(
            n_znf=int(x.size), n_nonznf=int(y.size), median_znf=float(np.median(x)),
            median_nonznf=float(np.median(y)),
            share_below_floor_znf=float(np.mean(x <= 0.02)),
            share_below_floor_nonznf=float(np.mean(y <= 0.02)),
            mwu_p_one_sided=float(mannwhitneyu(x, y, alternative="less").pvalue),
            auc=A.auc_lower(x, y))
        d = out["reliability_by_source"][k]
        print(f"  reliability {k:16s} ZNF {d['median_znf']:+.4f} vs non-ZNF "
              f"{d['median_nonznf']:+.4f}  AUC {d['auc']:.3f}  p {d['mwu_p_one_sided']:.1e}  "
              f"(<=0.02: {d['share_below_floor_znf']:.2f} vs {d['share_below_floor_nonznf']:.2f})",
              flush=True)

    # ---- 2. the reliability-corrected test, and its uncorrected comparator on the same genes ---
    for nm, s in (("C_full4_uncorrected", c_full4), ("C_rel", c_rel),
                  ("C_rel_no_hek293t", c_rel3)):
        f = np.isfinite(s)
        r = A.two_sample(s[is_z & f], s[is_n & f], b=2000)
        r["matched"] = A.matched_effect(s, cell, is_z & f, is_n & f, n_draw=20, b=1000)
        out[nm] = r
        print(f"{nm:22s} n {r['n_x']:4d}/{r['n_y']:6d} ZNF {r['median_x']:+.4f} vs "
              f"{r['median_y']:+.4f}  Δ {r['d_median']:+.4f} [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}] "
              f"AUC {r['auc']:.3f}; matched {r['matched']['d_median']:+.4f} "
              f"[{r['matched']['ci'][0]:+.4f}, {r['matched']['ci'][1]:+.4f}]", flush=True)

    # ---- 3. the descriptor: terciles of measured reliability -----------------------------------
    f = np.isfinite(rbar) & np.isfinite(c)
    q = np.quantile(rbar[f], [1 / 3, 2 / 3])
    terc = np.where(f, np.searchsorted(q, rbar, side="right"), -1)
    out["terciles"] = dict(cuts=[float(v) for v in q], rows=[])
    for t in range(3):
        m = f & (terc == t)
        x, y = c[is_z & m], c[is_n & m]
        row = dict(tercile=t, n_znf=int(x.size), n_nonznf=int(y.size),
                   r_median=float(np.median(rbar[m])),
                   median_znf=float(np.median(x)), median_nonznf=float(np.median(y)),
                   d_median=float(np.median(x) - np.median(y)), auc=A.auc_lower(x, y),
                   mwu_p_one_sided=float(mannwhitneyu(x, y, alternative="less").pvalue))
        row["boot"] = A.two_sample(x, y, b=2000)["ci"]
        row["matched"] = A.matched_effect(c, cell, is_z & m, is_n & m, n_draw=20, b=1000)
        out["terciles"]["rows"].append(row)
        print(f"  tercile {t} (r~{row['r_median']:+.3f}, n {row['n_znf']}/{row['n_nonznf']}): "
              f"Δ {row['d_median']:+.5f} {[round(v, 5) for v in row['boot']]} AUC {row['auc']:.3f}; "
              f"matched {row['matched']['d_median']:+.5f} "
              f"[{row['matched']['ci'][0]:+.5f}, {row['matched']['ci'][1]:+.5f}]", flush=True)

    # ---- 4. how much of the raw ZNF deficit does reliability explain? ---------------------------
    #        the same contrast on the genes where BOTH C_full4 and C_rel exist
    both = np.isfinite(c_full4) & np.isfinite(c_rel)
    d_unc = float(np.median(c_full4[is_z & both]) - np.median(c_full4[is_n & both]))
    d_cor = float(np.median(c_rel[is_z & both]) - np.median(c_rel[is_n & both]))
    scale = float(np.median(c_rel[is_n & both]) / np.median(c_full4[is_n & both]))
    out["explained"] = dict(n=int(both.sum()), d_uncorrected=d_unc, d_corrected=d_cor,
                            nonznf_scale_factor=scale,
                            d_uncorrected_scaled=d_unc * scale,
                            residual_after_correction=d_cor - d_unc * scale)
    print(f"EXPLAINED: on the {int(both.sum())} genes with both scores, uncorrected Δ {d_unc:+.5f}; "
          f"the correction multiplies the non-ZNF median by {scale:.2f}×, so a pure-noise "
          f"explanation predicts Δ {d_unc*scale:+.5f}; measured corrected Δ {d_cor:+.5f}", flush=True)

    # ---- append the reliability columns to the per-gene table -----------------------------------
    src = list(csv.reader(open(WORK / "znf279.csv")))
    hdr = src[0] + ["r_bar", "C_full4", "C_rel", "rel_tercile"]
    rows = [hdr]
    for i, r in enumerate(src[1:]):
        rows.append(r + [f"{rbar[i]:.6f}" if np.isfinite(rbar[i]) else "",
                         f"{c_full4[i]:.6f}" if np.isfinite(c_full4[i]) else "",
                         f"{c_rel[i]:.6f}" if np.isfinite(c_rel[i]) else "",
                         str(int(terc[i]))])
    with open(WORK / "znf279.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(rows)

    out["seconds"] = time.time() - t0
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
    (WORK / "desc279.json").write_text(json.dumps(out, indent=1))
    print(f"done in {out['seconds']:.0f}s, peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
