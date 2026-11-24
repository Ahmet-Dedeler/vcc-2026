"""279 — POST-HOC: is the held-set effect a ZINC-FINGER effect or a TRANSCRIPTION-FACTOR effect,
and does it survive a bootstrap that respects the paralogue clustering?

Three questions, all written after the preregistered numbers were read and all labelled POST-HOC:

  1. **Cytoband block bootstrap.** The preregistered interval resamples GENES; KRAB-ZNF paralogues
     sit in tandem arrays (chiefly 19p12/19q13), so genes are not independent.  Resample HGNC
     cytobands instead (`location` in the HGNC table), taking every gene of a sampled band.
  2. **TF class vs ZNF class.** CollecTRI's curated TF catalogue on disk
     (`data/external/regulon/collectri_omnipath.tsv`, 1,211 TFs) splits the held set into TF and
     non-TF.  Is conservation lower for TFs in general, and — inside the TF class — still lower
     for the zinc fingers?
  3. **The stated mechanism.** The claim names KRAB-ZNFs (lineage-restricted TE repressors).
     Inside the HGNC C2H2 zinc-finger family, are the KRAB-domain members lower than the rest?

    PYTHONDONTWRITEBYTECODE=1 uv run python tf279.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/279-znf-state-dependence"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(HERE))
import analyse279 as A  # noqa: E402
import inv279 as I  # noqa: E402

SEED = 20260906


def hgnc_bands() -> dict[str, str]:
    out = {}
    with open(I.HGNC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["symbol"]] = (r["location"] or "").strip()
    return out


def block_boot(x_band, x, y_band, y, b=5000, seed=SEED):
    """Δ = median(x) − median(y), resampling cytobands (not genes) inside each group."""
    rng = np.random.default_rng(seed)

    def prep(band, v):
        u, inv = np.unique(band, return_inverse=True)
        groups = [np.flatnonzero(inv == k) for k in range(u.size)]
        return groups
    gx, gy = prep(x_band, x), prep(y_band, y)
    d = np.empty(b)
    for i in range(b):
        sx = np.concatenate([gx[k] for k in rng.integers(0, len(gx), len(gx))])
        sy = np.concatenate([gy[k] for k in rng.integers(0, len(gy), len(gy))])
        d[i] = np.median(x[sx]) - np.median(y[sy])
    lo, hi = np.percentile(d, [2.5, 97.5])
    return dict(n_bands_x=len(gx), n_bands_y=len(gy), d_median=float(np.median(x) - np.median(y)),
                ci=[float(lo), float(hi)], excludes_zero_negative=bool(hi < 0))


def main() -> None:
    zh = np.load(WORK / "held279.npz", allow_pickle=True)
    genes = [str(g) for g in zh["genes"]]
    c = np.asarray(zh["c_held"], np.float64)
    expr = np.asarray(zh["expr"], np.float64)
    mag = np.asarray(zh["mag"], np.float64)
    cls, groups = A.classify(genes)
    is_z, is_n = cls == "ZNF", cls == "NONZNF"
    fin = np.isfinite(c)
    de, dm = A.deciles(expr), A.deciles(mag)
    cell = np.where((de >= 0) & (dm >= 0), de * 10 + dm, -1)
    out: dict = {}

    # ---- 1. cytoband block bootstrap -----------------------------------------------------------
    bands = hgnc_bands()
    bd = np.asarray([bands.get(g, "") or f"__{g}" for g in genes], dtype=object)
    mz, mn = is_z & fin, is_n & fin
    out["cytoband_block_bootstrap"] = block_boot(bd[mz], c[mz], bd[mn], c[mn])
    r = out["cytoband_block_bootstrap"]
    print(f"1. cytoband block bootstrap: Δ {r['d_median']:+.5f} [{r['ci'][0]:+.5f}, "
          f"{r['ci'][1]:+.5f}] over {r['n_bands_x']} ZNF bands / {r['n_bands_y']} non-ZNF bands",
          flush=True)
    zbands = {b for b in bd[mz]}
    top = sorted(((int((bd[mz] == b).sum()), b) for b in zbands), reverse=True)[:6]
    out["top_znf_bands"] = [[n, b] for n, b in top]
    print(f"   the six largest ZNF cytobands: {top}", flush=True)

    # ---- 2. TF class (CollecTRI) ---------------------------------------------------------------
    tfs = set()
    with open(ROOT / "data/external/regulon/collectri_omnipath.tsv") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            tfs.add(row["source_genesymbol"])
    is_tf = np.asarray([g in tfs for g in genes])
    out["collectri"] = dict(n_tf_catalogue=len(tfs), n_tf_in_held=int((is_tf & fin).sum()),
                            n_tf_znf=int((is_tf & is_z & fin).sum()))
    print(f"2. CollecTRI TFs: {len(tfs)} in the catalogue, {int((is_tf & fin).sum())} in the held "
          f"set, of which {int((is_tf & is_z & fin).sum())} are HGNC zinc fingers", flush=True)
    tests = {
        "TF_vs_nonTF": (is_tf & fin & ~is_z, (~is_tf) & fin & ~is_z),          # TFs, ZNFs removed
        "ZNF_vs_TF_nonZNF": (is_z & fin, is_tf & fin & ~is_z),                 # inside the TF class
        "ZNF_vs_nonTF_nonZNF": (is_z & fin, (~is_tf) & fin & ~is_z),
    }
    out["tf_tests"] = {}
    for nm, (mx, my) in tests.items():
        t = A.two_sample(c[mx], c[my], b=2000)
        t["matched"] = A.matched_effect(c, cell, mx, my, n_draw=20, b=1000)
        t["block"] = block_boot(bd[mx], c[mx], bd[my], c[my], b=2000)
        out["tf_tests"][nm] = t
        print(f"   {nm:22s} n {t['n_x']:5d}/{t['n_y']:6d}  {t['median_x']:+.5f} vs "
              f"{t['median_y']:+.5f}  Δ {t['d_median']:+.5f} [{t['ci'][0]:+.5f}, {t['ci'][1]:+.5f}]"
              f"  AUC {t['auc']:.3f}; matched {t['matched']['d_median']:+.5f} "
              f"[{t['matched']['ci'][0]:+.5f}, {t['matched']['ci'][1]:+.5f}]", flush=True)

    # ---- 3. KRAB vs non-KRAB C2H2 zinc fingers --------------------------------------------------
    krab = np.asarray([("KRAB domain containing" in groups.get(g, set())) for g in genes])
    c2h2 = np.asarray([("Zinc fingers C2H2-type" in groups.get(g, set())) for g in genes])
    mk, mnk = krab & fin, c2h2 & (~krab) & fin
    t = A.two_sample(c[mk], c[mnk], b=2000)
    t["matched"] = A.matched_effect(c, cell, mk, mnk, n_draw=20, b=1000)
    t["block"] = block_boot(bd[mk], c[mk], bd[mnk], c[mnk], b=2000)
    out["krab_vs_other_c2h2"] = t
    print(f"3. KRAB ({t['n_x']}) vs non-KRAB C2H2 ({t['n_y']}): {t['median_x']:+.5f} vs "
          f"{t['median_y']:+.5f}, Δ {t['d_median']:+.5f} [{t['ci'][0]:+.5f}, {t['ci'][1]:+.5f}], "
          f"AUC {t['auc']:.3f}, p {t['mwu_p_one_sided']:.3f}", flush=True)

    (WORK / "tf279.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
