#!/usr/bin/env python3
"""708 floor diagnostic (science [d80a10], 10-16; the lab's 801 finding, relayed by the coordinator): how many fold targets
sit below Arc's basal-expression floor IN THE HELD-OUT LINE, per panel and fold, and how their fold truth differs.

Arc's panel rule (AGENTS 2.1): every target has high basal expression in the screened lines. In a 708 / 438 fold the held-out
public line plays the context being predicted, so an Arc-like fold scores only targets expressed in THAT line. Today
lib438.build_fold438 keeps every panel target the held bank perturbed, with no expression filter.

Floors (both reported; neither is fitted to any score):
  F_abs : held-line control CPM >= 13.4122, 801's A/B/C threshold (the A/B/C panel's 5th percentile of pseudobulk control
          CPM, mean over A/B/C). Held CPM is on the official axis, 1e6 x expm1(b_ctrl) / sum, as build_fold438's gate.
  F_rank: the same floor by rank. q = the share of official genes below 13.4122 CPM in A/B/C; a held target passes if its
          held CPM >= the q-quantile of held CPM over all official genes (genes the held bank does not measure count as 0).
  unmeasured: the target gene is not on the held bank's measured axis (b_ctrl 0 by mapping; GWPS measures 8,248 genes).
KA: A/B/C's panel 5th percentile reproduces 801's 13.412203349132618 CPM (|d| < 1e-6).
Reads 438's fold caches (work/<panel>/fold-<held>-held.npz) and the banks' b_ctrl (memory-mapped); streams A/B/C's control
files once (column totals only). Writes work/floor708.json and prints one table. Light (< 2 GB).
usage: uv run python experiments/708-def-panel-repricing/floor708.py

DAY MODE (708 AMENDMENT 4, 10-22): the floor F from the bundle's OWN controls and panel by 801's rule (the panel's 5th
percentile of pseudobulk control CPM, mean over the bundle's contexts), then each fold's floor-kept targets among those in
dump708-<tag>-<fold>.npz: held CPM >= F; a target the held bank does not measure counts as below.
    uv run python experiments/708-def-panel-repricing/floor708.py --bundle <B> --tag <T> [--dumps-dir DIR] [--out JSON]
        [--expect-thr X]      # known answer: 180's relabelled A/B/C bundle must give 801's 13.412203349132618
Writes work/floormask708-<tag>.json (default). ~1 min, < 1 GB.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import h5py  # noqa: E402
import numpy as np  # noqa: E402

import lib438 as L  # noqa: E402
from lib438 import E062  # noqa: E402

W438 = ROOT / "experiments/438-october-stack-repricing/work"
CH = ROOT / "data/challenge"
THR801 = 13.412203349132618
PANELS = ("abc", "abcdef", "real", "efp", "efr")
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")


def official() -> list[str]:
    with open(CH / "gene_names.csv") as f:
        return [r["gene_name"] for r in csv.DictReader(f)]


def context_cpm(path: Path, n: int) -> np.ndarray:
    """Pseudobulk CPM of one control file: 1e6 x per-gene total / all counts, streamed over X's CSR arrays (as draw801)."""
    tot = np.zeros(n, np.float64)
    with h5py.File(path, "r") as f:
        X = f["X"]
        assert X.attrs.get("encoding-type") in ("csr_matrix", b"csr_matrix"), X.attrs.get("encoding-type")
        idx, dat = X["indices"], X["data"]
        step = 1 << 23
        for a in range(0, idx.shape[0], step):
            tot += np.bincount(idx[a:a + step], weights=dat[a:a + step].astype(np.float64), minlength=n)
    return 1e6 * tot / tot.sum()


KD = {}


def own_kd(h: str, targets: list[str]) -> np.ndarray:
    """The held bank's delta at each target's OWN gene (its knockdown in h; NaN when h does not measure that gene). The
    fold truth zeroes the target's own column, so this reads the bank (memory-mapped) directly."""
    if h not in KD:
        b = L.bank(h)
        delta = np.load(b.delta_path, mmap_mode="r") if hasattr(b, "delta_path") else b.delta
        KD[h] = (b, delta, {g: i for i, g in enumerate(b.genes)}, {q: i for i, q in enumerate(b.perts)})
    b, delta, gpos, ppos = KD[h]
    return np.array([float(delta[ppos[t], gpos[t]]) if (t in ppos and t in gpos) else np.nan for t in targets])


def held_profile(h: str):
    """(CPM on the official axis, measured mask) of the held bank's controls: build_fold438's b_ctrl and gate arithmetic."""
    b, cols = L.plolo._b_ctrl_on_axis(L.bank(h))
    e = np.expm1(b)
    m = np.zeros(len(b), bool)
    m[cols] = True
    return 1e6 * e / e.sum(), m


def day(bundle: Path, tag: str, dumps: Path, out_path: Path, expect: float | None) -> None:
    """708 AMENDMENT 4 on a day bundle; see the module docstring."""
    t0 = time.time()
    genes = official()
    with open(bundle / "gene_names.csv") as f:
        assert [r["gene_name"] for r in csv.DictReader(f)] == genes, "the bundle's gene order differs from the official one"
    ctxs = list(json.loads((bundle / "manifest.json").read_text())["contexts"])
    files = [bundle / f"context_{c}.h5ad" for c in ctxs]
    assert ctxs and all(f.exists() for f in files), files
    with open(bundle / "pert_counts.csv") as f:
        panel = [r["target_gene"] for r in csv.DictReader(f) if r["target_gene"] != "non-targeting"]
    assert len(panel) == len(set(panel)) == 300, len(panel)
    pos = {g: i for i, g in enumerate(genes)}
    cpm = np.mean([context_cpm(f, len(genes)) for f in files], axis=0)
    thr = float(np.percentile(cpm[[pos[t] for t in panel]], 5))
    print(f"floor F = {thr:.6f} CPM (the panel's 5th percentile of pseudobulk control CPM, mean over {ctxs})", flush=True)
    if expect is not None:
        assert abs(thr - expect) < 1e-6, f"floor {thr!r} != the known answer {expect!r}"
        print(f"KA: F == {expect!r} (|d| {abs(thr - expect):.1e})", flush=True)
    rec = dict(tag=tag, bundle=str(bundle), contexts=ctxs, threshold_cpm=thr, known_answer=expect,
               rule="801's rule on the bundle: the panel's 5th percentile of pseudobulk control CPM (mean over its contexts); "
                    "a fold target is kept iff the held bank measures it and its held control CPM >= F", folds={})
    for h in FOLDS:
        tg = [str(t) for t in np.load(dumps / f"dump708-{tag}-{h}.npz")["targets"]]
        hc, hm = held_profile(h)
        rec["folds"][h] = dict(n=len(tg), kept=[t for t in tg if hm[pos[t]] and hc[pos[t]] >= thr],
                               below=[t for t in tg if hm[pos[t]] and hc[pos[t]] < thr],
                               unmeasured=[t for t in tg if not hm[pos[t]]])
        r = rec["folds"][h]
        print(f"{h}: {r['n']} targets, kept {len(r['kept'])}, below {len(r['below'])}, unmeasured {len(r['unmeasured'])}",
              flush=True)
    rec["seconds"] = round(time.time() - t0, 1)
    out_path.write_text(json.dumps(rec, indent=1))
    print(f"wrote {out_path}", flush=True)


def main() -> None:
    genes = official()
    pos = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    cpm_abc = np.mean([context_cpm(CH / f"context_{c}.h5ad", n) for c in "ABC"], axis=0)
    with open(CH / "pert_counts.csv") as f:
        abc_targets = [r["target_gene"] for r in csv.DictReader(f)]
    thr = float(np.percentile(cpm_abc[[pos[t] for t in abc_targets]], 5))
    assert abs(thr - THR801) < 1e-6, (thr, THR801)
    q = float(np.mean(cpm_abc < thr))
    print(f"KA: A/B/C panel 5th percentile {thr:.6f} CPM == 801's {THR801:.6f}; share of official genes below it in A/B/C "
          f"q = {q:.4f}", flush=True)
    held_cpm, held_meas = {}, {}
    for h in FOLDS:
        b, cols = L.plolo._b_ctrl_on_axis(L.bank(h))
        e = np.expm1(b)
        held_cpm[h] = 1e6 * e / e.sum()
        m = np.zeros(n, bool)
        m[cols] = True
        held_meas[h] = m
    out = dict(threshold_cpm=thr, q_abc=q, folds={})
    print(f"{'panel':8s} {'fold':14s} {'n':>4s} {'unmeas':>6s} {'<abs':>5s} {'<rank':>5s} {'rank thr':>9s} "
          f"{'med CPM':>8s} {'|truth| below/above (rank)':>27s}  | own-gene kd / cos(pre,truth) medians: unmeasured, "
          f"measured < abs, >= abs", flush=True)
    for p in PANELS:
        for h in FOLDS:
            fp = W438 / p / f"fold-{h}-held.npz"
            if not fp.exists():
                continue
            z = np.load(fp, allow_pickle=True)
            tg = [str(t) for t in z["targets"]]
            truth_full = np.asarray(z["truth"], np.float64)
            sc = np.asarray(z["scored"])
            truth = truth_full[:, sc]
            pre = np.asarray(z["pre"], np.float64)[:, sc]
            norm = np.linalg.norm(truth - truth.mean(axis=0), axis=1)
            kd = own_kd(h, tg)                                        # the target's own log-effect in h (bank delta)
            with np.errstate(all="ignore"):
                cos = np.einsum("ij,ij->i", pre, truth) / (np.linalg.norm(pre, axis=1) * np.linalg.norm(truth, axis=1))
            tc = np.array([held_cpm[h][pos[t]] for t in tg])
            meas = np.array([held_meas[h][pos[t]] for t in tg])
            rthr = float(np.quantile(held_cpm[h], q))
            below_abs = tc < thr
            below_rank = tc < rthr
            nb = norm[below_rank]
            na = norm[~below_rank]
            row = dict(n=len(tg), unmeasured=int((~meas).sum()), below_abs=int(below_abs.sum()),
                       below_rank=int(below_rank.sum()), rank_thr_cpm=rthr, median_cpm=float(np.median(tc)),
                       truth_norm_median_below_rank=float(np.median(nb)) if len(nb) else None,
                       truth_norm_median_above_rank=float(np.median(na)) if len(na) else None,
                       below_rank_targets=[t for t, x in zip(tg, below_rank) if x],
                       below_abs_targets=[t for t, x in zip(tg, below_abs) if x])
            groups = {"unmeasured": ~meas, "measured_below_abs": meas & below_abs, "above_abs": ~below_abs}
            row["groups"] = {k: dict(n=int(g.sum()),
                                     kd_median=float(np.nanmedian(kd[g])) if np.isfinite(kd[g]).any() else None,
                                     kd_n=int(np.isfinite(kd[g]).sum()),
                                     kd_share_below_m0p5=float(np.mean(kd[g][np.isfinite(kd[g])] < -0.5)) if np.isfinite(kd[g]).any() else None,
                                     cos_median=float(np.nanmedian(cos[g])) if g.any() else None)
                             for k, g in groups.items()}
            out["folds"][f"{p}/{h}"] = row
            nbs = f"{np.median(nb):.3f}" if len(nb) else "-"
            gtxt = "  ".join((f"{v['n']}:" + (f"{v['kd_median']:+.2f}({v['kd_share_below_m0p5']:.2f})" if v["kd_median"] is not None else "na")
                              + f"/{v['cos_median']:+.3f}") if v["n"] else "0:-" for v in row["groups"].values())
            print(f"{p:8s} {h:14s} {len(tg):4d} {row['unmeasured']:6d} {row['below_abs']:5d} {row['below_rank']:5d} "
                  f"{rthr:9.3f} {row['median_cpm']:8.2f} {nbs:>13s} / {np.median(na):.3f}  | {gtxt}", flush=True)
    (HERE / "work" / "floor708.json").write_text(json.dumps(out, indent=1))
    print("wrote work/floor708.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=Path)
    ap.add_argument("--tag")
    ap.add_argument("--dumps-dir", type=Path, default=HERE / "work")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--expect-thr", type=float)
    a = ap.parse_args()
    if a.bundle is None:
        main()
    else:
        assert a.tag, "--tag is required with --bundle"
        day(a.bundle, a.tag, a.dumps_dir, a.out or HERE / "work" / f"floormask708-{a.tag}.json", a.expect_thr)
