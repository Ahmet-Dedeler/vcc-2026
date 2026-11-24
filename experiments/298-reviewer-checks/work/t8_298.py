"""298 item 3 (T8) — F14's study-vs-construct contrast in BOTH statistics on the same K562 pair.

F14 quotes a RELIABILITY RATIO (across-study median cosine / within-study median cosine = 0.79 at
n = 25 cells per side).  Replogle et al. 2022's same-line cross-experiment number is a RAW
CORRELATION (~0.50).  The two are not the same statistic, so this reports both, side by side, on:

  A. the STUDY pair   — K562-GWPS vs the K562 essential screen, 2,053 shared perturbations,
                        from 285's SAVED full-depth pseudobulks (`k562_*_study_full.npy`).
  B. the CONSTRUCT pair — the two guide constructs of each of the 775 two-construct perturbations
                        of K562-GWPS (285's 748-target cohort at n = 25), at FULL depth.  These
                        per-construct full-depth sums were not saved by 285, so this makes ONE
                        streaming accumulation pass over the public K562-GWPS bank (group sums
                        only; the cells are never held, unlike 285's `k562_guide.py` at 11.8 GiB).

  C. both, restricted to the SAME perturbations (study-shared ∩ two-construct), which is where the
     inversion has to be read.

Response arithmetic COPIED by path from `experiments/285-reliability-benchmark/work/lib285.py`
(prof / delta / resid / row_cos / boot_median / boot_median_diff), so the residual convention is
literally 285's.  Pearson is `row_cos` after removing each row's own across-gene mean.

    PYTHONDONTWRITEBYTECODE=1 uv run python t8_298.py

Writes work/t8_298.json, work/t8_298.log.
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E285 = ROOT / "experiments/285-reliability-benchmark/work"
WORK = ROOT / "experiments/298-reviewer-checks/work"
sys.path.insert(0, str(E285))
import lib285 as L  # noqa: E402

GWPS = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
BLOCK = 8192
P = L.Logger(WORK / "t8_298.log")
T0 = time.time()


def pearson_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Per-row Pearson r across genes = cosine of the row-centred vectors."""
    return L.row_cos(A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True))


def summarise(x: np.ndarray) -> dict:
    return dict(**L.stats(x), boot=L.boot_median(x), mean_boot=None)


def line(tag: str, d: dict) -> str:
    return (f"  {tag:<34} n {d['n']:4d}  median {L.fmt(d['boot'])}  mean {d['mean']:+.4f}  "
            f"<=0 {d['share_le0']:.3f}")


def main() -> None:
    out: dict = {}

    # ================================================================== KA rows
    P(f"298 t8_298.py — {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    k = json.loads((E285 / "k562_285.json").read_text())
    g = json.loads((E285 / "k562_guide.json").read_text())
    ka = []
    s25 = k["study"]["25"]
    ok_a = (abs(s25["within_mean"]["median"] - 0.0605) < 5e-5
            and abs(s25["cross"]["median"] - 0.0479) < 5e-5
            and abs(s25["ratio"] - 0.791) < 5e-4)
    ka.append(dict(row="KA-1 285's K562 study row at n=25 (within / across / ratio)",
                   want=[0.0605, 0.0479, 0.791],
                   got=[s25["within_mean"]["median"], s25["cross"]["median"], s25["ratio"]],
                   ok=bool(ok_a)))
    P(f"KA-1  285 study n=25: within {s25['within_mean']['median']:+.4f} (want +0.0605)  "
      f"across {s25['cross']['median']:+.4f} (want +0.0479)  ratio {s25['ratio']:.4f} "
      f"(want 0.791) -> {'PASS' if ok_a else 'FAIL'}")
    g25 = g["rows"]["25"]
    ok_b = abs(g25["rho"]["point"] - 0.0896) < 5e-4 and g25["n_targets"] == 748
    ka.append(dict(row="KA-2 285's K562 construct row at n=25 (rho_guide, cohort size)",
                   want=[0.0896, 748], got=[g25["rho"]["point"], g25["n_targets"]], ok=bool(ok_b)))
    P(f"KA-2  285 construct n=25: rho_guide {g25['rho']['point']:+.4f} (want +0.0896) on "
      f"{g25['n_targets']} perturbations (want 748) -> {'PASS' if ok_b else 'FAIL'}")
    if not (ok_a and ok_b):
        raise SystemExit("KA FAILED — abort (AGENTS.md:1375-1382)")

    out["ka"] = ka
    out["ratio_285"] = {n: dict(n_targets=k["study"][n]["n_targets"],
                                within=k["study"][n]["within_mean"]["median"],
                                within_boot=k["study"][n]["within_mean"]["boot"],
                                across=k["study"][n]["cross"]["median"],
                                across_boot=k["study"][n]["cross"]["boot"],
                                ratio=k["study"][n]["ratio"]) for n in ("25", "50")}
    out["ratio_285_construct"] = {n: dict(n_targets=g["rows"][n]["n_targets"],
                                          within=g["rows"][n]["cells"]["median"],
                                          across=g["rows"][n]["guide"]["median"],
                                          ratio=g["rows"][n]["rho"]) for n in ("25", "50")}

    # ================================================================== A. the STUDY pair
    P("\n=== A. STUDY pair at FULL depth (K562-GWPS vs K562-essential, Replogle 2022) ===")
    ma = np.load(E285 / "k562_gwps_meta.npz", allow_pickle=True)
    mb = np.load(E285 / "k562_essential_meta.npz", allow_pickle=True)
    ga_, gb_ = ma["genes"].astype(str), mb["genes"].astype(str)
    shared = ma["shared"].astype(str)
    assert list(shared) == list(mb["shared"].astype(str))
    common_g = sorted(set(ga_.tolist()) & set(gb_.tolist()))          # k562_analyse.py:234-239
    pa = {x: i for i, x in enumerate(ga_)}
    pb = {x: i for i, x in enumerate(gb_)}
    ia = np.asarray([pa[x] for x in common_g], dtype=np.int64)
    ib = np.asarray([pb[x] for x in common_g], dtype=np.int64)
    tset = set(shared.tolist())
    keep_col = np.flatnonzero([x not in tset for x in common_g])
    nA, nB = np.asarray(ma["ncell_shared"]), np.asarray(mb["ncell_shared"])
    cA = np.asarray(ma["ctrl"], dtype=np.float64)[0][ia]
    cB = np.asarray(mb["ctrl"], dtype=np.float64)[0][ib]
    FA = np.asarray(np.load(E285 / "k562_gwps_study_full.npy", mmap_mode="r"),
                    dtype=np.float64)[:, ia]
    FB = np.asarray(np.load(E285 / "k562_essential_study_full.npy", mmap_mode="r"),
                    dtype=np.float64)[:, ib]
    P(f"  {len(shared)} shared perturbations; common axis {len(common_g)} genes "
      f"({keep_col.size} scored); cells/perturbation median GWPS {np.median(nA):.0f} / "
      f"essential {np.median(nB):.0f}; UMI/cell median {float(ma['umi_median']):.0f} / "
      f"{float(mb['umi_median']):.0f}")

    DA = L.delta(FA, cA)[:, keep_col]
    DB = L.delta(FB, cB)[:, keep_col]
    study = dict(raw_cos=L.row_cos(DA, DB), raw_pearson=pearson_rows(DA, DB),
                 resid_cos=L.row_cos(L.resid(DA), L.resid(DB)))
    mag_study = np.sqrt(np.linalg.norm(L.resid(DA), axis=1) * np.linalg.norm(L.resid(DB), axis=1))
    del FA, FB
    out["study_full_depth"] = {kk: summarise(v) for kk, v in study.items()}
    out["study_full_depth"]["n_cells"] = dict(gwps_median=float(np.median(nA)),
                                              essential_median=float(np.median(nB)))
    for kk, v in study.items():
        P(line(f"study {kk}", out["study_full_depth"][kk]))

    # the n=50 cohort 285 used, at full depth (comparability with the ratio row)
    coh50 = np.flatnonzero((nA >= 100) & (nB >= 100))
    out["study_full_depth_cohort50"] = {kk: summarise(v[coh50]) for kk, v in study.items()}
    P(f"  restricted to 285's n=50 cohort ({coh50.size} perturbations):")
    for kk in study:
        P(line(f"  study {kk}", out["study_full_depth_cohort50"][kk]))

    # ================================================================== B. the CONSTRUCT pair
    P("\n=== B. CONSTRUCT pair at FULL depth (K562-GWPS, two-construct perturbations) ===")
    with h5py.File(GWPS, "r") as h:
        o = h["obs/perturbation"]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in o["categories"][:]])
        pert = cats[o["codes"][:]]
        o = h["obs/guide_id"]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in o["categories"][:]])
        guide = cats[o["codes"][:]]
        n_cells, n_genes = h["X"].shape
    pg = defaultdict(set)
    for t, gg in zip(pert, guide):
        pg[t].add(gg)
    multi = sorted([t for t in pg if t != "control" and len(pg[t]) == 2])   # k562_guide.py:47-49
    assert len(multi) == 775, f"two-construct cohort is {len(multi)}, 285 saw 775"
    tpos = {t: i for i, t in enumerate(multi)}
    # group id: 2*i for the LARGER construct, 2*i+1 for the smaller; the last group is the control
    cnt = defaultdict(int)
    for t, gg in zip(pert, guide):
        if t in tpos:
            cnt[(t, gg)] += 1
    rank: dict[tuple[str, str], int] = {}
    for t in multi:
        gs = sorted(pg[t], key=lambda x: (-cnt[(t, x)], x))
        rank[(t, gs[0])] = 0
        rank[(t, gs[1])] = 1
    NG = 2 * len(multi) + 1
    grp = np.full(n_cells, -1, dtype=np.int32)
    for i in range(n_cells):
        t = pert[i]
        if t in tpos:
            grp[i] = 2 * tpos[t] + rank[(t, guide[i])]
        elif t == "control":
            grp[i] = NG - 1
    P(f"  {len(multi)} two-construct perturbations, {int((grp >= 0).sum()):,} cells assigned "
      f"({int((grp == NG - 1).sum()):,} controls); {n_cells:,} x {n_genes} matrix")

    acc = np.zeros((NG, n_genes), dtype=np.float64)
    with h5py.File(GWPS, "r") as h:
        X = h["X"]
        for r0 in range(0, n_cells, BLOCK):
            r1 = min(r0 + BLOCK, n_cells)
            gb = grp[r0:r1]
            m = gb >= 0
            if not m.any():
                continue
            Xb = np.asarray(X[r0:r1][m], dtype=np.float64)
            np.add.at(acc, gb[m], Xb)
            if (r0 // BLOCK) % 40 == 0:
                P(f"    rows {r0:,}/{n_cells:,}  {time.time() - T0:.0f}s  "
                  f"peak {L.peak_gib():.2f} GiB")
    ncc = np.bincount(grp[grp >= 0], minlength=NG)
    P(f"  pass done {time.time() - T0:.0f}s, peak {L.peak_gib():.2f} GiB")

    ctrl_c = acc[NG - 1]
    C1 = acc[0:2 * len(multi):2]          # the larger construct
    C2 = acc[1:2 * len(multi):2]          # the smaller construct
    n1 = ncc[0:2 * len(multi):2]
    n2 = ncc[1:2 * len(multi):2]
    gpos = {x: i for i, x in enumerate(ga_)}
    own = {gpos[t] for t in multi if t in gpos}
    scored_c = np.asarray(sorted(set(range(n_genes)) - own), dtype=np.int64)
    DC1 = L.delta(C1, ctrl_c)[:, scored_c]
    DC2 = L.delta(C2, ctrl_c)[:, scored_c]
    cons = dict(raw_cos=L.row_cos(DC1, DC2), raw_pearson=pearson_rows(DC1, DC2),
                resid_cos=L.row_cos(L.resid(DC1), L.resid(DC2)))
    mag_cons = np.sqrt(np.linalg.norm(L.resid(DC1), axis=1) * np.linalg.norm(L.resid(DC2), axis=1))
    out["construct_full_depth"] = {kk: summarise(v) for kk, v in cons.items()}
    out["construct_full_depth"]["n_cells"] = dict(
        larger_median=float(np.median(n1)), smaller_median=float(np.median(n2)),
        larger_min=int(n1.min()), smaller_min=int(n2.min()), n_scored_genes=int(scored_c.size))
    P(f"  cells/construct median {np.median(n1):.0f} (larger) / {np.median(n2):.0f} (smaller); "
      f"{scored_c.size} scored genes")
    for kk in cons:
        P(line(f"construct {kk}", out["construct_full_depth"][kk]))

    # ================================================================== C. matched perturbations
    P("\n=== C. the SAME perturbations (study-shared ∩ two-construct) ===")
    inter = sorted(tset & set(multi))
    si = np.asarray([list(shared).index(t) for t in inter], dtype=np.int64)
    ci = np.asarray([tpos[t] for t in inter], dtype=np.int64)
    P(f"  {len(inter)} perturbations carry both contrasts")
    out["matched"] = dict(n=len(inter), perturbations=inter[:20])
    for tag, arrs, ix in (("study", study, si), ("construct", cons, ci)):
        out["matched"][tag] = {kk: summarise(v[ix]) for kk, v in arrs.items()}
        for kk in arrs:
            P(line(f"{tag} {kk}", out["matched"][tag][kk]))
    for kk in ("raw_cos", "raw_pearson", "resid_cos"):
        d = L.boot_median_diff(study[kk][si], cons[kk][ci])
        out["matched"][f"study_minus_construct_{kk}"] = d
        P(f"  STUDY − CONSTRUCT {kk:12s} {L.fmt(d)}  P(<=0) {d['p_le0']:.4f}")

    # ---- D. response-magnitude strata --------------------------------------------------------
    # Replogle et al. 2022's ~0.50 same-line cross-experiment Pearson is quoted on STRONGLY
    # RESPONDING perturbations; the rows above use every shared perturbation and every scored
    # gene.  This is the same statistic by tertile of the pair's response magnitude.
    P("\n=== D. by response-magnitude tertile (the pair's geometric-mean residual norm) ===")
    out["strata"] = {}
    for tag, arrs, mg, ix in (("study_all", study, mag_study, np.arange(len(shared))),
                              ("construct_all", cons, mag_cons, np.arange(len(multi))),
                              ("study_matched", study, mag_study, si),
                              ("construct_matched", cons, mag_cons, ci)):
        m = mg[ix]
        q = np.nanpercentile(m, [100 / 3, 200 / 3])
        ter = np.digitize(m, q)
        out["strata"][tag] = {}
        for t in (0, 1, 2):
            sel = ix[ter == t]
            out["strata"][tag][f"T{t + 1}"] = dict(
                n=int(sel.size), mag_median=float(np.median(mg[sel])),
                **{kk: dict(median=float(np.nanmedian(v[sel])),
                            boot=L.boot_median(v[sel])) for kk, v in arrs.items()})
            r = out["strata"][tag][f"T{t + 1}"]
            P(f"  {tag:<18} T{t + 1} n {sel.size:4d} mag {r['mag_median']:7.3f}  "
              f"raw cos {L.fmt(r['raw_cos']['boot'])}  raw r {L.fmt(r['raw_pearson']['boot'])}  "
              f"resid cos {L.fmt(r['resid_cos']['boot'])}")

    # the headline pair, in one common unit
    out["headline"] = dict(
        study_raw_cos=float(np.nanmedian(study["raw_cos"][si])),
        construct_raw_cos=float(np.nanmedian(cons["raw_cos"][ci])),
        study_raw_pearson=float(np.nanmedian(study["raw_pearson"][si])),
        construct_raw_pearson=float(np.nanmedian(cons["raw_pearson"][ci])),
        study_ratio_n25=k["study"]["25"]["ratio"], construct_ratio_n25=g["rows"]["25"]["rho"]["point"])
    P(f"\nHEADLINE (matched perturbations, full depth): raw cosine STUDY "
      f"{out['headline']['study_raw_cos']:+.4f} vs CONSTRUCT "
      f"{out['headline']['construct_raw_cos']:+.4f}; raw Pearson "
      f"{out['headline']['study_raw_pearson']:+.4f} vs "
      f"{out['headline']['construct_raw_pearson']:+.4f}; 285's ratios "
      f"{out['headline']['study_ratio_n25']:.3f} vs {out['headline']['construct_ratio_n25']:.3f}")

    out["seconds"] = round(time.time() - T0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(WORK / "t8_298.npz", shared=shared, multi=np.asarray(multi),
                        inter=np.asarray(inter), n_c1=n1, n_c2=n2, ncell_gwps=nA, ncell_ess=nB,
                        mag_study=mag_study, mag_construct=mag_cons,
                        **{f"study_{kk}": v for kk, v in study.items()},
                        **{f"construct_{kk}": v for kk, v in cons.items()})
    L.Logger.dump(P, out, WORK / "t8_298.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
