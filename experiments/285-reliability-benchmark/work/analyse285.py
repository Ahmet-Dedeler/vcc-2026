"""285 -- assembly: the benchmark table, the depth law (measured ladder vs Spearman-Brown), the response-
magnitude strata, the n at which the removable term reaches each irreducible one, and `bench285.csv`.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyse285.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

P = L.Logger(L.WORK / "analyse285.log")
TAB = L.WORK / "tables285.md"


def j(name):
    p = L.WORK / name
    return json.load(open(p)) if p.exists() else None


def main() -> None:
    t0 = time.time()
    P(f"285 analyse285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    h1, k5, cd, br, na_j, na_h, orn, cs, pd_ = (j(x) for x in
                                                ("h1_285.json", "k562_285.json", "cd4_285.json", "broad_285.json",
                                                 "nadig_jurkat.json", "nadig_hepg2.json", "orion_285.json",
                                                 "csplit_285.json", "pds285.json"))
    zh1 = np.load(L.WORK / "h1_285.npz", allow_pickle=True)
    zk5 = np.load(L.WORK / "k562_285.npz", allow_pickle=True)
    znj = np.load(L.WORK / "nadig_jurkat.npz", allow_pickle=True)
    znh = np.load(L.WORK / "nadig_hepg2.npz", allow_pickle=True)
    out: dict = {}
    rows: list[dict] = []

    # ------------------------------------------------------------------ 1. the depth law, per target
    P("\n=== 1. THE DEPTH LAW: measured ladder vs Spearman-Brown, per target on a FIXED cohort ===")
    out["depth_law"] = {}
    for src, z, key, ns in (("h1-test", zh1, "ladder_a|{}", (25, 50, 100, 200)),
                            ("k562-gwps", zk5, "a|{}", (25, 50, 100)),
                            ("nadig-jurkat", znj, "a|{}", (25, 50)),
                            ("nadig-hepg2", znh, "a|{}", (25, 50))):
        arrs = {n: np.asarray(z[key.format(n)], dtype=np.float64) for n in ns if key.format(n) in z}
        if len(arrs) < 2:
            continue
        fixed = np.all([np.isfinite(arrs[n]) for n in arrs], axis=0)
        lo, hi = min(arrs), max(arrs)
        d = dict(n_fixed=int(fixed.sum()), rungs={int(n): float(np.median(arrs[n][fixed])) for n in arrs})
        # SB from the smallest rung and from the largest, per target
        for base, label in ((lo, "from_smallest"), (hi, "from_largest")):
            pred = {}
            for n in arrs:
                if n == base:
                    continue
                k = n / base
                r0 = arrs[base][fixed]
                pr = k * r0 / (1 + (k - 1) * r0)
                meas = arrs[n][fixed]
                pred[int(n)] = dict(sb_median=float(np.median(pr)), measured_median=float(np.median(meas)),
                                    diff=L.boot_median_diff(meas, pr))
            d[label] = pred
        out["depth_law"][src] = d
        P(f"  {src}: fixed cohort {int(fixed.sum())} targets; measured "
          + ", ".join(f"n={n}: {d['rungs'][n]:+.4f}" for n in sorted(d["rungs"])))
        for label in ("from_smallest", "from_largest"):
            for n, v in sorted(d[label].items()):
                P(f"    SB {label:14s} -> n={n:4d}: predicted {v['sb_median']:+.4f}  measured "
                  f"{v['measured_median']:+.4f}  measured − SB {L.fmt(v['diff'])}")

    # ------------------------------------------------------------------ 2. cells vs UMIs at a fixed budget
    P("\n=== 2. CELLS vs UMI DEPTH at a FIXED TOTAL-UMI BUDGET (does the cell count matter beyond the UMIs?) ===")
    out["cells_vs_depth"] = {}
    for src, z, pref, pairs in (("h1-test", zh1, "thin|{}|{}", ((50, 0.5, 25), (100, 0.5, 50), (200, 0.5, 100))),
                                ("k562-gwps", zk5, "thin|{}|{}", ((50, 0.5, 25), (100, 0.5, 50)))):
        for n2, f, n1 in pairs:
            ka, kb = pref.format(n2, f), pref.format(n1, 1.0)
            if ka not in z or kb not in z:
                continue
            a, b = np.asarray(z[ka], dtype=np.float64), np.asarray(z[kb], dtype=np.float64)
            d = L.boot_median_diff(a, b)
            out["cells_vs_depth"][f"{src}|{n2}cells_half_depth_vs_{n1}cells_full"] = d
            P(f"  {src}: {n2} cells at 1/2 depth ({np.nanmedian(a):+.4f}) − {n1} cells at full depth "
              f"({np.nanmedian(b):+.4f}) = {L.fmt(d)}   [same total UMIs]")

    # ------------------------------------------------------------------ 3. the count split vs the cell split
    P("\n=== 3. COUNT SPLIT vs CELL SPLIT at matched UMIs (the calibration the CD4 bank needs) ===")
    out["csplit"] = {}
    for src in ("h1_test", "k562_gwps"):
        for n, v in sorted(cs[src].items(), key=lambda kv: int(kv[0])):
            out["csplit"][f"{src}|{n}"] = dict(cell=v["cell"]["median"], count=v["count"]["median"],
                                               ratio=v["ratio_of_medians"], offset=v["diff"]["point"])
            P(f"  {src:10s} n={int(n):4d}: cell {v['cell']['median']:+.4f}  count {v['count']['median']:+.4f}  "
              f"ratio {v['ratio_of_medians']:.3f}  count − cell {L.fmt(v['diff'])}")
    offs = [v["offset"] for v in out["csplit"].values()]
    out["csplit_offset_mean"] = float(np.mean(offs))
    P(f"  the count-minus-cell OFFSET is nearly constant across n and source: mean {np.mean(offs):+.4f} "
      f"(range {min(offs):+.4f} .. {max(offs):+.4f})")

    # ------------------------------------------------------------------ 4. the response-magnitude strata
    P("\n=== 4. RESPONSE-MAGNITUDE STRATA (terciles of the full-depth residual row norm) ===")
    out["strata"] = {}
    for src, z, key, mag in (("h1-test n=100", zh1, "ladder_a|100", np.asarray(zh1["magnitude"])),
                             ("k562-gwps n=50", zk5, "a|50", np.asarray(zk5["magnitude"]))):
        x = np.asarray(z[key], dtype=np.float64)
        ok = np.isfinite(x)
        q = np.percentile(mag[ok], [100 / 3, 200 / 3])
        ter = {"T1 weak": mag <= q[0], "T2 mid": (mag > q[0]) & (mag <= q[1]), "T3 strong": mag > q[1]}
        out["strata"][src] = {k: dict(**L.stats(x[v & ok]), boot=L.boot_median(x[v & ok])) for k, v in ter.items()}
        P(f"  {src}: " + "  ".join(f"{k} {out['strata'][src][k]['median']:+.4f} (n {out['strata'][src][k]['n']})"
                                   for k in ter))

    # ------------------------------------------------------------------ 5. the benchmark table + csv
    P("\n=== 5. THE BENCHMARK TABLE (rows = source x unit x n) ===")

    def add(source, unit, n, cells, val, n_t, extra=""):
        if not val or "median" not in val:
            return
        rows.append(dict(source=source, unit=unit, n_per_side=n, cells_per_side=cells,
                         median_cos_resid=round(float(val["median"]), 5),
                         lo=round(float(val["boot"]["lo"]), 5) if "boot" in val else "",
                         hi=round(float(val["boot"]["hi"]), 5) if "boot" in val else "",
                         share_le0=round(float(val["share_le0"]), 4),
                         n_targets=n_t, note=extra))

    for n, v in sorted(h1["ladder_a"].items(), key=lambda kv: int(kv[0])):
        add("h1-test (H1 hESC, Flex)", "a cells (same guide)", int(n), int(n), v, v["n"])
    for n, v in sorted(h1["lanes"].items(), key=lambda kv: int(kv[0])):
        if v.get("n"):
            add("h1-test (H1 hESC, Flex)", "d lane", int(n), int(n), v, v["n"])
    for gn, v in h1["guide"].items():
        add("h1-test (H1 hESC, Flex)", "b guide construct", int(gn), int(gn),
            dict(**v["guide"], boot=L.boot_median(np.asarray([x["guide"] for x in v["per_target"].values()]))),
            v["guide"]["n"])
        add("h1-test (H1 hESC, Flex)", "a cells, one construct", int(gn), int(gn), dict(**v["cells"], boot=L.boot_median(np.asarray([x["cells"] for x in v["per_target"].values()]))), v["cells"]["n"])
    for k, v in sorted(h1["thin"].items()):
        add("h1-test (H1 hESC, Flex)", f"a' thinned f={v['frac']}", v["n_cells"], v["n_cells"], v, v["n"])
    for n, v in sorted(k5["ladder_a"].items(), key=lambda kv: int(kv[0])):
        add("K562-GWPS (panel targets)", "a cells (same guide)", int(n), int(n), v["fixed"], v["fixed"]["n"], "fixed cohort")
        add("K562-GWPS (panel targets)", "a cells (same guide)", int(n), int(n), v["all"], v["all"]["n"], "all available")
    for n, v in sorted(k5["gemgroup"].items(), key=lambda kv: int(kv[0])):
        if v.get("n"):
            add("K562-GWPS (panel targets)", "d gemgroup", int(n), int(n), v, v["n"])
    for gn, v in k5["guide"].items():
        add("K562-GWPS (panel targets)", "b guide construct", int(gn), int(gn),
            dict(**v["guide"], boot=L.boot_median(np.asarray([x["guide"] for x in v["per_target"].values()]))), v["guide"]["n"])
        add("K562-GWPS (panel targets)", "a cells, one construct", int(gn), int(gn),
            dict(**v["cells"], boot=L.boot_median(np.asarray([x["cells"] for x in v["per_target"].values()]))), v["cells"]["n"])
    for k, v in sorted(k5["thin"].items()):
        add("K562-GWPS (panel targets)", f"a' thinned f={v['frac']}", v["n_cells"], v["n_cells"], v, v["n"])
    for n, v in sorted(k5["study"].items(), key=lambda kv: int(kv[0])):
        add("K562 GWPS vs essential (2,053 shared)", "a cells within study", int(n), int(n), v["within_mean"], v["n_targets"])
        add("K562 GWPS vs essential (2,053 shared)", "e different STUDY, same line", int(n), int(n), v["cross"], v["n_targets"])
    for src, na in (("Nadig Jurkat", na_j), ("Nadig HepG2", na_h)):
        if not na:
            continue
        for n, v in sorted(na["budgets"].items(), key=lambda kv: int(kv[0])):
            add(src, "a cells (all guides)", int(n), int(n), v["cells"], v["cells"]["n"])
            add(src, "d gemgroup", int(n), int(n), v["gem"], v["gem"]["n"])
            add(src, "b guide construct", int(n), int(n), v["guide"], v["guide"]["n"])
            add(src, "a cells, one construct", int(n), int(n), v["guide_cellctrl"], v["guide_cellctrl"]["n"])
    for T, v in br["budgets"].items():
        if "pooled_cell" not in v:
            continue
        add("Broad-16 (12 lines)", "a cells", f"T={T}", "", v["pooled_cell"], v["pooled_cell"]["n"])
        add("Broad-16 (12 lines)", "b guide construct", f"T={T}", "", v["pooled_guide"], v["pooled_guide"]["n"])
        if "lines_matched" in v:
            add("Broad-16 (12 lines)", "f different LINE", f"T={T}", "", v["lines_matched"], v["lines_matched"]["n"])
    for T, blk in cd["budgets"].items():
        for key, v in sorted(blk.items()):
            if "|vs_" in key or "|" not in key:
                continue
            cond, unit = key.split("|")
            nm = {"countsplit": "a* count split (Poisson only)", "guide": "b guide construct",
                  "donor": "c different DONOR", "state": "s different STATE"}.get(unit, unit)
            add(f"CD4 GWCD4i {cond}", nm, f"T={T}", round(v["cells_per_side"] or 0, 1), v, v["n"])
    for line, v in orn.items():
        add(f"267 cache {line}", "a cells", int(v["cells_rand"]), int(v["cells_rand"]), v["rand"], v["n"])
        add(f"267 cache {line}", "d batch", int(v["cells_batch"]), int(v["cells_batch"]), v["batch"], v["n"])
        add(f"267 cache {line}", "b guide construct", int(v["guide"]["cells"]), int(v["guide"]["cells"]),
            v["guide"], v["guide"]["n_guide"])
    with open(L.WORK / "bench285.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    P(f"  wrote {len(rows)} rows to work/bench285.csv")

    # ------------------------------------------------------------------ 6. n* -- where the ladder meets a ceiling
    P("\n=== 6. THE n AT WHICH MORE CELLS OF ONE CONSTRUCT REACH AN IRREDUCIBLE TERM (Spearman-Brown) ===")
    out["n_star"] = {}

    def nstar(r0, n0, target):
        if not (0 < r0 < 1) or not (0 < target < 1):
            return None
        k = target * (1 - r0) / (r0 * (1 - target))
        return float(k * n0)

    ceilings = [
        ("h1-test", float(h1["ladder_a"]["25"]["median"]), 25,
         [("b guide construct (n=37, 5 targets)", h1["guide"]["37"]["guide"]["median"]),
          ("d lane (n=25)", h1["lanes"]["25"]["median"])]),
        ("K562-GWPS panel", float(k5["ladder_a"]["25"]["fixed"]["median"]), 25,
         [("b guide construct (n=35, 8 targets)", k5["guide"]["35"]["guide"]["median"]),
          ("d gemgroup (n=25)", k5["gemgroup"]["25"]["median"])]),
        ("K562 shared (essential) targets", float(k5["study"]["25"]["within_mean"]["median"]), 25,
         [("e different STUDY, same line (n=25)", k5["study"]["25"]["cross"]["median"])]),
        ("Nadig Jurkat", float(na_j["budgets"]["25"]["guide_cellctrl"]["median"]), 25,
         [("b guide construct (n=25)", na_j["budgets"]["25"]["guide"]["median"]),
          ("d gemgroup (n=25)", na_j["budgets"]["25"]["gem"]["median"])]),
        ("Nadig HepG2", float(na_h["budgets"]["25"]["guide_cellctrl"]["median"]), 25,
         [("b guide construct (n=25)", na_h["budgets"]["25"]["guide"]["median"])]),
    ]
    for src, r0, n0, cl in ceilings:
        for name, tv in cl:
            ns = nstar(r0, n0, tv)
            out["n_star"][f"{src} | {name}"] = dict(r_cells_at_n0=r0, n0=n0, ceiling=tv, n_star=ns)
            P(f"  {src:32s} vs {name:42s} ceiling {tv:+.4f}: reached at n = "
              + (f"{ns:.0f} cells/side" if ns and ns > 0 else "n/a (ceiling <= 0)")
              + (f"  [already exceeded at n={n0}]" if ns and ns < n0 else ""))

    # ------------------------------------------------------------------ 7. DISATTENUATION: the unit's own ceiling
    P("\n=== 7. THE UNIT'S OWN CEILING (rho = the unit's agreement with the cell-sampling attenuation divided out) ===")
    P("   observed cosine of a unit at n cells/side ~= rho_unit * a(n), and the same-construct CELL split has")
    P("   rho = 1, so a(n) IS the cell-split row: rho_unit = cos_unit(n) / cos_cells(n) at the SAME n.")
    out["rho"] = {}

    def rho(num, den, label, n, cells_note=""):
        num = np.asarray(num, dtype=np.float64); den = np.asarray(den, dtype=np.float64)
        m = np.isfinite(num) & np.isfinite(den)
        num, den = num[m], den[m]
        if num.size < 5:
            return
        rng = np.random.default_rng(L.SEED)
        idx = rng.integers(0, num.size, size=(L.NBOOT, num.size))
        b = np.median(num[idx], 1) / np.median(den[idx], 1)
        v = dict(rho=float(np.median(num) / np.median(den)), lo=float(np.percentile(b, 2.5)),
                 hi=float(np.percentile(b, 97.5)), n=int(num.size), n_per_side=n, note=cells_note,
                 cos_unit=float(np.median(num)), cos_cells=float(np.median(den)))
        out["rho"][label] = v
        P(f"  {label:56s} n={n:>7}  rho {v['rho']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}]  "
          f"(unit {v['cos_unit']:+.4f} / cells {v['cos_cells']:+.4f}, {v['n']} targets) {cells_note}")

    zk = zk5
    for gn in ("37", "62"):
        if gn in h1["guide"]:
            pt = h1["guide"][gn]["per_target"]
            rho([x["guide"] for x in pt.values()], [x["cells"] for x in pt.values()],
                "h1-test  b GUIDE construct / same-construct cells", gn)
    for gn in ("35", "75"):
        if gn in k5["guide"]:
            pt = k5["guide"][gn]["per_target"]
            rho([x["guide"] for x in pt.values()], [x["cells"] for x in pt.values()],
                "K562-GWPS  b GUIDE construct / same-construct cells", gn)
    for n in ("25", "50", "100"):
        if n in h1["lanes"] and f"lane|{n}" in zh1 and f"cells_gc|{n}" in zh1:
            rho(zh1[f"lane|{n}"], zh1[f"cells_gc|{n}"], "h1-test  d LANE / random cells", n)
        if n in k5["gemgroup"] and f"gem|{n}" in zk and f"a|{n}" in zk:
            rho(zk[f"gem|{n}"], zk[f"a|{n}"], "K562-GWPS  d GEMGROUP / random cells", n)
    for src, z, na in (("Nadig Jurkat", znj, na_j), ("Nadig HepG2", znh, na_h)):
        for n in ("25", "50"):
            if f"guide|{n}" in z and f"guide_cellctrl|{n}" in z:
                rho(z[f"guide|{n}"], z[f"guide_cellctrl|{n}"], f"{src}  b GUIDE construct / same-construct cells", n)
            if f"gem|{n}" in z and f"a|{n}" in z:
                rho(z[f"gem|{n}"], z[f"a|{n}"], f"{src}  d GEMGROUP / random cells", n)
    for n in ("25", "50"):
        if f"study_cross|{n}" in zk:
            rho(zk[f"study_cross|{n}"], zk[f"study_within|{n}"],
                "K562 GWPS-vs-essential  e different STUDY / within-study cells", n)
    for line, v in orn.items():
        out["rho"][f"267 {line}  d BATCH / random cells"] = dict(
            rho=float(v["batch"]["median"] / v["rand"]["median"]), n_per_side=v["cells_batch"],
            cos_unit=v["batch"]["median"], cos_cells=v["rand"]["median"], n=v["n"], note="ratio of medians (cached)")
        P(f"  {'267 ' + line + '  d BATCH / random cells':56s} n={v['cells_batch']:>7.0f}  "
          f"rho {v['batch']['median'] / v['rand']['median']:+.3f} (ratio of medians, cached arrays)")
    for T, v in br["budgets"].items():
        if "pooled_guide" in v:
            r_ = v["pooled_guide"]["median"] / v["pooled_cell"]["median"]
            out["rho"][f"Broad-16  b GUIDE construct / cells (T={T})"] = dict(
                rho=float(r_), n_per_side=T, cos_unit=v["pooled_guide"]["median"], cos_cells=v["pooled_cell"]["median"],
                n=v["pooled_guide"]["n"], note="pooled over 12 lines, ratio of medians")
            P(f"  {'Broad-16  b GUIDE construct / cells':56s} n={T:>7}  rho {r_:+.3f} (pooled 12 lines)")
        if "lines_matched" in v:
            r_ = v["lines_matched"]["median"] / v["pooled_cell"]["median"]
            out["rho"][f"Broad-16  f different LINE / cells (T={T})"] = dict(
                rho=float(r_), n_per_side=T, cos_unit=v["lines_matched"]["median"], cos_cells=v["pooled_cell"]["median"],
                n=v["lines_matched"]["n"], note="pooled over 66 line pairs, ratio of medians")
            P(f"  {'Broad-16  f different LINE / cells':56s} n={T:>7}  rho {r_:+.3f} (66 line pairs)")
    # CD4: the count split corrected to a cell split by the +0.0438 offset measured in section 3
    off = out["csplit_offset_mean"]
    for T, blk in cd["budgets"].items():
        for cond in ("Rest", "Stim8hr", "Stim48hr"):
            base = blk.get(f"{cond}|countsplit")
            if not base:
                continue
            cellsplit = base["median"] - off
            for unit, nm in (("guide", "b GUIDE construct"), ("donor", "c different DONOR"), ("state", "s different STATE")):
                v = blk.get(f"{cond}|{unit}")
                if not v or cellsplit <= 0:
                    continue
                lab = f"CD4 {cond}  {nm} / cells (count split − {off:.3f})"
                out["rho"][f"{lab} (T={T})"] = dict(rho=float(v["median"] / cellsplit), n_per_side=T,
                                                    cos_unit=v["median"], cos_cells=float(cellsplit), n=v["n"],
                                                    note="count split corrected to a cell split by the section-3 offset")
                P(f"  {lab:56s} n={T:>7}  rho {v['median'] / cellsplit:+.3f}  (unit {v['median']:+.4f} / "
                  f"corrected cells {cellsplit:+.4f}, {v['n']} targets, ~{(v['cells_per_side'] or 0):.0f} cells/side)")

    # ------------------------------------------------------------------ 8. where the ladder meets each ceiling
    P("\n=== 8. THE n AT WHICH MORE CELLS OF ONE CONSTRUCT STOP HELPING (a(n) = rho_unit, Spearman-Brown) ===")
    out["n_stop"] = {}
    ladders = {"h1-test": (25, float(h1["ladder_a"]["25"]["median"])),
               "K562-GWPS": (25, float(k5["ladder_a"]["25"]["fixed"]["median"])),
               "Nadig Jurkat": (25, float(na_j["budgets"]["25"]["cells"]["median"])),
               "Nadig HepG2": (25, float(na_h["budgets"]["25"]["cells"]["median"]))}
    for lab, v in out["rho"].items():
        src = lab.split("  ")[0].strip()
        if src not in ladders or not (0 < v["rho"] < 1):
            continue
        n0, a0 = ladders[src]
        k = v["rho"] * (1 - a0) / (a0 * (1 - v["rho"]))
        out["n_stop"][lab] = float(k * n0)
        P(f"  {lab:60s}: a(n) reaches rho = {v['rho']:+.3f} at n = {k * n0:8.0f} cells/side")

    out["seconds"] = round(time.time() - t0, 1)
    L.Logger.dump(P, out, L.WORK / "analyse285.json")
    P(f"\ndone {out['seconds']}s")
    P.save()


if __name__ == "__main__":
    main()
