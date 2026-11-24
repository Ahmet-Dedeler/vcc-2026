"""290 step 5 — REPLICATION on the declared partitions #2 and #3 (PREREG §2, §8), plus a labelled
POST-HOC sweep over all 20 admissible partitions.

Partitions #2 and #3 mix the primary partition's source halves, so their feature blocks must be
re-aggregated per partition — `work/persrc290.npz` holds the per-SOURCE per-gene values, so no bank
is re-read.  KA-5 first: the aggregation reproduces `work/feat290.csv` for the primary partition at
0.0e+00.

    PYTHONDONTWRITEBYTECODE=1 uv run python rep290.py

Writes work/rep290.json, work/rep290.log.
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import pandas as pd

import lib290 as L
from fit290 import draw_partitions

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def aggregate(Q: dict, rows: np.ndarray, on_axis: np.ndarray) -> np.ndarray:
    """PREREG §5's aggregation over an arbitrary set of source rows -> the 17 frozen columns."""
    with np.errstate(invalid="ignore"):
        mag = np.nanmean(Q["mag"][rows], axis=0)
        cols = [np.log(np.where(mag > 0, mag, np.nan)),
                np.nanmean(Q["breadth"][rows], axis=0),
                np.nanmean(Q["gen_cos"][rows], axis=0),
                np.nanmean(np.abs(Q["gen_cos"][rows]), axis=0),
                np.nanmean(Q["pc5"][rows], axis=0),
                np.nanmean(Q["pc20"][rows], axis=0)]
        cols += [np.nanmean(Q[c][rows], axis=0) for c in L.PE_COLS]
        cols += [np.nanmean(np.abs(Q["proj_prolif"][rows]), axis=0),
                 np.nanmean(Q["ontarget"][rows], axis=0),
                 on_axis.astype(float)]
    return np.column_stack(cols)


def main() -> None:
    h = L.check_prereg()
    P(f"prereg {h} MATCH")

    z = np.load(L.WORK / "persrc290.npz", allow_pickle=True)
    keys = np.asarray(z["keys"]).astype(str)
    lines = np.asarray(z["lines"]).astype(str)
    on_axis = np.asarray(z["on_axis"])
    Q = {k: np.asarray(z[k]) for k in
         (["mag", "breadth", "gen_cos", "pc5", "pc20"] + L.PE_COLS
          + ["proj_prolif", "proj_stress", "ontarget"])}
    h1 = np.asarray(z["half1"])
    h2 = np.asarray(z["half2"])

    # ---- KA-5: the aggregation reproduces feat290.csv on the primary halves --------------------
    feat = pd.read_csv(L.WORK / "feat290.csv")
    ka5 = {}
    for tag, rows in (("h1", h1), ("h2", h2)):
        A = aggregate(Q, rows, on_axis)
        B = feat[[c if c == "on_axis" else f"{c}_{tag}" for c in L.FEATS290]].to_numpy(float)
        d = float(np.nanmax(np.abs(A - B)))
        same_nan = bool(np.array_equal(np.isnan(A), np.isnan(B)))
        # the frozen artifact is a CSV of float32 aggregates, so the round-trip carries ~1e-7;
        # 279's own KA-0 / KA-1 use the same < 1e-6 bar for the same reason (build279.py:131,58)
        ka5[tag] = dict(max_abs_diff=d, same_nan_pattern=same_nan,
                        pass_=bool(d < 1e-6 and same_nan))
        P(f"KA-5  re-aggregated {tag} vs feat290.csv: max|diff| {d:.1e} (bar < 1e-6, the float32 "
          f"CSV round-trip), NaN pattern identical {same_nan} -> "
          f"{'PASS' if ka5[tag]['pass_'] else 'FAIL'}")
    if not all(v["pass_"] for v in ka5.values()):
        raise SystemExit("KA-5 FAILED — abort")

    # ---- the shared machinery ------------------------------------------------------------------
    held = L.load_held()
    X21 = held[L.FEATS].to_numpy(float)
    comp = np.load(L.SPLITS284, allow_pickle=True)["comp"]
    fold = L.component_folds(comp, 5)
    zz = np.load(L.HELD279, allow_pickle=True)
    Ccross = np.asarray(zz["Ccross"])
    pairs = np.asarray(zz["pairs_cross"])
    parts, _ = draw_partitions(lines, pairs)
    ref284 = json.loads(L.REM284J.read_text())

    def one(i: int, declared: bool, full: bool = True) -> dict:
        p1, p2, ha, hb = parts[i]
        s1 = np.array([q for q, ln in enumerate(lines) if ln in ha])
        s2 = np.array([q for q, ln in enumerate(lines) if ln in hb])
        C1 = L.cmap(Ccross, p1).astype(float)
        C2 = L.cmap(Ccross, p2).astype(float)
        n1 = np.isfinite(Ccross[p1]).sum(axis=0)
        n2 = np.isfinite(Ccross[p2]).sum(axis=0)
        both = (n1 >= 1) & (n2 >= 1)
        pf1, _ = L.oof_ridge(X21, C1, fold)
        pf2, _ = L.oof_ridge(X21, C2, fold)
        r1, r2 = C1 - pf1, C2 - pf2
        ok = both & np.isfinite(r1) & np.isfinite(r2)
        u = np.flatnonzero(ok)
        w = ref284["partitions"][i]
        ceil_, rem_ = L.spearman(C1[u], C2[u]), L.spearman(r1[u], r2[u])
        d = dict(tag=f"partition#{i + 1}", declared=declared, n=len(u),
                 n_sources_half1=len(s1), n_sources_half2=len(s2),
                 half1_lines=ha, half2_lines=hb, ceiling=ceil_, remainder=rem_,
                 reproduces_284=bool(abs(ceil_ - w["ceiling"]) < 5e-5
                                     and abs(rem_ - w["rho_remainder"]) < 5e-5))
        res = {}
        for direction, (srows, y, ytag) in (("h1->r2", (s1, r2, "remainder_half2")),
                                            ("h2->r1", (s2, r1, "remainder_half1"))):
            X = aggregate(Q, srows, on_axis)[u]
            yy = y[u]
            pr, _ = L.oof_ridge(X, yy, fold[u])
            pa, _ = L.oof_ridge(X21[u], yy, fold[u])
            j_pa = L.FEATS290.index("pe_prolif_abs")
            j_pc = L.FEATS290.index("gen_pc5")
            j_lm = L.FEATS290.index("logmag")
            if full:
                bb = L.block_boot(yy, dict(response=pr, annotation=pa), comp[u],
                                  np.random.default_rng(L.SEED + 121 + i))
                cells = dict(response=bb["response"]["rho"], annotation=bb["annotation"]["rho"],
                             paired=bb["_paired_rho"]["response-annotation"],
                             pe_prolif_abs=L.boot_rho(X[:, j_pa], yy, comp[u],
                                                      np.random.default_rng(L.SEED + 131 + i)),
                             gen_pc5=L.boot_rho(X[:, j_pc], yy, comp[u],
                                                np.random.default_rng(L.SEED + 141 + i)),
                             logmag=L.boot_rho(X[:, j_lm], yy, comp[u],
                                               np.random.default_rng(L.SEED + 151 + i)))
            else:
                cells = dict(response=dict(point=L.spearman(yy, pr)),
                             annotation=dict(point=L.spearman(yy, pa)),
                             paired=dict(mean=L.spearman(yy, pr) - L.spearman(yy, pa)),
                             pe_prolif_abs=dict(point=L.spearman(X[:, j_pa], yy)),
                             gen_pc5=dict(point=L.spearman(X[:, j_pc], yy)),
                             logmag=dict(point=L.spearman(X[:, j_lm], yy)))
            res[direction] = dict(target=ytag, **cells,
                                  other_half_remainder=L.spearman(
                                      yy, (r1 if ytag.endswith("2") else r2)[u]))
        d["directions"] = res
        return d

    out = {"prereg_sha256": h, "ka5": ka5, "declared": {}, "sweep": []}
    P("\nDECLARED replication partitions (#2 and #3):")
    for i in (1, 2):
        d = one(i, True)
        out["declared"][d["tag"]] = d
        P(f"  {d['tag']}: n {d['n']}, {d['n_sources_half1']}/{d['n_sources_half2']} sources, "
          f"ceiling {d['ceiling']:+.4f} remainder {d['remainder']:+.4f} "
          f"(reproduces 284: {d['reproduces_284']})")
        for k, v in d["directions"].items():
            P(f"    {k}: response {v['response']['point']:+.4f} [{v['response']['lo']:+.4f}, "
              f"{v['response']['hi']:+.4f}]  annotation {v['annotation']['point']:+.4f}  "
              f"paired D {v['paired']['mean']:+.4f} [{v['paired']['lo']:+.4f}, "
              f"{v['paired']['hi']:+.4f}]  other-half remainder {v['other_half_remainder']:+.4f}")
            P(f"       pe_prolif_abs {v['pe_prolif_abs']['point']:+.4f} "
              f"[{v['pe_prolif_abs']['lo']:+.4f}, {v['pe_prolif_abs']['hi']:+.4f}]  gen_pc5 "
              f"{v['gen_pc5']['point']:+.4f} [{v['gen_pc5']['lo']:+.4f}, "
              f"{v['gen_pc5']['hi']:+.4f}]  logmag {v['logmag']['point']:+.4f}")

    P("\nPOST-HOC sweep over all 20 admissible partitions (both directions each):")
    for i in range(len(parts)):
        d = one(i, i in (0, 1, 2), full=False)
        out["sweep"].append(d)
        a, b = d["directions"]["h1->r2"], d["directions"]["h2->r1"]
        P(f"  #{i + 1:2d} n {d['n']:5d} ceiling {d['ceiling']:+.4f} rem {d['remainder']:+.4f} | "
          f"h1->r2 resp {a['response']['point']:+.4f} ann {a['annotation']['point']:+.4f} "
          f"pe {a['pe_prolif_abs']['point']:+.4f} pc5 {a['gen_pc5']['point']:+.4f} | "
          f"h2->r1 resp {b['response']['point']:+.4f} ann {b['annotation']['point']:+.4f} "
          f"pe {b['pe_prolif_abs']['point']:+.4f} pc5 {b['gen_pc5']['point']:+.4f}")

    resp = np.array([d["directions"][k]["response"]["point"] for d in out["sweep"]
                     for k in ("h1->r2", "h2->r1")])
    ann = np.array([d["directions"][k]["annotation"]["point"] for d in out["sweep"]
                    for k in ("h1->r2", "h2->r1")])
    pe = np.array([d["directions"][k]["pe_prolif_abs"]["point"] for d in out["sweep"]
                   for k in ("h1->r2", "h2->r1")])
    pc = np.array([d["directions"][k]["gen_pc5"]["point"] for d in out["sweep"]
                   for k in ("h1->r2", "h2->r1")])
    out["sweep_summary"] = dict(
        n_cells=len(resp),
        response=dict(median=float(np.median(resp)), lo=float(resp.min()), hi=float(resp.max()),
                      n_ge_0p10=int((resp >= 0.10).sum()), n_positive=int((resp > 0).sum())),
        annotation=dict(median=float(np.median(ann)), lo=float(ann.min()), hi=float(ann.max())),
        pe_prolif_abs=dict(median=float(np.median(pe)), n_negative=int((pe < 0).sum())),
        gen_pc5=dict(median=float(np.median(pc)), n_negative=int((pc < 0).sum())))
    s = out["sweep_summary"]
    P(f"\nsweep: response median {s['response']['median']:+.4f} "
      f"[{s['response']['lo']:+.4f}, {s['response']['hi']:+.4f}], "
      f"{s['response']['n_ge_0p10']}/{s['n_cells']} >= +0.10; annotation median "
      f"{s['annotation']['median']:+.4f}; pe_prolif_abs {s['pe_prolif_abs']['n_negative']}"
      f"/{s['n_cells']} negative (median {s['pe_prolif_abs']['median']:+.4f}); gen_pc5 "
      f"{s['gen_pc5']['n_negative']}/{s['n_cells']} negative "
      f"(median {s['gen_pc5']['median']:+.4f})")

    out["seconds"] = time.time() - T0
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    (L.WORK / "rep290.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "rep290.log").write_text("\n".join(LOG))
    P(f"peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
