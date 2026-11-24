"""279 — the preregistered test on the held genes.

Reads work/held279.npz (built by build279.py) and applies PREREGISTRATION.md §4-§6 exactly:
the raw effect, the matched effect, the AUC, the random-family null, the HGNC family ranking,
the two TF-family controls, the three single-dataset replications and the in-sample panel
reference.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyse279.py

Writes work/analyse279.json, work/znf279.csv.  Read-only outside work/.
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
import lib276 as L  # noqa: E402
import inv279 as I  # noqa: E402
import build279 as B  # noqa: E402

SEED = 20260906
B_BOOT = 10000
B_MATCH = 2000
N_DRAW = 50
N_RANDFAM = 100
ZNF_GROUPS = I.ZNF_GROUPS


# ------------------------------------------------------------------ group membership (frozen §2)
def classify(genes: list[str]) -> tuple[np.ndarray, dict[str, set[str]]]:
    groups = I.hgnc_groups()
    cls = np.empty(len(genes), dtype=object)
    for i, g in enumerate(genes):
        gg = groups.get(g, set())
        if any(x in gg for x in ZNF_GROUPS):
            cls[i] = "ZNF"
        elif any(("inc finger" in x) or x.startswith("Zinc") for x in gg):
            cls[i] = "OTHER_ZF"
        else:
            cls[i] = "NONZNF"
    return cls, groups


def symbol_znf(g: str) -> bool:
    """276 §5's post-hoc symbol rule (declared secondary)."""
    return g.startswith("ZNF") or g.startswith("ZKSCAN") or g.startswith("ZSCAN")


# ------------------------------------------------------------------ statistics
def auc_lower(x: np.ndarray, y: np.ndarray) -> float:
    """P(x < y) + 0.5 P(x = y): >0.5 means group x sits LOWER than group y."""
    u1 = mannwhitneyu(x, y, alternative="two-sided").statistic     # #(x>y) + 0.5 ties
    return float(1.0 - u1 / (x.size * y.size))


def two_sample(x: np.ndarray, y: np.ndarray, b: int = B_BOOT, seed: int = SEED) -> dict:
    """Δ = median(x) − median(y) with a two-sample bootstrap over genes."""
    rng = np.random.default_rng(seed)
    d0 = float(np.median(x) - np.median(y))
    dm = float(np.mean(x) - np.mean(y))
    a0 = auc_lower(x, y)
    draws = np.empty(b); adr = np.empty(b)
    for i in range(b):
        xs = x[rng.integers(0, x.size, x.size)]
        ys = y[rng.integers(0, y.size, y.size)]
        draws[i] = np.median(xs) - np.median(ys)
        adr[i] = np.nan
    lo, hi = np.percentile(draws, [2.5, 97.5])
    p = float(mannwhitneyu(x, y, alternative="less").pvalue)
    return dict(n_x=int(x.size), n_y=int(y.size), median_x=float(np.median(x)),
                median_y=float(np.median(y)), d_median=d0, ci=[float(lo), float(hi)],
                d_mean=dm, auc=a0, mwu_p_one_sided=p,
                excludes_zero_negative=bool(hi < 0), excludes_zero_positive=bool(lo > 0))


def auc_ci(x: np.ndarray, y: np.ndarray, b: int = 2000, seed: int = SEED) -> list[float]:
    rng = np.random.default_rng(seed + 1)
    dr = np.empty(b)
    for i in range(b):
        dr[i] = auc_lower(x[rng.integers(0, x.size, x.size)], y[rng.integers(0, y.size, y.size)])
    return [float(v) for v in np.percentile(dr, [2.5, 97.5])]


def deciles(v: np.ndarray) -> np.ndarray:
    """decile index 0..9 over the finite values; -1 where not finite."""
    out = np.full(v.size, -1, np.int64)
    ok = np.isfinite(v)
    q = np.quantile(v[ok], np.linspace(0, 1, 11)[1:-1])
    out[ok] = np.searchsorted(q, v[ok], side="right")
    return out


def match_draw(cells_z: np.ndarray, cells_n: np.ndarray, idx_z: np.ndarray, idx_n: np.ndarray,
               rng) -> tuple[np.ndarray, np.ndarray]:
    """1:1 without replacement inside the same (expr-decile, mag-decile) cell."""
    pool: dict[int, list[int]] = {}
    for c, i in zip(cells_n, idx_n):
        pool.setdefault(int(c), []).append(int(i))
    for c in pool:
        arr = np.asarray(pool[c]); rng.shuffle(arr); pool[c] = list(arr)
    zi, ni = [], []
    order = rng.permutation(idx_z.size)
    for q in order:
        c = int(cells_z[q])
        lst = pool.get(c)
        if lst:
            ni.append(lst.pop())
            zi.append(int(idx_z[q]))
    return np.asarray(zi, np.int64), np.asarray(ni, np.int64)


def matched_effect(score: np.ndarray, cell: np.ndarray, is_z: np.ndarray, is_n: np.ndarray,
                   n_draw: int = N_DRAW, b: int = B_MATCH, seed: int = SEED) -> dict:
    ok = (cell >= 0) & np.isfinite(score)
    idx_z = np.flatnonzero(is_z & ok)
    idx_n = np.flatnonzero(is_n & ok)
    rng = np.random.default_rng(seed)
    draws, pairs = [], []
    for _ in range(n_draw):
        zi, ni = match_draw(cell[idx_z], cell[idx_n], idx_z, idx_n, rng)
        pairs.append((zi, ni))
        draws.append(float(np.median(score[zi]) - np.median(score[ni])))
    d0 = float(np.mean(draws))
    rng2 = np.random.default_rng(seed + 2)
    bt = np.empty(b)
    for i in range(b):
        zi, ni = pairs[int(rng2.integers(0, len(pairs)))]
        s = rng2.integers(0, zi.size, zi.size)
        bt[i] = np.median(score[zi[s]]) - np.median(score[ni[s]])
    lo, hi = np.percentile(bt, [2.5, 97.5])
    zi, ni = pairs[0]
    return dict(n_matched=int(zi.size), n_znf_eligible=int(idx_z.size),
                n_dropped=int(idx_z.size - zi.size), n_draws=n_draw,
                d_median=d0, ci=[float(lo), float(hi)],
                per_draw_min=float(np.min(draws)), per_draw_max=float(np.max(draws)),
                auc_first_draw=auc_lower(score[zi], score[ni]),
                excludes_zero_negative=bool(hi < 0), excludes_zero_positive=bool(lo > 0))


def submap(Ccross: np.ndarray, pairs: np.ndarray, keys: list[str], want: list[tuple[str, str]],
           require_all: bool = False) -> np.ndarray:
    """mean cosine over a named subset of source pairs."""
    ki = {k: i for i, k in enumerate(keys)}
    sel = []
    for a, b in want:
        ia, ib = ki[a], ki[b]
        for p, (i, j) in enumerate(pairs):
            if (i, j) == (ia, ib) or (i, j) == (ib, ia):
                sel.append(p)
    S = Ccross[sel]
    with np.errstate(invalid="ignore"):
        out = np.nanmean(S, axis=0)
    if require_all:
        out = np.where(np.isfinite(S).all(axis=0), out, np.nan)
    return out


def main() -> None:
    t0 = time.time()
    z = np.load(WORK / "held279.npz", allow_pickle=True)
    genes = [str(g) for g in z["genes"]]
    keys = [str(k) for k in z["keys"]]
    c = np.asarray(z["c_held"], dtype=np.float64)
    Ccross = np.asarray(z["Ccross"], dtype=np.float64)
    pairs = np.asarray(z["pairs_cross"])
    mag = np.asarray(z["mag"], dtype=np.float64)
    expr = np.asarray(z["expr"], dtype=np.float64)
    n_pairs = np.asarray(z["n_pairs"])
    n_lines = np.asarray(z["n_lines"])
    cls, groups = classify(genes)
    is_z = cls == "ZNF"; is_n = cls == "NONZNF"; is_o = cls == "OTHER_ZF"
    print(f"held {len(genes)}: ZNF {int(is_z.sum())}, non-ZNF {int(is_n.sum())}, "
          f"other-ZF {int(is_o.sum())}", flush=True)

    out: dict = dict(n_held=len(genes), n_znf=int(is_z.sum()), n_nonznf=int(is_n.sum()),
                     n_otherzf=int(is_o.sum()))

    # ---------- primary ------------------------------------------------------------------------
    fin = np.isfinite(c)
    x, y = c[is_z & fin], c[is_n & fin]
    prim = two_sample(x, y)
    prim["auc_ci"] = auc_ci(x, y)
    out["primary"] = prim
    print(f"PRIMARY  ZNF median {prim['median_x']:+.5f} (n {prim['n_x']}) vs non-ZNF "
          f"{prim['median_y']:+.5f} (n {prim['n_y']}); Δ {prim['d_median']:+.5f} "
          f"[{prim['ci'][0]:+.5f}, {prim['ci'][1]:+.5f}]  AUC {prim['auc']:.3f} "
          f"{prim['auc_ci']}  MWU p {prim['mwu_p_one_sided']:.2e}", flush=True)

    # confound descriptives
    out["confounds"] = {}
    for nm, v in (("expr_ccle_median", expr), ("mag", mag), ("n_pairs", n_pairs.astype(float)),
                  ("n_lines", n_lines.astype(float))):
        vz = v[is_z & np.isfinite(v)]; vn = v[is_n & np.isfinite(v)]
        out["confounds"][nm] = dict(znf_median=float(np.median(vz)), nonznf_median=float(np.median(vn)),
                                    mwu_p=float(mannwhitneyu(vz, vn, alternative="two-sided").pvalue))
        print(f"  confound {nm:18s} ZNF {np.median(vz):+.4f} vs non-ZNF {np.median(vn):+.4f}")

    # ---------- matched ------------------------------------------------------------------------
    de, dm = deciles(expr), deciles(mag)
    cell = np.where((de >= 0) & (dm >= 0), de * 10 + dm, -1)
    mt = matched_effect(c, cell, is_z & fin, is_n & fin)
    out["matched"] = mt
    print(f"MATCHED  n {mt['n_matched']} pairs (dropped {mt['n_dropped']}), Δ {mt['d_median']:+.5f} "
          f"[{mt['ci'][0]:+.5f}, {mt['ci'][1]:+.5f}], per-draw [{mt['per_draw_min']:+.5f}, "
          f"{mt['per_draw_max']:+.5f}], AUC {mt['auc_first_draw']:.3f}", flush=True)

    # ---------- verdict ------------------------------------------------------------------------
    confirmed = bool(prim["d_median"] < 0 and prim["excludes_zero_negative"]
                     and mt["d_median"] < 0 and mt["excludes_zero_negative"])
    out["verdict"] = dict(confirmed=confirmed,
                          rule="raw Δmed<0 with 95% CI<0 AND matched Δmed<0 with 95% CI<0")
    print(f"VERDICT: {'CONFIRMED' if confirmed else 'NOT CONFIRMED'}", flush=True)

    # ---------- control (a): random size-matched families + the HGNC family ranking -------------
    rng = np.random.default_rng(SEED)
    pool = np.flatnonzero(is_n & fin)
    nz = int((is_z & fin).sum())
    rnd = np.asarray([np.median(c[rng.choice(pool, nz, replace=False)]) - np.median(c[pool])
                      for _ in range(N_RANDFAM)])
    out["control_random_family"] = dict(n_draws=N_RANDFAM, mean=float(rnd.mean()),
                                        band=[float(v) for v in np.percentile(rnd, [2.5, 97.5])],
                                        observed=prim["d_median"],
                                        share_below_observed=float(np.mean(rnd <= prim["d_median"])))
    print(f"CONTROL (a) random families: band {out['control_random_family']['band']}, "
          f"observed {prim['d_median']:+.5f}, share of draws below it "
          f"{out['control_random_family']['share_below_observed']:.3f}", flush=True)

    fam_members: dict[str, list[int]] = {}
    for i, g in enumerate(genes):
        if not fin[i]:
            continue
        for f in groups.get(g, set()):
            fam_members.setdefault(f, []).append(i)
    fam_rows = []
    base_nonzf = c[is_n & fin]
    for f, mem in fam_members.items():
        if len(mem) < 30:
            continue
        mem = np.asarray(mem)
        rest = c[np.setdiff1d(np.flatnonzero(is_n & fin), mem, assume_unique=False)]
        fam_rows.append(dict(family=f, n=int(mem.size),
                             d_median=float(np.median(c[mem]) - np.median(rest)),
                             median=float(np.median(c[mem]))))
    fam_rows.sort(key=lambda r: r["d_median"])
    znf_like = [r for r in fam_rows if r["family"] in ZNF_GROUPS]
    out["family_ranking"] = dict(n_families=len(fam_rows), rows=fam_rows[:15],
                                 rows_top=fam_rows[-10:],
                                 znf_rows=[dict(rank=fam_rows.index(r) + 1, **r) for r in znf_like])
    print(f"CONTROL (a2) {len(fam_rows)} HGNC families with >=30 held genes; ZNF families rank "
          f"{[ (fam_rows.index(r)+1) for r in znf_like ]} of {len(fam_rows)} (1 = lowest)", flush=True)

    # ---------- control (b): the other TF families ---------------------------------------------
    def fam_mask(pred) -> np.ndarray:
        m = np.zeros(len(genes), bool)
        for i, g in enumerate(genes):
            if any(pred(f) for f in groups.get(g, set())):
                m[i] = True
        return m

    homeo = fam_mask(lambda f: "homeobox" in f.lower())
    bhlh = fam_mask(lambda f: f == "Basic helix-loop-helix proteins")
    out["control_tf_families"] = {}
    for nm, m in (("homeobox", homeo), ("bHLH", bhlh)):
        mm = m & is_n & fin           # they are non-ZNF by construction; keep the comparison clean
        r = two_sample(c[mm], c[is_n & fin & ~m], b=2000)
        r["matched"] = matched_effect(c, cell, mm, is_n & fin & ~m, n_draw=20, b=1000)
        out["control_tf_families"][nm] = r
        print(f"CONTROL (b) {nm:9s} n {int(mm.sum()):4d} Δ {r['d_median']:+.5f} "
              f"[{r['ci'][0]:+.5f}, {r['ci'][1]:+.5f}] AUC {r['auc']:.3f}; matched "
              f"{r['matched']['d_median']:+.5f} {r['matched']['ci']}", flush=True)

    # ---------- control (c): single-dataset replications ----------------------------------------
    out["control_single_dataset"] = {}
    variants = {
        "c1_orion_hct116_x_hek293t": ([("orion_hct116", "orion_hek293t")], False),
        "c2_replogle_k562_x_rpe1": ([("k562_gwps", "rpe1")], False),
        "c3_balanced_core": ([("k562_gwps", "orion_hct116"), ("k562_gwps", "orion_hek293t"),
                              ("orion_hct116", "orion_hek293t")], True),
    }
    for nm, (want, req) in variants.items():
        s = submap(Ccross, pairs, keys, want, require_all=req)
        f = np.isfinite(s)
        r = two_sample(s[is_z & f], s[is_n & f], b=2000)
        r["matched"] = matched_effect(s, cell, is_z & f, is_n & f, n_draw=20, b=1000)
        out["control_single_dataset"][nm] = r
        print(f"CONTROL (c) {nm:28s} n {r['n_x']:4d}/{r['n_y']:6d} Δ {r['d_median']:+.5f} "
              f"[{r['ci'][0]:+.5f}, {r['ci'][1]:+.5f}] AUC {r['auc']:.3f}; matched "
              f"{r['matched']['d_median']:+.5f} {r['matched']['ci']}", flush=True)

    # ---------- declared secondaries ------------------------------------------------------------
    krab = fam_mask(lambda f: f == "KRAB domain containing") & fin
    symb = np.asarray([symbol_znf(g) for g in genes]) & fin & ~is_o
    out["secondary"] = {}
    for nm, m in (("KRAB_only", krab), ("symbol_rule_276", symb)):
        r = two_sample(c[m], c[is_n & fin & ~m], b=2000)
        r["matched"] = matched_effect(c, cell, m, is_n & fin & ~m, n_draw=20, b=1000)
        out["secondary"][nm] = r
        print(f"SECONDARY {nm:16s} n {int(m.sum()):4d} Δ {r['d_median']:+.5f} "
              f"[{r['ci'][0]:+.5f}, {r['ci'][1]:+.5f}] AUC {r['auc']:.3f}; matched "
              f"{r['matched']['d_median']:+.5f} {r['matched']['ci']}", flush=True)
    # without nadig_jurkat (challenge line A)
    nj = keys.index("nadig_jurkat")
    keep = np.asarray([p for p, (i, j) in enumerate(pairs) if i != nj and j != nj])
    with np.errstate(invalid="ignore"):
        c_nj = np.nanmean(Ccross[keep], axis=0)
    f = np.isfinite(c_nj)
    out["secondary"]["no_nadig_jurkat"] = two_sample(c_nj[is_z & f], c_nj[is_n & f], b=2000)
    print(f"SECONDARY no_nadig_jurkat Δ {out['secondary']['no_nadig_jurkat']['d_median']:+.5f} "
          f"{out['secondary']['no_nadig_jurkat']['ci']}", flush=True)

    # ---------- in-sample reference (the 300 panel targets, S279's sources and axis) -------------
    axis = [str(g) for g in z["axis"]]
    per_src = {k: I.CANDIDATES[k][0] for k in keys}
    idxs, Us = {}, {}
    for k in keys:
        i, u, _ = B.load_unit(per_src[k], axis, L.PANEL)
        idxs[k], Us[k] = i, u
    lines = [str(x) for x in z["lines"]]
    Cp = np.full((len(pairs), len(L.PANEL)), np.nan)
    for p, (i, j) in enumerate(pairs):
        g, cc = B.pair_cos(idxs[keys[i]], Us[keys[i]], idxs[keys[j]], Us[keys[j]])
        Cp[p, g] = cc
    with np.errstate(invalid="ignore"):
        cp = np.nanmean(Cp, axis=0)
    clsp, _ = classify(L.PANEL)
    zp = (clsp == "ZNF") & np.isfinite(cp)
    npz_ = (clsp == "NONZNF") & np.isfinite(cp)
    ins = two_sample(cp[zp], cp[npz_], b=2000)
    ins["auc_ci"] = auc_ci(cp[zp], cp[npz_], b=1000)
    symbol_p = np.asarray([symbol_znf(t) for t in L.PANEL]) & np.isfinite(cp)
    ins_symbol = two_sample(cp[symbol_p], cp[np.isfinite(cp) & ~symbol_p], b=2000)
    out["in_sample_panel"] = dict(hgnc_rule=ins, symbol_rule=ins_symbol,
                                  n_znf_hgnc=int(zp.sum()), n_znf_symbol=int(symbol_p.sum()),
                                  median_panel=float(np.nanmedian(cp)))
    print(f"IN-SAMPLE panel (S279 sources/axis): HGNC ZNF n {int(zp.sum())} Δ {ins['d_median']:+.5f} "
          f"{ins['ci']} AUC {ins['auc']:.3f}; symbol rule n {int(symbol_p.sum())} "
          f"Δ {ins_symbol['d_median']:+.5f} {ins_symbol['ci']}", flush=True)

    # ---------- the per-gene table ---------------------------------------------------------------
    with open(WORK / "znf279.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "class", "krab", "symbol_rule", "C_held", "n_pairs", "n_lines",
                    "C_within", "mag", "ccle_median", "expr_decile", "mag_decile",
                    "null_mean", "c1_orion_pair", "c3_balanced_core"])
        s1 = submap(Ccross, pairs, keys, [("orion_hct116", "orion_hek293t")])
        s3 = submap(Ccross, pairs, keys, variants["c3_balanced_core"][0], require_all=True)
        cw = np.asarray(z["c_within"]); nm_ = np.asarray(z["null_mean"])
        for i, g in enumerate(genes):
            w.writerow([g, cls[i], int(krab[i]), int(symb[i]), f"{c[i]:.6f}", int(n_pairs[i]),
                        int(n_lines[i]), "" if not np.isfinite(cw[i]) else f"{cw[i]:.6f}",
                        f"{mag[i]:.6f}" if np.isfinite(mag[i]) else "",
                        f"{expr[i]:.4f}" if np.isfinite(expr[i]) else "",
                        int(de[i]), int(dm[i]),
                        f"{nm_[i]:.6f}" if np.isfinite(nm_[i]) else "",
                        f"{s1[i]:.6f}" if np.isfinite(s1[i]) else "",
                        f"{s3[i]:.6f}" if np.isfinite(s3[i]) else ""])

    out["seconds"] = time.time() - t0
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
    (WORK / "analyse279.json").write_text(json.dumps(out, indent=1))
    print(f"\ndone in {out['seconds']:.0f}s, peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
