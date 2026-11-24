"""284 — THE REMAINDER TEST (PREREG §7): is the gene-level remainder reproducible across two source
sets that share no source and no cell line?

  PYTHONDONTWRITEBYTECODE=1 uv run python rem284.py

Outputs: work/rem284.json, work/rem284.log, work/rem284.npz
"""
from __future__ import annotations

import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd

import lib284 as L

T0 = time.time()
LOG: list[str] = []
N_PART = 20
N_NULL = 200


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


# ---- the two out-of-fold predictors, both fitted inside each fold's training rows --------------
def oof_ridge(X: np.ndarray, y: np.ndarray, fold: np.ndarray) -> np.ndarray:
    """Grouped out-of-fold ridge: folds are whole components, so a test gene shares no family with
    its training rows (PREREG §6). alpha re-chosen by 5-fold inner CV inside the training rows."""
    p = np.full(len(y), np.nan)
    ok = np.isfinite(y)
    for k in range(int(fold.max()) + 1):
        te = ok & (fold == k)
        tr = ok & (fold != k)
        if te.sum() == 0 or tr.sum() < 50:
            continue
        p[te], _ = L.fit_ridge(X[tr], y[tr], X[te])
    return p


def oof_family_median(mem: list[list[str]], y: np.ndarray, fold: np.ndarray) -> np.ndarray:
    """Out-of-fold family-median-only predictor: the mean over the gene's HGNC families of the
    median y of that family's TRAINING genes; the training global median where none exists."""
    p = np.full(len(y), np.nan)
    ok = np.isfinite(y)
    for k in range(int(fold.max()) + 1):
        te = np.flatnonzero(ok & (fold == k))
        tri = np.flatnonzero(ok & (fold != k))
        if len(te) == 0 or len(tri) == 0:
            continue
        acc: dict[str, list[float]] = {}
        for i in tri:
            for f in mem[i]:
                acc.setdefault(f, []).append(float(y[i]))
        med = {f: float(np.median(v)) for f, v in acc.items()}
        gmed = float(np.median(y[tri]))
        for i in te:
            vals = [med[f] for f in mem[i] if f in med]
            p[i] = float(np.mean(vals)) if vals else gmed
    return p


def spearman_brown(r: float) -> float:
    return float(2.0 * r / (1.0 + r)) if np.isfinite(r) and (1.0 + r) != 0 else float("nan")


def main() -> None:
    h = hashlib.sha256((L.HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    if h != (L.WORK / "PREREGISTRATION.sha256").read_text().split()[0]:
        sys.exit("preregistration digest mismatch — aborting")
    P(f"prereg {h} MATCH")

    held = L.load_held()
    mem = L.memberships(held)
    sp = np.load(L.WORK / "splits284.npz", allow_pickle=True)
    comp = sp["comp"]
    fold = L.component_folds(comp, 5)
    X = held[L.FEATS].to_numpy(float)
    rbar = held["r_bar"].to_numpy(float)

    z = np.load(L.HELD279, allow_pickle=True)
    Ccross = np.asarray(z["Ccross"])
    keys = np.asarray(z["keys"]).astype(str)
    lines = np.asarray(z["lines"]).astype(str)
    pairs = np.asarray(z["pairs_cross"])
    c_held = np.asarray(z["c_held"])
    P(f"{len(keys)} sources, {len(set(lines))} lines, {len(pairs)} cross-line pairs; "
      f"fold sizes {np.bincount(fold).tolist()}")

    # how lines repeat across sources
    rep: dict[str, list[str]] = {}
    for k, ln in zip(keys, lines):
        rep.setdefault(ln, []).append(k)
    multi = {ln: v for ln, v in rep.items() if len(v) > 1}
    P("  lines with >1 source: " + "; ".join(f"{ln}({len(v)})" for ln, v in sorted(multi.items())))

    out: dict = {"prereg_sha256": h, "sources": len(keys), "lines": int(len(set(lines))),
                 "lines_multi_source": {k: v for k, v in multi.items()}, "partitions": []}

    line_groups = sorted(set(lines.tolist()))
    src_of_line = {ln: np.flatnonzero(lines == ln) for ln in line_groups}

    def make_partition(rng) -> tuple[np.ndarray, np.ndarray, list[str], list[str]] | None:
        assign = rng.integers(0, 2, size=len(line_groups))
        h1 = [ln for ln, a in zip(line_groups, assign) if a == 0]
        h2 = [ln for ln, a in zip(line_groups, assign) if a == 1]
        if len(h1) < 8 or len(h2) < 8:
            return None
        s1 = np.concatenate([src_of_line[ln] for ln in h1])
        s2 = np.concatenate([src_of_line[ln] for ln in h2])
        p1 = np.array([q for q, (i, j) in enumerate(pairs) if i in set(s1) and j in set(s1)])
        p2 = np.array([q for q, (i, j) in enumerate(pairs) if i in set(s2) and j in set(s2)])
        if len(p1) < 20 or len(p2) < 20:
            return None
        return p1, p2, h1, h2

    rng = np.random.default_rng(L.SEED)
    parts = []
    tries = 0
    while len(parts) < N_PART and tries < 10000:
        tries += 1
        r = make_partition(rng)
        if r is not None:
            parts.append(r)
    P(f"drew {len(parts)} admissible partitions in {tries} tries")

    def analyse(p1, p2, h1, h2, primary: bool, tag: str) -> dict:
        C1 = L.cmap(Ccross, p1).astype(float)
        C2 = L.cmap(Ccross, p2).astype(float)
        n1 = np.isfinite(Ccross[p1]).sum(axis=0)
        n2 = np.isfinite(Ccross[p2]).sum(axis=0)
        both = (n1 >= 1) & (n2 >= 1)
        pf1, pf2 = oof_ridge(X, C1, fold), oof_ridge(X, C2, fold)
        pm1, pm2 = oof_family_median(mem, C1, fold), oof_family_median(mem, C2, fold)
        r1, r2 = C1 - pf1, C2 - pf2
        m1, m2 = C1 - pm1, C2 - pm2
        ok = both & np.isfinite(r1) & np.isfinite(r2)
        okm = both & np.isfinite(m1) & np.isfinite(m2)
        d = dict(tag=tag, half1_lines=h1, half2_lines=h2, n_pairs_1=int(len(p1)),
                 n_pairs_2=int(len(p2)), n_both=int(both.sum()), n_ok=int(ok.sum()),
                 ceiling=L.spearman(C1[both], C2[both]),
                 rho_remainder=L.spearman(r1[ok], r2[ok]),
                 rho_remainder_familymedian=L.spearman(m1[okm], m2[okm]),
                 oof_rho_half1=L.spearman(C1[np.isfinite(pf1)], pf1[np.isfinite(pf1)]),
                 oof_rho_half2=L.spearman(C2[np.isfinite(pf2)], pf2[np.isfinite(pf2)]))
        if not primary:
            return d

        rb = np.random.default_rng(L.SEED + 301)
        d["boot_remainder"] = L.boot_rho(r1[ok], r2[ok], comp[ok], rb)
        d["boot_ceiling"] = L.boot_rho(C1[both], C2[both], comp[both],
                                       np.random.default_rng(L.SEED + 302))
        d["boot_remainder_familymedian"] = L.boot_rho(m1[okm], m2[okm], comp[okm],
                                                      np.random.default_rng(L.SEED + 303))
        # attenuation: sub-halves of each half's PAIRS, remainder against the same OOF prediction
        rr = np.random.default_rng(L.SEED + 401)
        att = {}
        for nm, pp, pf, CC in (("half1", p1, pf1, C1), ("half2", p2, pf2, C2)):
            q = rr.permutation(len(pp))
            a, b = pp[q[: len(pp) // 2]], pp[q[len(pp) // 2:]]
            Ca, Cb = L.cmap(Ccross, a).astype(float), L.cmap(Ccross, b).astype(float)
            ra, rbb = Ca - pf, Cb - pf
            m = np.isfinite(ra) & np.isfinite(rbb) & ok
            rho_half = L.spearman(ra[m], rbb[m])
            att[nm] = dict(n_pairs_a=int(len(a)), n_pairs_b=int(len(b)), n=int(m.sum()),
                           rho_subhalf=rho_half, rho_full=spearman_brown(rho_half))
        denom = float(np.sqrt(att["half1"]["rho_full"] * att["half2"]["rho_full"]))
        d["attenuation"] = dict(**att, denominator=denom,
                                corrected=float(d["rho_remainder"] / denom) if denom > 0 else
                                float("nan"))
        # null: permute remainder2 WITHIN components
        rn = np.random.default_rng(L.SEED + 501)
        idx = np.flatnonzero(ok)
        cc = comp[idx]
        order = np.argsort(cc, kind="stable")
        bounds = np.flatnonzero(np.diff(cc[order])) + 1
        groups = np.split(order, bounds)
        nulls = []
        base2 = r2[idx]
        for _ in range(N_NULL):
            perm = np.arange(len(idx))
            for g in groups:
                if len(g) > 1:
                    perm[g] = g[rn.permutation(len(g))]
            nulls.append(L.spearman(r1[idx], base2[perm]))
        nulls = np.asarray(nulls)
        d["null"] = dict(mean=float(np.nanmean(nulls)), sd=float(np.nanstd(nulls)),
                         n_reaching=int((nulls >= d["rho_remainder"]).sum()), n_draws=N_NULL,
                         z=float((d["rho_remainder"] - np.nanmean(nulls)) / np.nanstd(nulls)))
        # descriptive: by r_bar tercile and by min pair count
        strat = {}
        rb_ok = ok & np.isfinite(rbar)
        if rb_ok.sum() > 100:
            q = np.nanquantile(rbar[rb_ok], [1 / 3, 2 / 3])
            for lab, m in (("r_bar low", rb_ok & (rbar <= q[0])),
                           ("r_bar mid", rb_ok & (rbar > q[0]) & (rbar <= q[1])),
                           ("r_bar high", rb_ok & (rbar > q[1]))):
                strat[lab] = dict(n=int(m.sum()), rho=L.spearman(r1[m], r2[m]),
                                  ceiling=L.spearman(C1[m], C2[m]))
        npmin = np.minimum(n1, n2)
        for lab, m in (("min pairs 1", ok & (npmin == 1)), ("min pairs 2-4", ok & (npmin >= 2) &
                                                            (npmin <= 4)),
                       ("min pairs >=5", ok & (npmin >= 5))):
            strat[lab] = dict(n=int(m.sum()), rho=L.spearman(r1[m], r2[m]),
                              ceiling=L.spearman(C1[m], C2[m]))
        d["strata"] = strat
        d["_arrays"] = dict(C1=C1, C2=C2, r1=r1, r2=r2, m1=m1, m2=m2, ok=ok, n1=n1, n2=n2)
        return d

    # -------- primary partition --------
    p1, p2, h1, h2 = parts[0]
    P(f"\n=== PRIMARY partition #1: half1 {len(h1)} lines / {len(p1)} pairs, half2 {len(h2)} "
      f"lines / {len(p2)} pairs ===")
    P(f"  half1 lines: {h1}")
    P(f"  half2 lines: {h2}")
    prim = analyse(p1, p2, h1, h2, True, "partition#1")
    arrays = prim.pop("_arrays")
    P(f"  genes with >=1 pair in both halves: {prim['n_both']}; usable {prim['n_ok']}")
    P(f"  CEILING   rho(C1,C2)   {prim['boot_ceiling']['point']:+.4f} "
      f"[{prim['boot_ceiling']['lo']:+.4f}, {prim['boot_ceiling']['hi']:+.4f}]")
    P(f"  REMAINDER (features)  {prim['boot_remainder']['point']:+.4f} "
      f"[{prim['boot_remainder']['lo']:+.4f}, {prim['boot_remainder']['hi']:+.4f}] "
      f"(sd {prim['boot_remainder']['sd']:.4f}, n {prim['boot_remainder']['n']})")
    P(f"  REMAINDER (family median only) {prim['boot_remainder_familymedian']['point']:+.4f} "
      f"[{prim['boot_remainder_familymedian']['lo']:+.4f}, "
      f"{prim['boot_remainder_familymedian']['hi']:+.4f}]")
    P(f"  attenuation: half1 sub-half rho {prim['attenuation']['half1']['rho_subhalf']:+.4f} -> SB "
      f"{prim['attenuation']['half1']['rho_full']:+.4f}; half2 "
      f"{prim['attenuation']['half2']['rho_subhalf']:+.4f} -> "
      f"{prim['attenuation']['half2']['rho_full']:+.4f}; corrected "
      f"{prim['attenuation']['corrected']:+.4f}")
    P(f"  null (permute within components): {prim['null']['mean']:+.4f} ± {prim['null']['sd']:.4f}, "
      f"{prim['null']['n_reaching']}/{N_NULL} reach the observed, z {prim['null']['z']:+.1f}")
    for k, v in prim["strata"].items():
        P(f"    {k:16s} n {v['n']:6d}  remainder rho {v['rho']:+.4f}  ceiling {v['ceiling']:+.4f}")
    out["primary"] = prim

    # -------- 20 partitions --------
    P(f"\n=== {len(parts)} admissible partitions ===")
    for i, (a, b, ha, hb) in enumerate(parts):
        d = analyse(a, b, ha, hb, False, f"partition#{i + 1}")
        out["partitions"].append(d)
        P(f"  #{i + 1:2d} lines {len(ha):2d}/{len(hb):2d} pairs {len(a):3d}/{len(b):3d} "
          f"n {d['n_ok']:6d}  ceiling {d['ceiling']:+.4f}  remainder {d['rho_remainder']:+.4f}  "
          f"remainder(fam-median) {d['rho_remainder_familymedian']:+.4f}")
    rem = np.array([d["rho_remainder"] for d in out["partitions"]])
    cei = np.array([d["ceiling"] for d in out["partitions"]])
    remm = np.array([d["rho_remainder_familymedian"] for d in out["partitions"]])
    out["distribution"] = dict(
        remainder=dict(median=float(np.median(rem)), mean=float(rem.mean()), sd=float(rem.std()),
                       lo=float(rem.min()), hi=float(rem.max()),
                       n_positive=int((rem > 0).sum()), n_ge_0p05=int((rem >= 0.05).sum())),
        remainder_familymedian=dict(median=float(np.median(remm)), lo=float(remm.min()),
                                    hi=float(remm.max()), n_positive=int((remm > 0).sum())),
        ceiling=dict(median=float(np.median(cei)), lo=float(cei.min()), hi=float(cei.max())))
    P(f"  remainder over {len(parts)} partitions: median {np.median(rem):+.4f} "
      f"[{rem.min():+.4f}, {rem.max():+.4f}], {int((rem > 0).sum())}/{len(parts)} > 0, "
      f"{int((rem >= 0.05).sum())}/{len(parts)} >= 0.05")
    P(f"  ceiling   over {len(parts)} partitions: median {np.median(cei):+.4f} "
      f"[{cei.min():+.4f}, {cei.max():+.4f}]")

    # -------- the frozen pass rule --------
    lo, hi = prim["boot_remainder"]["lo"], prim["boot_remainder"]["hi"]
    pt = prim["boot_remainder"]["point"]
    excludes0 = bool(lo > 0 or hi < 0)
    verdict = ("GENE PROPERTY the features do not predict" if (excludes0 and pt >= 0.05)
               else "not distinguishable from measurement noise at this depth" if not excludes0
               else "reproducible but below the 0.05 threshold")
    med = float(np.median(rem))
    equivocal = bool(np.sign(med) != np.sign(pt) or ((med >= 0.05) != (pt >= 0.05)))
    out["verdict"] = dict(point=pt, lo=lo, hi=hi, excludes_zero=excludes0,
                          partition_median=med, equivocal=equivocal, verdict=verdict)
    P(f"\nPASS RULE: point {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]; excludes 0 {excludes0}; "
      f"20-partition median {med:+.4f}; EQUIVOCAL {equivocal} -> {verdict}")

    np.savez_compressed(L.WORK / "rem284.npz", **{k: v for k, v in arrays.items()},
                        gene=held["gene"].to_numpy().astype(str), fold=fold, comp=comp)
    out["log"] = LOG
    (L.WORK / "rem284.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "rem284.log").write_text("\n".join(LOG))
    P("wrote work/rem284.json, work/rem284.npz")


if __name__ == "__main__":
    main()
