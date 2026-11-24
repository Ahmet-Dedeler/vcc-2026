"""218 Phase 0c — the LOAD-BEARING regression: does INDEPENDENCE predict CONTRIBUTION?

Joins `resid218.json` (what each source is worth as an un-whitened post-whitening residual,
in relative aggregate-cosine gain at a matched dose) to `independence218.json` (how
independent that source is from the production mixture, measured in the same frame), and
regresses one on the other.

Two rival predictors are fitted on the SAME points, which is the whole design:

  H1 (the mandate's)  contribution ~ INDEPENDENCE
                      predictors: orth_frac, 1 - sub_10, -excess_pre, and the categorical
                      platform/lineage/transformation/lab novelty score vs Replogle+Orion.
  H0 (the boring one) contribution ~ TRANSFER ALIGNMENT to the recipient
                      predictors: cos_truth, excess_truth.

A predictor "wins" only if it (a) has the right sign, (b) explains more variance, and
(c) still orders the three MAXIMALLY-independent sources (CD4, KOLF2, Feng iPSC), which
all score identically on every categorical independence axis.  That last test is the one
that can kill H1 on its own: if three sources that are equally independent contribute very
differently, independence cannot be the mechanism.

    uv run python regress218.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from sources218 import SOURCES  # noqa: E402

PRIMARY_F = 0.10          # the shipped CD4 dose: every source is compared at matched `f`
DOSES = (0.05, 0.10, 0.20, 0.34, 0.55, 0.80)
SEED = 20260904

# Categorical novelty against the production mixture {Replogle K562, Orion HCT116,
# Orion HEK293T}.  One point each for a platform, a lineage, a transformation state and a
# lab that mixture does not already contain.
NOVELTY = {
    "cd4":            dict(platform=1, lineage=1, transformation=1, lab=1),
    "kolf2":          dict(platform=1, lineage=1, transformation=1, lab=1),
    "feng_ipsc":      dict(platform=1, lineage=1, transformation=1, lab=1),
    "telohaec":       dict(platform=1, lineage=1, transformation=1, lab=1),
    "gse295214_h1":   dict(platform=1, lineage=1, transformation=1, lab=1),
    "k562_gwps":      dict(platform=0, lineage=0, transformation=0, lab=0),
    "orion_hct116":   dict(platform=0, lineage=0, transformation=0, lab=0),
    "orion_hek293t":  dict(platform=0, lineage=0, transformation=0, lab=0),
    "perturbmatch":   dict(platform=0, lineage=0, transformation=0, lab=0),
    "self_pre":       dict(platform=0, lineage=0, transformation=0, lab=0),
}
# A source that IS the held line (or the same cells as it) is an oracle, not a donor.
LEAK = {("k562_gwps", "perturbmatch")}          # (fold, source): same K562 cells


def arm_name(src: str, mode: str, f: float) -> str:
    if src == "cd4":
        return f"cd4_{'gn' if mode == 'gene_norm' else 'raw'}{f:g}"
    tag = "gn" if mode == "gene_norm" else "rw"
    return f"{src}_{tag}{f:g}"


def ols(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4 or x.std() == 0:
        return dict(n=int(len(x)), slope=None, r=None, r2=None, p_sign=None)
    A = np.vstack([x, np.ones_like(x)]).T
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r = float(np.corrcoef(x, y)[0, 1])
    # bootstrap over points for the sign of the slope
    rng = np.random.default_rng(SEED)
    slopes = []
    for _ in range(4000):
        idx = rng.integers(0, len(x), len(x))
        if x[idx].std() == 0:
            continue
        b = np.polyfit(x[idx], y[idx], 1)[0]
        slopes.append(b)
    slopes = np.asarray(slopes)
    return dict(n=int(len(x)), slope=float(beta[0]), intercept=float(beta[1]),
                r=r, r2=1.0 - ss_res / ss_tot if ss_tot > 0 else None,
                p_slope_gt0=float((slopes > 0).mean()),
                slope_ci=[float(np.percentile(slopes, 2.5)),
                          float(np.percentile(slopes, 97.5))])


def main() -> None:
    resid = json.loads((WORK / "resid218.json").read_text())
    ind = json.loads((WORK / "independence218.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in resid}
    indep = {(r["fold"], r["source"]): r for r in ind["rows"]}
    folds = sorted({r["fold"] for r in resid})

    # ---------------------------------------------------------------- contribution table
    rows = []
    srcs = ["cd4"] + [s for s in SOURCES if s in {r["source"] for r in ind["rows"]}]
    srcs = list(dict.fromkeys(srcs)) + ["self_pre"]
    for fold in folds:
        base = by.get((fold, "f0"))
        if base is None:
            continue
        b_cos, b_pds, b_q = base["cos_agg"], base["pds"], base["a15_sgn_w"]
        for src in srcs:
            for mode in ("gene_norm", "raw"):
                for f in DOSES:
                    if src == "self_pre":
                        arm = f"self{f:g}"
                        if mode != "gene_norm":
                            continue
                    else:
                        arm = arm_name(src, mode, f)
                    r = by.get((fold, arm))
                    if r is None:
                        continue
                    a = r["audit"]
                    rows.append(dict(
                        fold=fold, source=src, mode=mode, f=f, arm=arm,
                        vacuous=bool(a.get("VACUOUS")),
                        leak=bool((fold, src) in LEAK or src == fold),
                        carrying=a["carrying_rows"], coverage=a["coverage_frac"],
                        fraction=a["fraction"],
                        d_cos=r["cos_agg"] - b_cos,
                        rel_cos=(r["cos_agg"] - b_cos) / b_cos if b_cos else None,
                        d_pds=r["pds"] - b_pds, d_qw=r["a15_sgn_w"] - b_q,
                        cos_agg=r["cos_agg"], pds=r["pds"]))
    (WORK / "contrib218.json").write_text(json.dumps(rows, indent=2))

    # -------------------------------------------------------- the headline: per source
    print("=" * 108)
    print(f"CONTRIBUTION at the matched shipped dose f = {PRIMARY_F} "
          f"(un-whitened post-whitening residual, gene_norm)")
    print("=" * 108)
    print(f"{'source':16s} {'folds':>5s} {'rel dcos %':>11s} {'d cos_agg':>11s} "
          f"{'d pds':>10s} {'d q_w':>9s} {'vs CD4':>8s}   note")
    hi = {}
    for src in srcs:
        sel = [r for r in rows if r["source"] == src and r["f"] == PRIMARY_F
               and r["mode"] == "gene_norm" and not r["leak"]]
        if not sel:
            continue
        vac = all(r["vacuous"] for r in sel)
        rel = float(np.mean([r["rel_cos"] for r in sel]))
        hi[src] = dict(rel=rel, n=len(sel),
                       d_cos=float(np.mean([r["d_cos"] for r in sel])),
                       d_pds=float(np.mean([r["d_pds"] for r in sel])),
                       d_qw=float(np.mean([r["d_qw"] for r in sel])),
                       vacuous=vac,
                       coverage=float(np.mean([r["coverage"] for r in sel])))
    cd4rel = hi.get("cd4", {}).get("rel", np.nan)
    for src, h in sorted(hi.items(), key=lambda kv: -kv[1]["rel"]):
        note = "VACUOUS (median sizing rule -> fraction 0)" if h["vacuous"] else \
            f"coverage {h['coverage']:.2f}"
        print(f"{src:16s} {h['n']:5d} {100*h['rel']:+11.3f} {h['d_cos']:+11.6f} "
              f"{h['d_pds']:+10.6f} {h['d_qw']:+9.6f} "
              f"{h['rel']/cd4rel if cd4rel else float('nan'):8.3f}   {note}")

    # ---------------------------------------------------------------- the regression
    print()
    print("=" * 108)
    print("H1 (independence) vs H0 (transfer alignment) -- same points, same endpoint")
    print("=" * 108)
    pts = []
    for r in rows:
        if r["f"] != PRIMARY_F or r["mode"] != "gene_norm" or r["leak"] or r["vacuous"]:
            continue
        i = indep.get((r["fold"], r["source"]))
        if i is None:
            continue
        nv = NOVELTY.get(r["source"], {})
        pts.append(dict(**r, novelty=sum(nv.values()),
                        orth_frac=i["orth_frac"], sub_10=i["sub_10"], sub_50=i["sub_50"],
                        excess_pre=i["excess_pre"], cos_pre=i["cos_pre"],
                        cos_truth=i["cos_truth"], excess_truth=i["excess_truth"]))
    y = [p["rel_cos"] for p in pts]
    fits = {}
    for label, key, sign in (
            ("H1  novelty score (0-4)", "novelty", +1),
            ("H1  orth_frac (energy outside consensus row space)", "orth_frac", +1),
            ("H1  1 - sub_10 (outside consensus top-10 PCs)", "sub_10", -1),
            ("H1  -excess_pre (target-specific redundancy, negated)", "excess_pre", -1),
            ("H0  cos_truth (alignment to the recipient's truth)", "cos_truth", +1),
            ("H0  excess_truth (target-SPECIFIC transfer)", "excess_truth", +1),
            ("--  coverage (share of targets the source carries)", "coverage", +1)):
        x = [sign * p[key] for p in pts]
        f = ols(x, y)
        fits[label] = f
        if f["slope"] is None:
            print(f"{label:52s} n={f['n']:3d}  (degenerate)")
            continue
        print(f"{label:52s} n={f['n']:3d}  r={f['r']:+.3f}  r2={f['r2']:+.3f}  "
              f"slope>0 in {100*f['p_slope_gt0']:5.1f}% of bootstraps")

    # ------------------------------------------- the test that can kill H1 on its own
    print()
    print("=" * 108)
    print("THE DECISIVE ROW: three sources with IDENTICAL categorical independence (4/4)")
    print("=" * 108)
    for src in ("cd4", "kolf2", "feng_ipsc"):
        h = hi.get(src)
        i = [p for p in pts if p["source"] == src]
        if not h or not i:
            continue
        print(f"  {src:12s} novelty 4/4  rel dcos {100*h['rel']:+7.3f}%  "
              f"({h['rel']/cd4rel if cd4rel else float('nan'):.2f}x CD4)  "
              f"orth {np.mean([p['orth_frac'] for p in i]):.4f}  "
              f"cos_truth {np.mean([p['cos_truth'] for p in i]):+.5f}  "
              f"excess_truth {np.mean([p['excess_truth'] for p in i]):+.5f}")

    # ------------------------------------------------- the relaxation floor, and matched
    # 207 section 5.3: 099's self-residual -- adding back OUR OWN un-whitened effect -- already
    # buys part of this channel for free.  It is the correct null for an EXTERNAL source, and
    # it is not zero.  A source below it is worse than using nothing new at all.
    floor = hi.get("self_pre", {}).get("rel", 0.0)
    print()
    print("=" * 108)
    print(f"AGAINST THE RIGHT NULL: 099's self-residual (whitening relaxation, free) "
          f"= {100*floor:+.3f}%")
    print("=" * 108)
    for src, h in sorted(hi.items(), key=lambda kv: -kv[1]["rel"]):
        if h["vacuous"]:
            continue
        print(f"  {src:16s} {100*h['rel']:+8.3f}%   excess over floor "
              f"{100*(h['rel']-floor):+8.3f} pp"
              f"{'   <- clears it' if h['rel'] > floor else ''}")

    # -------------------------------------- matched-target alignment as the predictor
    mp = WORK / "matched218.json"
    matched_fit = None
    if mp.exists():
        m = {(r["fold"], r["source"]): r for r in json.loads(mp.read_text())}
        xs, ys, labels = [], [], []
        for p in pts:
            k = m.get((p["fold"], p["source"]))
            if k is None or not k.get("ratio_vs_cd4"):
                continue
            xs.append(k["ratio_vs_cd4"])
            ys.append(p["rel_cos"])
            labels.append((p["fold"], p["source"]))
        matched_fit = ols(xs, ys)
        print()
        print("=" * 108)
        print("H0, CONFOUND-CONTROLLED: contribution vs MATCHED-target alignment "
              "(166 section 4.2's control)")
        print("=" * 108)
        if matched_fit["slope"] is not None:
            print(f"  ratio_vs_cd4 (same rows)                        n={matched_fit['n']:3d}  "
                  f"r={matched_fit['r']:+.3f}  r2={matched_fit['r2']:+.3f}  "
                  f"slope>0 in {100*matched_fit['p_slope_gt0']:5.1f}% of bootstraps")
        # within-fold rank agreement, which removes the fold-level truth-quality term
        rhos = []
        for fold in folds:
            xf = [x for x, l in zip(xs, labels) if l[0] == fold]
            yf = [y for y, l in zip(ys, labels) if l[0] == fold]
            if len(xf) < 4:
                continue
            rx = np.argsort(np.argsort(xf))
            ry = np.argsort(np.argsort(yf))
            rhos.append(float(np.corrcoef(rx, ry)[0, 1]))
        if rhos:
            print(f"  within-fold Spearman rho of contribution on matched alignment: "
                  f"{np.round(rhos, 3).tolist()}  mean {np.mean(rhos):+.3f}")
        # the same, for the independence predictor
        rhos_i = []
        for fold in folds:
            sel = [p for p in pts if p["fold"] == fold]
            if len(sel) < 4:
                continue
            rx = np.argsort(np.argsort([p["orth_frac"] for p in sel]))
            ry = np.argsort(np.argsort([p["rel_cos"] for p in sel]))
            rhos_i.append(float(np.corrcoef(rx, ry)[0, 1]))
        if rhos_i:
            print(f"  within-fold Spearman rho of contribution on INDEPENDENCE (orth_frac): "
                  f"{np.round(rhos_i, 3).tolist()}  mean {np.mean(rhos_i):+.3f}")

    (WORK / "regress218.json").write_text(json.dumps(
        dict(primary_f=PRIMARY_F, headline=hi, fits=fits, points=pts,
             relaxation_floor=floor, matched_fit=matched_fit), indent=2))
    print(f"\nwrote {WORK/'regress218.json'} and {WORK/'contrib218.json'}")


if __name__ == "__main__":
    main()
