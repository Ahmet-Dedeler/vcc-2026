"""308 -- the analysis: paired target bootstraps of every arm and control against the incumbent on the three PLOLO
cancer folds, the arm-minus-control contrasts, the `nmae` leg, the disagreeing-set diagnostics, the |kd| covariate
structure (is |kd| just row norm?), the depth conditioning, and `work/kd308.csv`.

ARITHMETIC ON SAVED ARRAYS ONLY -- no fold is rebuilt, no operator loaded.
    PYTHONDONTWRITEBYTECODE=1 uv run python analyse308.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
B = 20_000
SEED = 20260907
rng = np.random.default_rng(SEED)
LINES: list[str] = []
OUT: dict = dict(seed=SEED, B=B, folds=list(CANCER))


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def boot(diff: np.ndarray) -> dict:
    d = diff[np.isfinite(diff)]
    idx = rng.integers(0, d.size, size=(B, d.size))
    means = d[idx].mean(1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(d.mean()), lo=float(lo), hi=float(hi), sd=float(means.std(ddof=1)),
                p_le0=float((means <= 0).mean()), n=int(d.size))


def boot_pool(diffs: list[np.ndarray]) -> dict:
    acc = np.zeros(B)
    for d in diffs:
        d = d[np.isfinite(d)]
        acc += d[rng.integers(0, d.size, size=(B, d.size))].mean(1)
    means = acc / len(diffs)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return dict(point=float(np.mean([float(np.nanmean(d)) for d in diffs])), lo=float(lo), hi=float(hi),
                sd=float(means.std(ddof=1)), p_le0=float((means <= 0).mean()),
                per_fold=[float(np.nanmean(d)) for d in diffs])


def fb(b: dict) -> str:
    return f"{b['point']:+.5f} [{b['lo']:+.5f}, {b['hi']:+.5f}]"


def main() -> None:
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in CANCER:
        d = json.loads((WORK / f"sweep308-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertarget308-{held}.npz", allow_pickle=True)
        pt[held] = {k: np.asarray(z[k]) for k in z.files}
        zi = np.load(WORK / f"identity308-{held}.npz")
        idn[held] = {k: np.asarray(zi[k]) for k in zi.files}
    arms = [r["arm"] for r in json.loads((WORK / "sweep308-k562_gwps.json").read_text())["rows"]]

    def nmae(arm: str) -> float:
        U = np.concatenate([idn[f][f"{arm}|U"] for f in CANCER])
        M = np.concatenate([idn[f][f"{arm}|M"] for f in CANCER])
        D = np.concatenate([idn[f][f"{arm}|D"] for f in CANCER])
        return float(np.mean(1 + (U - 2 * M) / D))

    P("# 308 -- |kd|-weighted donors: every arm and control vs the incumbent `cd4svd5`, 3 PLOLO cancer folds")
    P(f"# paired target bootstrap, B = {B}, seed {SEED}; `nmae` = 188's exact identity concatenated over the folds")
    P("")
    base_nmae = nmae("cd4svd5")
    P(f"incumbent: pds " + " / ".join(f"{rows[(f, 'cd4svd5')]['pds']:.6f}" for f in CANCER)
      + f"; nmae identity {base_nmae:.6f}")
    P("")

    # -------------------------------------------------------------- every arm vs the incumbent
    table: dict[str, dict] = {}
    P(f"{'arm':32s} {'pooled dpds [95% CI]':34s} {'per-fold dpds':30s} {'folds>0':>7s} {'dcos%':>7s} {'dnmae':>9s}")
    for a in arms:
        if a == "cd4svd5":
            continue
        diffs = [pt[f][f"{a}|pds"] - pt[f]["cd4svd5|pds"] for f in CANCER]
        b = boot_pool(diffs)
        dc = float(np.mean([rows[(f, a)]["cos_agg"] / rows[(f, 'cd4svd5')]["cos_agg"] - 1 for f in CANCER])) * 100
        dn = nmae(a) - base_nmae
        nf = sum(1 for x in b["per_fold"] if x > 0)
        table[a] = dict(boot=b, dcos_pct=dc, dnmae=dn, folds_pos=nf,
                        pds={f: rows[(f, a)]["pds"] for f in CANCER},
                        n_weighted={f: rows[(f, a)]["n_weighted"] for f in CANCER},
                        audit={f: rows[(f, a)]["audit"] for f in CANCER})
        P(f"{a:32s} {fb(b):34s} " + " ".join(f"{x:+.4f}" for x in b["per_fold"]) + f"   {nf}/3  {dc:+6.2f}  {dn:+9.5f}")
    OUT["arms"] = {a: dict(pooled=v["boot"], dcos_pct=v["dcos_pct"], dnmae=v["dnmae"], folds_pos=v["folds_pos"],
                           pds=v["pds"], n_weighted=v["n_weighted"]) for a, v in table.items()}

    # -------------------------------------------------------------- arm minus control contrasts
    P("")
    P("## arm - control, paired over targets (pooled over the three folds)")
    CONTRASTS = []
    for cname in ("abs", "rel"):
        for a in [f"gate:{cname}", f"gate:{cname}:dis"] + [f"soft:{cname}:p{p}" for p in (0.5, 1.0, 2.0, 4.0)] \
                 + [f"soft:{cname}:p{p}:dis" for p in (1.0, 2.0)]:
            if a not in table:
                continue
            ctrls = ["relw", "relw_glob", "normrel1.5", f"gate:{cname}:rev" if a.startswith("gate") else None]
            if a.startswith("soft"):
                p = a.split(":p")[1].split(":")[0]
                ctrls.append(f"soft:{cname}:p{p}:rev" if f"soft:{cname}:p{p}:rev" in table else None)
                ctrls += [f"soft:{cname}:p1.0:shuf:s{s}" for s in (0, 1, 2)]
            else:
                ctrls += [f"gate:{cname}:shuf:s{s}" for s in (0, 1, 2)]
                ctrls += [f"gate:{cname}:rand:s{s}" for s in (0, 1, 2) if f"gate:{cname}:rand:s{s}" in table]
            if a.startswith("soft") and f"soft:{cname}:p1.0:rand:s0" in table:
                ctrls += [f"soft:{cname}:p1.0:rand:s{s}" for s in (0, 1, 2)]
            for c in [x for x in ctrls if x and x in table]:
                CONTRASTS.append((a, c))
    seen = set()
    for a, c in CONTRASTS:
        if (a, c) in seen:
            continue
        seen.add((a, c))
        diffs = [pt[f][f"{a}|pds"] - pt[f][f"{c}|pds"] for f in CANCER]
        b = boot_pool(diffs)
        OUT.setdefault("contrasts", {})[f"{a} - {c}"] = b
        P(f"  {a:30s} - {c:26s} {fb(b):34s} folds>0 {sum(1 for x in b['per_fold'] if x > 0)}/3")

    # -------------------------------------------------------------- the mechanism: is |kd| just row norm?
    P("")
    P("## the covariate structure of |kd| (per fold, over the two-donor targets)")
    cov = {}
    for f in CANCER:
        z = pt[f]
        A = np.asarray(z["abs_kd_filled"], float)
        avail = np.asarray(z["avail"], bool)
        R = np.asarray(z["R_raw"], float)
        cells = np.asarray(z["cells"], float)
        don = meta[f]["donors"]
        for j, dn_ in enumerate(don):
            m = avail[j] & np.isfinite(A[j])
            rr = float(spearmanr(A[j][m & np.isfinite(R[j])], R[j][m & np.isfinite(R[j])]).statistic)
            mc = m & np.isfinite(cells[j])
            rc = float(spearmanr(A[j][mc], cells[j][mc]).statistic)
            cov[f"{f}|{dn_}"] = dict(rho_kd_reliability=rr, rho_kd_cells=rc, n=int(m.sum()))
            P(f"  [{f:14s}] {dn_:14s} n {int(m.sum()):3d}  Spearman(|kd|, SB reliability) {rr:+.3f}   "
              f"Spearman(|kd|, cells/target) {rc:+.3f}")
    # is |kd| just the row NORM in disguise?  The `normrel1.5` weights encode the norm ratio exactly:
    # for J = 2, w_0/w_1 = (K_0/K_1)^-1.5, so log(K_0/K_1) = -(2/3) log(w_0/w_1).
    for f in CANCER:
        W = np.load(WORK / f"weights308-{f}.npz")["normrel1.5"].astype(float)
        A = np.asarray(pt[f]["abs_kd_filled"], float)
        avail = np.asarray(pt[f]["avail"], bool)
        ok = np.isfinite(W).all(1) & avail.all(0) & np.isfinite(A).all(0) & (A > 0).all(0)
        lognorm = -(2.0 / 3.0) * (np.log(W[ok, 0]) - np.log(W[ok, 1]))
        logkd = np.log(A[0, ok]) - np.log(A[1, ok])
        r = float(spearmanr(logkd, lognorm).statistic)
        cov[f"{f}|kd_vs_norm"] = dict(rho_log_ratio=r, n=int(ok.sum()))
        P(f"  [{f:14s}] Spearman(log |kd| ratio, log row-norm ratio) over the two donors = {r:+.3f} (n {int(ok.sum())})")
    OUT["covariates"] = cov

    # -------------------------------------------------------------- does |kd| pick the better donor? (the oracle test)
    P("")
    P("## does the higher-|kd| donor carry the row closer to the held truth?  (the arbiter question, 263's frame)")
    P("   -- read on the ORACLE arrays 263 measured; here approximated by the per-target `pds` of the two")
    P("      always-drop-one arms is not available, so this is the direct agreement of the |kd| pick with the")
    P("      per-target advantage of the gate arm over the incumbent.")
    for cname in ("abs", "rel"):
        for f in CANCER:
            z = pt[f]
            g = np.asarray(z[f"gate|{cname}"], bool)
            d = pt[f][f"gate:{cname}|pds"] - pt[f]["cd4svd5|pds"]
            if g.sum() == 0:
                continue
            up = int((d[g] > 0).sum())
            dn_ = int((d[g] < 0).sum())
            P(f"  [{cname}] {f:14s} disagreeing {int(g.sum()):3d}: gate moves {up} targets UP, {dn_} DOWN, "
              f"{int(g.sum()) - up - dn_} unchanged; mean d pds on the set {float(np.mean(d[g])):+.5f}")
            OUT.setdefault("disagree_moves", {})[f"{cname}|{f}"] = dict(n=int(g.sum()), up=up, down=dn_,
                                                                        mean_dpds=float(np.mean(d[g])))

    # -------------------------------------------------------------- depth conditioning (programme rule 9)
    P("")
    P("## depth conditioning (rule 9): the arm re-read inside terciles of the donors' cells/target")
    for cname in ("abs", "rel"):
        a = f"gate:{cname}"
        rows_t = []
        for f in CANCER:
            z = pt[f]
            cells = np.asarray(z["cells"], float)
            with np.errstate(all="ignore"):
                cw = np.where(np.asarray(z["avail"], bool), cells, np.nan)
                mc = np.where(np.isfinite(cw).any(0), np.nanmin(np.where(np.isfinite(cw), cw, np.inf), axis=0), np.nan)
            d = pt[f][f"{a}|pds"] - pt[f]["cd4svd5|pds"]
            g = np.asarray(z[f"gate|{cname}"], bool)
            ok = g & np.isfinite(mc)
            if ok.sum() < 6:
                continue
            q = np.quantile(mc[ok], [1 / 3, 2 / 3])
            for lab, m in (("T1 shallow", ok & (mc <= q[0])), ("T2", ok & (mc > q[0]) & (mc <= q[1])),
                           ("T3 deep", ok & (mc > q[1]))):
                rows_t.append((cname, f, lab, int(m.sum()), float(np.mean(d[m])) if m.sum() else float("nan")))
        for r in rows_t:
            P(f"  [{r[0]}] {r[1]:14s} {r[2]:10s} n {r[3]:3d}  mean d pds {r[4]:+.5f}")
        OUT.setdefault("depth", {})[cname] = [dict(fold=r[1], tercile=r[2], n=r[3], mean_dpds=r[4]) for r in rows_t]

    # -------------------------------------------------------------- the h1-test instrument (273's reference)
    P("")
    P("## h1-test (273's Flex-400 reference, union OFF, 85 targets); KA against 273's own published arms")
    fx = json.loads((WORK / "flex308.json").read_text())
    ref = json.loads((ROOT / "experiments/273-flex-reference-instrument/work/instrument273.json").read_text())
    fa, ra = fx["arms"], ref.get("arms", {})
    ka_ok = True
    for nm in ("mean:rebuilt", "relw", "relw:rev", "ebrow"):
        d = max(abs(float(fa[nm][k]) - float(ra[nm][k])) for k in ("pds", "cos", "a23_cos_agg"))
        ka_ok &= d == 0.0
        P(f"  KA {'PASS' if d == 0.0 else 'FAIL'} h1-test {nm:14s} == 273's own row on pds / cos / cos_agg: |d| {d:.1e}")
        OUT.setdefault("flex_ka", {})[nm] = dict(max_abs_diff=d, ok=bool(d == 0.0))
    fbase = fa["mean:rebuilt"]["pds"]
    zf = np.load(WORK / "flex308_pertarget.npz", allow_pickle=True)
    P(f"  incumbent (mean:rebuilt) pds {fbase:.6f}, cos_agg {fa['mean:rebuilt']['a23_cos_agg']:.6f}, 85 targets")
    P(f"  {'arm':28s} {'d pds':>10s} {'paired 95% CI':>26s} {'d cos_agg%':>11s}")
    for nm in sorted(fa):
        if nm == "mean:rebuilt":
            continue
        d = np.asarray(zf[f"{nm}|pds"], float) - np.asarray(zf["mean:rebuilt|pds"], float)
        b = boot(d)
        OUT.setdefault("flex", {})[nm] = dict(d_pds=fa[nm]["pds"] - fbase, boot=b,
                                              dcos_pct=100 * (fa[nm]["a23_cos_agg"] / fa["mean:rebuilt"]["a23_cos_agg"] - 1))
        P(f"  {nm:28s} {fa[nm]['pds']-fbase:+10.5f} [{b['lo']:+.5f}, {b['hi']:+.5f}]"
          f" {100*(fa[nm]['a23_cos_agg']/fa['mean:rebuilt']['a23_cos_agg']-1):+10.2f}")
    OUT["flex_meta"] = dict(n_disagree=fx["n_disagree"], n_targets=fx["n_targets"], ka_all_pass=bool(ka_ok))

    # -------------------------------------------------------------- the CSV
    feat = np.load(WORK / "kd308_feature.npz", allow_pickle=True)
    panel = [str(t) for t in feat["panel"]]
    donors = [str(d) for d in feat["donors"]]
    kd = np.asarray(feat["kd"], float)
    measured = np.asarray(feat["measured"], bool)
    cells_all = np.asarray(feat["cells"], float)
    cut_rel = np.asarray(feat["cut_rel"], float)
    cut_abs = float(feat["cut_abs"])
    per_fold_pds: dict[tuple[str, str], dict[str, float]] = {}
    for f in CANCER:
        tz = [str(t) for t in pt[f]["targets"]]
        for a in ("cd4svd5", "gate:abs", "gate:rel", "soft:rel:p1.0", "soft:abs:p1.0", "relw"):
            key = f"{a}|pds"
            if key not in pt[f]:
                continue
            for i, t in enumerate(tz):
                per_fold_pds.setdefault((f, t), {})[a] = float(pt[f][key][i])
    with (WORK / "kd308.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target", "donor", "kd_lfc", "abs_kd", "measured", "cells", "cut_abs", "cut_rel",
                    "pass_abs", "pass_rel", "w_gate_abs", "w_gate_rel", "w_soft_rel_p1",
                    "fold", "disagree_abs", "disagree_rel",
                    "pds_cd4svd5", "pds_gate_abs", "pds_gate_rel", "pds_soft_rel_p1", "pds_relw"])
        wcache = {f: np.load(WORK / f"weights308-{f}.npz") for f in CANCER}
        for f in CANCER:
            z = pt[f]
            tz = [str(t) for t in z["targets"]]
            don = meta[f]["donors"]
            avail = np.asarray(z["avail"], bool)
            for i, t in enumerate(tz):
                pi = panel.index(t)
                for j, dn_ in enumerate(don):
                    dj = donors.index(dn_)
                    def wv(key):
                        return (float(wcache[f][key][i, j]) if key in wcache[f].files
                                and np.isfinite(wcache[f][key][i, j]) else "")
                    w.writerow([t, dn_, f"{kd[pi, dj]:.6f}" if measured[pi, dj] else "",
                                f"{abs(kd[pi, dj]):.6f}" if measured[pi, dj] else "", int(measured[pi, dj]),
                                "" if not np.isfinite(cells_all[pi, dj]) else int(cells_all[pi, dj]),
                                f"{cut_abs:.4f}", f"{cut_rel[dj]:.4f}",
                                int(measured[pi, dj] and abs(kd[pi, dj]) >= cut_abs),
                                int(measured[pi, dj] and abs(kd[pi, dj]) >= cut_rel[dj]),
                                wv("gate:abs"), wv("gate:rel"), wv("soft:rel:p1.0"), f,
                                int(np.asarray(z["gate|abs"], bool)[i]), int(np.asarray(z["gate|rel"], bool)[i]),
                                *[f"{per_fold_pds[(f, t)].get(a, float('nan')):.6f}" for a in
                                  ("cd4svd5", "gate:abs", "gate:rel", "soft:rel:p1.0", "relw")]])
    P("")
    P("wrote work/kd308.csv, work/analyse308.json, work/analyse308.log")
    (WORK / "analyse308.json").write_text(json.dumps(OUT, indent=2))
    (WORK / "analyse308.log").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
