"""290 step 2 — KA-2 / KA-3 and the preregistered fits (PREREG §3, §6, §7, §8).

    PYTHONDONTWRITEBYTECODE=1 uv run python fit290.py

Writes work/fit290.json, work/fit290.log, work/pred290.csv.
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import pandas as pd

import lib290 as L

T0 = time.time()
LOG: list[str] = []
N_PART = 20
N_SHUF = 200


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


# ---- rem284.py:109-131 verbatim (the admissible-partition draw) --------------------------------
def draw_partitions(lines: np.ndarray, pairs: np.ndarray):
    line_groups = sorted(set(lines.tolist()))
    src_of_line = {ln: np.flatnonzero(lines == ln) for ln in line_groups}

    def make_partition(rng):
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
    parts, tries = [], 0
    while len(parts) < N_PART and tries < 10000:
        tries += 1
        r = make_partition(rng)
        if r is not None:
            parts.append(r)
    return parts, tries


def main() -> None:
    h = L.check_prereg()
    P(f"prereg {h} MATCH")

    held = L.load_held()
    X21 = held[L.FEATS].to_numpy(float)
    rbar = held["r_bar"].to_numpy(float)
    gene = held["gene"].to_numpy().astype(str)
    comp = np.load(L.SPLITS284, allow_pickle=True)["comp"]
    fold = L.component_folds(comp, 5)
    blk281 = L.block281(held)

    z = np.load(L.HELD279, allow_pickle=True)
    Ccross = np.asarray(z["Ccross"])
    lines = np.asarray(z["lines"]).astype(str)
    pairs = np.asarray(z["pairs_cross"])

    feat = pd.read_csv(L.WORK / "feat290.csv")
    assert list(feat["gene"].astype(str)) == list(gene), "feat290 gene order differs from feat281"

    ref = np.load(L.REM284, allow_pickle=True)
    ref284 = json.loads(L.REM284J.read_text())

    def maps(p1, p2):
        """rem284.analyse's construction, verbatim."""
        C1 = L.cmap(Ccross, p1).astype(float)
        C2 = L.cmap(Ccross, p2).astype(float)
        n1 = np.isfinite(Ccross[p1]).sum(axis=0)
        n2 = np.isfinite(Ccross[p2]).sum(axis=0)
        both = (n1 >= 1) & (n2 >= 1)
        pf1, _ = L.oof_ridge(X21, C1, fold)
        pf2, _ = L.oof_ridge(X21, C2, fold)
        r1, r2 = C1 - pf1, C2 - pf2
        ok = both & np.isfinite(r1) & np.isfinite(r2)
        return dict(C1=C1, C2=C2, r1=r1, r2=r2, pf1=pf1, pf2=pf2, n1=n1, n2=n2, ok=ok)

    parts, tries = draw_partitions(lines, pairs)
    P(f"drew {len(parts)} admissible partitions in {tries} tries "
      f"(rem284: 20 in 22 tries -> {'MATCH' if tries == 22 else 'DIFFER'})")

    # ================= KA-2: the primary partition, rebuilt =====================================
    p1, p2, h1l, h2l = parts[0]
    assert h1l == L.HALF1_LINES and h2l == L.HALF2_LINES, "partition #1 is not the declared one"
    M = maps(p1, p2)
    d = {k: float(np.nanmax(np.abs(M[k] - np.asarray(ref[k])))) for k in ("C1", "C2", "r1", "r2")}
    same_ok = bool(np.array_equal(M["ok"], np.asarray(ref["ok"])))
    ok = M["ok"]
    ceil_r = L.spearman(M["C1"][ok], M["C2"][ok])
    rem_r = L.spearman(M["r1"][ok], M["r2"][ok])
    ka2 = dict(max_abs_diff=d, same_ok=same_ok, n=int(ok.sum()), ceiling=ceil_r, remainder=rem_r,
               pass_=bool(all(v == 0.0 for v in d.values()) and same_ok and int(ok.sum()) == 2848
                          and abs(ceil_r - 0.3816) < 5e-5 and abs(rem_r - 0.3277) < 5e-5))
    P(f"KA-2  primary partition rebuilt: max|diff| C1 {d['C1']:.1e} C2 {d['C2']:.1e} "
      f"r1 {d['r1']:.1e} r2 {d['r2']:.1e}; ok identical {same_ok}; n {int(ok.sum())} (want 2848); "
      f"ceiling {ceil_r:+.4f} (want +0.3816); remainder {rem_r:+.4f} (want +0.3277) "
      f"-> {'PASS' if ka2['pass_'] else 'FAIL'}")
    if not ka2["pass_"]:
        raise SystemExit("KA-2 FAILED — abort (AGENTS.md:1375-1382)")

    # ================= KA-3: replication partitions #2 and #3 ===================================
    rep = {}
    ka3ok = True
    for i in (1, 2):
        a, b, ha, hb = parts[i]
        Mi = maps(a, b)
        oki = Mi["ok"]
        c_, r_ = L.spearman(Mi["C1"][oki], Mi["C2"][oki]), L.spearman(Mi["r1"][oki], Mi["r2"][oki])
        w = ref284["partitions"][i]
        good = abs(c_ - w["ceiling"]) < 5e-5 and abs(r_ - w["rho_remainder"]) < 5e-5
        ka3ok &= bool(good)
        rep[f"partition#{i + 1}"] = dict(half1_lines=ha, half2_lines=hb, n=int(oki.sum()),
                                         ceiling=c_, remainder=r_, want_ceiling=w["ceiling"],
                                         want_remainder=w["rho_remainder"], reproduces=bool(good),
                                         _M=Mi)
        P(f"KA-3  partition#{i + 1}: n {int(oki.sum())} ceiling {c_:+.4f} (want "
          f"{w['ceiling']:+.4f}) remainder {r_:+.4f} (want {w['rho_remainder']:+.4f}) "
          f"-> {'PASS' if good else 'FAIL'}")
    if not ka3ok:
        raise SystemExit("KA-3 FAILED — abort")

    # ================= the unit and the feature matrices ========================================
    u = np.flatnonzero(ok)
    compu, foldu, blku = comp[u], fold[u], blk281[u]
    P(f"\nunit: {len(u)} genes; fold sizes {np.bincount(foldu, minlength=5).tolist()}; "
      f"{len(set(compu.tolist()))} components; {len(set(blku.tolist()))} first-family blocks")

    def Xof(tag: str) -> np.ndarray:
        cols = [c if c == "on_axis" else f"{c}_{tag}" for c in L.FEATS290]
        return feat.loc[u, cols].to_numpy(float)

    Xh1, Xh2 = Xof("h1"), Xof("h2")
    cov1 = np.isfinite(Xh1).all(axis=1).mean()
    cov2 = np.isfinite(Xh2).all(axis=1).mean()
    P(f"complete feature rows on the unit: half1 {cov1:.3f}, half2 {cov2:.3f} "
      f"(on-target is NaN off the axis by design)")
    covtab = {c: dict(h1=float(np.isfinite(feat.loc[u, c if c == 'on_axis' else c + '_h1']).mean()),
                      h2=float(np.isfinite(feat.loc[u, c if c == 'on_axis' else c + '_h2']).mean()))
              for c in L.FEATS290}

    out = {"prereg_sha256": h, "ka2": ka2, "ka3": {k: {kk: vv for kk, vv in v.items()
                                                       if kk != "_M"} for k, v in rep.items()},
           "unit": dict(n=len(u), fold_sizes=np.bincount(foldu, minlength=5).tolist(),
                        n_components=len(set(compu.tolist())), n_blocks281=len(set(blku.tolist())),
                        feature_coverage=covtab)}

    # ================= the preregistered fits ===================================================
    def run(Xf: np.ndarray, y: np.ndarray, tag: str, with_gbm: bool = True) -> dict:
        pr, ir = L.oof_ridge(Xf, y, foldu)
        res = dict(ridge=pr, ridge_info={str(k): v for k, v in ir.items()})
        if with_gbm:
            pg, ig = L.oof_gbm(Xf, y, foldu)
            res["gbm"] = pg
            res["gbm_info"] = {str(k): v for k, v in ig.items()}
        P(f"  {tag:34s} ridge rho {L.spearman(y, pr):+.4f}"
          + (f"   gbm rho {L.spearman(y, res['gbm']):+.4f}" if with_gbm else ""))
        return res

    y_r2 = M["r2"][u]        # PRIMARY target
    y_C2 = M["C2"][u]        # secondary target
    y_r1 = M["r1"][u]        # symmetric target
    y_C1 = M["C1"][u]

    P("\nout-of-fold fits (grouped 5-fold over components):")
    fit_pri = run(Xh1, y_r2, "half1 features -> remainder_half2")
    fit_ann = run(X21[u], y_r2, "annotation 21 -> remainder_half2")
    fit_sec = run(Xh1, y_C2, "half1 features -> C_half2 (raw)")
    fit_ann_sec = run(X21[u], y_C2, "annotation 21 -> C_half2 (raw)")
    fit_sym = run(Xh2, y_r1, "half2 features -> remainder_half1")
    fit_sym_ann = run(X21[u], y_r1, "annotation 21 -> remainder_half1")
    fit_sym_sec = run(Xh2, y_C1, "half2 features -> C_half1 (raw)")
    fit_leak = run(Xh2, y_r2, "CONTROL C2 half2 features -> remainder_half2 (LEAK)", False)

    # ---- baselines (single columns) ------------------------------------------------------------
    base = {
        "annotation_ridge": fit_ann["ridge"],
        "annotation_pred_C2": M["pf2"][u],
        "r_bar": rbar[u],
        "logmag_h1": feat.loc[u, "logmag_h1"].to_numpy(float),
        "C_half1": y_C1,
        "remainder_half1": y_r1,
        "coverage_h1": feat.loc[u, "log_n_src_h1"].to_numpy(float),
    }

    preds = dict(response_ridge=fit_pri["ridge"], response_gbm=fit_pri["gbm"], **base)
    rng = np.random.default_rng(L.SEED + 11)
    bb_comp = L.block_boot(y_r2, preds, compu, rng)
    bb_fam = L.block_boot(y_r2, preds, blku, np.random.default_rng(L.SEED + 12))

    P("\nPRIMARY TARGET remainder_half2 — out-of-fold Spearman, component / first-family bootstrap:")
    for k in preds:
        c, f = bb_comp[k]["rho"], bb_fam[k]["rho"]
        P(f"  {k:22s} {c['point']:+.4f} [{c['lo']:+.4f}, {c['hi']:+.4f}]   "
          f"fam [{f['lo']:+.4f}, {f['hi']:+.4f}]   R2 {bb_comp[k]['r2']['point']:+.4f}")
    P(f"  {'mean predictor':22s}  rho undefined (constant); R2 0 by construction")

    # ---- per-feature Spearman -------------------------------------------------------------------
    P("\nper-feature Spearman with the primary target (component bootstrap):")
    perfeat = {}
    for c in L.FEATS290 + L.EXTRA_COLS:
        col = feat.loc[u, c if c == "on_axis" else f"{c}_h1"].to_numpy(float)
        perfeat[c] = L.boot_rho(col, y_r2, compu, np.random.default_rng(L.SEED + 21))
        v = perfeat[c]
        P(f"  {c:34s} {v['point']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]  n {v['n']}")

    # ---- CONTROL C1: shuffled target --------------------------------------------------------------
    P(f"\nCONTROL C1 — {N_SHUF} shuffles of the target (ridge only, 281 KA-4' form):")
    rs = np.random.default_rng(L.SEED + 31)
    nulls = []
    for _ in range(N_SHUF):
        ysh = y_r2[rs.permutation(len(y_r2))]
        p, _ = L.oof_ridge(Xh1, ysh, foldu)
        nulls.append(L.spearman(ysh, p))
    nulls = np.asarray(nulls)
    obs = L.spearman(y_r2, fit_pri["ridge"])
    c1 = dict(mean=float(np.nanmean(nulls)), sd=float(np.nanstd(nulls)),
              n_reaching=int((nulls >= obs).sum()), n_draws=N_SHUF,
              z=float((obs - np.nanmean(nulls)) / np.nanstd(nulls)))
    P(f"  null {c1['mean']:+.4f} +- {c1['sd']:.4f}; {c1['n_reaching']}/{N_SHUF} reach "
      f"{obs:+.4f}; z {c1['z']:+.1f}")

    # ---- POST-HOC: coverage added, and the strongest-baseline combinations ------------------------
    Xcov = np.column_stack([Xh1, feat.loc[u, "log_n_src_h1"].to_numpy(float)])
    fit_cov = run(Xcov, y_r2, "POST-HOC half1 features + coverage", False)
    Xc1 = np.column_stack([Xh1, y_C1])
    fit_withC1 = run(Xc1, y_r2, "POST-HOC half1 features + C_half1", False)

    # ================= the frozen pass rules ====================================================
    rho_resp = bb_comp["response_ridge"]["rho"]
    p1_pass = bool(rho_resp["point"] >= 0.10 and (rho_resp["lo"] > 0 or rho_resp["hi"] < 0))
    dkey = "response_ridge-annotation_ridge"
    dpair = bb_comp["_paired_rho"][dkey]
    p2_pass = bool(dpair["lo"] > 0)
    p3 = perfeat["pe_prolif_abs"]
    p3_pass = bool(p3["point"] < 0 and p3["hi"] < 0)
    p4 = perfeat["gen_pc5"]
    p4_pass = bool(p4["point"] < 0 and p4["hi"] < 0)

    rules = dict(
        P1=dict(rule="OOF Spearman >= +0.10 AND component interval excludes 0",
                value=rho_resp, passed=p1_pass),
        P2=dict(rule="paired d rho (response - annotation) interval entirely above 0",
                value=dpair, passed=p2_pass),
        P3=dict(rule="Spearman(pe_prolif_abs_h1, remainder_half2) < 0, interval excludes 0",
                value=p3, passed=p3_pass),
        P4=dict(rule="Spearman(gen_pc5_h1, remainder_half2) < 0, interval excludes 0",
                value=p4, passed=p4_pass))
    P("\nFROZEN PASS RULES")
    for k, v in rules.items():
        val = v["value"]
        pt = val.get("point", val.get("mean"))
        P(f"  {k}  {pt:+.4f} [{val['lo']:+.4f}, {val['hi']:+.4f}] -> "
          f"{'PASS' if v['passed'] else 'FAIL'}")

    # ================= replication on partitions #2 and #3 ======================================
    P("\nreplication on the declared partitions #2 and #3:")
    reps = {}
    for i in (1, 2):
        Mi = rep[f"partition#{i + 1}"]["_M"]
        oki = Mi["ok"]
        ui = np.flatnonzero(oki)
        # which half is which: parts[i][2] are half-1 LINES; features follow the SOURCE half, so
        # the feature tag is h1 when the partition's half 1 equals the primary's half-1 line set,
        # otherwise the features must be rebuilt.  Instead: use both stored feature blocks and
        # report the direction whose feature half is a SUBSET of the partition's half-1 lines.
        ha = set(parts[i][2])
        tag = "h1" if ha.issubset(set(L.HALF1_LINES)) else ("h2" if ha.issubset(set(L.HALF2_LINES))
                                                            else None)
        if tag is None:
            reps[f"partition#{i + 1}"] = dict(
                note="the partition's halves mix the primary partition's source halves; the frozen "
                     "feature blocks are not defined for it (features are built per SOURCE half)",
                n=int(oki.sum()))
            P(f"  partition#{i + 1}: half-1 lines are not a subset of either frozen feature half — "
              f"no frozen feature block applies")
            continue
        cols = [c if c == "on_axis" else f"{c}_{tag}" for c in L.FEATS290]
        Xi = feat.loc[ui, cols].to_numpy(float)
        yi = Mi["r2"][ui]
        pi, _ = L.oof_ridge(Xi, yi, fold[ui])
        pai, _ = L.oof_ridge(X21[ui], yi, fold[ui])
        bbi = L.block_boot(yi, dict(response=pi, annotation=pai), comp[ui],
                           np.random.default_rng(L.SEED + 41 + i))
        reps[f"partition#{i + 1}"] = dict(
            n=int(oki.sum()), feature_half=tag,
            response=bbi["response"]["rho"], annotation=bbi["annotation"]["rho"],
            paired=bbi["_paired_rho"]["response-annotation"],
            ceiling=rep[f"partition#{i + 1}"]["ceiling"],
            remainder=rep[f"partition#{i + 1}"]["remainder"],
            pe_prolif_abs=L.boot_rho(feat.loc[ui, f"pe_prolif_abs_{tag}"].to_numpy(float), yi,
                                     comp[ui], np.random.default_rng(L.SEED + 51 + i)),
            gen_pc5=L.boot_rho(feat.loc[ui, f"gen_pc5_{tag}"].to_numpy(float), yi, comp[ui],
                               np.random.default_rng(L.SEED + 61 + i)))
        r = reps[f"partition#{i + 1}"]
        P(f"  partition#{i + 1} (features from the {tag} source half, n {r['n']}): response "
          f"{r['response']['point']:+.4f} [{r['response']['lo']:+.4f}, {r['response']['hi']:+.4f}]"
          f"  annotation {r['annotation']['point']:+.4f}  paired D "
          f"{r['paired']['mean']:+.4f} [{r['paired']['lo']:+.4f}, {r['paired']['hi']:+.4f}]")
        P(f"      pe_prolif_abs {r['pe_prolif_abs']['point']:+.4f} "
          f"[{r['pe_prolif_abs']['lo']:+.4f}, {r['pe_prolif_abs']['hi']:+.4f}]   gen_pc5 "
          f"{r['gen_pc5']['point']:+.4f} [{r['gen_pc5']['lo']:+.4f}, {r['gen_pc5']['hi']:+.4f}]")

    # ================= symmetric direction and secondary target =================================
    def summarise(y, pdict, tag, groups=None):
        b = L.block_boot(y, pdict, compu if groups is None else groups,
                         np.random.default_rng(L.SEED + 71))
        P(f"\n{tag}:")
        for k in pdict:
            v = b[k]["rho"]
            P(f"  {k:22s} {v['point']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]  "
              f"R2 {b[k]['r2']['point']:+.4f}")
        return {k: {m: b[k][m] for m in ("rho", "r2", "auc")} for k in pdict} | {
            "_paired_rho": b["_paired_rho"]}

    sec = summarise(y_C2, dict(response_ridge=fit_sec["ridge"], response_gbm=fit_sec["gbm"],
                               annotation_ridge=fit_ann_sec["ridge"], C_half1=y_C1,
                               r_bar=rbar[u], logmag_h1=base["logmag_h1"]),
                    "SECONDARY TARGET C_half2 (raw)")
    sym = summarise(y_r1, dict(response_ridge=fit_sym["ridge"], response_gbm=fit_sym["gbm"],
                               annotation_ridge=fit_sym_ann["ridge"], C_half2=y_C2,
                               remainder_half2=y_r2, r_bar=rbar[u],
                               logmag_h2=feat.loc[u, "logmag_h2"].to_numpy(float)),
                    "SYMMETRIC half2 features -> remainder_half1")
    sym_sec = summarise(y_C1, dict(response_ridge=fit_sym_sec["ridge"],
                                   response_gbm=fit_sym_sec["gbm"], C_half2=y_C2),
                        "SYMMETRIC secondary C_half1 (raw)")
    leak = summarise(y_r2, dict(response_ridge_SAMEHALF=fit_leak["ridge"],
                                response_ridge_otherhalf=fit_pri["ridge"]),
                     "CONTROL C2 — features on the SAME half as the target (LEAKAGE reading)")
    posthoc = summarise(y_r2, dict(plus_coverage=fit_cov["ridge"], plus_C1=fit_withC1["ridge"],
                                   frozen=fit_pri["ridge"]),
                        "POST-HOC additions")

    # symmetric P3/P4
    symp = {c: L.boot_rho(feat.loc[u, f"{c}_h2"].to_numpy(float), y_r1, compu,
                          np.random.default_rng(L.SEED + 81))
            for c in ("pe_prolif_abs", "gen_pc5", "gen_pc20", "gen_cos", "gen_abscos",
                      "pe_stress_abs", "logmag", "breadth")}
    P("\nsymmetric per-feature (half2 features -> remainder_half1):")
    for c, v in symp.items():
        P(f"  {c:20s} {v['point']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]")

    out.update(dict(primary=dict(component={k: bb_comp[k] for k in preds},
                                 family={k: bb_fam[k] for k in preds},
                                 paired_component=bb_comp["_paired_rho"],
                                 paired_family=bb_fam["_paired_rho"],
                                 ridge_info=fit_pri["ridge_info"], gbm_info=fit_pri["gbm_info"]),
                    per_feature=perfeat, control_shuffle=c1, control_samehalf=leak,
                    secondary=sec, symmetric=sym, symmetric_secondary=sym_sec,
                    symmetric_per_feature=symp, posthoc=posthoc, replication=reps, rules=rules,
                    verdict=dict(P1=p1_pass, P2=p2_pass, P3=p3_pass, P4=p4_pass)))

    pred = pd.DataFrame(dict(
        gene=gene[u], component=compu, fold=foldu, block281=blku,
        n_pairs_half1=M["n1"][u], n_pairs_half2=M["n2"][u],
        C_half1=y_C1, C_half2=y_C2, remainder_half1=y_r1, remainder_half2=y_r2,
        annot_pred_C_half1=M["pf1"][u], annot_pred_C_half2=M["pf2"][u], r_bar=rbar[u],
        oof_response_ridge=fit_pri["ridge"], oof_response_gbm=fit_pri["gbm"],
        oof_annotation_ridge=fit_ann["ridge"],
        oof_response_ridge_C2=fit_sec["ridge"], oof_response_ridge_sym=fit_sym["ridge"],
        oof_response_ridge_samehalf=fit_leak["ridge"]))
    pred.to_csv(L.WORK / "pred290.csv", index=False)

    out["seconds"] = time.time() - T0
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    (L.WORK / "fit290.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "fit290.log").write_text("\n".join(LOG))
    P(f"\npeak RSS {out['peak_rss_gib']:.2f} GiB; wrote work/fit290.json, work/pred290.csv")


if __name__ == "__main__":
    main()
