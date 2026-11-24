"""291 surface B — claim C4: does 290's response geometry (measured on the 23 K562-FREE banks)
predict the GWCD4i authors' own published CD4-vs-K562 conservation (`logfc_pearson_r`)?
KA-3 (287's random-control row) first; nothing else is computed if it fails.

    PYTHONDONTWRITEBYTECODE=1 uv run python surfB291.py
    -> work/pred291_gwcd4i.csv, work/surfB291.json, work/surfB291.log
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import pandas as pd

import lib291 as L

T0 = time.time()
LOG: list[str] = []


def P(m=""):
    line = f"[{time.time() - T0:6.1f}s] {m}" if m else ""
    print(line, flush=True)
    LOG.append(line)


def main():
    h = L.check_prereg()
    P(f"prereg {h} MATCH")
    k5 = pd.read_csv(L.K562TAB, low_memory=False)

    # ---- KA-3 (287's KA-1): the random controls read at chance ---------------------------------
    ka3 = {}
    okall = True
    for cond in L.CONDS:
        s = k5[k5["condition"] == cond]
        r = s["logfc_pearson_r"].dropna()
        rand = pd.concat([s["random_r1"], s["random_r2"], s["random_r3"]]).dropna()
        row = dict(n_real=int(len(r)), real_median=float(r.median()),
                   random_median=float(rand.median()), random_p95=float(rand.quantile(0.95)),
                   pass_=bool(abs(float(rand.median())) < 0.01 and 0.25 < float(r.median()) < 0.36))
        ka3[cond] = row
        okall &= row["pass_"]
        P(f"KA-3 {cond:9s} real median {row['real_median']:+.4f} (n {row['n_real']}); random median "
          f"{row['random_median']:+.4f}, p95 {row['random_p95']:+.4f} "
          f"-> {'PASS' if row['pass_'] else 'FAIL'}")
    if not okall:
        raise SystemExit("KA-3 FAILED — the parse is wrong, abort (AGENTS.md:1375-1382)")

    # ---- models (PREREG §2), features from the 23 K562-free banks -------------------------------
    store = L.Store()
    rem = np.load(L.REM284, allow_pickle=True)
    fold = np.asarray(rem["fold"], int)
    comp = np.asarray(rem["comp"], int)
    z = np.load(L.HELD279, allow_pickle=True)
    c_held = np.asarray(z["c_held"], float)

    rows_free = store.subset(exclude=L.K562_SOURCES)
    rows_all = store.subset()
    P(f"K562 sources excluded: {L.K562_SOURCES} -> {len(rows_free)} of {len(rows_all)} banks")
    agg_free = store.aggregate(rows_free)
    agg_all = store.aggregate(rows_all)

    ann = L.load_annotation()
    assert list(ann["gene"]) == store.genes
    Xa = ann[L.FEATS].to_numpy(float)
    Xf = agg_free[L.FEATS290].to_numpy(float)
    Xall = agg_all[L.FEATS290].to_numpy(float)

    # a K562-free training TARGET, labelled sensitivity
    Cc = np.asarray(z["Ccross"], float)
    pairs = np.asarray(z["pairs_cross"])
    keys = [str(k) for k in z["keys"]]
    kfree = [i for i, (a, b) in enumerate(pairs)
             if keys[a] not in L.K562_SOURCES and keys[b] not in L.K562_SOURCES]
    with np.errstate(invalid="ignore"):
        c_held_free = np.nanmean(Cc[kfree], axis=0)
    P(f"K562-free training target: {len(kfree)} of {len(pairs)} cross-line pairs, "
      f"n {int(np.isfinite(c_held_free).sum())}, median {np.nanmedian(c_held_free):+.6f}")
    del Cc

    pred_f, info_f = L.fit_apply_grouped(Xf, c_held, fold, Xf)
    pred_all, info_a28 = L.fit_apply_grouped(Xall, c_held, fold, Xall)
    pred_ann, info_ann = L.fit_apply_grouped(Xa, c_held, fold, Xa)
    pred_ff, info_ff = L.fit_apply_grouped(Xf, c_held_free, fold, Xf)
    oof_f, _ = L.oof_ridge(Xf, c_held, fold)
    P(f"M_resp(23 banks) alpha {info_f['alpha']:g} cv {info_f['grouped_cv_spearman']:+.4f}; "
      f"M_resp(28) alpha {info_a28['alpha']:g} cv {info_a28['grouped_cv_spearman']:+.4f}; "
      f"M_ann alpha {info_ann['alpha']:g} cv {info_ann['grouped_cv_spearman']:+.4f}; "
      f"M_resp(23, K562-free target) alpha {info_ff['alpha']:g} cv {info_ff['grouped_cv_spearman']:+.4f}")

    gpos = {g: i for i, g in enumerate(store.genes)}
    held = set(store.genes)
    znf = pd.read_csv(L.ZNF279, low_memory=False)
    cls = dict(zip(znf["gene"], znf["class"]))

    results, tables = {}, []
    for cond in L.CONDS:
        s = k5[k5["condition"] == cond].dropna(subset=["logfc_pearson_r"])
        s = s[s["target_contrast_gene_name"].isin(held)].copy()
        g = s["target_contrast_gene_name"].to_numpy()
        assert len(set(g)) == len(g)
        gi = np.array([gpos[x] for x in g])
        y = s["logfc_pearson_r"].to_numpy(float)
        nd_k = s["n_degs_MASH_K562"].to_numpy(float)
        nd_c = s[f"n_degs_MASH_{cond}"].to_numpy(float)
        dcm = s["donor_correlation_mean"].to_numpy(float)
        blk = comp[gi]
        preds = {
            "M_resp": pred_f[gi],
            "M_ann": pred_ann[gi],
            "C_held": c_held[gi],
            "n_degs_K562": nd_k,
            "n_degs_cond": nd_c,
            "gen_pc20": agg_free["gen_pc20"].to_numpy(float)[gi],
            "breadth": agg_free["breadth"].to_numpy(float)[gi],
            "logmag": agg_free["logmag"].to_numpy(float)[gi],
            "M_resp_28banks": pred_all[gi],
            "M_resp_K562free_target": pred_ff[gi],
            "M_resp_oof": oof_f[gi],
        }
        bb = L.block_boot(y, preds, blk, np.random.default_rng(L.SEED))
        P(f"\n=== {cond} (n {len(y)}, {len(set(blk.tolist()))} component blocks) ===")
        for k in preds:
            P(f"  {k:24s} rho {bb[k]['point']:+.4f} [{bb[k]['lo']:+.4f}, {bb[k]['hi']:+.4f}]")
        pr = bb["_paired_rho"]
        for k in ("M_ann", "C_held", "n_degs_K562", "gen_pc20"):
            d = pr[f"M_resp-{k}"]
            P(f"  paired M_resp - {k:16s} {d['mean']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}] "
              f"P(<=0) {d['p_le_0']:.3f}")

        # donor_correlation_mean, on its own subset
        md = np.isfinite(dcm)
        bd = L.block_boot(y[md], {"donor_corr": dcm[md], "M_resp": preds["M_resp"][md]},
                          blk[md], np.random.default_rng(L.SEED))
        P(f"  donor_correlation_mean (n {int(md.sum())}) rho {bd['donor_corr']['point']:+.4f} "
          f"[{bd['donor_corr']['lo']:+.4f}, {bd['donor_corr']['hi']:+.4f}]; M_resp on the same subset "
          f"{bd['M_resp']['point']:+.4f} [{bd['M_resp']['lo']:+.4f}, {bd['M_resp']['hi']:+.4f}]")

        # partials given the authors' own magnitude
        Z = np.column_stack([nd_k, nd_c])
        part = {k: L.boot_partial(preds[k], y, Z, blk, np.random.default_rng(L.SEED))
                for k in ("M_resp", "M_ann", "gen_pc20", "breadth", "C_held")}
        for k, v in part.items():
            P(f"  partial | (n_degs_K562, n_degs_{cond}): {k:10s} {v['point']:+.4f} "
              f"[{v['lo']:+.4f}, {v['hi']:+.4f}]")

        # null: 200 gene shuffles
        rng = np.random.default_rng(L.SEED)
        nv = np.array([L.spearman(rng.permutation(y), preds["M_resp"]) for _ in range(200)])
        obs = bb["M_resp"]["point"]
        nullrow = dict(mean=float(nv.mean()), sd=float(nv.std()),
                       n_reaching=int((nv >= obs).sum()),
                       z=float((obs - nv.mean()) / max(nv.std(), 1e-12)))
        P(f"  null (200 gene shuffles) {nullrow['mean']:+.4f} +- {nullrow['sd']:.4f}; "
          f"{nullrow['n_reaching']}/200 reach {obs:+.4f}; z {nullrow['z']:+.1f}")

        rule = dict(
            cond1=bool(bb["M_resp"]["point"] >= 0.15 and bb["M_resp"]["lo"] > 0),
            cond2=bool(pr["M_resp-M_ann"]["lo"] > 0),
            cond3=bool(part["M_resp"]["point"] > 0 and part["M_resp"]["lo"] > 0),
            cond4=bool(bb["gen_pc20"]["point"] > 0 and bb["gen_pc20"]["lo"] > 0
                       and bb["breadth"]["point"] > 0 and bb["breadth"]["lo"] > 0))
        rule["PASS"] = bool(all(rule[c] for c in ("cond1", "cond2", "cond3", "cond4")))
        P(f"  frozen rule in {cond}: {rule}")

        results[cond] = dict(n=len(y), n_blocks=len(set(blk.tolist())), boot=bb,
                             donor=bd, n_donor=int(md.sum()), partial=part, null=nullrow,
                             rule=rule)
        tables.append(pd.DataFrame(dict(
            gene=g, condition=cond, logfc_pearson_r=y,
            n_degs_MASH_K562=nd_k, n_degs_MASH_cond=nd_c, donor_correlation_mean=dcm,
            pred_M_resp=preds["M_resp"], pred_M_ann=preds["M_ann"],
            pred_M_resp_oof=preds["M_resp_oof"], pred_M_resp_28banks=preds["M_resp_28banks"],
            C_held=preds["C_held"], gen_pc20=preds["gen_pc20"], breadth=preds["breadth"],
            logmag=preds["logmag"], component=blk,
            znf_class=[cls.get(x, "") for x in g])))

    pd.concat(tables).to_csv(L.WORK / "pred291_gwcd4i.csv", index=False)
    prim = results[L.PRIMARY_COND]["rule"]
    P(f"\nC4 ({L.PRIMARY_COND}) {'PASS' if prim['PASS'] else 'FAIL'} by the frozen rule: {prim}")

    res = dict(prereg_sha256=h, ka3=ka3, k562_sources=L.K562_SOURCES,
               n_banks_free=len(rows_free),
               model_info=dict(M_resp23=info_f, M_resp28=info_a28, M_ann=info_ann,
                               M_resp23_k562free_target=info_ff,
                               n_k562free_pairs=len(kfree)),
               conditions=results, primary_condition=L.PRIMARY_COND,
               C4_PASS=bool(prim["PASS"]),
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30)
    (L.WORK / "surfB291.json").write_text(json.dumps(res, indent=1, default=float))
    (L.WORK / "surfB291.log").write_text("\n".join(LOG))
    P(f"peak RSS {res['peak_rss_gib']:.2f} GiB; wrote work/pred291_gwcd4i.csv")


if __name__ == "__main__":
    main()
