"""291 surface A — claim C3: does 290's response geometry (measured on the 28 banks) predict the
cross-line conservation of the SAME genes' Cas9-knockout responses over the 12 Broad-16 working
lines?  KA-4 first (268's per-target decomposition recomputed from work/halves/).

    PYTHONDONTWRITEBYTECODE=1 uv run python surfA291.py
    -> work/pred291_broad16.csv, work/surfA291.json, work/surfA291.log
"""
from __future__ import annotations

import json
import resource
import sys
import time

import numpy as np
import pandas as pd

import lib291 as L

sys.path.insert(0, str(L.W268))
import load268 as B                                             # noqa: E402  (asserts the reserve)

T0 = time.time()
LOG: list[str] = []


def P(m=""):
    line = f"[{time.time() - T0:6.1f}s] {m}" if m else ""
    print(line, flush=True)
    LOG.append(line)


def decompose_shared(RA, RB):
    """decomp268.decompose, the per-target energy part only (verbatim arithmetic)."""
    mA, mB = RA.mean(1), RB.mean(1)
    rA, rB = RA - mA[:, None], RB - mB[:, None]
    nL = RA.shape[1]
    E = 0.5 * ((RA ** 2).sum((1, 2)) + (RB ** 2).sum((1, 2)))
    X_sh = nL * (mA * mB).sum(1)
    X_ln = (rA * rB).sum((1, 2))
    tot_full = X_sh + X_ln + 0.5 * (E - X_sh - X_ln)
    return dict(shared_full=X_sh / tot_full, line_full=X_ln / tot_full,
                noise_full=0.5 * (E - X_sh - X_ln) / tot_full,
                E_tot_half=E, X_shared=X_sh, X_line=X_ln)


def main():
    h = L.check_prereg()
    P(f"prereg {h} MATCH")
    rng = np.random.default_rng(L.SEED)

    # ---- 268's responses, working lines only ---------------------------------------------------
    H, genes = B.load_halves(B.WORKING)
    lines = B.WORKING
    P(f"lines {len(lines)}: {lines}")
    tg = B.ALL_TARGETS
    ctrl = list(B.CONTROLS)
    R, ok, avail = B.build_responses(H, genes, lines, tg, keys=("full", "A", "B"), min_cpm=5.0)
    G = int(ok.sum())
    assert G == 12199, f"gene filter drifted: {G}"
    P(f"genes kept {G} (min mean control CPM 5.0), targets {len(tg)}, lines {len(lines)}")

    # ---- KA-4: 268's per-target decomposition --------------------------------------------------
    usable = B.usable_targets(H, lines, 6)
    tpos = {t: i for i, t in enumerate(tg)}
    ui = np.array([tpos[t] for t in usable])
    dd = decompose_shared(R["A"][ui].astype(np.float64), R["B"][ui].astype(np.float64))
    ref = pd.read_csv(L.W268 / "per_target_decomp_cpm5.csv").set_index("target")
    assert sorted(ref.index) == sorted(usable), "usable target set drifted"
    ka4 = {}
    for c in ("shared_full", "line_full", "noise_full", "E_tot_half", "X_shared", "X_line"):
        a = dd[c]
        b = ref.loc[usable, c].to_numpy(float)
        ka4[c] = float(np.max(np.abs(a - b)))
    tot = dd["E_tot_half"].sum()
    tf = dd["X_shared"].sum() + dd["X_line"].sum() + 0.5 * (tot - dd["X_shared"].sum() - dd["X_line"].sum())
    pooled = dict(shared=float(dd["X_shared"].sum() / tf), line=float(dd["X_line"].sum() / tf),
                  noise=float(0.5 * (tot - dd["X_shared"].sum() - dd["X_line"].sum()) / tf))
    want = json.loads((L.W268 / "decomp268_cpm5.json").read_text())["cell_halves"]["pooled_full_depth"]
    ka4["pooled_full_depth"] = pooled
    ka4["pooled_full_depth_required"] = want
    ka4["pooled_max_abs_diff"] = float(max(abs(pooled[k] - want[k]) for k in want))
    # per_target_decomp_cpm5.csv is float32; the reproduction bar is float32 round-trip
    ka4["pass_"] = bool(max(ka4[c] for c in ("shared_full", "line_full", "noise_full")) < 1e-6
                        and ka4["pooled_max_abs_diff"] < 1e-6)
    P(f"KA-4 268's decomposition on {len(usable)} usable targets: max|diff| shared_full "
      f"{ka4['shared_full']:.1e}, line_full {ka4['line_full']:.1e}, noise_full {ka4['noise_full']:.1e}; "
      f"pooled_full_depth {pooled['shared']:.6f}/{pooled['line']:.6f}/{pooled['noise']:.6f} vs required "
      f"{want['shared']:.6f}/{want['line']:.6f}/{want['noise']:.6f} (max|diff| {ka4['pooled_max_abs_diff']:.1e}) "
      f"-> {'PASS' if ka4['pass_'] else 'FAIL'}")
    if not ka4["pass_"]:
        raise SystemExit("KA-4 FAILED — abort (AGENTS.md:1375-1382)")

    # ---- the unit (PREREG §3) ------------------------------------------------------------------
    store = L.Store()
    held = set(store.genes)
    allok = avail.all(1)
    avail_t = [t for t, a in zip(tg, allok) if a]
    unit = [t for t in avail_t if t not in ctrl and t in held]
    or_unit = [t for t in avail_t if t in ctrl and t in held]
    assert len(unit) == 82 and len(or_unit) == 7, (len(unit), len(or_unit))
    P(f"unit: {len(unit)} non-OR targets (available on all 12 lines, in the held map); "
      f"{len(or_unit)} OR controls in the held map; unavailable {sorted(set(tg) - set(avail_t))}")

    # ---- the conservation statistic (DEFINITION-276 §3) ----------------------------------------
    uidx = np.array([tpos[t] for t in unit])
    Ru = R["full"][uidx].astype(np.float64)                      # (82, 12, G)
    mu_line = Ru.mean(0)                                         # per line, mean over the 82 rows
    U = Ru - mu_line[None, :, :]
    nrm = np.linalg.norm(U, axis=2)
    U /= np.maximum(nrm, 1e-12)[:, :, None]
    pairs = [(i, j) for i in range(len(lines)) for j in range(i + 1, len(lines))]
    assert len(pairs) == 66
    Cp = np.stack([(U[:, i] * U[:, j]).sum(1) for i, j in pairs], axis=1)   # (82, 66)
    C_b16 = Cp.mean(1)
    P(f"C_broad16 over {len(pairs)} line pairs: median {np.median(C_b16):+.5f}, "
      f"IQR [{np.percentile(C_b16, 25):+.5f}, {np.percentile(C_b16, 75):+.5f}], "
      f"min {C_b16.min():+.5f}, max {C_b16.max():+.5f}, <=0 {int((C_b16 <= 0).sum())}")

    # OR controls, projected on the SAME line centring, rebuilt against the other 9 controls
    def ctrl_excl(l, c):
        tl = H[l]["targets"]
        ci = [tl.index(x) for x in ctrl if x != c]
        return B.logprof(np.nansum(H[l]["full"][ci], 0), ok)
    Ro = np.zeros((len(or_unit), len(lines), G))
    for a, c in enumerate(or_unit):
        for j, l in enumerate(lines):
            tl = H[l]["targets"]
            Ro[a, j] = B.logprof(H[l]["full"][tl.index(c)], ok) - ctrl_excl(l, c)
    Uo = Ro - mu_line[None, :, :]
    Uo /= np.maximum(np.linalg.norm(Uo, axis=2), 1e-12)[:, :, None]
    C_or = np.stack([(Uo[:, i] * Uo[:, j]).sum(1) for i, j in pairs], axis=1).mean(1)
    P(f"OR-control negative row (n {len(or_unit)}, control = the other 9): median "
      f"{np.median(C_or):+.5f}, range [{C_or.min():+.5f}, {C_or.max():+.5f}]")

    del Ru, U, Uo, Ro, R, H
    # secondary unit: the 82 + the 7 OR controls
    unit89 = unit + or_unit
    C89 = np.concatenate([C_b16, C_or])

    # ---- the models (PREREG §2) ----------------------------------------------------------------
    rem = np.load(L.REM284, allow_pickle=True)
    fold = np.asarray(rem["fold"], int)
    comp = np.asarray(rem["comp"], int)
    z = np.load(L.HELD279, allow_pickle=True)
    c_held = np.asarray(z["c_held"], float)

    agg = store.aggregate(store.subset())                        # ALL 28 banks
    Xr = agg[L.FEATS290].to_numpy(float)
    ann = L.load_annotation()
    assert list(ann["gene"]) == store.genes
    Xa = ann[L.FEATS].to_numpy(float)

    gpos = {g: i for i, g in enumerate(store.genes)}
    ui2 = np.array([gpos[t] for t in unit])
    ui89 = np.array([gpos[t] for t in unit89])

    pred_r, info_r = L.fit_apply_grouped(Xr, c_held, fold, Xr)
    pred_a, info_a = L.fit_apply_grouped(Xa, c_held, fold, Xa)
    P(f"M_resp alpha {info_r['alpha']:g} (grouped CV rho {info_r['grouped_cv_spearman']:+.4f}); "
      f"M_ann alpha {info_a['alpha']:g} (grouped CV rho {info_a['grouped_cv_spearman']:+.4f})")
    oof_r, _ = L.oof_ridge(Xr, c_held, fold)
    oof_a, _ = L.oof_ridge(Xa, c_held, fold)

    preds = {
        "M_resp": pred_r[ui2],
        "M_ann": pred_a[ui2],
        "C_held": c_held[ui2],
        "logmag": agg["logmag"].to_numpy(float)[ui2],
        "gen_pc20": agg["gen_pc20"].to_numpy(float)[ui2],
        "breadth": agg["breadth"].to_numpy(float)[ui2],
        "M_resp_oof": oof_r[ui2],
        "M_ann_oof": oof_a[ui2],
    }

    # ---- C3, primary outcome -------------------------------------------------------------------
    tgt = np.arange(len(unit))                                   # paired TARGET bootstrap
    bb = L.block_boot(C_b16, preds, tgt, np.random.default_rng(L.SEED))
    for k in preds:
        P(f"  C3 {k:12s} rho {bb[k]['point']:+.4f} [{bb[k]['lo']:+.4f}, {bb[k]['hi']:+.4f}]")
    pr = bb["_paired_rho"]
    for k in ("M_ann", "C_held", "logmag"):
        d = pr[f"M_resp-{k}"]
        P(f"  C3 paired M_resp - {k:8s} {d['mean']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}] "
          f"P(<=0) {d['p_le_0']:.3f}")

    # second outcome: 268's shared-energy fraction
    sh = ref["shared_full"].to_dict()
    m2 = np.array([t in sh for t in unit])
    preds2 = {k: v[m2] for k, v in preds.items()}
    y2 = np.array([sh[t] for t in unit if t in sh])
    bb2 = L.block_boot(y2, preds2, np.arange(len(y2)), np.random.default_rng(L.SEED))
    P(f"  C3 second outcome (268 shared_full, n {len(y2)}): "
      + "; ".join(f"{k} {bb2[k]['point']:+.4f} [{bb2[k]['lo']:+.4f}, {bb2[k]['hi']:+.4f}]"
                  for k in ("M_resp", "M_ann", "C_held", "gen_pc20")))

    # secondary unit with the OR controls
    preds89 = {"M_resp": pred_r[ui89], "M_ann": pred_a[ui89], "C_held": c_held[ui89],
               "gen_pc20": agg["gen_pc20"].to_numpy(float)[ui89]}
    bb89 = L.block_boot(C89, preds89, np.arange(len(unit89)), np.random.default_rng(L.SEED))
    P("  C3 secondary unit (82 + 7 OR): "
      + "; ".join(f"{k} {bb89[k]['point']:+.4f} [{bb89[k]['lo']:+.4f}, {bb89[k]['hi']:+.4f}]"
                  for k in preds89))

    # ---- null: 200 target shuffles --------------------------------------------------------------
    nullv = []
    for _ in range(200):
        nullv.append(L.spearman(rng.permutation(C_b16), preds["M_resp"]))
    nullv = np.asarray(nullv)
    obs = bb["M_resp"]["point"]
    nullrow = dict(mean=float(nullv.mean()), sd=float(nullv.std()),
                   n_reaching=int((nullv >= obs).sum()),
                   z=float((obs - nullv.mean()) / max(nullv.std(), 1e-12)))
    P(f"  C3 null (200 target shuffles): {nullrow['mean']:+.4f} +- {nullrow['sd']:.4f}; "
      f"{nullrow['n_reaching']}/200 reach {obs:+.4f}; z {nullrow['z']:+.1f}")

    # ---- depth (programme rule 9) ---------------------------------------------------------------
    Hn, _ = B.load_halves(lines)
    ncell = np.array([[int(Hn[l]["ncell"][Hn[l]["targets"].index(t), 0]) for l in lines]
                      for t in unit], float)
    del Hn
    med_cells, min_cells = np.median(ncell, 1), ncell.min(1)
    cq = pd.read_csv(L.CQCOV).set_index("cell_line")
    umi = np.array([cq.loc[l, "avg_umis_per_ko"] for l in lines], float)
    pair_cos = Cp.mean(0)
    pair_umi = np.array([np.sqrt(umi[i] * umi[j]) for i, j in pairs])
    depth = dict(
        rho_conservation_median_cells=L.spearman(C_b16, med_cells),
        rho_conservation_min_cells=L.spearman(C_b16, min_cells),
        rho_pair_cos_pair_umi=L.spearman(pair_cos, pair_umi),
        n_pairs=len(pairs),
        rho_C3_partial_log_med_cells=L.partial_spearman(preds["M_resp"], C_b16, np.log(med_cells)),
        rho_gen_pc20_partial_log_med_cells=L.partial_spearman(preds["gen_pc20"], C_b16,
                                                              np.log(med_cells)),
        umi_per_line={l: float(umi[i]) for i, l in enumerate(lines)},
    )
    P(f"  depth: rho(C_broad16, median cells/KO) {depth['rho_conservation_median_cells']:+.3f}; "
      f"min cells {depth['rho_conservation_min_cells']:+.3f}; "
      f"rho(pair cosine, pair geo-mean UMI/KO) over 66 pairs "
      f"{depth['rho_pair_cos_pair_umi']:+.3f}; C3 partial given log median cells/KO "
      f"{depth['rho_C3_partial_log_med_cells']:+.4f}")

    # ---- outputs ---------------------------------------------------------------------------------
    out = pd.DataFrame(dict(
        target=unit, C_broad16=C_b16, n_pairs=len(pairs),
        shared_full=[sh.get(t, np.nan) for t in unit],
        pred_M_resp=preds["M_resp"], pred_M_ann=preds["M_ann"],
        pred_M_resp_oof=preds["M_resp_oof"], pred_M_ann_oof=preds["M_ann_oof"],
        C_held=preds["C_held"], logmag=preds["logmag"], gen_pc20=preds["gen_pc20"],
        breadth=preds["breadth"], median_cells_per_ko=med_cells, min_cells_per_ko=min_cells,
        component=comp[ui2], is_OR_control=False))
    out_or = pd.DataFrame(dict(
        target=or_unit, C_broad16=C_or, n_pairs=len(pairs), shared_full=np.nan,
        pred_M_resp=pred_r[np.array([gpos[t] for t in or_unit])],
        pred_M_ann=pred_a[np.array([gpos[t] for t in or_unit])],
        pred_M_resp_oof=oof_r[np.array([gpos[t] for t in or_unit])],
        pred_M_ann_oof=oof_a[np.array([gpos[t] for t in or_unit])],
        C_held=c_held[np.array([gpos[t] for t in or_unit])],
        logmag=agg["logmag"].to_numpy(float)[np.array([gpos[t] for t in or_unit])],
        gen_pc20=agg["gen_pc20"].to_numpy(float)[np.array([gpos[t] for t in or_unit])],
        breadth=agg["breadth"].to_numpy(float)[np.array([gpos[t] for t in or_unit])],
        median_cells_per_ko=np.nan, min_cells_per_ko=np.nan,
        component=comp[np.array([gpos[t] for t in or_unit])], is_OR_control=True))
    pd.concat([out, out_or]).to_csv(L.WORK / "pred291_broad16.csv", index=False)

    # frozen rule
    rule = dict(
        cond1_M_resp_interval_above_0=bool(bb["M_resp"]["lo"] > 0),
        cond2_gen_pc20_positive_interval_above_0=bool(bb["gen_pc20"]["point"] > 0
                                                      and bb["gen_pc20"]["lo"] > 0))
    rule["C3_PASS"] = bool(rule["cond1_M_resp_interval_above_0"]
                           and rule["cond2_gen_pc20_positive_interval_above_0"])
    P(f"\nC3 {'PASS' if rule['C3_PASS'] else 'FAIL'} by the frozen rule: "
      f"M_resp {bb['M_resp']['point']:+.4f} [{bb['M_resp']['lo']:+.4f}, {bb['M_resp']['hi']:+.4f}]; "
      f"gen_pc20 {bb['gen_pc20']['point']:+.4f} [{bb['gen_pc20']['lo']:+.4f}, "
      f"{bb['gen_pc20']['hi']:+.4f}]")

    res = dict(prereg_sha256=h, ka4=ka4, n_unit=len(unit), n_or=len(or_unit),
               unit=unit, or_unit=or_unit, genes_kept=G, n_pairs=len(pairs),
               C_broad16=dict(median=float(np.median(C_b16)), q25=float(np.percentile(C_b16, 25)),
                              q75=float(np.percentile(C_b16, 75)), min=float(C_b16.min()),
                              max=float(C_b16.max()), n_le_0=int((C_b16 <= 0).sum())),
               C_or=dict(median=float(np.median(C_or)), min=float(C_or.min()),
                         max=float(C_or.max())),
               model_info=dict(M_resp=info_r, M_ann=info_a),
               primary=bb, second_outcome=bb2, secondary_unit_89=bb89,
               null=nullrow, depth=depth, rule=rule,
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30)
    (L.WORK / "surfA291.json").write_text(json.dumps(res, indent=1, default=float))
    (L.WORK / "surfA291.log").write_text("\n".join(LOG))
    P(f"peak RSS {res['peak_rss_gib']:.2f} GiB; wrote work/pred291_broad16.csv")


if __name__ == "__main__":
    main()
