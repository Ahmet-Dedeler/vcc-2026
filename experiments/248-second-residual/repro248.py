"""248 KNOWN-ANSWER ROW 3 and the PLOLO price -- the second union through the PRODUCTION BUILDER.

222's `repro222.py` / 249's `repro249.py` pattern: 162's PLOLO folds, 218's scoring, and the
builder's OWN objects substituted in -- `prepare_h1_union(sources=...)` for both unions and the
WHOLE of `apply_frozen_cd4_residual` (CD4 by its rule, union 1 and union 2 by the carrying-row
rule, ONE merged add, ONE preserve) for 225's `merged()` chain.  The builder is imported with the
arm's env (VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_UNION2_SOURCES=dld1_lowmoi
VCC_UNION2_G=0.35); `UNION2_G` is rebound per dose and `UNION2_SOURCES` is passed explicitly.

    KA-A  the base (cd4@0.34 + h1pair@0.28) through the builder == 218 compose218's
          `cd40.34+h1PAIRCARRY0.28`  (cos_agg, a23_cos, pds to 1e-9)
    KA-B  CUM1 at 225's six doses, real / SHUF / CD4ROWS variants, through the builder ==
          225's `stack225d.json` rows (d_cos, d_pds to 1e-12; d_qw reported) -- the mandate's
          row 3 at g = 0.35: rel macro cos +1.405%, d pds +0.00256 over the 3 cancer folds
    KA-C  CUM4 likewise (225's four as union 2)
    KA-D  the builder's output at g = 0.35 == compose218.merged() chained by hand, bit for bit
    KA-E  225's leave-top-3-out on dld1 (pds 25.2%, cos 54.5%) and CUM4 reproduces (conc225)

    FALSIFIERS (all through the builder's function, marginal over the incumbent h1union):
      target-shuffle at every dose (two seeds), row-matched CD4 control (CD4ROWS), the
      non-carrying-row split (rowsplit), leave-top-k-out on cosine and pds SEPARATELY,
      `kolf2` apart.  Every arm also reports rel cos_agg (243's functional) and d q_w.

One fold per process (memory); `--summary` aggregates.

    uv run python repro248.py --fold k562_gwps      (orion_hct116, orion_hek293t, kolf2)
    uv run python repro248.py --summary
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
os.environ["VCC_UNION2_SOURCES"] = "dld1_lowmoi"
os.environ["VCC_UNION2_G"] = "0.35"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_CD4_DEFLATE_K"):
    os.environ.pop(_k, None)

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E218 = ROOT / "experiments/218-independent-donors"
E225 = ROOT / "experiments/225-source-sweep"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"          # this experiment's copies of 229's fold caches (sha-identical)

from plolo import (FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, solve_fraction)
from resid218 import ALPHA, GAMMA, LAM, score_effects, shuffled  # noqa: E402
from compose218 import merged, restrict_rows  # noqa: E402

H1_PAIR = ("gse295214_h1", "h1_2025")
CUM = {"CUM1": ("dld1_lowmoi",),
       "CUM4": ("dld1_lowmoi", "jiang_5line_noK562", "igvf_cardiomyocyte", "myllia_a549")}
G = (0.10, 0.17, 0.28, 0.35, 0.45, 0.60)
BUILT_G = 0.35
SEEDS = (20726, 771103)             # 225's F1 seeds; stack225d carries the first
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KS = (1, 3, 5, 10)
TOL = 1e-9
TOL_KA = 1e-12


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("repro248_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repro248_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def public(e: dict) -> dict:
    return {k: v for k, v in e.items() if not k.startswith("_")}


def sources_of(MOD, names):
    return tuple((n, MOD.UNION_SOURCE_REGISTRY[n]["path"]) for n in names)


def carrying(res) -> np.ndarray:
    return np.asarray(res.target_condition_count) > 0


def run_fold(held: str) -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34 and MOD.UNION2_G == 0.35
    assert [n for n, _ in MOD.UNION2_SOURCES] == ["dld1_lowmoi"]
    print(f"builder {builder_path}  sha256 {MOD.file_sha256(builder_path)}", flush=True)
    ref218 = {(r["fold"], r["arm"]): r
              for r in json.loads((E218 / "work/compose218.json").read_text())}
    ref225 = {(r["fold"], r["tag"], r.get("variant", ""), r["g"]): r
              for r in json.loads((E225 / "work/stack225d.json").read_text())}

    fold = build_fold(held)
    targets = list(fold["targets"])
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    pp = np.asarray(fold["panel_positions"], dtype=np.int64)
    cd4 = cd4_for(targets)
    op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
               else PROD_OPS / "whitener_union_pruned.npz")
    whitened = Operator(op_path).apply(pre * ALPHA, wcols, GAMMA, LAM)
    frac_cd4_218, _ = solve_fraction(whitened, cd4, 0.34)

    out: dict = dict(fold=held, n_targets=len(targets), builder_sha256=MOD.file_sha256(builder_path),
                     ka={}, arms=[], rowsplit=[], conc={}, failures=[])

    def fail(msg):
        print(f"  *** {msg}", flush=True)
        out["failures"].append(msg)

    # ---- the unions the BUILD makes -------------------------------------------------------
    h1, h1_details = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0,
                                          sources=sources_of(MOD, H1_PAIR))
    unions = {}
    for tag, names in CUM.items():
        unions[tag], _d = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0,
                                               sources=sources_of(MOD, names))
        print(f"[{held}] {tag} union {names}: carries {int(carrying(unions[tag]).sum())}/{len(targets)}",
              flush=True)
    print(f"[{held}] H1 pair carries {int(carrying(h1).sum())}/{len(targets)}; staged "
          f"({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)

    def through_builder(u2, g2):
        """The production entry point: CD4 by its rule, union 1 at 0.28, union 2 at g2."""
        MOD.UNION2_G = g2
        audits = {}
        final, cd4_audit, u1_audit = MOD.apply_frozen_cd4_residual(
            whitened, cd4, wcols, pp, extra=h1, extra_audit={},
            extra2=u2, extra2_audit=None if u2 is None else {}, audits_out=audits)
        final[:, pp] = 0.0
        return final, cd4_audit, u1_audit, audits.get("union2")

    # ---- KA-A: the incumbent through the builder == 218 ------------------------------------
    base_final, base_cd4, base_u1, _ = through_builder(None, None)
    base = score_effects(base_final, fold, truth)
    ref = ref218[(held, "cd40.34+h1PAIRCARRY0.28")]
    d = max(abs(ref[k] - base[k]) for k in ("cos_agg", "a23_cos", "pds"))
    out["ka"]["A_base_vs_218"] = dict(max_abs_diff=float(d), ok=bool(d < TOL),
                                      frac_cd4=base_cd4["fraction"],
                                      frac_h1=base_u1["residual_frame"]["solved_fraction"])
    print(f"  {'PASS' if d < TOL else 'FAIL'} KA-A base through the builder vs 218 compose218: "
          f"max|d| {d:.2e}  (cd4 frac {base_cd4['fraction']:.6f} == 218 {frac_cd4_218:.6f}: "
          f"{abs(base_cd4['fraction'] - frac_cd4_218):.1e})", flush=True)
    if d >= TOL:
        fail(f"KA-A: base differs from 218 by {d:.3e}")
    b_pds = np.asarray(base["_pds_per_target_a23"], dtype=np.float64)
    b_cos = np.asarray(base["_cos_per_target_a23"], dtype=np.float64)
    del base_final

    # ---- the arms and their controls, every dose ------------------------------------------
    for tag, u2 in unions.items():
        variants = {"": u2, "CD4ROWS": restrict_rows(cd4, u2)}
        for seed in SEEDS:
            variants[f"SHUF{seed}"] = shuffled(u2, seed)
        carry = carrying(u2)
        for g in G:
            for vtag, res in variants.items():
                final, cd4_a, u1_a, u2_a = through_builder(res, g)
                e = score_effects(final, fold, truth)
                fr = u2_a["residual_frame"]
                row = dict(
                    fold=held, arm=tag, variant=vtag, g=g,
                    fraction=fr["solved_fraction"], per_unit=fr["g_per_unit_fraction"],
                    carrying_rows=fr["carrying_rows"], n_rows=len(targets),
                    coverage=fr["coverage_frac"],
                    shipped_median_rule_would_solve=fr["shipped_median_rule_would_solve"],
                    frac_cd4=cd4_a["fraction"], frac_h1=u1_a["residual_frame"]["solved_fraction"],
                    preserve_scale=u2_a["merged_preserve"]["post_add_global_preserve_scale"],
                    d_cos=e["cos"] - base["cos"],
                    rel_cos_pct=100.0 * (e["cos"] - base["cos"]) / max(abs(base["cos"]), 1e-300),
                    d_cos_agg=e["cos_agg"] - base["cos_agg"],
                    rel_cos_agg_pct=100.0 * (e["cos_agg"] - base["cos_agg"]) / max(abs(base["cos_agg"]), 1e-300),
                    d_pds=e["pds"] - base["pds"],
                    d_qw=e["a15_sgn_w"] - base["a15_sgn_w"],
                    d_a15_cos_agg=e["a15_cos_agg"] - base["a15_cos_agg"],
                    rel_a15_cos_agg_pct=100.0 * (e["a15_cos_agg"] - base["a15_cos_agg"]) / max(abs(base["a15_cos_agg"]), 1e-300),
                    **{f"arm_{k}": v for k, v in public(e).items()},
                    **{f"base_{k}": v for k, v in public(base).items()},
                )
                # CD4 and union 1 must be untouched by union 2 -- the single-variable constraint
                if abs(row["frac_cd4"] - base_cd4["fraction"]) != 0.0 or \
                        row["frac_h1"] != base_u1["residual_frame"]["solved_fraction"]:
                    fail(f"{tag}{vtag}@{g}: CD4 or union-1 fraction moved")
                # 225's row, where it exists (stack225d: variants '', 'SHUF' (seed 20726), 'CD4ROWS')
                key225 = (held, tag, {"": "", "SHUF20726": "SHUF", "CD4ROWS": "CD4ROWS"}.get(vtag, "?"), g)
                r = ref225.get(key225)
                if r is not None:
                    dc = row["d_cos"] - r["d_cos"]
                    dp = row["d_pds"] - r["d_pds"]
                    dq = row["d_qw"] - r["d_qw"]
                    ok = abs(dc) < TOL_KA and abs(dp) < TOL_KA
                    out["ka"][f"B_{tag}{vtag or 'real'}@{g:g}_vs_225"] = dict(
                        d_cos_delta=dc, d_pds_delta=dp, d_qw_delta=dq, ok=bool(ok),
                        ref_d_cos=r["d_cos"], ref_d_pds=r["d_pds"],
                        ref_fraction=r["audit"]["fraction"], our_fraction=row["fraction"])
                    row["vs225_ok"] = bool(ok)
                    if not ok:
                        fail(f"KA-B {tag}{vtag}@{g}: 225's row did not reproduce (d_cos {dc:+.2e}, d_pds {dp:+.2e})")
                if vtag == "":
                    row["rowsplit"] = dict(
                        n_carry=int(carry.sum()),
                        rel_in=float(np.mean(np.asarray(e["_cos_per_target_a23"])[carry]) / b_cos[carry].mean() - 1),
                        rel_out=float(np.mean(np.asarray(e["_cos_per_target_a23"])[~carry]) / b_cos[~carry].mean() - 1),
                        d_pds_in=float((np.asarray(e["_pds_per_target_a23"])[carry] - b_pds[carry]).mean()),
                        d_pds_out=float((np.asarray(e["_pds_per_target_a23"])[~carry] - b_pds[~carry]).mean()),
                    )
                    if g == BUILT_G:
                        out["conc"][tag] = dict(
                            targets=targets, carry=carry.tolist(),
                            d_pds=(np.asarray(e["_pds_per_target_a23"]) - b_pds).tolist(),
                            d_cos=(np.asarray(e["_cos_per_target_a23"]) - b_cos).tolist(),
                            fold_d_pds=float(e["pds"] - base["pds"]),
                            fold_d_cos=float(e["cos"] - base["cos"]))
                out["arms"].append(row)
                tagv = f"{tag}{('/' + vtag) if vtag else ''}"
                print(f"    {tagv:16s} g={g:<5g} fx {row['fraction']:.5f} carry {row['carrying_rows']:3d}  "
                      f"macro {row['rel_cos_pct']:+7.3f}%  agg {row['rel_cos_agg_pct']:+7.3f}%  "
                      f"dpds {row['d_pds']:+.6f}  dqw {row['d_qw']:+.5f}"
                      + (f"  vs225 {'OK' if row.get('vs225_ok') else 'FAIL'}" if r is not None else ""),
                      flush=True)
                del final

        # ---- KA-D: the builder's arm at 0.35 == compose218.merged() chained by hand --------
        final_b, _, _, u2_a = through_builder(u2, BUILT_G)
        f1, _ = MOD.solve_carrying_row_fraction(whitened, h1, 0.28)
        f2, _ = MOD.solve_carrying_row_fraction(whitened, u2, BUILT_G)
        chained = merged(merged(cd4, frac_cd4_218, h1, f1), 1.0, u2, f2)
        ref_final, _p = plolo.E062.add_post_whiten_residual(whitened, chained, 1.0, wcols)
        ref_final[:, pp] = 0.0
        same = bool(np.array_equal(final_b, ref_final))
        out["ka"][f"D_{tag}@{BUILT_G:g}_builder_vs_merged_chain"] = dict(
            bit_identical=same, max_abs_diff=float(np.abs(final_b - ref_final).max()))
        print(f"  {'PASS' if same else 'FAIL'} KA-D {tag}@{BUILT_G}: builder output == compose218.merged() "
              f"chained + one add + one preserve: max|d| {np.abs(final_b - ref_final).max():.1e}", flush=True)
        if not same:
            fail(f"KA-D {tag}: builder != merged chain")
        del final_b, ref_final, variants

    out["peak_rss_gib"] = peak_gib()
    out["seconds"] = round(time.time() - t0, 1)
    (WORK / f"repro248-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"[{held}] {len(out['failures'])} failures; {out['seconds']}s; peak {peak_gib():.2f} GiB", flush=True)
    if out["failures"]:
        raise SystemExit(1)


def summary() -> None:
    folds = {}
    for held in FOLDS:
        p = WORK / f"repro248-{held}.json"
        if p.exists():
            folds[held] = json.loads(p.read_text())
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    fails = {h: f["failures"] for h, f in folds.items()}
    P(f"folds present: {sorted(folds)}; failures: { {h: len(v) for h, v in fails.items()} }")
    ka_all = all(v["ok"] for f in folds.values() for v in f["ka"].values()
                 if "ok" in v) and all(v.get("bit_identical", True) for f in folds.values()
                                       for v in f["ka"].values())
    n_ka = sum(len(f["ka"]) for f in folds.values())
    P(f"known-answer rows: {n_ka} across folds, all pass: {ka_all}")
    conc225 = json.loads((E225 / "work/conc225.json").read_text())["summary"]

    res = {"folds": sorted(folds), "ka_all_pass": bool(ka_all), "n_ka": n_ka, "arms": {}, "conc": {}}
    for tag in CUM:
        P(f"\n=== {tag}: union 2 = {CUM[tag]} on top of the shipped h1union, three cancer folds "
          f"(kolf2 apart) ===")
        P(f"  {'variant':12s}{'g':>6s}{'rel macro':>11s}{'rel agg':>10s}{'rel a15agg':>11s}{'d pds':>11s}"
          f"{'pds+':>6s}{'d q_w':>10s}{'carry':>18s}{'vs225':>7s}{'kolf2 macro':>13s}{'kolf2 dpds':>12s}")
        for vtag in ("", "SHUF20726", "SHUF771103", "CD4ROWS"):
            for g in G:
                rows = [r for h in CANCER if h in folds for r in folds[h]["arms"]
                        if r["arm"] == tag and r["variant"] == vtag and r["g"] == g]
                if len(rows) != 3:
                    continue
                kr = [r for r in folds.get("kolf2", {"arms": []})["arms"]
                      if r["arm"] == tag and r["variant"] == vtag and r["g"] == g]
                rec = dict(
                    rel_macro=float(np.mean([r["rel_cos_pct"] for r in rows])) / 100,
                    rel_agg=float(np.mean([r["rel_cos_agg_pct"] for r in rows])) / 100,
                    rel_a15_agg=float(np.mean([r["rel_a15_cos_agg_pct"] for r in rows])) / 100,
                    d_pds=float(np.mean([r["d_pds"] for r in rows])),
                    pds_pos=int(sum(r["d_pds"] > 0 for r in rows)),
                    d_qw=float(np.mean([r["d_qw"] for r in rows])),
                    per_fold={r["fold"]: dict(rel_macro=r["rel_cos_pct"] / 100, rel_agg=r["rel_cos_agg_pct"] / 100,
                                              d_pds=r["d_pds"], d_qw=r["d_qw"], fraction=r["fraction"],
                                              carrying=r["carrying_rows"], n_rows=r["n_rows"],
                                              preserve_scale=r["preserve_scale"]) for r in rows},
                    vs225=[r.get("vs225_ok") for r in rows],
                    kolf2=None if not kr else dict(rel_macro=kr[0]["rel_cos_pct"] / 100, rel_agg=kr[0]["rel_cos_agg_pct"] / 100,
                                                   d_pds=kr[0]["d_pds"], d_qw=kr[0]["d_qw"]),
                )
                res["arms"][f"{tag}{vtag}@{g:g}"] = rec
                carry = "/".join(str(r["carrying_rows"]) for r in rows)
                v225 = "".join("Y" if x else ("n" if x is False else "-") for x in rec["vs225"])
                k2 = rec["kolf2"]
                P(f"  {vtag or 'real':12s}{g:>6g}{100*rec['rel_macro']:>+10.3f}%{100*rec['rel_agg']:>+9.3f}%"
                  f"{100*rec['rel_a15_agg']:>+10.3f}%{rec['d_pds']:>+11.6f}{rec['pds_pos']:>4d}/3{rec['d_qw']:>+10.5f}"
                  f"{carry:>18s}{v225:>7s}"
                  + (f"{100*k2['rel_macro']:>+12.3f}%{k2['d_pds']:>+12.6f}" if k2 else ""))
        # rowsplit at the built dose
        rows = [r for h in CANCER if h in folds for r in folds[h]["arms"]
                if r["arm"] == tag and r["variant"] == "" and r["g"] == BUILT_G]
        if len(rows) == 3:
            rs = dict(rel_in=float(np.mean([r["rowsplit"]["rel_in"] for r in rows])),
                      rel_out=float(np.mean([r["rowsplit"]["rel_out"] for r in rows])),
                      d_pds_in=float(np.mean([r["rowsplit"]["d_pds_in"] for r in rows])),
                      d_pds_out=float(np.mean([r["rowsplit"]["d_pds_out"] for r in rows])),
                      max_abs_rel_out_pp=float(100 * max(abs(r["rowsplit"]["rel_out"]) for r in rows)),
                      per_fold={r["fold"]: r["rowsplit"] for r in rows})
            res["arms"][f"{tag}@{BUILT_G:g}"]["rowsplit"] = rs
            P(f"  rowsplit @{BUILT_G}: carrying rows rel cos {100*rs['rel_in']:+.3f}%, NOT-carrying "
              f"{100*rs['rel_out']:+.5f}% (largest |excursion| {rs['max_abs_rel_out_pp']:.5f} pp; 218's bar 0.018 pp); "
              f"d pds carrying {rs['d_pds_in']:+.6f}, NOT-carrying {rs['d_pds_out']:+.6f}")
        # leave-top-k-out, 222/225's definition, at the built dose
        pf = {h: folds[h]["conc"][tag] for h in CANCER if h in folds and tag in folds[h]["conc"]}
        if len(pf) == 3:
            full_p = float(np.mean([pf[h]["fold_d_pds"] for h in CANCER]))
            full_c = float(np.mean([pf[h]["fold_d_cos"] for h in CANCER]))
            rec = dict(full_d_pds=full_p, full_d_cos=full_c, k={})
            P(f"  leave-top-k-out @{BUILT_G} (each channel ranked by its OWN per-target delta; mean of per-fold LOO means):")
            P(f"    {'k':>3s}{'d pds full':>12s}{'d pds LOO':>11s}{'pds ret':>9s}{'top share':>11s}"
              f"{'d cos full':>12s}{'d cos LOO':>11s}{'cos ret':>9s}")
            for k in KS:
                lp, lc, sp, sc, top = [], [], [], [], {}
                for h in CANCER:
                    dp = np.asarray(pf[h]["d_pds"]); dc = np.asarray(pf[h]["d_cos"])
                    for arr, L, S in ((dp, lp, sp), (dc, lc, sc)):
                        o = np.argsort(-arr)[:k]
                        keep = np.ones(arr.size, bool); keep[o] = False
                        L.append(float(arr[keep].mean()))
                        S.append(float(arr[o].sum() / arr.sum()) if arr.sum() else float("nan"))
                    if k == 3:
                        o = np.argsort(-dp)[:3]
                        top[h] = [pf[h]["targets"][i] for i in o]
                mp, mc = float(np.mean(lp)), float(np.mean(lc))
                rec["k"][str(k)] = dict(loo_d_pds=mp, loo_d_cos=mc,
                                        pds_retained=mp / full_p if full_p else float("nan"),
                                        cos_retained=mc / full_c if full_c else float("nan"),
                                        top_share_pds=float(np.mean(sp)), top_share_cos=float(np.mean(sc)),
                                        **({"top3_targets_by_pds": top} if k == 3 else {}))
                r = rec["k"][str(k)]
                P(f"    {k:>3d}{full_p:>+12.5f}{mp:>+11.5f}{100*r['pds_retained']:>8.1f}%{100*r['top_share_pds']:>10.1f}%"
                  f"{full_c:>+12.6f}{mc:>+11.6f}{100*r['cos_retained']:>8.1f}%")
            P(f"    top-3 d pds targets: " + " | ".join(f"{h[:8]}: {','.join(v)}" for h, v in rec['k']['3']['top3_targets_by_pds'].items()))
            # KA-E: 225's conc row for the same construction (dld1 @0.35; CUM4 was 0.28 in conc225 -- not comparable)
            ref_tag = {"CUM1": "dld1_lowmoi"}.get(tag)
            if ref_tag is not None:
                r3 = conc225[ref_tag]["k"]["3"]
                d_ret = (abs(rec["k"]["3"]["pds_retained"] - r3["pds_retained"]),
                         abs(rec["k"]["3"]["cos_retained"] - r3["cos_retained"]))
                ok = max(d_ret) < 1e-6 and abs(full_p - conc225[ref_tag]["full_d_pds"]) < 1e-9
                rec["ka_E_vs_conc225"] = dict(ok=bool(ok), ref_pds_retained=r3["pds_retained"],
                                              ref_cos_retained=r3["cos_retained"], d=d_ret)
                P(f"    {'PASS' if ok else 'FAIL'} KA-E leave-top-3-out vs 225 conc225 ({ref_tag}@0.35): pds retained "
                  f"{100*rec['k']['3']['pds_retained']:.1f}% (225: {100*r3['pds_retained']:.1f}%), cos retained "
                  f"{100*rec['k']['3']['cos_retained']:.1f}% (225: {100*r3['cos_retained']:.1f}%)")
            res["conc"][tag] = rec
    (WORK / "repro248-summary.json").write_text(json.dumps(res, indent=1, default=float))
    (WORK / "repro248-summary.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'repro248-summary.json'}")
    if not ka_all or any(fails.values()):
        raise SystemExit(1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    if a.summary:
        summary()
    else:
        run_fold(a.fold)


if __name__ == "__main__":
    main()
