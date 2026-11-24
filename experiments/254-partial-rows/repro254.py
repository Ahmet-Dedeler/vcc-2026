"""254 KNOWN-ANSWER ROWS 2 and 3 -- the builder's partial-row residual through 162's PLOLO.

246's rows re-derived with the PRODUCTION builder's code substituted in (222's `repro222.py` /
249's `repro249.py` pattern): `extract_cd4_partial` for `lib246.cd4_partial`, and the WHOLE of
`apply_frozen_cd4_residual` ([deflate ->] size -> merge -> one preserve) for 246's `compose`.
Env is read at the builder's import, so ONE PROCESS = ONE CONFIGURATION:

  --mode partial   VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_PARTIAL_G=0.34
      KA-A  STACK (lib246)                       == 246 sweep246 `STACK`           (cos_agg, a23_cos, pds; 1e-9)
      KA-B  builder extract_cd4_partial(prep)    == lib246.cd4_partial(targets, equalise=False), bit for bit;
            builder freeze_complete_mean3(prep)  == plolo.cd4_for(targets), bit for bit (the same prep)
      KA-C  the partial arm through the builder  == 246 `C2_partial@0.34` (1e-9); its fraction == 246's;
            CD4's and the union's fractions untouched; effect matrix == 246's compose, max|d| 0
  --mode composed  ... + VCC_CD4_DEFLATE_K=5
      KA-K  K=5 alone through the builder        == 246 `B2_cd4_svd5_resolved` (249's row, re-verified)
      KA-R  the COMPOSED arm, the builder's real path (deflated CD4 at the RE-SOLVED fraction + union +
            partial): CD4's fraction/deflation == KA-K's, the partial's fraction == KA-C's, effect matrix
            == lib246 compose at the re-solved fraction, max|d| 0 -- the object the package carries
      KA-E  the same with `solve_cd4_dose` pinned to size the UN-deflated residual (= the stack's FIXED
            fraction, 246's construction of `E_cd4_svd5+partial@0.34`) == 246's E row (1e-9);
            effect matrix == 246's compose, max|d| 0

Per-target `pds` / cosine arrays are saved for the leave-top-k-out statistic (`--summary`).
One fold per process (memory; 246 measured 4.1-5.0 GiB).

    uv run python repro254.py --mode partial  --fold k562_gwps   (orion_hct116, orion_hek293t, kolf2)
    uv run python repro254.py --mode composed --fold k562_gwps   (...)
    uv run python repro254.py --summary
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

ap = argparse.ArgumentParser()
ap.add_argument("--fold")
ap.add_argument("--mode", choices=("partial", "composed"))
ap.add_argument("--summary", action="store_true")
ARGS = ap.parse_args()

os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
os.environ["VCC_CD4_PARTIAL_G"] = "0.34"
if ARGS.mode == "composed":
    os.environ["VCC_CD4_DEFLATE_K"] = "5"
else:
    os.environ.pop("VCC_CD4_DEFLATE_K", None)
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_CD4_DEFLATE_BASIS",
           "VCC_UNION2_SOURCES", "VCC_UNION2_G"):
    os.environ.pop(_k, None)

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E246 = ROOT / "experiments/246-direction-hypotheses"
sys.path.insert(0, str(E246))
import lib246 as L  # noqa: E402  (imports plolo / resid218 / compose218; 246's construction)
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"          # this experiment's copies of the fold caches (sha-identical)
L.WORK = WORK

TOL = 1e-9
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
KEYS = ("cos_agg", "a23_cos", "pds")
G = 0.34
K = 5


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("repro254_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repro254_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def public(e: dict) -> dict:
    return {k: v for k, v in e.items() if not k.startswith("_")}


def same_prepared(a, b) -> bool:
    return (np.array_equal(np.asarray(a.mean_effect), np.asarray(b.mean_effect))
            and np.array_equal(np.asarray(a.gene_votes), np.asarray(b.gene_votes))
            and np.array_equal(np.asarray(a.target_condition_count), np.asarray(b.target_condition_count)))


def run_fold(held: str, mode: str) -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.CD4_PARTIAL_G == G and MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    assert MOD.CD4_DEFLATE_K == (K if mode == "composed" else None)
    print(f"builder {builder_path}  sha256 {MOD.file_sha256(builder_path)}  mode {mode}", flush=True)
    ref_a = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    ref_b = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}

    st = L.stage(held)
    fold, targets = st["fold"], st["targets"]
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    pp = np.asarray(fold["panel_positions"], dtype=np.int64)
    print(f"[{held}] staged {len(targets)} targets; stack frac_cd4 {st['frac_cd4']:.6f} "
          f"frac_x {st['frac_x']:.6f}  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    out: dict = dict(fold=held, mode=mode, n_targets=len(targets), builder_sha256=MOD.file_sha256(builder_path),
                     frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"], ka={}, failures=[])
    per_target: dict[str, np.ndarray] = {}

    def fail(msg):
        print(f"  *** {msg}", flush=True)
        out["failures"].append(msg)

    def ka(name, e, ref):
        d = max(abs(ref[k] - e[k]) for k in KEYS)
        ok = d < TOL
        out["ka"][name] = dict(max_abs_diff=float(d), ok=bool(ok), ref_arm=ref["arm"])
        print(f"  {'PASS' if ok else 'FAIL'} {name:60s} vs 246 {ref['arm']:26s} max|d| {d:.2e}", flush=True)
        if not ok:
            fail(f"{name}: differs from 246 {ref['arm']} by {d:.3e}")

    def ka_bool(name, ok, detail=""):
        out["ka"][name] = dict(ok=bool(ok), detail=detail)
        print(f"  {'PASS' if ok else 'FAIL'} {name:60s} {detail}", flush=True)
        if not ok:
            fail(name)

    def score_final(final):
        final = np.asarray(final, dtype=np.float64).copy()
        final[:, pp] = 0.0
        return L.score(st, final)

    def keep(name, e):
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])

    # ---------------------------------------------------------------- KA-A: the stack
    e_stack = L.score(st, L.compose(st, L.stack_parts(st)))
    ka("KA-A STACK (lib246)", e_stack, ref_a["STACK"])
    out["stack"] = public(e_stack)
    keep("STACK", e_stack)

    # ---------------------------------------------------------------- KA-B: the partial object itself
    banks = [L.E062.Bank.load(n, L.CD4_ROOT / n, require_cd4_scale=True) for n in L.CD4_NAMES]
    prep = L.E062.prepare_cd4_mean(banks, targets, L.OFFICIAL, L.PANEL_SET, "gene_norm",
                                   chunk_rows=256, max_rss_gib=40.0)
    del banks
    part246 = L.cd4_partial(targets, equalise=False)
    partB = MOD.extract_cd4_partial(prep)
    frozenB = MOD.freeze_complete_mean3(prep)
    del prep
    ka_bool("KA-B1 extract_cd4_partial(prep) == lib246.cd4_partial (arrays + details)",
            same_prepared(partB, part246) and partB.condition_details == part246.condition_details,
            f"n_partial {partB.condition_details[0]['n_partial']} (one {partB.condition_details[0]['n_one']}, "
            f"two {partB.condition_details[0]['n_two']}), median complete norm "
            f"{partB.condition_details[0]['median_complete_norm']:.6f}")
    ka_bool("KA-B2 freeze_complete_mean3(prep) == plolo.cd4_for(targets) (the complete-row object)",
            same_prepared(frozenB, st["cd4"]))
    del frozenB
    paudit = MOD.cd4_partial_support_audit(partB, target_g=G)
    fr246, per_unit246 = L.solve_fraction_carrying(st["whitened"], part246, G)
    out["partial_info"] = partB.condition_details[0]
    out["partial_audit"] = paudit
    out["partial_fraction_246"] = float(fr246)

    def builder_arm(partial, *, fixed=False):
        ao: dict = {}
        original = MOD.solve_cd4_dose
        if fixed:
            # 246's E construction: the deflated CD4 residual at the STACK's fraction, i.e. the
            # all-rows rule solved on the UN-deflated residual (exactly `st["frac_cd4"]`).
            MOD.solve_cd4_dose = lambda w, c: MOD.solve_cd4_fraction(w, st["cd4"])
        try:
            o, ca, ua = MOD.apply_frozen_cd4_residual(
                st["whitened"], st["cd4"], wcols, pp, extra=st["union"], extra_audit={},
                partial=partial, partial_audit=(paudit if partial is not None else None),
                audits_out=ao)
        finally:
            MOD.solve_cd4_dose = original
        return o, ca, ua, ao

    def chk_common(tag, ca, ua, ao, want_cd4_fraction):
        d_cd4 = abs(float(ca["fraction"]) - want_cd4_fraction)
        ka_bool(f"{tag}: CD4's fraction {ca['fraction']:.12f} == expected {want_cd4_fraction:.12f}", d_cd4 < 1e-12,
                f"|d| {d_cd4:.1e}")
        d_union = abs(float(ua["residual_frame"]["solved_fraction"]) - st["frac_x"])
        ka_bool(f"{tag}: the union's solved fraction == the stack's", d_union < 1e-12, f"|d| {d_union:.1e}")
        if ao.get("cd4_partial") is not None:
            fp = float(ao["cd4_partial"]["residual_frame"]["solved_fraction"])
            ka_bool(f"{tag}: the partial's solved fraction {fp:.12f} == lib246 solve_fraction_carrying {fr246:.12f}",
                    fp == fr246, f"|d| {abs(fp - fr246):.1e}; shipped all-rows rule would solve "
                                 f"{ao['cd4_partial']['residual_frame']['shipped_median_rule_would_solve']!r}")

    def chk_effects(tag, o, ref_final):
        chk = np.asarray(o, dtype=np.float64).copy()
        chk[:, pp] = 0.0
        d = float(np.abs(chk - ref_final).max())
        ka_bool(f"{tag}: effect matrix, builder path vs lib246 compose, max|d| == 0", d == 0.0, f"max|d| {d:.2e}")
        return d

    def report(tag, e, base=e_stack):
        print(f"  {tag}: macro {100*(e['a23_cos']/base['a23_cos']-1):+.4f}%  agg "
              f"{100*(e['cos_agg']/base['cos_agg']-1):+.4f}%  d pds {e['pds']-base['pds']:+.6f}  "
              f"d q_w {e['a15_sgn_w']-base['a15_sgn_w']:+.6f}", flush=True)

    if mode == "partial":
        # ------------------------------------------------------------ KA-C: the partial arm
        o, ca, ua, ao = builder_arm(partB)
        e = score_final(o)
        ka("KA-C the partial arm through apply_frozen_cd4_residual", e, ref_b["C2_partial@0.34"])
        fr_ref = float(ref_b["C2_partial@0.34"]["audit"]["fraction"])
        ka_bool("KA-C the partial's fraction == 246's json fraction", abs(float(ao["cd4_partial"]["residual_frame"]["solved_fraction"]) - fr_ref) < 1e-12,
                f"{ao['cd4_partial']['residual_frame']['solved_fraction']:.12f} vs {fr_ref:.12f}")
        chk_common("KA-C", ca, ua, ao, st["frac_cd4"])
        ref_final = L.compose(st, L.stack_parts(st) + [(part246, fr246)])
        out["partial_vs_246_compose_max_abs_diff"] = chk_effects("KA-C", o, ref_final)
        report("PARTIAL", e)
        out["partial"] = public(e)
        out["partial_frame"] = ao["cd4_partial"]["residual_frame"]
        out["partial_block"] = {k: v for k, v in ao["cd4_partial"].items() if k != "residual_frame"}
        out["cd4_audit"] = {k: v for k, v in ca.items() if k not in ("residual_frame",)}
        keep("partial", e)
    else:
        d246, a246 = L.deflate(st["cd4"], "svd", K)
        fr_res, _ = L.solve_fraction(st["whitened"], d246, L.CD4_F)
        # ------------------------------------------------------------ KA-K: K=5 alone (249's arm)
        o_k, ca_k, ua_k, ao_k = builder_arm(None)
        e_k = score_final(o_k)
        ka("KA-K cd4svd5 (K=5 alone) through the builder", e_k, ref_b["B2_cd4_svd5_resolved"])
        chk_common("KA-K", ca_k, ua_k, ao_k, fr_res)
        ka_bool("KA-K no `cd4_partial` block when no partial is handed in", "cd4_partial" not in ao_k)
        report("SVD5 (resolved, 249's arm)", e_k)
        out["svd5"] = public(e_k)
        out["svd5_fraction"] = float(ca_k["fraction"])
        out["deflation_audit"] = ca_k["deflation"]
        keep("svd5", e_k)
        # ------------------------------------------------------------ KA-R: the composed arm, the builder's real path
        o_r, ca_r, ua_r, ao_r = builder_arm(partB)
        e_r = score_final(o_r)
        chk_common("KA-R", ca_r, ua_r, ao_r, fr_res)
        ka_bool("KA-R CD4's deflation audit identical to KA-K's (the partial never touches the deflation)",
                json.dumps(ca_r["deflation"], sort_keys=True) == json.dumps(ca_k["deflation"], sort_keys=True))
        ref_r = L.compose(st, [(d246, fr_res), (st["union"], st["frac_x"]), (part246, fr246)])
        out["composed_resolved_vs_lib246_compose_max_abs_diff"] = chk_effects("KA-R", o_r, ref_r)
        report("COMPOSED (resolved -- the package's object)", e_r)
        report("COMPOSED (resolved) vs SVD5", e_r, base=e_k)
        out["composed_resolved"] = public(e_r)
        out["composed_resolved_fraction"] = float(ca_r["fraction"])
        out["composed_resolved_partial_frame"] = ao_r["cd4_partial"]["residual_frame"]
        out["composed_resolved_preserve"] = ca_r["post_add_preserve"]
        keep("composed_resolved", e_r)
        # ------------------------------------------------------------ KA-E: 246's E construction (fixed fraction)
        o_f, ca_f, ua_f, ao_f = builder_arm(partB, fixed=True)
        e_f = score_final(o_f)
        ka("KA-E composed at the FIXED stack fraction == 246 E_cd4_svd5+partial@0.34", e_f, ref_b["E_cd4_svd5+partial@0.34"])
        chk_common("KA-E", ca_f, ua_f, ao_f, st["frac_cd4"])
        ref_f = L.compose(st, [(d246, st["frac_cd4"]), (st["union"], st["frac_x"]), (part246, fr246)])
        out["composed_fixed_vs_246_compose_max_abs_diff"] = chk_effects("KA-E", o_f, ref_f)
        report("COMPOSED (fixed -- 246's E row)", e_f)
        out["composed_fixed"] = public(e_f)
        out["composed_fixed_fraction"] = float(ca_f["fraction"])
        keep("composed_fixed", e_f)

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = peak_gib()
    (WORK / f"repro254-{mode}-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / f"pertarget254-{mode}-{held}.npz", **per_target)
    print(f"[{held}/{mode}] {len(out['failures'])} failures; {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; "
          f"wrote {WORK/f'repro254-{mode}-{held}.json'}", flush=True)
    if out["failures"]:
        raise SystemExit(1)


# --------------------------------------------------------------------------- summary
def retention(pt, arm, base, folds, k=3):
    """analyse246.concentration's leave-top-k-out statistic, verbatim in form: per fold the mean
    per-target delta with the top-k targets (on the delta's own ordering) deleted, averaged over
    folds, divided by the averaged full mean.  Cosine and pds each on their own ordering."""
    rp, rc = [], []
    for f in folds:
        dp = pt[(f, f"{arm}|pds")] - pt[(f, f"{base}|pds")]
        dc = pt[(f, f"{arm}|cos")] - pt[(f, f"{base}|cos")]
        op = np.argsort(-dp)[:k]
        oc = np.argsort(-dc)[:k]
        keep_p = np.ones(dp.size, bool); keep_p[op] = False
        keep_c = np.ones(dc.size, bool); keep_c[oc] = False
        rp.append((dp[keep_p].mean(), dp.mean()))
        rc.append((dc[keep_c].mean(), dc.mean()))
    full_p = float(np.mean([t for _, t in rp])); loo_p = float(np.mean([l for l, _ in rp]))
    full_c = float(np.mean([t for _, t in rc])); loo_c = float(np.mean([l for l, _ in rc]))
    return dict(pds_full=full_p, pds_loo=loo_p, pds_retained=(loo_p / full_p) if full_p else float("nan"),
                cos_full=full_c, cos_loo=loo_c, cos_retained=(loo_c / full_c) if full_c else float("nan"))


def summary() -> None:
    per = {(m, h): json.loads((WORK / f"repro254-{m}-{h}.json").read_text())
           for m in ("partial", "composed") for h in FOLDS}
    pt = {}
    for (m, h), _ in per.items():
        z = np.load(WORK / f"pertarget254-{m}-{h}.npz")
        for key in z.files:
            pt[(h, key)] = z[key]
    ref_b = {h: {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{h}.json").read_text())["rows"]} for h in FOLDS}
    ref_a = {h: {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{h}.json").read_text())["rows"]} for h in FOLDS}
    an = json.loads((E246 / "work/analyse246.json").read_text())
    # the statistic's own known answer: 246's per-target arrays through `retention` reproduce analyse246.json
    pt246 = {}
    for h in FOLDS:
        for tag in ("", "b"):
            z = np.load(E246 / f"work/pertarget246{tag}-{h}.npz")
            for key in z.files:
                pt246[(h, key)] = z[key]

    def agg(mode, key, folds, base_key="stack"):
        rm, ra, dp, dq = [], [], [], []
        for h in folds:
            s = per[(mode, h)][base_key]
            a = per[(mode, h)][key]
            rm.append(a["a23_cos"] / s["a23_cos"] - 1)
            ra.append(a["cos_agg"] / s["cos_agg"] - 1)
            dp.append(a["pds"] - s["pds"])
            dq.append(a["a15_sgn_w"] - s["a15_sgn_w"])
        return dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)),
                    d_qw=float(np.mean(dq)), per_fold_d_pds=[float(x) for x in dp],
                    per_fold_rel_macro=[float(x) for x in rm], pds_positive_folds=int(sum(x > 0 for x in dp)))

    def agg_ref(arm, folds):
        rm, ra, dp = [], [], []
        for h in folds:
            s, a = ref_a[h]["STACK"], ref_b[h][arm]
            rm.append(a["a23_cos"] / s["a23_cos"] - 1)
            ra.append(a["cos_agg"] / s["cos_agg"] - 1)
            dp.append(a["pds"] - s["pds"])
        return dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)))

    all_ka = all(v["ok"] for p in per.values() for v in p["ka"].values()) and not any(p["failures"] for p in per.values())
    n_ka = sum(len(p["ka"]) for p in per.values())
    lines = []
    P = lambda s: (print(s), lines.append(s))  # noqa: E731

    partial = agg("partial", "partial", CANCER)
    ref_c2 = agg_ref("C2_partial@0.34", CANCER)
    svd5 = agg("composed", "svd5", CANCER)
    ref_svd5 = agg_ref("B2_cd4_svd5_resolved", CANCER)
    fixed = agg("composed", "composed_fixed", CANCER)
    ref_e = agg_ref("E_cd4_svd5+partial@0.34", CANCER)
    resolved = agg("composed", "composed_resolved", CANCER)
    incr = agg("composed", "composed_resolved", CANCER, base_key="svd5")
    conc_partial = retention(pt, "partial", "STACK", CANCER)
    conc_fixed = retention(pt, "composed_fixed", "STACK", CANCER)
    conc_resolved = retention(pt, "composed_resolved", "STACK", CANCER)
    conc_svd5 = retention(pt, "svd5", "STACK", CANCER)
    conc_incr = retention(pt, "composed_resolved", "svd5", CANCER)
    ka_stat_e = retention(pt246, "E_cd4_svd5+partial@0.34", "STACK", CANCER)
    ka_stat_c = retention(pt246, "C2_partial@0.34", "STACK", CANCER)
    d_stat = max(abs(ka_stat_e["pds_retained"] - an["concentration"]["E_cd4_svd5+partial@0.34"]["3"]["pds_retained"]),
                 abs(ka_stat_c["pds_retained"] - an["concentration"]["C2_partial@0.34"]["3"]["pds_retained"]))
    kolf_partial = agg("partial", "partial", ("kolf2",))
    kolf_resolved = agg("composed", "composed_resolved", ("kolf2",))
    kolf_fixed = agg("composed", "composed_fixed", ("kolf2",))
    kolf_incr = agg("composed", "composed_resolved", ("kolf2",), base_key="svd5")

    P("=" * 104)
    P(f"per-fold known-answer rows: {'ALL PASS' if all_ka else 'FAILURES'} ({n_ka} rows over {len(per)} fold x mode cells)")
    P(f"leave-top-3-out statistic's own known answer on 246's arrays: E {ka_stat_e['pds_retained']:.6f} "
      f"(analyse246 {an['concentration']['E_cd4_svd5+partial@0.34']['3']['pds_retained']:.6f}), C2 "
      f"{ka_stat_c['pds_retained']:.6f} (analyse246 {an['concentration']['C2_partial@0.34']['3']['pds_retained']:.6f}); "
      f"max|d| {d_stat:.1e}")
    P("=" * 104)
    P("KNOWN-ANSWER ROW 2 -- VCC_CD4_PARTIAL_G=0.34 alone reproduces 246's C2_partial@0.34 (3 cancer folds)")
    P("=" * 104)
    P(f"  {'':30s}{'this run (builder)':>22s}{'246 (json)':>14s}{'|d|':>10s}{'mandate':>10s}")
    for label, mine, want, m in (("rel MACRO cosine", partial["rel_macro"], ref_c2["rel_macro"], +0.00587),
                                 ("rel AGGREGATE cosine", partial["rel_agg"], ref_c2["rel_agg"], float("nan")),
                                 ("d pds", partial["d_pds"], ref_c2["d_pds"], 0.00178)):
        P(f"  {label:30s}{mine:>+22.9f}{want:>+14.9f}{abs(mine-want):>10.1e}{m:>+10.5f}")
    P(f"  rounded to the mandate's digits: macro {100*partial['rel_macro']:+.3f}%  d pds {partial['d_pds']:+.5f}   "
      f"(mandate: +0.587% / +0.00178); pds positive on {partial['pds_positive_folds']}/3 cancer folds "
      f"{['%+.5f' % x for x in partial['per_fold_d_pds']]}")
    P(f"  leave-top-3-out on pds: {100*conc_partial['pds_retained']:.1f}% retained (246: "
      f"{100*an['concentration']['C2_partial@0.34']['3']['pds_retained']:.1f}%); kolf2 apart: macro "
      f"{100*kolf_partial['rel_macro']:+.3f}% d pds {kolf_partial['d_pds']:+.5f}")
    P(f"  partial rows per fold: " + ", ".join(f"{h} {per[('partial', h)]['partial_info']['n_partial']} "
      f"({per[('partial', h)]['partial_info']['n_one']}+{per[('partial', h)]['partial_info']['n_two']})" for h in FOLDS)
      + f"; fractions {[round(per[('partial', h)]['partial_frame']['solved_fraction'], 6) for h in FOLDS]}")
    row2 = (all_ka and abs(partial["rel_macro"] - ref_c2["rel_macro"]) < TOL
            and abs(partial["rel_agg"] - ref_c2["rel_agg"]) < TOL and abs(partial["d_pds"] - ref_c2["d_pds"]) < TOL)
    P("")
    P("=" * 104)
    P("KNOWN-ANSWER ROW 3 -- VCC_CD4_DEFLATE_K=5 + VCC_CD4_PARTIAL_G=0.34 and 246's E_cd4_svd5+partial@0.34")
    P("=" * 104)
    P("  (a) 246's E construction = the deflated CD4 residual at the STACK's FIXED fraction + union + partial,")
    P("      through the builder with `solve_cd4_dose` pinned to the un-deflated residual (the fixed fraction):")
    P(f"  {'':30s}{'this run (builder)':>22s}{'246 (json)':>14s}{'|d|':>10s}{'mandate':>10s}")
    for label, mine, want, m in (("rel MACRO cosine", fixed["rel_macro"], ref_e["rel_macro"], -0.00699),
                                 ("rel AGGREGATE cosine", fixed["rel_agg"], ref_e["rel_agg"], float("nan")),
                                 ("d pds", fixed["d_pds"], ref_e["d_pds"], 0.01405)):
        P(f"  {label:30s}{mine:>+22.9f}{want:>+14.9f}{abs(mine-want):>10.1e}{m:>+10.5f}")
    P(f"  leave-top-3-out on pds: {100*conc_fixed['pds_retained']:.2f}% retained (246: "
      f"{100*an['concentration']['E_cd4_svd5+partial@0.34']['3']['pds_retained']:.2f}%, mandate 61.5%)")
    P(f"  rounded: macro {100*fixed['rel_macro']:+.3f}%  d pds {fixed['d_pds']:+.5f}  LT3O {100*conc_fixed['pds_retained']:.1f}%; "
      f"pds positive {fixed['pds_positive_folds']}/3 {['%+.5f' % x for x in fixed['per_fold_d_pds']]}; kolf2 apart macro "
      f"{100*kolf_fixed['rel_macro']:+.3f}% d pds {kolf_fixed['d_pds']:+.5f}")
    row3 = (all_ka and abs(fixed["rel_macro"] - ref_e["rel_macro"]) < TOL and abs(fixed["rel_agg"] - ref_e["rel_agg"]) < TOL
            and abs(fixed["d_pds"] - ref_e["d_pds"]) < TOL
            and abs(conc_fixed["pds_retained"] - an["concentration"]["E_cd4_svd5+partial@0.34"]["3"]["pds_retained"]) < 1e-6)
    P("")
    P("  (b) the builder's REAL path with both flags = the shipped cd4svd5 construction (CD4's fraction RE-SOLVED on")
    P("      the deflated object, 249's resolved variant) + union + partial -- the object the package carries:")
    P(f"  {'':30s}{'composed resolved':>22s}{'svd5 resolved':>14s}{'246 svd5_res':>14s}")
    for label, mine, want, w246 in (("rel MACRO cosine", resolved["rel_macro"], svd5["rel_macro"], ref_svd5["rel_macro"]),
                                    ("rel AGGREGATE cosine", resolved["rel_agg"], svd5["rel_agg"], ref_svd5["rel_agg"]),
                                    ("d pds", resolved["d_pds"], svd5["d_pds"], ref_svd5["d_pds"])):
        P(f"  {label:30s}{mine:>+22.9f}{want:>+14.9f}{w246:>+14.9f}")
    P(f"  composed resolved: macro {100*resolved['rel_macro']:+.3f}%  agg {100*resolved['rel_agg']:+.3f}%  d pds "
      f"{resolved['d_pds']:+.5f} (positive {resolved['pds_positive_folds']}/3 {['%+.5f' % x for x in resolved['per_fold_d_pds']]})  "
      f"d q_w {resolved['d_qw']:+.5f}  LT3O {100*conc_resolved['pds_retained']:.1f}%; kolf2 apart macro "
      f"{100*kolf_resolved['rel_macro']:+.3f}% d pds {kolf_resolved['d_pds']:+.5f}")
    P(f"  INCREMENT over cd4svd5 (the slot's one variable): macro {100*incr['rel_macro']:+.3f}% "
      f"{['%+.3f%%' % (100*x) for x in incr['per_fold_rel_macro']]}, agg {100*incr['rel_agg']:+.3f}%, d pds {incr['d_pds']:+.5f} "
      f"(positive {incr['pds_positive_folds']}/3 {['%+.5f' % x for x in incr['per_fold_d_pds']]}), d q_w {incr['d_qw']:+.5f}, "
      f"LT3O of the increment {100*conc_incr['pds_retained']:.1f}%; kolf2 apart macro {100*kolf_incr['rel_macro']:+.3f}% "
      f"d pds {kolf_incr['d_pds']:+.5f}")
    P(f"  svd5 alone through this builder vs 246 B2_cd4_svd5_resolved: |d| macro {abs(svd5['rel_macro']-ref_svd5['rel_macro']):.1e} "
      f"d pds {abs(svd5['d_pds']-ref_svd5['d_pds']):.1e}; LT3O svd5 {100*conc_svd5['pds_retained']:.1f}% (246 "
      f"{100*an['concentration']['B2_cd4_svd5_resolved']['3']['pds_retained']:.1f}%)")
    P(f"  fractions (PLOLO frame): stack {np.mean([per[('composed', h)]['frac_cd4_stack'] for h in CANCER]):.4f} -> svd5 "
      f"{np.mean([per[('composed', h)]['svd5_fraction'] for h in CANCER]):.4f}; partial "
      f"{np.mean([per[('composed', h)]['composed_resolved_partial_frame']['solved_fraction'] for h in CANCER]):.4f}")
    P(f"\nROW 2 {'PASS' if row2 else 'FAIL'}   ROW 3 (246's construction reproduced) {'PASS' if row3 else 'FAIL'}   "
      f"statistic KA {'PASS' if d_stat < 1e-9 else 'FAIL'}")
    (WORK / "repro254-summary.json").write_text(json.dumps(
        dict(all_ka=bool(all_ka), n_ka=n_ka, row2_pass=bool(row2), row3_pass=bool(row3), statistic_ka_max_abs_diff=d_stat,
             partial=partial, ref_c2=ref_c2, svd5=svd5, ref_svd5=ref_svd5, composed_fixed=fixed, ref_e=ref_e,
             composed_resolved=resolved, increment_over_svd5=incr,
             concentration=dict(partial=conc_partial, composed_fixed=conc_fixed, composed_resolved=conc_resolved,
                                svd5=conc_svd5, increment=conc_incr, ka_e_246=ka_stat_e, ka_c2_246=ka_stat_c),
             kolf2=dict(partial=kolf_partial, composed_resolved=kolf_resolved, composed_fixed=kolf_fixed, increment=kolf_incr),
             per_fold={f"{m}-{h}": {k: v for k, v in p.items() if k in ("stack", "partial", "svd5", "composed_resolved",
                                                                           "composed_fixed", "partial_fraction_246",
                                                                           "svd5_fraction", "composed_resolved_fraction",
                                                                           "frac_cd4_stack", "frac_x", "peak_rss_gib", "seconds")}
                       for (m, h), p in per.items()}),
        indent=1, default=float))
    (WORK / "repro254-summary.log").write_text("\n".join(lines) + "\n")
    if not (row2 and row3 and d_stat < 1e-9):
        raise SystemExit(1)


if __name__ == "__main__":
    if ARGS.summary:
        summary()
    else:
        run_fold(ARGS.fold, ARGS.mode)
