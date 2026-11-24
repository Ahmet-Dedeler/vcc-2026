"""255 KNOWN-ANSWER ROW 1 -- the builder's k = 10 deflation through 162's PLOLO instrument.

249's `repro249.py` (222's `repro222.py` pattern) with k = 10: 246's rows re-derived with the
PRODUCTION builder's code substituted in -- `deflate_cd4_residual` for `lib246.deflate`, and the
WHOLE of `apply_frozen_cd4_residual` (deflate -> size -> merge -> one preserve) for 246's
`compose`.  The env the builder reads at import is the arm's: VCC_ALPHA=0.982
VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=10.

    KA-A   STACK (lib246)                        == 246 sweep246 `STACK`              (cos_agg, a23_cos, pds; 1e-9)
    KA-B   builder deflation k=10                == lib246.deflate(cd4, "svd", 10), bit for bit, + its audit
    KA-C   ARM (k=10) through the builder's entry == 246 `B2_cd4_svd10_resolved`     (1e-9), and its solved fraction
    KA-C'  the k = 5 arm through the SAME builder (module global switched) == 246 `B2_cd4_svd5_resolved`
           and 249's own per-fold numbers -- the INCREMENT k 5 -> 10 is read off these two
    KA-D   energy-matched random basis (246's seed and bisection, matched to svd10's energy) through
           the same builder path: subspace bit-identical to lib246.deflate(cd4, "rand", k', seed); fixed
           sizing == 246 `B2_cd4_randE10` (1e-9); re-solved sizing reported (d pds <= +0.001 expected)

One fold per process (memory; 249 measured 4.9-5.2 GiB); `--summary` aggregates.

    uv run python repro255.py --fold k562_gwps      (orion_hct116, orion_hek293t, kolf2)
    uv run python repro255.py --summary
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
os.environ["VCC_CD4_DEFLATE_K"] = "10"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_CD4_DEFLATE_BASIS",
           "VCC_CD4_PARTIAL_G", "VCC_UNION2_SOURCES", "VCC_UNION2_G"):
    os.environ.pop(_k, None)

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E246 = ROOT / "experiments/246-direction-hypotheses"
E249 = ROOT / "experiments/249-cd4-deflate"
sys.path.insert(0, str(E246))
import lib246 as L  # noqa: E402  (imports plolo / resid218 / compose218; 246's construction)
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"          # this experiment's copies of the fold caches (sha-identical to 246's)
L.WORK = WORK

K = 10
K_INC = 5                            # the incumbent's k (cd4svd5)
TOL = 1e-9
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KEYS = ("cos_agg", "a23_cos", "pds")


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("repro255_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repro255_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def public(e: dict) -> dict:
    return {k: v for k, v in e.items() if not k.startswith("_")}


def lt_out(delta: np.ndarray, n: int) -> float:
    """leave-top-n-out retention of a per-target delta: the sum over all but the n largest
    contributions, divided by the full sum (246 section 3.4's statistic, on per-target `pds`)."""
    total = float(delta.sum())
    if total == 0.0:
        return float("nan")
    order = np.argsort(-delta)
    return float(delta[order[n:]].sum() / total)


def run_fold(held: str) -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.CD4_DEFLATE_K == K and MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    assert getattr(MOD, "CD4_PARTIAL_G", None) is None and getattr(MOD, "CD4_DEFLATE_BASIS", None) is None
    assert getattr(MOD, "UNION2_G", None) is None and getattr(MOD, "CD4_SIZING", None) is None
    print(f"builder {builder_path}  sha256 {MOD.file_sha256(builder_path)}", flush=True)
    ref_a = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    ref_b = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}
    ref249 = json.loads((E249 / f"work/repro249-{held}.json").read_text())

    st = L.stage(held)
    fold, targets = st["fold"], st["targets"]
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    pp = np.asarray(fold["panel_positions"], dtype=np.int64)
    print(f"[{held}] staged {len(targets)} targets; stack frac_cd4 {st['frac_cd4']:.6f} "
          f"frac_x {st['frac_x']:.6f}  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    out: dict = dict(fold=held, n_targets=len(targets), builder_sha256=MOD.file_sha256(builder_path),
                     frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"], ka={}, failures=[])

    def fail(msg):
        print(f"  *** {msg}", flush=True)
        out["failures"].append(msg)

    def ka(name, e, ref):
        d = max(abs(ref[k] - e[k]) for k in KEYS)
        ok = d < TOL
        out["ka"][name] = dict(max_abs_diff=float(d), ok=bool(ok), ref_arm=ref["arm"])
        print(f"  {'PASS' if ok else 'FAIL'} {name:56s} vs 246 {ref['arm']:24s} max|d| {d:.2e}", flush=True)
        if not ok:
            fail(f"{name}: differs from 246 {ref['arm']} by {d:.3e}")

    def score_final(final):
        final = np.asarray(final, dtype=np.float64).copy()
        final[:, pp] = 0.0
        return L.score(st, final)

    # ---------------------------------------------------------------- KA-A: the stack
    e_stack = L.score(st, L.compose(st, L.stack_parts(st)))
    ka("KA-A STACK (lib246)", e_stack, ref_a["STACK"])
    out["stack"] = public(e_stack)

    # ---------------------------------------------------------------- KA-B: the deflation itself, k = 10
    d246, a246 = L.deflate(st["cd4"], "svd", K)
    dB, aB = MOD.deflate_cd4_residual(st["cd4"], K)
    same = (np.array_equal(np.asarray(dB.mean_effect), np.asarray(d246.mean_effect))
            and np.array_equal(np.asarray(dB.gene_votes), np.asarray(d246.gene_votes))
            and np.array_equal(np.asarray(dB.target_condition_count), np.asarray(d246.target_condition_count)))
    same_audit = (aB["energy_removed_frac"] == a246["energy_removed_frac"]
                  and aB["top_k_energy_frac_of_centred"] == a246["top_k_energy_frac_of_centred"]
                  and aB["n_carried"] == a246["n_carried"])
    ref_aud = ref_b["B2_cd4_svd10_resolved"]["audit"]
    same_ref = (aB["energy_removed_frac"] == ref_aud["energy_removed_frac"]
                and aB["n_carried"] == ref_aud["n_carried"] and aB["k"] == ref_aud["k"] == K)
    out["ka"]["KA-B deflation k=10 bit-identical to lib246"] = dict(
        ok=bool(same and same_audit and same_ref),
        max_abs_diff=float(np.abs(np.asarray(dB.mean_effect) - np.asarray(d246.mean_effect)).max()))
    print(f"  {'PASS' if same and same_audit and same_ref else 'FAIL'} KA-B builder deflation k=10 == lib246.deflate: arrays "
          f"{'bit-identical' if same else 'DIFFER'}, audit {'identical' if same_audit else 'DIFFERS'} "
          f"(n_carried {aB['n_carried']}, energy removed {aB['energy_removed_frac']:.9f}, top-k share "
          f"{aB['top_k_energy_frac_of_centred']:.9f}; 246 json {ref_aud['energy_removed_frac']:.9f})", flush=True)
    if not (same and same_audit and same_ref):
        fail("KA-B: the builder's k=10 deflation is not 246's")
    out["deflation_audit"] = aB

    # ---------------------------------------------------------------- KA-C: the arm, builder path (module K = 10)
    arm_out, cd4_audit, union_audit = MOD.apply_frozen_cd4_residual(
        st["whitened"], st["cd4"], wcols, pp, extra=st["union"], extra_audit={})
    assert cd4_audit["deflation"]["k"] == K
    e_arm = score_final(arm_out)
    ka("KA-C ARM k=10 through apply_frozen_cd4_residual", e_arm, ref_b["B2_cd4_svd10_resolved"])
    fr_ref = float(ref_aud["fraction"])
    fr_arm = float(cd4_audit["fraction"])
    d_fr = abs(fr_arm - fr_ref)
    print(f"  {'PASS' if d_fr < 1e-12 else 'FAIL'} KA-C solved fraction {fr_arm:.12f} vs 246 {fr_ref:.12f} "
          f"(|d| {d_fr:.1e}); stack {st['frac_cd4']:.6f} -> {fr_arm:.6f} (x{fr_arm/st['frac_cd4']:.4f})", flush=True)
    if d_fr >= 1e-12:
        fail("KA-C: the re-solved fraction differs from 246's")
    fr246, _ = L.solve_fraction(st["whitened"], d246, L.CD4_F)
    ref_final = L.compose(st, [(d246, fr246), (st["union"], st["frac_x"])])
    chk = np.asarray(arm_out, dtype=np.float64).copy()
    chk[:, pp] = 0.0
    d_eff = float(np.abs(chk - ref_final).max())
    print(f"  effect matrix: builder path vs 246 compose max|d| {d_eff:.2e}", flush=True)
    d_union = abs(float(union_audit["residual_frame"]["solved_fraction"]) - st["frac_x"])
    print(f"  {'PASS' if d_union < 1e-12 else 'FAIL'} the union's solved fraction through the builder "
          f"{union_audit['residual_frame']['solved_fraction']:.12f} == stack's {st['frac_x']:.12f} (|d| {d_union:.1e})", flush=True)
    if d_union >= 1e-12:
        fail("the union's fraction moved under the CD4 flag")
    out["arm"] = public(e_arm)
    out["arm_fraction"] = fr_arm
    out["arm_frame"] = cd4_audit["residual_frame"]
    out["arm_preserve"] = cd4_audit["post_add_preserve"]
    out["arm_deflation"] = cd4_audit["deflation"]
    out["arm_vs_246_compose_max_abs_diff"] = d_eff
    out["union_frame_solved_fraction"] = float(union_audit["residual_frame"]["solved_fraction"])
    print(f"  ARM k=10: macro {100*(e_arm['a23_cos']/e_stack['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_arm['cos_agg']/e_stack['cos_agg']-1):+.4f}%  d pds {e_arm['pds']-e_stack['pds']:+.6f}  "
          f"d q_w {e_arm['a15_sgn_w']-e_stack['a15_sgn_w']:+.6f}", flush=True)

    # ---------------------------------------------------------------- KA-C': the k = 5 arm (cd4svd5) through the SAME builder
    MOD.CD4_DEFLATE_K = K_INC
    try:
        arm5_out, cd4_audit5, union_audit5 = MOD.apply_frozen_cd4_residual(
            st["whitened"], st["cd4"], wcols, pp, extra=st["union"], extra_audit={})
    finally:
        MOD.CD4_DEFLATE_K = K
    assert cd4_audit5["deflation"]["k"] == K_INC
    e_arm5 = score_final(arm5_out)
    ka("KA-C' k=5 (cd4svd5) through the same builder", e_arm5, ref_b["B2_cd4_svd5_resolved"])
    d249 = max(abs(ref249["arm"][k] - e_arm5[k]) for k in KEYS)
    ok249 = d249 < 1e-12 and abs(ref249["arm_fraction"] - float(cd4_audit5["fraction"])) < 1e-12
    out["ka"]["KA-C' k=5 == 249's repro249 row"] = dict(max_abs_diff=float(d249), ok=bool(ok249))
    print(f"  {'PASS' if ok249 else 'FAIL'} KA-C' k=5 through this builder == 249's own repro249-{held} arm row "
          f"(max|d| {d249:.1e}; fraction {float(cd4_audit5['fraction']):.12f} vs 249 {ref249['arm_fraction']:.12f})", flush=True)
    if not ok249:
        fail("KA-C': the k=5 arm through the current builder is not 249's")
    out["arm5"] = public(e_arm5)
    out["arm5_fraction"] = float(cd4_audit5["fraction"])
    out["arm5_deflation"] = cd4_audit5["deflation"]
    print(f"  ARM k=5:  macro {100*(e_arm5['a23_cos']/e_stack['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_arm5['cos_agg']/e_stack['cos_agg']-1):+.4f}%  d pds {e_arm5['pds']-e_stack['pds']:+.6f}", flush=True)
    print(f"  INCREMENT k 5->10: macro {100*(e_arm['a23_cos']/e_arm5['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_arm['cos_agg']/e_arm5['cos_agg']-1):+.4f}%  d pds {e_arm['pds']-e_arm5['pds']:+.6f}  "
          f"d q_w {e_arm['a15_sgn_w']-e_arm5['a15_sgn_w']:+.6f}  fraction {fr_arm:.6f} vs {float(cd4_audit5['fraction']):.6f} "
          f"(x{fr_arm/float(cd4_audit5['fraction']):.4f})", flush=True)

    # per-target `pds` for the concentration statistics (246 section 3.4's leave-top-n-out)
    p_stack = np.asarray(e_stack["_pds_per_target_a23"], dtype=np.float64)
    p10 = np.asarray(e_arm["_pds_per_target_a23"], dtype=np.float64)
    p5 = np.asarray(e_arm5["_pds_per_target_a23"], dtype=np.float64)
    mean_chk = abs(float(p10.mean()) - float(e_arm["pds"]))
    out["concentration"] = dict(
        per_target_mean_equals_pds=bool(mean_chk < 1e-9), per_target_mean_abs_diff=float(mean_chk),
        arm10_vs_stack=dict(lt1o=lt_out(p10 - p_stack, 1), lt3o=lt_out(p10 - p_stack, 3)),
        arm5_vs_stack=dict(lt1o=lt_out(p5 - p_stack, 1), lt3o=lt_out(p5 - p_stack, 3)),
        increment=dict(lt1o=lt_out(p10 - p5, 1), lt3o=lt_out(p10 - p5, 3),
                       n_targets_up=int(np.sum(p10 - p5 > 0)), n_targets_down=int(np.sum(p10 - p5 < 0))),
    )
    c = out["concentration"]
    print(f"  concentration (LT3O of per-target d pds): k=10 vs stack {100*c['arm10_vs_stack']['lt3o']:.1f}%  "
          f"k=5 vs stack {100*c['arm5_vs_stack']['lt3o']:.1f}%  increment {100*c['increment']['lt3o']:+.1f}% "
          f"(targets up {c['increment']['n_targets_up']} / down {c['increment']['n_targets_down']}; "
          f"per-target mean == pds: {c['per_target_mean_equals_pds']})", flush=True)
    np.savez_compressed(WORK / f"pertarget255-{held}.npz", stack=p_stack, arm10=p10, arm5=p5,
                        cos_stack=np.asarray(e_stack["_cos_per_target_a23"]),
                        cos_arm10=np.asarray(e_arm["_cos_per_target_a23"]),
                        cos_arm5=np.asarray(e_arm5["_cos_per_target_a23"]))

    # ---------------------------------------------------------------- KA-D: energy-matched random control (matched to svd10)
    target = aB["energy_removed_frac"]
    seed = L.RAND_SEED + 7 * 0 + K                      # sweep246b: seed_i = 0, k = 10
    rows = np.flatnonzero(np.asarray(st["cd4"].target_condition_count) > 0)
    lo, hi, best = 1, 260, None
    n_eval = 0
    while lo <= hi:                                     # sweep246b.rand_energy_matched, verbatim
        mid = (lo + hi) // 2
        rng = np.random.default_rng(seed + 1000 * mid)  # lib246.deflate's 'rand' rng
        Uk, _ = np.linalg.qr(rng.standard_normal((rows.size, mid)))
        dR, aR = MOD.deflate_cd4_residual(st["cd4"], mid, basis=Uk)
        n_eval += 1
        if aR["energy_removed_frac"] >= target:
            best = (dR, aR, mid, Uk)
            hi = mid - 1
        else:
            lo = mid + 1
    assert best is not None
    dR, aR, kk, Uk = best
    d246r, a246r = L.deflate(st["cd4"], "rand", kk, seed=seed)
    same_r = np.array_equal(np.asarray(dR.mean_effect), np.asarray(d246r.mean_effect))
    ref_r = ref_b["B2_cd4_randE10"]
    print(f"  {'PASS' if same_r and kk == ref_r['audit']['k'] else 'FAIL'} KA-D random subspace: k'={kk} "
          f"(246: {ref_r['audit']['k']}), energy {aR['energy_removed_frac']:.6f} (246 {ref_r['audit']['energy_removed_frac']:.6f}); "
          f"builder-path array {'bit-identical' if same_r else 'DIFFERS'} to lib246.deflate('rand'); "
          f"{n_eval} bisection evaluations", flush=True)
    if not (same_r and kk == ref_r["audit"]["k"]):
        fail("KA-D: the random control subspace is not 246's")
    e_ctl_fixed = L.score(st, L.compose(st, [(dR, st["frac_cd4"]), (st["union"], st["frac_x"])]))
    ka("KA-D control, fixed sizing (246's row)", e_ctl_fixed, ref_r)
    original = MOD.deflate_cd4_residual

    def control_deflation(cd4, k, basis=None):        # the SAME builder path, the control basis
        return original(cd4, kk, basis=Uk)

    MOD.deflate_cd4_residual = control_deflation
    try:
        ctl_out, ctl_audit, _ = MOD.apply_frozen_cd4_residual(
            st["whitened"], st["cd4"], wcols, pp, extra=st["union"], extra_audit={})
    finally:
        MOD.deflate_cd4_residual = original
    e_ctl_res = score_final(ctl_out)
    print(f"  CONTROL fixed:    macro {100*(e_ctl_fixed['a23_cos']/e_stack['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_ctl_fixed['cos_agg']/e_stack['cos_agg']-1):+.4f}%  d pds {e_ctl_fixed['pds']-e_stack['pds']:+.6f}", flush=True)
    print(f"  CONTROL resolved: macro {100*(e_ctl_res['a23_cos']/e_stack['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_ctl_res['cos_agg']/e_stack['cos_agg']-1):+.4f}%  d pds {e_ctl_res['pds']-e_stack['pds']:+.6f}  "
          f"fraction {ctl_audit['fraction']:.6f}", flush=True)
    out["control"] = dict(k_prime=int(kk), seed=int(seed), energy_removed_frac=aR["energy_removed_frac"],
                          fixed=public(e_ctl_fixed), resolved=public(e_ctl_res),
                          resolved_fraction=float(ctl_audit["fraction"]),
                          deflation_audit=ctl_audit["deflation"])

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = peak_gib()
    (WORK / f"repro255-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"[{held}] {len(out['failures'])} failures; {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; "
          f"wrote {WORK/f'repro255-{held}.json'}", flush=True)
    if out["failures"]:
        raise SystemExit(1)


def summary() -> None:
    per = {h: json.loads((WORK / f"repro255-{h}.json").read_text()) for h in CANCER + ("kolf2",)}
    ref_b = {h: {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{h}.json").read_text())["rows"]}
             for h in per}
    ref_a = {h: {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{h}.json").read_text())["rows"]}
             for h in per}
    p246 = json.loads((E246 / "work/price246.json").read_text())
    s249 = json.loads((E249 / "work/repro249-summary.json").read_text())

    def agg(folds, arm_key, base_key="stack", ctl=None):
        rm, ra, dp, dq = [], [], [], []
        for h in folds:
            s = per[h][base_key]
            a = per[h][arm_key] if ctl is None else per[h]["control"][ctl]
            rm.append(a["a23_cos"] / s["a23_cos"] - 1)
            ra.append(a["cos_agg"] / s["cos_agg"] - 1)
            dp.append(a["pds"] - s["pds"])
            dq.append(a["a15_sgn_w"] - s["a15_sgn_w"])
        return dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)),
                    d_qw=float(np.mean(dq)), per_fold_d_pds=[float(x) for x in dp],
                    per_fold_rel_macro=[float(x) for x in rm], per_fold_rel_agg=[float(x) for x in ra],
                    pds_positive_folds=int(sum(x > 0 for x in dp)))

    def agg_ref(folds, arm, base="STACK"):
        rm, ra, dp = [], [], []
        for h in folds:
            s = ref_a[h][base] if base == "STACK" else ref_b[h][base]
            a = ref_b[h][arm]
            rm.append(a["a23_cos"] / s["a23_cos"] - 1)
            ra.append(a["cos_agg"] / s["cos_agg"] - 1)
            dp.append(a["pds"] - s["pds"])
        return dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)))

    arm = agg(CANCER, "arm")
    ref = agg_ref(CANCER, "B2_cd4_svd10_resolved")
    arm5 = agg(CANCER, "arm5")
    ref5 = agg_ref(CANCER, "B2_cd4_svd5_resolved")
    inc = agg(CANCER, "arm", base_key="arm5")
    ref_inc = agg_ref(CANCER, "B2_cd4_svd10_resolved", base="B2_cd4_svd5_resolved")
    ctl_fixed = agg(CANCER, None, ctl="fixed")
    ctl_res = agg(CANCER, None, ctl="resolved")
    ref_ctl = agg_ref(CANCER, "B2_cd4_randE10")
    kolf = agg(("kolf2",), "arm")
    kolf_inc = agg(("kolf2",), "arm", base_key="arm5")
    e3 = float(np.mean([per[h]["deflation_audit"]["energy_removed_frac"] for h in CANCER]))
    e4 = float(np.mean([per[h]["deflation_audit"]["energy_removed_frac"] for h in per]))
    fr_stack = float(np.mean([per[h]["frac_cd4_stack"] for h in CANCER]))
    fr_arm = float(np.mean([per[h]["arm_fraction"] for h in CANCER]))
    fr_arm5 = float(np.mean([per[h]["arm5_fraction"] for h in CANCER]))
    conc = {key: dict(lt1o=float(np.mean([per[h]["concentration"][key]["lt1o"] for h in CANCER])),
                      lt3o=float(np.mean([per[h]["concentration"][key]["lt3o"] for h in CANCER])),
                      per_fold_lt3o=[per[h]["concentration"][key]["lt3o"] for h in CANCER])
            for key in ("arm10_vs_stack", "arm5_vs_stack", "increment")}
    all_ka = all(v["ok"] for h in per for v in per[h]["ka"].values()) and not any(per[h]["failures"] for h in per)
    lines = []
    P = lambda s: (print(s), lines.append(s))  # noqa: E731
    P("=" * 100)
    P("KNOWN-ANSWER ROW 1 -- the staged k = 10 arm reproduces 246's B2_cd4_svd10_resolved (three cancer folds, kolf2 apart)")
    P("=" * 100)
    P(f"  per-fold KA rows: {'ALL PASS' if all_ka else 'FAILURES'} "
      f"({sum(len(per[h]['ka']) for h in per)} rows over {len(per)} folds)")
    P(f"  {'':34s}{'this run (builder)':>22s}{'246 (json rows)':>16s}{'|d|':>10s}{'price246.json':>15s}{'mandate':>10s}")
    for label, mine, want, p246v, m in (("rel MACRO cosine", arm["rel_macro"], ref["rel_macro"], p246["B2_cd4_svd10_resolved"]["rel_macro"], -0.0167),
                                        ("rel AGGREGATE cosine", arm["rel_agg"], ref["rel_agg"], float("nan"), -0.0154),
                                        ("d pds", arm["d_pds"], ref["d_pds"], p246["B2_cd4_svd10_resolved"]["d_pds"], 0.01281)):
        P(f"  {label:34s}{mine:>+22.9f}{want:>+16.9f}{abs(mine-want):>10.1e}{p246v:>+15.9f}{m:>+10.5f}")
    P(f"  rounded to the mandate's digits: macro {100*arm['rel_macro']:+.2f}%  agg {100*arm['rel_agg']:+.2f}%  "
      f"d pds {arm['d_pds']:+.5f}   (mandate: -1.67% / -1.54% / +0.01281)")
    P(f"  pds positive on {arm['pds_positive_folds']}/3 cancer folds: {['%+.5f' % x for x in arm['per_fold_d_pds']]}; "
      f"per-fold macro {['%+.3f%%' % (100*x) for x in arm['per_fold_rel_macro']]}, agg {['%+.3f%%' % (100*x) for x in arm['per_fold_rel_agg']]}")
    P(f"  energy_removed_frac: 3-cancer-fold mean {e3:.4f}, 4-fold mean {e4:.4f} (mandate ~0.286); per fold "
      f"{[round(per[h]['deflation_audit']['energy_removed_frac'], 5) for h in per]}")
    P(f"  CD4 solved fraction, PLOLO's frame: stack {fr_stack:.4f} -> k=5 {fr_arm5:.4f} -> k=10 {fr_arm:.4f} "
      f"(x{fr_arm/fr_arm5:.4f} over k=5)")
    P(f"  kolf2 apart: macro {100*kolf['rel_macro']:+.3f}%  agg {100*kolf['rel_agg']:+.3f}%  d pds {kolf['d_pds']:+.5f}")
    P(f"  effect matrix, builder path vs 246's compose, max|d| per fold: "
      f"{[per[h]['arm_vs_246_compose_max_abs_diff'] for h in per]}")
    P("")
    P("=" * 100)
    P("KA-C' -- the k = 5 arm (cd4svd5) through the SAME builder, and the INCREMENT k 5 -> 10")
    P("=" * 100)
    P(f"  k=5 through this builder: macro {100*arm5['rel_macro']:+.6f}%  agg {100*arm5['rel_agg']:+.6f}%  d pds {arm5['d_pds']:+.9f}")
    P(f"  246 B2_cd4_svd5_resolved:  macro {100*ref5['rel_macro']:+.6f}%  agg {100*ref5['rel_agg']:+.6f}%  d pds {ref5['d_pds']:+.9f}")
    P(f"  249 repro249-summary:      macro {100*s249['arm']['rel_macro']:+.6f}%  agg {100*s249['arm']['rel_agg']:+.6f}%  d pds {s249['arm']['d_pds']:+.9f}  "
      f"max|d| {max(abs(arm5[k]-s249['arm'][k]) for k in ('rel_macro','rel_agg','d_pds')):.1e}")
    P(f"  INCREMENT (k=10 vs k=5, per-fold ratios then mean): macro {100*inc['rel_macro']:+.4f}%  agg {100*inc['rel_agg']:+.4f}%  "
      f"d pds {inc['d_pds']:+.6f} ({inc['pds_positive_folds']}/3 folds positive: {['%+.5f' % x for x in inc['per_fold_d_pds']]})  d q_w {inc['d_qw']:+.6f}")
    P(f"  the same increment off 246's own json rows: macro {100*ref_inc['rel_macro']:+.4f}%  agg {100*ref_inc['rel_agg']:+.4f}%  d pds {ref_inc['d_pds']:+.6f}")
    P(f"  kolf2 apart, increment: macro {100*kolf_inc['rel_macro']:+.3f}%  agg {100*kolf_inc['rel_agg']:+.3f}%  d pds {kolf_inc['d_pds']:+.5f}")
    P(f"  leave-top-3-out retention of per-target d pds (mean over cancer folds): k=10 vs stack {100*conc['arm10_vs_stack']['lt3o']:.1f}%  "
      f"k=5 vs stack {100*conc['arm5_vs_stack']['lt3o']:.1f}% (246 section 3.4 quotes 56% for the fixed-sizing k=5)  "
      f"increment {100*conc['increment']['lt3o']:+.1f}% (per fold {['%+.1f%%' % (100*x) for x in conc['increment']['per_fold_lt3o']]})")
    P("")
    P("=" * 100)
    P("KA-D -- energy-matched random control (matched to svd10's energy) through the same builder path (one seed)")
    P("=" * 100)
    P(f"  random dims per fold {[per[h]['control']['k_prime'] for h in per]} (246: 73/79/78/72), energy "
      f"{[round(per[h]['control']['energy_removed_frac'], 4) for h in per]}")
    P(f"  fixed sizing (246's row):     macro {100*ctl_fixed['rel_macro']:+.3f}%  agg {100*ctl_fixed['rel_agg']:+.3f}%  "
      f"d pds {ctl_fixed['d_pds']:+.6f}   [246 json: {100*ref_ctl['rel_macro']:+.3f}% / {100*ref_ctl['rel_agg']:+.3f}% / "
      f"{ref_ctl['d_pds']:+.6f}]")
    P(f"  RE-SOLVED through the builder: macro {100*ctl_res['rel_macro']:+.3f}%  agg {100*ctl_res['rel_agg']:+.3f}%  "
      f"d pds {ctl_res['d_pds']:+.6f}   per fold {['%+.5f' % x for x in ctl_res['per_fold_d_pds']]}")
    ctl_ok = ctl_res["d_pds"] <= 0.001 and ctl_fixed["d_pds"] <= 0.001
    P(f"  control {'PASS' if ctl_ok else 'FAIL'}: d pds <= +0.001 (the same ~28% energy in random directions buys "
      f"{ctl_res['d_pds']:+.4f} against the arm's {arm['d_pds']:+.4f})")
    row1 = all_ka and abs(arm["rel_macro"] - ref["rel_macro"]) < TOL and abs(arm["rel_agg"] - ref["rel_agg"]) < TOL \
        and abs(arm["d_pds"] - ref["d_pds"]) < TOL \
        and abs(arm["rel_macro"] - p246["B2_cd4_svd10_resolved"]["rel_macro"]) < TOL \
        and abs(arm["d_pds"] - p246["B2_cd4_svd10_resolved"]["d_pds"]) < TOL
    P(f"\nROW 1 {'PASS' if row1 else 'FAIL'} (k=10 arm == 246 to 1e-9 on every fold and in aggregate, incl. price246.json's digits)   "
      f"CONTROL {'PASS' if ctl_ok else 'FAIL'}")
    (WORK / "repro255-summary.json").write_text(json.dumps(
        dict(arm=arm, ref246=ref, arm5=arm5, ref246_svd5=ref5, ref249_svd5=s249["arm"], increment=inc, ref246_increment=ref_inc,
             control_fixed=ctl_fixed, control_resolved=ctl_res, ref246_control=ref_ctl,
             kolf2=kolf, kolf2_increment=kolf_inc, energy_removed_3fold=e3, energy_removed_4fold=e4,
             fraction_stack=fr_stack, fraction_arm5=fr_arm5, fraction_arm=fr_arm, concentration=conc,
             row1_pass=bool(row1), control_pass=bool(ctl_ok), all_ka=bool(all_ka)),
        indent=1))
    (WORK / "repro255-summary.log").write_text("\n".join(lines) + "\n")
    if not (row1 and ctl_ok):
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
