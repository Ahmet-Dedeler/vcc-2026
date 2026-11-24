"""249 KNOWN-ANSWER ROWS 2 and 3 -- the builder's deflation through 162's PLOLO instrument.

246's rows re-derived with the PRODUCTION builder's code substituted in (222's `repro222.py`
pattern): `deflate_cd4_residual` for `lib246.deflate`, and the WHOLE of
`apply_frozen_cd4_residual` (deflate -> size -> merge -> one preserve) for 246's `compose`.
The env the builder reads at import is the arm's: VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34
VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5.

    KA-A  STACK (lib246)                    == 246 sweep246 `STACK`           (cos_agg, a23_cos, pds; 1e-9)
    KA-B  builder deflation                 == lib246.deflate(cd4, "svd", 5), bit for bit, + its audit
    KA-C  ARM through the builder's entry   == 246 `B2_cd4_svd5_resolved`     (1e-9), and its solved fraction
    KA-D  energy-matched random basis through the SAME builder path, 246's seed and bisection:
          the subspace == lib246.deflate(cd4, "rand", k', seed) bit for bit; fixed sizing ==
          246 `B2_cd4_randE5` (1e-9); re-solved sizing through `apply_frozen_cd4_residual` is
          the mandate's row 3 (mean d pds over the cancer folds <= +0.001)

One fold per process (memory; 246 measured 4.1-5.0 GiB); `--summary` aggregates.

    uv run python repro249.py --fold k562_gwps      (orion_hct116, orion_hek293t, kolf2)
    uv run python repro249.py --summary
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
os.environ["VCC_CD4_DEFLATE_K"] = "5"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES"):
    os.environ.pop(_k, None)

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E246 = ROOT / "experiments/246-direction-hypotheses"
sys.path.insert(0, str(E246))
import lib246 as L  # noqa: E402  (imports plolo / resid218 / compose218; 246's construction)
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"          # this experiment's copies of the fold caches (sha-identical)
L.WORK = WORK

K = 5
TOL = 1e-9
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KEYS = ("cos_agg", "a23_cos", "pds")


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("repro249_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repro249_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def public(e: dict) -> dict:
    return {k: v for k, v in e.items() if not k.startswith("_")}


def run_fold(held: str) -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.CD4_DEFLATE_K == K and MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    print(f"builder {builder_path}  sha256 {MOD.file_sha256(builder_path)}", flush=True)
    ref_a = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    ref_b = {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}

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
        print(f"  {'PASS' if ok else 'FAIL'} {name:52s} vs 246 {ref['arm']:24s} max|d| {d:.2e}", flush=True)
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

    # ---------------------------------------------------------------- KA-B: the deflation itself
    d246, a246 = L.deflate(st["cd4"], "svd", K)
    dB, aB = MOD.deflate_cd4_residual(st["cd4"], K)
    same = (np.array_equal(np.asarray(dB.mean_effect), np.asarray(d246.mean_effect))
            and np.array_equal(np.asarray(dB.gene_votes), np.asarray(d246.gene_votes))
            and np.array_equal(np.asarray(dB.target_condition_count), np.asarray(d246.target_condition_count)))
    same_audit = (aB["energy_removed_frac"] == a246["energy_removed_frac"]
                  and aB["top_k_energy_frac_of_centred"] == a246["top_k_energy_frac_of_centred"]
                  and aB["n_carried"] == a246["n_carried"])
    ref_aud = ref_b["B2_cd4_svd5_resolved"]["audit"]
    same_ref = (aB["energy_removed_frac"] == ref_aud["energy_removed_frac"]
                and aB["n_carried"] == ref_aud["n_carried"])
    out["ka"]["KA-B deflation bit-identical to lib246"] = dict(ok=bool(same and same_audit and same_ref),
                                                              max_abs_diff=float(np.abs(np.asarray(dB.mean_effect) - np.asarray(d246.mean_effect)).max()))
    print(f"  {'PASS' if same and same_audit and same_ref else 'FAIL'} KA-B builder deflation == lib246.deflate: arrays "
          f"{'bit-identical' if same else 'DIFFER'}, audit {'identical' if same_audit else 'DIFFERS'} "
          f"(n_carried {aB['n_carried']}, energy removed {aB['energy_removed_frac']:.9f}, top-k share "
          f"{aB['top_k_energy_frac_of_centred']:.9f}; 246 json {ref_aud['energy_removed_frac']:.9f})", flush=True)
    if not (same and same_audit and same_ref):
        fail("KA-B: the builder's deflation is not 246's")
    out["deflation_audit"] = aB

    # ---------------------------------------------------------------- KA-C: the arm, builder path
    arm_out, cd4_audit, union_audit = MOD.apply_frozen_cd4_residual(
        st["whitened"], st["cd4"], wcols, pp, extra=st["union"], extra_audit={})
    e_arm = score_final(arm_out)
    ka("KA-C ARM through apply_frozen_cd4_residual", e_arm, ref_b["B2_cd4_svd5_resolved"])
    fr_ref = float(ref_aud["fraction"])
    fr_arm = float(cd4_audit["fraction"])
    d_fr = abs(fr_arm - fr_ref)
    print(f"  {'PASS' if d_fr < 1e-12 else 'FAIL'} KA-C solved fraction {fr_arm:.12f} vs 246 {fr_ref:.12f} "
          f"(|d| {d_fr:.1e}); stack {st['frac_cd4']:.6f} -> {fr_arm:.6f} (x{fr_arm/st['frac_cd4']:.4f})", flush=True)
    if d_fr >= 1e-12:
        fail("KA-C: the re-solved fraction differs from 246's")
    # 246's own composition of the same object, for the max|d| on the effect matrix
    fr246, _ = L.solve_fraction(st["whitened"], d246, L.CD4_F)
    ref_final = L.compose(st, [(d246, fr246), (st["union"], st["frac_x"])])
    chk = np.asarray(arm_out, dtype=np.float64).copy()
    chk[:, pp] = 0.0
    d_eff = float(np.abs(chk - ref_final).max())
    print(f"  effect matrix: builder path vs 246 compose max|d| {d_eff:.2e}", flush=True)
    # the union never sees the flag
    d_union = abs(float(union_audit["residual_frame"]["solved_fraction"]) - st["frac_x"])
    print(f"  {'PASS' if d_union < 1e-12 else 'FAIL'} the union's solved fraction through the builder "
          f"{union_audit['residual_frame']['solved_fraction']:.12f} == stack's {st['frac_x']:.12f} (|d| {d_union:.1e})", flush=True)
    if d_union >= 1e-12:
        fail("the union's fraction moved under the CD4 flag")
    out["arm"] = public(e_arm)
    out["arm_fraction"] = fr_arm
    out["arm_frame"] = cd4_audit["residual_frame"]
    out["arm_preserve"] = cd4_audit["post_add_preserve"]
    out["arm_vs_246_compose_max_abs_diff"] = d_eff
    out["union_frame_solved_fraction"] = float(union_audit["residual_frame"]["solved_fraction"])
    print(f"  ARM: macro {100*(e_arm['a23_cos']/e_stack['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_arm['cos_agg']/e_stack['cos_agg']-1):+.4f}%  d pds {e_arm['pds']-e_stack['pds']:+.6f}  "
          f"d q_w {e_arm['a15_sgn_w']-e_stack['a15_sgn_w']:+.6f}", flush=True)

    # ---------------------------------------------------------------- KA-D: energy-matched random control
    target = aB["energy_removed_frac"]
    seed = L.RAND_SEED + 7 * 0 + K                      # sweep246b: seed_i = 0, k = 5
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
    ref_r = ref_b["B2_cd4_randE5"]
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
    (WORK / f"repro249-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"[{held}] {len(out['failures'])} failures; {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; "
          f"wrote {WORK/f'repro249-{held}.json'}", flush=True)
    if out["failures"]:
        raise SystemExit(1)


def summary() -> None:
    per = {h: json.loads((WORK / f"repro249-{h}.json").read_text()) for h in CANCER + ("kolf2",)}
    ref_b = {h: {r["arm"]: r for r in json.loads((E246 / f"work/sweep246b-{h}.json").read_text())["rows"]}
             for h in per}
    ref_a = {h: {r["arm"]: r for r in json.loads((E246 / f"work/sweep246-{h}.json").read_text())["rows"]}
             for h in per}

    def agg(folds, arm_key, ctl=None):
        rm, ra, dp, dq = [], [], [], []
        for h in folds:
            s = per[h]["stack"]
            a = per[h][arm_key] if ctl is None else per[h]["control"][ctl]
            rm.append(a["a23_cos"] / s["a23_cos"] - 1)
            ra.append(a["cos_agg"] / s["cos_agg"] - 1)
            dp.append(a["pds"] - s["pds"])
            dq.append(a["a15_sgn_w"] - s["a15_sgn_w"])
        return dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)),
                    d_qw=float(np.mean(dq)), per_fold_d_pds=[float(x) for x in dp],
                    pds_positive_folds=int(sum(x > 0 for x in dp)))

    def agg_ref(folds, arm):
        rm, ra, dp = [], [], []
        for h in folds:
            s, a = ref_a[h]["STACK"], ref_b[h][arm]
            rm.append(a["a23_cos"] / s["a23_cos"] - 1)
            ra.append(a["cos_agg"] / s["cos_agg"] - 1)
            dp.append(a["pds"] - s["pds"])
        return dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)))

    arm = agg(CANCER, "arm")
    ref = agg_ref(CANCER, "B2_cd4_svd5_resolved")
    ctl_fixed = agg(CANCER, None, "fixed")
    ctl_res = agg(CANCER, None, "resolved")
    ref_ctl = agg_ref(CANCER, "B2_cd4_randE5")
    kolf = agg(("kolf2",), "arm")
    e3 = float(np.mean([per[h]["deflation_audit"]["energy_removed_frac"] for h in CANCER]))
    e4 = float(np.mean([per[h]["deflation_audit"]["energy_removed_frac"] for h in per]))
    fr_stack = float(np.mean([per[h]["frac_cd4_stack"] for h in CANCER]))
    fr_arm = float(np.mean([per[h]["arm_fraction"] for h in CANCER]))
    all_ka = all(v["ok"] for h in per for v in per[h]["ka"].values()) and not any(per[h]["failures"] for h in per)
    lines = []
    P = lambda s: (print(s), lines.append(s))  # noqa: E731
    P("=" * 100)
    P("KNOWN-ANSWER ROW 2 -- the staged arm reproduces 246 (three cancer folds, kolf2 apart)")
    P("=" * 100)
    P(f"  per-fold KA rows: {'ALL PASS' if all_ka else 'FAILURES'} "
      f"({sum(len(per[h]['ka']) for h in per)} rows over {len(per)} folds, every one within 1e-9)")
    P(f"  {'':34s}{'this run (builder)':>22s}{'246 (json)':>14s}{'|d|':>10s}{'mandate':>10s}")
    for label, mine, want, m in (("rel MACRO cosine", arm["rel_macro"], ref["rel_macro"], -0.0081),
                                 ("rel AGGREGATE cosine", arm["rel_agg"], ref["rel_agg"], -0.0046),
                                 ("d pds", arm["d_pds"], ref["d_pds"], 0.01226)):
        P(f"  {label:34s}{mine:>+22.9f}{want:>+14.9f}{abs(mine-want):>10.1e}{m:>+10.5f}")
    P(f"  rounded to the mandate's digits: macro {100*arm['rel_macro']:+.2f}%  agg {100*arm['rel_agg']:+.2f}%  "
      f"d pds {arm['d_pds']:+.5f}   (mandate: -0.81% / -0.46% / +0.01226)")
    P(f"  pds positive on {arm['pds_positive_folds']}/3 cancer folds: {['%+.5f' % x for x in arm['per_fold_d_pds']]}")
    P(f"  energy_removed_frac: 3-cancer-fold mean {e3:.4f}, 4-fold mean {e4:.4f} (mandate 0.202); per fold "
      f"{[round(per[h]['deflation_audit']['energy_removed_frac'], 5) for h in per]}")
    P(f"  CD4 solved fraction, PLOLO's frame: stack {fr_stack:.4f} -> arm {fr_arm:.4f} (246 section 3.3: 0.374 -> 0.396)")
    P(f"  kolf2 apart: macro {100*kolf['rel_macro']:+.3f}%  agg {100*kolf['rel_agg']:+.3f}%  d pds {kolf['d_pds']:+.5f}")
    P(f"  effect matrix, builder path vs 246's compose, max|d| per fold: "
      f"{[per[h]['arm_vs_246_compose_max_abs_diff'] for h in per]}")
    P("")
    P("=" * 100)
    P("KNOWN-ANSWER ROW 3 -- energy-matched random control through the same builder path (one seed)")
    P("=" * 100)
    P(f"  random dims per fold {[per[h]['control']['k_prime'] for h in per]} (246: 52/54/55/52), energy "
      f"{[round(per[h]['control']['energy_removed_frac'], 4) for h in per]}")
    P(f"  fixed sizing (246's row):     macro {100*ctl_fixed['rel_macro']:+.3f}%  agg {100*ctl_fixed['rel_agg']:+.3f}%  "
      f"d pds {ctl_fixed['d_pds']:+.6f}   [246 json: {100*ref_ctl['rel_macro']:+.3f}% / {100*ref_ctl['rel_agg']:+.3f}% / "
      f"{ref_ctl['d_pds']:+.6f}]")
    P(f"  RE-SOLVED through the builder: macro {100*ctl_res['rel_macro']:+.3f}%  agg {100*ctl_res['rel_agg']:+.3f}%  "
      f"d pds {ctl_res['d_pds']:+.6f}   per fold {['%+.5f' % x for x in ctl_res['per_fold_d_pds']]}")
    row3 = ctl_res["d_pds"] <= 0.001 and ctl_fixed["d_pds"] <= 0.001
    P(f"  ROW 3 {'PASS' if row3 else 'FAIL'}: control d pds <= +0.001 (the builder applies a projector, "
      f"not a norm cut: the same energy in random directions buys {ctl_res['d_pds']:+.4f} against the arm's "
      f"{arm['d_pds']:+.4f})")
    row2 = all_ka and abs(arm["rel_macro"] - ref["rel_macro"]) < TOL and abs(arm["rel_agg"] - ref["rel_agg"]) < TOL \
        and abs(arm["d_pds"] - ref["d_pds"]) < TOL
    P(f"\nROW 2 {'PASS' if row2 else 'FAIL'}   ROW 3 {'PASS' if row3 else 'FAIL'}")
    (WORK / "repro249-summary.json").write_text(json.dumps(
        dict(arm=arm, ref246=ref, control_fixed=ctl_fixed, control_resolved=ctl_res, ref246_control=ref_ctl,
             kolf2=kolf, energy_removed_3fold=e3, energy_removed_4fold=e4, fraction_stack=fr_stack,
             fraction_arm=fr_arm, row2_pass=bool(row2), row3_pass=bool(row3), all_ka=bool(all_ka)),
        indent=1))
    (WORK / "repro249-summary.log").write_text("\n".join(lines) + "\n")
    if not (row2 and row3):
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
