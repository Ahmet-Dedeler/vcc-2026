"""251 KNOWN-ANSWER ROW 3 -- the builder's PINNED path on 233's disjoint panels, through 162's
PLOLO instrument as 250 ran it (233's per-panel re-fit operator and CD4 re-extraction, 247's
carrying-row sizing at g_cd4), against 250's `pinABC5` rows to the digit.

249's `repro249.py` pattern with 250's `lib250` (byte copy in this directory; its caches land
under 251's `work/plolo/`): the WHOLE of the builder's `apply_frozen_cd4_residual` (pinned deflate
-> carrying-row size -> merge with the H1 union -> one preserve) is substituted for 250's
`deflate_gene_basis` + `solve_cd4` + `compose`.  The env the builder reads at import is October's
(`def_pipeline.DEF_ENV` on top of the shipped arms) plus the pinned basis:
    VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28
    VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641 VCC_CD4_DEFLATE_BASIS=V5_abc

    KA-A  BASE (250's construction, carrying sizing)   == 250 sweep250 `BASE`     (cos_agg, a23_cos, pds; 1e-12)
    KA-B  builder pinned deflation                     == lib250.deflate_gene_basis(cd4, V5_abc) BIT FOR BIT, audit equal
    KA-C  ARM through apply_frozen_cd4_residual        == 250 `pinABC5`             (1e-9), its solved fraction (1e-12),
          the union's fraction untouched, the effect matrix == 250's compose (max |d|)
    KA-D  the diagnostic: captured_frac_of_centred == 250's basis250.json `<panel>|by_abc|5` (1e-12)

One (panel, fold) per process (250 measured 4.2-5.1 GiB); `--summary` aggregates the three cancer
folds per panel and compares with 250 section 4.3's rows (recomputed from 250's own json).

    uv run python experiments/251-pinned-basis/repro251.py --panel real --fold k562_gwps
    uv run python experiments/251-pinned-basis/repro251.py --summary
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
os.environ["VCC_CD4_SIZING"] = "carrying"
os.environ["VCC_CD4_TARGET_G"] = "0.3443655278649641"
os.environ["VCC_CD4_DEFLATE_BASIS"] = "V5_abc"
for _k in ("VCC_CD4_DEFLATE_K", "VCC_UNION_SOURCES", "VCC_UNION2_SOURCES", "VCC_UNION2_G"):
    os.environ.pop(_k, None)

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E250 = ROOT / "experiments/250-deflation-october"
sys.path.insert(0, str(HERE))
import lib250 as L  # noqa: E402  (the byte copy; bind250 -> 233's instrument; caches under 251/work)

assert L.WORK == WORK, f"lib250 bound to {L.WORK}, not 251's work dir"
PANELS = ("real", "efr")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KEYS = ("cos_agg", "a23_cos", "pds")
TOL_KA = 1e-12
TOL_ARM = 1e-9


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("repro251_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repro251_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def public(e: dict) -> dict:
    return {k: v for k, v in e.items() if not k.startswith("_")}


def run_cell(panel: str, held: str) -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.CD4_DEFLATE_BASIS == "V5_abc" and MOD.CD4_DEFLATE_K is None
    assert MOD.CD4_SIZING == "carrying" and abs(MOD.CD4_TARGET_G - L.G_CD4) == 0.0
    assert MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    print(f"builder {builder_path}  sha256 {MOD.file_sha256(builder_path)}", flush=True)
    ref = {r["arm"]: r for r in json.loads((E250 / f"work/sweep250-{panel}-{held}.json").read_text())["rows"]}
    basis250 = json.loads((E250 / "work/basis250.json").read_text())
    V250 = np.load(E250 / "work/basis250-production.npz")["V5_abc"]

    st = L.stage(panel, held)
    fold, targets = st["fold"], st["targets"]
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    pp = np.asarray(fold["panel_positions"], dtype=np.int64)
    cd4, whitened = st["cd4"], st["whitened"]
    print(f"[{panel}/{held}] staged {len(targets)} targets; op {st['opinfo']['path']} axis "
          f"{st['opinfo']['axis_sha256'][:12]}; union frac {st['frac_u']:.6f} "
          f"({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    out: dict = dict(panel=panel, fold=held, n_targets=len(targets),
                     builder_sha256=MOD.file_sha256(builder_path), frac_u=st["frac_u"], ka={}, failures=[])

    def fail(msg):
        print(f"  *** {msg}", flush=True)
        out["failures"].append(msg)

    def ka(name, e, r, tol):
        d = max(abs(r[k] - e[k]) for k in KEYS)
        ok = d < tol
        out["ka"][name] = dict(max_abs_diff=float(d), ok=bool(ok), ref_arm=r["arm"], tol=tol)
        print(f"  {'PASS' if ok else 'FAIL'} {name:56s} vs 250 {r['arm']:10s} max|d| {d:.2e}", flush=True)
        if not ok:
            fail(f"{name}: differs from 250 {r['arm']} by {d:.3e}")

    def score_final(final):
        final = np.asarray(final, dtype=np.float64).copy()
        final[:, pp] = 0.0
        return L.score(st, final)

    # ---------------------------------------------------------------- KA-A: the panel's October base
    frac0, pu0 = L.solve_cd4(whitened, cd4, "carrying")
    base_final, base_pres = L.compose(st, cd4, frac0)
    e_base = L.score(st, base_final)
    ka("KA-A BASE (carrying @ g_cd4, 250's construction)", e_base, ref["BASE"], TOL_KA)
    if abs(frac0 - ref["BASE"]["audit"]["frac_cd4"]) >= 1e-12:
        fail("KA-A: the base fraction differs from 250's")
    out["base"] = public(e_base)
    out["base_fraction"] = float(frac0)

    # ---------------------------------------------------------------- KA-B: the pinned deflation itself
    V, info = MOD.load_cd4_deflation_basis("V5_abc", official_genes=L.OFFICIAL)
    d250, a250 = L.deflate_gene_basis(cd4, V250, "ABCprod")
    dB, aB = MOD.deflate_cd4_residual_pinned(cd4, V, basis_info=info)
    same = (np.array_equal(np.asarray(dB.mean_effect), np.asarray(d250.mean_effect))
            and np.array_equal(np.asarray(dB.gene_votes), np.asarray(d250.gene_votes))
            and np.array_equal(np.asarray(dB.target_condition_count), np.asarray(d250.target_condition_count)))
    same_audit = (aB["energy_removed_frac"] == a250["energy_removed_frac"]
                  and aB["captured_frac_of_centred"] == a250["captured_frac_of_centred"]
                  and aB["n_carried"] == a250["n_carried"])
    ref_aud = ref["pinABC5"]["audit"]
    same_ref = (aB["energy_removed_frac"] == ref_aud["energy_removed_frac"]
                and aB["captured_frac_of_centred"] == ref_aud["captured_frac_of_centred"]
                and aB["n_carried"] == ref_aud["n_carried"])
    okB = bool(np.array_equal(V, V250) and same and same_audit and same_ref)
    out["ka"]["KA-B builder pinned deflation == lib250.deflate_gene_basis"] = dict(
        ok=okB, v_equal=bool(np.array_equal(V, V250)),
        max_abs_diff=float(np.abs(np.asarray(dB.mean_effect) - np.asarray(d250.mean_effect)).max()))
    print(f"  {'PASS' if okB else 'FAIL'} KA-B builder pinned deflation == lib250.deflate_gene_basis: V "
          f"{'bit-identical' if np.array_equal(V, V250) else 'DIFFERS'} to 250's, arrays "
          f"{'bit-identical' if same else 'DIFFER'}, audit {'identical' if same_audit and same_ref else 'DIFFERS'} "
          f"(n_carried {aB['n_carried']}, energy removed {aB['energy_removed_frac']:.9f}, captured "
          f"{aB['captured_frac_of_centred']:.9f}; 250 json {ref_aud['energy_removed_frac']:.9f} / "
          f"{ref_aud['captured_frac_of_centred']:.9f})", flush=True)
    if not okB:
        fail("KA-B: the builder's pinned deflation is not 250's")
    out["deflation_audit"] = {k: v for k, v in aB.items() if k != "provenance"}

    # ---------------------------------------------------------------- KA-D: the census diagnostic vs basis250.json
    want_capt = basis250["cross"][f"{panel}|by_abc|5"]["captured"]
    d_capt = abs(aB["captured_frac_of_centred"] - want_capt)
    okD = d_capt < TOL_KA and aB["n_carried"] == basis250["audits"][f"{panel}|5"]["n_carried"]
    out["ka"]["KA-D captured_frac_of_centred == basis250.json <panel>|by_abc|5"] = dict(
        ok=bool(okD), max_abs_diff=float(d_capt), want=want_capt, got=aB["captured_frac_of_centred"])
    print(f"  {'PASS' if okD else 'FAIL'} KA-D diagnostic captured_frac_of_centred {aB['captured_frac_of_centred']!r} "
          f"vs 250 basis250.json {want_capt!r} (|d| {d_capt:.1e}); own top-5 share "
          f"{aB['own_top_k_energy_frac_of_centred']!r} (250: {basis250['audits'][f'{panel}|5']['top_k_energy_frac_of_centred']!r})",
          flush=True)
    if not okD:
        fail("KA-D: the diagnostic is not 250's")

    # ---------------------------------------------------------------- KA-C: the arm, the builder's entry point
    arm_out, cd4_audit, union_audit = MOD.apply_frozen_cd4_residual(
        whitened, cd4, wcols, pp, extra=st["h1"], extra_audit={})
    e_arm = score_final(arm_out)
    ka("KA-C ARM through apply_frozen_cd4_residual (pinned V5_abc)", e_arm, ref["pinABC5"], TOL_ARM)
    fr_ref = float(ref_aud["frac_cd4"])
    fr_arm = float(cd4_audit["fraction"])
    d_fr = abs(fr_arm - fr_ref)
    print(f"  {'PASS' if d_fr < 1e-12 else 'FAIL'} KA-C solved fraction {fr_arm:.12f} vs 250 {fr_ref:.12f} "
          f"(|d| {d_fr:.1e}); base {frac0:.6f} -> {fr_arm:.6f} (x{fr_arm/frac0:.4f}); sizing rule "
          f"{cd4_audit['residual_frame'].get('sizing_rule')} g={cd4_audit['residual_frame'].get('target_g_requested')}",
          flush=True)
    if d_fr >= 1e-12:
        fail("KA-C: the re-solved fraction differs from 250's")
    fr250, _ = L.solve_cd4(whitened, d250, "carrying")
    ref_final, _ = L.compose(st, d250, fr250)
    chk = np.asarray(arm_out, dtype=np.float64).copy()
    chk[:, pp] = 0.0
    d_eff = float(np.abs(chk - ref_final).max())
    print(f"  effect matrix: builder path vs 250 compose max|d| {d_eff:.2e}", flush=True)
    d_union = abs(float(union_audit["residual_frame"]["solved_fraction"]) - st["frac_u"])
    print(f"  {'PASS' if d_union < 1e-12 else 'FAIL'} the union's solved fraction through the builder "
          f"{union_audit['residual_frame']['solved_fraction']:.12f} == 250's {st['frac_u']:.12f} (|d| {d_union:.1e})",
          flush=True)
    if d_union >= 1e-12:
        fail("the union's fraction moved under the CD4 flag")
    assert cd4_audit["deflation"]["mode"] == "pinned" and cd4_audit["deflation"]["basis"] == "V5_abc"
    out["arm"] = public(e_arm)
    out["arm_fraction"] = fr_arm
    out["arm_frame"] = cd4_audit["residual_frame"]
    out["arm_deflation"] = {k: v for k, v in cd4_audit["deflation"].items() if k != "provenance"}
    out["arm_vs_250_compose_max_abs_diff"] = d_eff
    out["union_frame_solved_fraction"] = float(union_audit["residual_frame"]["solved_fraction"])
    print(f"  ARM: macro {100*(e_arm['a23_cos']/e_base['a23_cos']-1):+.4f}%  agg "
          f"{100*(e_arm['cos_agg']/e_base['cos_agg']-1):+.4f}%  d pds {e_arm['pds']-e_base['pds']:+.6f}  "
          f"d q_w {e_arm['a15_sgn_w']-e_base['a15_sgn_w']:+.6f}", flush=True)

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = peak_gib()
    (WORK / f"repro251-{panel}-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"[{panel}/{held}] {len(out['failures'])} failures; {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB; "
          f"wrote {WORK/f'repro251-{panel}-{held}.json'}", flush=True)
    if out["failures"]:
        raise SystemExit(1)


def summary() -> None:
    lines = []
    P = lambda s: (print(s), lines.append(s))  # noqa: E731
    P("=" * 100)
    P("KNOWN-ANSWER ROW 3 -- the pinned path on 233's `real` and 238's `efr` reproduces 250's pinABC5 rows")
    P("=" * 100)
    all_ok = True
    res = {}
    for panel in PANELS:
        per = {h: json.loads((WORK / f"repro251-{panel}-{h}.json").read_text()) for h in CANCER}
        ref = {h: {r["arm"]: r for r in json.loads((E250 / f"work/sweep250-{panel}-{h}.json").read_text())["rows"]}
               for h in CANCER}
        rm, ra, dp, dq = [], [], [], []
        rm_r, ra_r, dp_r = [], [], []
        for h in CANCER:
            b, a = per[h]["base"], per[h]["arm"]
            rm.append(a["a23_cos"] / b["a23_cos"] - 1); ra.append(a["cos_agg"] / b["cos_agg"] - 1)
            dp.append(a["pds"] - b["pds"]); dq.append(a["a15_sgn_w"] - b["a15_sgn_w"])
            rb, rp = ref[h]["BASE"], ref[h]["pinABC5"]
            rm_r.append(rp["a23_cos"] / rb["a23_cos"] - 1); ra_r.append(rp["cos_agg"] / rb["cos_agg"] - 1)
            dp_r.append(rp["pds"] - rb["pds"])
        ok_ka = all(v["ok"] for h in CANCER for v in per[h]["ka"].values()) and not any(per[h]["failures"] for h in CANCER)
        mine = dict(rel_macro=float(np.mean(rm)), rel_agg=float(np.mean(ra)), d_pds=float(np.mean(dp)),
                    d_qw=float(np.mean(dq)), per_fold_d_pds=[float(x) for x in dp],
                    pds_positive_folds=int(sum(x > 0 for x in dp)))
        want = dict(rel_macro=float(np.mean(rm_r)), rel_agg=float(np.mean(ra_r)), d_pds=float(np.mean(dp_r)))
        d = max(abs(mine[k] - want[k]) for k in want)
        ok = ok_ka and d < TOL_ARM
        all_ok &= ok
        capt = [per[h]["deflation_audit"]["captured_frac_of_centred"] for h in CANCER]
        res[panel] = dict(mine=mine, ref250=want, max_abs_diff=float(d), ka_all=bool(ok_ka), ok=bool(ok),
                          captured=capt, n_carried=[per[h]["deflation_audit"]["n_carried"] for h in CANCER],
                          energy_removed=[per[h]["deflation_audit"]["energy_removed_frac"] for h in CANCER],
                          fractions=[per[h]["arm_fraction"] for h in CANCER],
                          base_fractions=[per[h]["base_fraction"] for h in CANCER],
                          eff_max_abs_diff=[per[h]["arm_vs_250_compose_max_abs_diff"] for h in CANCER],
                          peak_rss_gib=[per[h]["peak_rss_gib"] for h in CANCER],
                          seconds=[per[h]["seconds"] for h in CANCER])
        P(f"\n  panel `{panel}`: per-fold KA rows {'ALL PASS' if ok_ka else 'FAILURES'} "
          f"({sum(len(per[h]['ka']) for h in CANCER)} rows over 3 folds)")
        P(f"  {'':30s}{'this run (builder)':>22s}{'250 (json)':>14s}{'|d|':>10s}{'250 sec. 4.3':>14s}")
        for label, key, digits in (("rel MACRO cosine", "rel_macro", "%"), ("rel AGGREGATE cosine", "rel_agg", "%"),
                                   ("d pds", "d_pds", "")):
            m, w = mine[key], want[key]
            P(f"  {label:30s}{m:>+22.9f}{w:>+14.9f}{abs(m-w):>10.1e}"
              f"{('%+.2f%%' % (100*w)) if digits else ('%+.5f' % w):>14s}")
        P(f"  rounded: macro {100*mine['rel_macro']:+.2f}%  agg {100*mine['rel_agg']:+.2f}%  d pds {mine['d_pds']:+.5f}  "
          f"d q_w {mine['d_qw']:+.5f}; pds positive on {mine['pds_positive_folds']}/3 folds "
          f"{['%+.5f' % x for x in mine['per_fold_d_pds']]}")
        P(f"  captured_frac_of_centred per fold {[round(c, 6) for c in capt]} (250: 0.10-0.14); n_carried "
          f"{res[panel]['n_carried']}; energy removed {[round(e, 6) for e in res[panel]['energy_removed']]}; "
          f"CD4 fraction base -> arm {[round(f, 5) for f in res[panel]['base_fractions']]} -> "
          f"{[round(f, 5) for f in res[panel]['fractions']]}")
        P(f"  effect matrix vs 250's compose, max|d| per fold {res[panel]['eff_max_abs_diff']}; peak RSS "
          f"{[round(x, 2) for x in res[panel]['peak_rss_gib']]} GiB")
        P(f"  ROW 3 [{panel}] {'PASS' if ok else 'FAIL'} (max |d| on the three-fold means {d:.2e} < {TOL_ARM})")
    P(f"\nROW 3 {'PASS' if all_ok else 'FAIL'}")
    (WORK / "repro251-summary.json").write_text(json.dumps(res, indent=1))
    (WORK / "repro251-summary.log").write_text("\n".join(lines) + "\n")
    if not all_ok:
        raise SystemExit(1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=PANELS)
    ap.add_argument("--fold", choices=CANCER)
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    if a.summary:
        summary()
    else:
        run_cell(a.panel, a.fold)


if __name__ == "__main__":
    main()
