"""261 -- the PRODUCTION frame without a 20 GiB staging (258/259's reconstruction off 249's staged artifacts).

    uv run python prod261.py

Recovers the incumbent's whitened base on the real 300-target panel from 249's `effects-off.npy` / `cd4-off.npz` /
`union-off.npz` (`wh = final/s - merged`; the preserve is one scalar), then:
  KA-P1  re-composed h1union  == 249 effects-off.npy  (abeb0749...)           builder path, max|d| <= 1e-9
  KA-P2  re-composed cd4svd5  == 249 effects-svd5.npy (a28d5be2...)           builder path (VCC_CD4_DEFLATE_K=5),
         fractions / preserve at 1e-12, array sha checked                      -- R1 "on the production frame"
  KA-P3  the SAME cd4svd5 through THIS experiment's path (lib246 deflate -> plolo.solve_fraction -> carrying-row
         sizing -> merged_many -> one preserve) == effects-svd5.npy, max|d| <= 1e-9    -- the two code paths agree
  KA-P4  the amplitude calibration factor (package A_c over effect-level A_c on ARM6 / cd4svd5 / cd4svd10)
         == 258's `amp_calibration` at 1e-12
then every real cell of the grid on the production panel, BOTH ways: the builder's own `apply_frozen_cd4_residual`
with the pre-deflated union handed in as `extra` (exactly what a `VCC_UNION_DEFLATE_K` flag would do: deflate the
union, then `solve_carrying_row_fraction` at `H1_UNION_G`) and this experiment's compose; asserted equal at 1e-9.
Each cell's delivered amplitude is read as 211/243's (A_c, B_c) terms on the effect-level render (`deliver`,
minnorm, the scorer's own row weights), plus 216's median-norm proxy.  -> work/prod261.json
"""
from __future__ import annotations

import contextlib
import gc
import importlib.util
import io
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
# the builder reads its knobs at import: the incumbent's configuration (258/259's env), nothing else live
os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
os.environ["VCC_CD4_DEFLATE_K"] = "5"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_UNION2_SOURCES", "VCC_UNION2_G",
           "VCC_CD4_DEFLATE_BASIS", "VCC_CD4_PARTIAL_G", "VCC_CD4_FRACTION", "VCC_SIGN_REPAIR"):
    os.environ.pop(_k, None)

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import lib261 as Y  # noqa: E402

L, X = Y.L, Y.X
ROOT, WORK = Y.ROOT, Y.WORK
E229 = ROOT / "experiments/229-stack4-build"
E243 = ROOT / "experiments/243-transfer-factor"
E255 = ROOT / "experiments/255-cd4-deflate-k10"
sys.path.insert(0, str(E243))
import plolo  # noqa: E402
from plolo import E062, OFFICIAL, PANEL, axis_sha256, deliver  # noqa: E402
from lib243 import CP46, KC  # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402

CTX = ("A", "B", "C")
N_GENES = 18533
K_GRID = (1, 2, 3, 5, 8)


def peak_gib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024**3


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod261_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["prod261_builder"] = module
    spec.loader.exec_module(module)
    return module, path


def builder_compose(MOD, wh, cd4, wcols, panel_positions, k, extra, g):
    """The production path: CD4 deflated at k (None = off) before `solve_cd4_dose`; `extra` (a union object, here
    possibly pre-deflated) sized by `solve_carrying_row_fraction` at g; ONE merged add, ONE preserve."""
    MOD.CD4_DEFLATE_K = k
    MOD.H1_UNION_G = g
    MOD.UNION2_G = None
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        final, cd4_audit, u_audit = MOD.apply_frozen_cd4_residual(
            wh, cd4, np.asarray(wcols, dtype=np.int64), np.asarray(panel_positions, dtype=np.int64),
            extra=extra, extra_audit=None if extra is None else {})
    info = dict(cd4_fraction=float(cd4_audit["fraction"]),
                preserve_scale=float(cd4_audit["post_add_preserve"]["post_add_global_preserve_scale"]),
                union_fraction=None if u_audit is None else float(u_audit["residual_frame"]["solved_fraction"]),
                union_carrying=None if u_audit is None else int(u_audit["residual_frame"]["carrying_rows"]),
                union_median_norm=None if u_audit is None else float(u_audit["residual_frame"]["residual_median_row_norm_over_carrying_rows"]),
                energy_removed=None if "deflation" not in cd4_audit else float(cd4_audit["deflation"]["energy_removed_frac"]))
    return final, info


def my_compose(wh, d_cd4, fr_cd4, u, g, wcols, panel_positions):
    """This experiment's path (lib246 `compose`): the same objects through `merged_many` and one preserve."""
    st = dict(whitened=wh, fold=dict(wcols=np.asarray(wcols, dtype=np.int64), panel_positions=np.asarray(panel_positions, dtype=np.int64)),
              frac_x=float("nan"))
    parts, fu, pu = Y.union_parts(st, d_cd4, fr_cd4, u, g)
    # lib246.compose == this add, then `final[:, panel_positions] = 0` (PLOLO's harness convention; the builder
    # keeps the panel columns at whitened x preserve and asserts it).  Compare the construction, not the zeroing.
    final, _ = E062.add_post_whiten_residual(wh, L.merged_many(parts), 1.0, np.asarray(wcols, dtype=np.int64))
    return final, fu, pu


def terms_from_final(final, targets, b, gate):
    """211/243's per-context (A_c, B_c, Kc_c, cap_c) on the EFFECT-LEVEL render (258/259, verbatim)."""
    gpos = {g: i for i, g in enumerate(OFFICIAL)}
    per = {}
    for c in CTX:
        dh = deliver(final, b[c], gate[c], "minnorm")
        mu = dh + b[c][None, :]
        cols = {j: gpos[t] for j, t in enumerate(targets) if t in gpos}
        w = _row_weights(len(targets), N_GENES, cols)
        a = np.einsum("ij,ij->i", dh, dh)
        for j, col in cols.items():
            a[j] -= dh[j, col] ** 2
        per[c] = dict(A=float((a * w).sum()), B=float(_across_pert_budget(mu, cols, w)), Kc=KC[c],
                      cap=float(w.sum()) * CP46[c], n=len(targets), resolved=len(cols), wsum=float(w.sum()))
    return per


def amp_ratio(per_arm, per_base):
    return float(np.mean([np.sqrt(per_arm[c]["A"] / per_base[c]["A"]) for c in CTX]))


def main() -> None:
    t0 = time.time()
    MOD, builder_path = load_builder()
    assert MOD.CD4_DEFLATE_K == 5 and MOD.H1_UNION_G == 0.28 and MOD.CD4_TARGET_F == 0.34
    for name in ("CD4_SIZING", "UNION2_G", "UNION2_SOURCES", "CD4_DEFLATE_BASIS", "CD4_PARTIAL_G"):
        assert getattr(MOD, name) is None, f"builder global {name} is live"
    bsha = MOD.file_sha256(builder_path)
    print(f"builder {builder_path.relative_to(ROOT)} sha256 {bsha[:16]}", flush=True)
    stage = Y.E249 / "work/stage"
    z = np.load(stage / "cd4-off.npz", allow_pickle=True)
    targets = tuple(str(t) for t in z["targets"])
    cd4 = E062.CD4Prepared(mode="complete_mean3", mean_effect=np.asarray(z["mean_effect"], dtype=np.float64),
                           gene_votes=np.asarray(z["gene_votes"], dtype=np.uint8),
                           target_condition_count=np.asarray(z["target_condition_count"], dtype=np.uint8), condition_details=[])
    u = np.load(stage / "union-off.npz", allow_pickle=True)
    union = E062.CD4Prepared(mode="union", mean_effect=np.asarray(u["mean_effect"], dtype=np.float64),
                             gene_votes=np.asarray(u["gene_votes"], dtype=np.uint8),
                             target_condition_count=np.asarray(u["target_condition_count"], dtype=np.uint8), condition_details=[])
    final_off = np.load(stage / "effects-off.npy")
    final_svd5 = np.load(stage / "effects-svd5.npy")
    aud_off = json.loads((stage / "stage-off.json").read_text())
    aud_svd5 = json.loads((stage / "stage-svd5.json").read_text())
    s_off = float(aud_off["audit"]["h1_union"]["merged_preserve"]["post_add_global_preserve_scale"])
    fc = float(aud_off["audit"]["cd4"]["fraction"])
    fu = float(aud_off["audit"]["h1_union"]["residual_frame"]["solved_fraction"])
    base_norm = float(aud_off["audit"]["h1_union"]["residual_frame"]["whitened_base_median_row_norm"])
    zk = np.load(WORK / "plolo/fold-kolf2.npz", allow_pickle=True)
    wcols = np.asarray(zk["wcols"], dtype=np.int64)
    axis = axis_sha256(OFFICIAL[i] for i in wcols)
    assert axis == plolo.OPERATOR_REGISTRY["union_pruned"]["axis_sha256"], axis
    gpos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([gpos[t] for t in targets], dtype=np.int64)
    vc = np.asarray(cd4.gene_votes) > 0
    vu = np.asarray(union.gene_votes) > 0
    merged = np.where(vc, np.asarray(cd4.mean_effect), 0.0) * fc + np.where(vu, np.asarray(union.mean_effect), 0.0) * fu
    wh = final_off / s_off - merged
    del merged
    ka = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        print(f"  {'PASS' if ok else 'FAIL'} {label:96s} got {got:.12g} want {want:.12g} |d| {d:.1e}", flush=True)

    print("\nKNOWN-ANSWER ROWS (production frame)", flush=True)
    check("recovered base: median row norm == stage-off whitened_base_median_row_norm", float(np.median(np.linalg.norm(wh, axis=1))), base_norm, 1e-9)
    f_off, i_off = builder_compose(MOD, wh, cd4, wcols, panel_positions, None, union, 0.28)
    check("KA-P1 re-composed h1union == 249 effects-off.npy (abeb0749), max|d|", float(np.abs(f_off - final_off).max()), 0.0, 1e-9)
    check("KA-P1 effects-off.npy array sha == abeb0749...", 1.0 if MOD.array_sha256(final_off).startswith("abeb0749") else 0.0, 1.0, 0.0)
    del f_off
    f_5, i_5 = builder_compose(MOD, wh, cd4, wcols, panel_positions, 5, union, 0.28)
    check("KA-P2 re-composed cd4svd5 (builder path) == 249 effects-svd5.npy (a28d5be2), max|d|", float(np.abs(f_5 - final_svd5).max()), 0.0, 1e-9)
    check("KA-P2 cd4 fraction == stage-svd5", i_5["cd4_fraction"], float(aud_svd5["audit"]["cd4"]["fraction"]), 1e-12)
    check("KA-P2 union fraction == stage-svd5", i_5["union_fraction"], float(aud_svd5["audit"]["h1_union"]["residual_frame"]["solved_fraction"]), 1e-12)
    check("KA-P2 preserve == stage-svd5", i_5["preserve_scale"], float(aud_svd5["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]), 1e-12)
    check("KA-P2 effects-svd5.npy array sha == a28d5be2...", 1.0 if MOD.array_sha256(final_svd5).startswith("a28d5be2") else 0.0, 1.0, 0.0)
    # this experiment's path on the production frame
    d5, aud5 = L.deflate(cd4, "svd", 5)
    fr5, pu5 = L.solve_fraction(wh, d5, Y.CD4_F)
    check("KA-P3 lib246 deflate('svd',5) energy removed == builder's", aud5["energy_removed_frac"], float(aud_svd5["audit"]["cd4"]["deflation"]["energy_removed_frac"]), 0.0)
    check("KA-P3 plolo.solve_fraction on the deflated CD4 == stage-svd5 fraction", float(fr5), float(aud_svd5["audit"]["cd4"]["fraction"]), 1e-12)
    m_5, fu_5, pu_5 = my_compose(wh, d5, float(fr5), union, 0.28, wcols, panel_positions)
    check("KA-P3 cd4svd5 through THIS experiment's compose == 249 effects-svd5.npy, max|d|", float(np.abs(m_5 - final_svd5).max()), 0.0, 1e-9)
    check("KA-P3 carrying-row union fraction (compose218) == builder's", fu_5, i_5["union_fraction"], 1e-12)
    del m_5

    # amplitude calibration (258/259's): package A_c over effect-level A_c on the three package pairs
    ctx = {}
    for c in CTX:
        b = plolo.challenge_bctrl(c)
        e = np.expm1(b)
        ctx[c] = (b, (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (b > 0))
    b = {c: ctx[c][0] for c in CTX}
    gate = {c: ctx[c][1] for c in CTX}
    bA, gA = ctx["A"]
    final_arm6 = np.load(E229 / "work/effects-stack6g035.npy")
    final_svd10 = np.load(E255 / "work/stage/effects-svd10.npy")
    inv243 = json.loads((E243 / "work/invert243.json").read_text())
    amp249 = json.loads((Y.E249 / "work/amp249.json").read_text())
    amp255 = json.loads((E255 / "work/amp255.json").read_text())
    arm6_243 = json.loads((E243 / "work/arm6_243.json").read_text())
    pkg = dict(h1union=inv243["arms"]["h1union"]["per"], cd4svd5=amp249["per_arm"], ARM6=arm6_243["per"], cd4svd10=amp255["per_arm"])
    terms: dict[str, dict] = {}
    info_all: dict[str, dict] = {}

    def measure(name, final):
        t1 = time.time()
        per = terms_from_final(final, targets, b, gate)
        H = deliver(final, bA, gA, None)
        terms[name] = dict(per=per, delivered_norm_A=float(np.median(np.linalg.norm(H[:, wcols], axis=1))),
                           median_row_norm_final=float(np.median(np.linalg.norm(final, axis=1))))
        print(f"  {name:22s} A {per['A']['A']:.6f}/{per['B']['A']:.6f}/{per['C']['A']:.6f}  dlvA {terms[name]['delivered_norm_A']:.5f} ({time.time()-t1:.0f}s)", flush=True)
        del H

    for name, final in (("h1union (file)", final_off), ("cd4svd5 (file)", final_svd5), ("ARM6 (file)", final_arm6), ("cd4svd10 (file)", final_svd10)):
        measure(name, final)
    del final_arm6, final_svd10
    calib = []
    for arm, base, fa, fb in (("ARM6", "h1union", "ARM6 (file)", "h1union (file)"),
                              ("cd4svd5", "h1union", "cd4svd5 (file)", "h1union (file)"),
                              ("cd4svd10", "cd4svd5", "cd4svd10 (file)", "cd4svd5 (file)")):
        r_pkg = amp_ratio(pkg[arm], pkg[base])
        r_eff = amp_ratio(terms[fa]["per"], terms[fb]["per"])
        calib.append(dict(arm=arm, base=base, package_pct=100 * (r_pkg - 1), effect_level_pct=100 * (r_eff - 1),
                          ratio_pkg_over_eff=(r_pkg - 1) / (r_eff - 1)))
    amp_cal = float(np.mean([c["ratio_pkg_over_eff"] for c in calib]))
    p258 = json.loads((Y.E258 / "work/price258.json").read_text())
    check("KA-P4 amplitude calibration factor == 258's amp_calibration", amp_cal, p258["factors"]["amp_calibration"], 1e-12)
    n_ok = sum(r["ok"] for r in ka)
    print(f"known-answer rows {n_ok}/{len(ka)} pass  ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    if n_ok != len(ka):
        (WORK / "prod261.json").write_text(json.dumps(dict(known_answer=ka, failed=True), indent=1, default=float))
        raise SystemExit(1)

    # ------------------------------------------------------------------ the grid on the production panel
    print("\nTHE GRID on the production panel (builder path with the pre-deflated union as `extra`; this experiment's compose asserted equal)", flush=True)
    spec = Y.spectrum(union)
    print(f"  union: {spec['n_carried']} carrying rows, centred energy {spec['centred_energy_frac_of_raw']:.4f} of raw, "
          f"top shares {[round(x, 4) for x in spec['top_shares_of_centred'][:8]]}", flush=True)
    cells = [("u_c0_g0.28", "centre", 0, 0.28)] + [(f"u_k{k}_g0.28", "svd", k, 0.28) for k in K_GRID]
    cells += [(f"u_k0_g{g:g}", "none", 0, g) for g in (0.20, 0.35)]
    cells += [(f"u_k{k}_g{g:g}", "svd", k, g) for g in (0.20, 0.35) for k in (3, 5)]
    cross = []
    for name, mode, k, g in cells:
        t1 = time.time()
        uk, audk = L.deflate(union, mode, k) if mode != "none" else L.deflate(union, "none")
        f_b, i_b = builder_compose(MOD, wh, cd4, wcols, panel_positions, 5, uk, g)
        f_m, fu_m, pu_m = my_compose(wh, d5, float(fr5), uk, g, wcols, panel_positions)
        dmax = float(np.abs(f_b - f_m).max())
        cross.append(dict(cell=name, max_abs_diff=dmax, ok=dmax <= 1e-9))
        check(f"cross-path {name}: builder(extra=deflated union) == this experiment's compose, max|d|", dmax, 0.0, 1e-9)
        check(f"cross-path {name}: union fraction", fu_m, i_b["union_fraction"], 1e-12)
        del f_m
        info_all[name] = dict(**i_b, g=g, k=k, mode=mode, deflation=audk, fraction_ratio=i_b["union_fraction"] / i_5["union_fraction"])
        measure(name, f_b)
        del f_b
        gc.collect()
    # controls at every k (amplitude only; they are not priced)
    for k in K_GRID:
        efrac = info_all[f"u_k{k}_g0.28"]["deflation"]["energy_removed_frac"]
        dr, audr, kk = X.rand_energy_matched(union, efrac, Y.RAND_SEED + 5)
        ds, auds = Y.shuffled_projector(union, k, Y.SHUF_SEED + k)
        for name, obj, aud in ((f"randE_k{k}", dr, dict(rand_dim=int(kk), **audr)), (f"shufP_k{k}", ds, auds)):
            f_b, i_b = builder_compose(MOD, wh, cd4, wcols, panel_positions, 5, obj, 0.28)
            info_all[name] = dict(**i_b, g=0.28, k=k, mode=aud["mode"], deflation=aud, fraction_ratio=i_b["union_fraction"] / i_5["union_fraction"])
            measure(name, f_b)
            del f_b
    base = terms["cd4svd5 (file)"]
    amps = {}
    for name in info_all:
        r_eff = amp_ratio(terms[name]["per"], base["per"])
        amps[name] = dict(effect_level_pct=100 * (r_eff - 1), calibrated_pct=100 * amp_cal * (r_eff - 1),
                          median_norm_A_pct=100 * (terms[name]["delivered_norm_A"] / base["delivered_norm_A"] - 1),
                          union_fraction=info_all[name]["union_fraction"], fraction_ratio=info_all[name]["fraction_ratio"],
                          carrying=info_all[name]["union_carrying"], energy_removed=info_all[name]["deflation"].get("energy_removed_frac"))
        print(f"  {name:14s} union fraction {info_all[name]['union_fraction']:.6f} (x{info_all[name]['fraction_ratio']:.4f}), carrying {info_all[name]['union_carrying']}, "
              f"E removed {amps[name]['energy_removed']}, amplitude effect-level {amps[name]['effect_level_pct']:+.4f}% calibrated {amps[name]['calibrated_pct']:+.4f}% "
              f"(216 median-norm A {amps[name]['median_norm_A_pct']:+.4f}%)", flush=True)
    n_ok = sum(r["ok"] for r in ka)
    out = dict(builder_sha256=bsha, known_answer=ka, ka_pass=n_ok, ka_total=len(ka), cross_path=cross, calibration=calib,
               amp_calibration=amp_cal, union_spectrum=spec, terms=terms, compose_info=info_all, amplitude=amps,
               package_terms=pkg, cd4svd5_prod=dict(cd4_fraction=i_5["cd4_fraction"], union_fraction=i_5["union_fraction"],
                                                     union_carrying=i_5["union_carrying"], preserve=i_5["preserve_scale"]),
               seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (WORK / "prod261.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"known-answer rows {n_ok}/{len(ka)}; wrote {WORK/'prod261.json'}  ({out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB)", flush=True)
    if n_ok != len(ka):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
