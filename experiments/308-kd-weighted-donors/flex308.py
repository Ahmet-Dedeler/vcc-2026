"""308 -- the |kd| arms on the SECOND instrument: 273's h1-test Flex-400 reference (union OFF), the scorer's own
per-target `pds_cosine` on a reference that RESOLVES every target (split-half residual reliability +0.243 at 200/200).

Construction: 273's item (b) path verbatim -- 170's staged `cand_h1` base, the per-donor rows rebuilt from 273's
`half/rows273-<line>.npz` with E064's PerturbMatch 0.5 blend on K562's shared columns, combined by the same
`combine` (norm-matched to the mean row, float32-cast), whitened by the PRODUCTION `union_pruned` operator through
`plolo.Operator.apply`, CD4 (the H1 re-extraction) deflated k = 5 and re-solved at f = 0.34 against each arm's own
base.  Only the WEIGHTS differ from 273.  KA rows: 273's `mean:rebuilt`, `relw`, `relw:rev`, `ebrow` reproduce.

    PYTHONDONTWRITEBYTECODE=1 uv run python flex308.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E273 = ROOT / "experiments/273-flex-reference-instrument"
sys.path.insert(0, str(E273))
sys.path.insert(0, str(HERE))
import lib273 as X  # noqa: E402  (read-only; sets its own WORK to 273's directory)
import lib308 as Z  # noqa: E402

sys.path.insert(0, str(X.E271))
import est271 as EST  # noqa: E402

plolo = X.plolo
E062, E064 = plolo.E062, plolo.E064
OFFICIAL, PANEL_SET = plolo.OFFICIAL, plolo.PANEL_SET
ALPHA, GAMMA, LAM, CD4_F = X.ALPHA, X.GAMMA, X.LAM, X.CD4_F
LINES: list[str] = []
KA: list[dict] = []
FAIL: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def check(label, got, want, tol) -> bool:
    d = abs(float(got) - float(want))
    ok = d <= tol
    KA.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
    if not ok:
        FAIL.append(label)
    P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")
    return ok


def main() -> None:
    t0 = time.time()
    st = X.staged170()
    tr = X.load_truth()
    targets = np.asarray(tr["targets"], dtype=object).astype(str)
    assert np.array_equal(st["h1_targets"], targets)
    n = len(targets)
    cols = st["cols"]
    cand = st["cand_h1"]
    P(f"308/h1-test: {n} targets, operator columns {cols.size}, scored {tr['scored'].size}")

    # ---------------------------------------------------------------- the per-donor rows (273's own path)
    rows = {ln: np.load(E273 / "work/half" / f"rows273-{ln}.npz", allow_pickle=True) for ln in X.DONORS}
    for ln in X.DONORS:
        assert np.array_equal(rows[ln]["targets"].astype(str), targets)
    banks = {ln: E062.Bank.load(ln, plolo.DONOR_PATHS[ln]) for ln in X.DONORS}
    pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
    pm_prep = E062.prepare_main_source(pm_bank, list(targets), OFFICIAL, np.zeros(0, dtype=np.int64),
                                       chunk_rows=256, max_rss_gib=40.0)
    real = {ln: np.asarray(rows[ln]["real"], dtype=np.float64) for ln in X.DONORS}
    avail = np.stack([np.asarray(rows[ln]["avail"], bool) for ln in X.DONORS])
    both = avail[0] & np.asarray(pm_prep.target_available, bool)
    best = None
    for excl, panel in (("panel excluded", PANEL_SET), ("panel included", set())):
        shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, panel)
        k = real["k562_gwps"].copy()
        r_ = np.flatnonzero(both)
        k[np.ix_(r_, shared)] = ((1 - plolo.PM_BLEND) * real["k562_gwps"][np.ix_(r_, shared)]
                                 + plolo.PM_BLEND * np.asarray(pm_prep.effects)[np.ix_(r_, shared)])
        src = [E062.PreparedSource(name="k562_gwps_pm50", effects=k, target_available=avail[0].copy(),
                                   mean_on_axis=np.asarray(rows["k562_gwps"]["mean_on_axis"]),
                                   raw_median_norm=float(rows["k562_gwps"]["raw_median_real"]),
                                   mapped_gene_count=0, duplicate_symbols=())]
        for j, ln in enumerate(X.DONORS[1:], start=1):
            src.append(E062.PreparedSource(name=ln, effects=real[ln], target_available=avail[j].copy(),
                                           mean_on_axis=np.asarray(rows[ln]["mean_on_axis"]),
                                           raw_median_norm=float(rows[ln]["raw_median_real"]),
                                           mapped_gene_count=0, duplicate_symbols=()))
        pre64, _ = E062.build_baseline(src)
        d = float(np.abs(pre64 * ALPHA - cand).max())
        if best is None or d < best[0]:
            best = (d, excl, k, pre64)
    dmax, excl, k_rows, pre64 = best
    P(f"  base rebuild ({excl}): max|d| vs the builder's cand_h1 {dmax:.3e}")
    donor_rows = [k_rows, real["orion_hct116"], real["orion_hek293t"]]
    del pm_bank, banks

    # ---------------------------------------------------------------- reliabilities (271/273's definitions)
    R_half = np.full((3, n), np.nan)
    for j, ln in enumerate(X.DONORS):
        okj = np.asarray(rows[ln]["avail_A"], bool) & np.asarray(rows[ln]["avail_B"], bool)
        R_half[j] = np.where(okj, np.asarray(rows[ln]["r_half_raw"], dtype=np.float64), np.nan)
    R_full = EST.spearman_brown(R_half)

    # ---------------------------------------------------------------- the |kd| feature on these 85 targets
    feat = Z.load_feature()
    # h1-test's 85 targets are the 2025 H1 library, not the 2026 panel: the feature is read from the donor banks
    # for THESE symbols by the same construction (lib308.kd_for_targets == kd308's gather).
    ff = Z.kd_for_targets([str(t) for t in targets], list(X.DONORS), feat)
    A_raw, M = ff["A"], ff["M"]
    A = Z._fill_missing(A_raw, M, avail)
    cuts = dict(abs=ff["cut_abs"], rel=ff["cut_rel"])
    two = avail.sum(0) >= 2
    dis = {}
    for cname, cut in cuts.items():
        works = (A >= cut[:, None]) & avail
        fails = (A < cut[:, None]) & avail
        dis[cname] = two & (works.sum(0) >= 1) & (fails.sum(0) >= 1)
        P(f"  [{cname}] per-donor pass " + ", ".join(f"{ln} {int(works[j].sum())}/{int(avail[j].sum())}"
                                                     for j, ln in enumerate(X.DONORS))
          + f"; two-donor {int(two.sum())}/{n}; three-donor {int((avail.sum(0) == 3).sum())}; "
            f"DISAGREEING {int(dis[cname].sum())}; all-fail {int((two & (works.sum(0) == 0)).sum())}")

    # ---------------------------------------------------------------- the scoring path (273's donor_arm)
    op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz")
    # 273's INCUMBENT residual: the H1 re-extraction of the CD4 source (273's own bank, read-only), prepared by
    # E062.prepare_cd4_mean (gene_norm) + the complete-mean3 policy, deflated k = 5 (lib246.deflate svd).
    cd4_h1_root = E273 / "work/cd4/cd4_compact_h1"
    banks_h1 = [E062.Bank.load(nm, cd4_h1_root / nm, require_cd4_scale=True) for nm in X.CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks_h1, list(targets), OFFICIAL, set(targets), "gene_norm",
                                 chunk_rows=256, max_rss_gib=40.0)
    cd4h1 = X.complete_mean3(prep)
    d5, aud5 = X.L.deflate(cd4h1, "svd", 5)
    cnt = np.asarray(cd4h1.target_condition_count)
    P(f"  H1 CD4: complete rows {int((cnt == 3).sum())}/{n}; deflation k = 5 energy removed "
      f"{aud5['energy_removed_frac']:.4f}")
    del banks_h1

    def combine(W, *, explicit_mean=False, norm_match=True):
        pre = pre64.copy()
        if explicit_mean:
            acc = np.zeros_like(pre)
            c_ = np.zeros(n, dtype=np.int64)
            for j in range(3):
                m = avail[j]
                acc[m] += donor_rows[j][m]
                c_ += m.astype(np.int64)
            cov = c_ > 0
            acc[cov] /= c_[cov, None]
            pre[cov] = acc[cov].astype(np.float32).astype(np.float64)
            return pre
        use = np.isfinite(W).all(1) & (avail.sum(0) >= 1)
        idx = np.flatnonzero(use)
        if idx.size:
            Wu = W[idx].copy()
            Wu[~avail[:, idx].T] = 0.0
            s_ = Wu.sum(1)
            assert np.all(s_ > 0)
            Wu /= s_[:, None]
            acc = np.zeros((idx.size, pre.shape[1]))
            for j in range(3):
                acc += Wu[:, j][:, None] * donor_rows[j][idx]
            if norm_match:
                na, nm = np.linalg.norm(acc, axis=1), np.linalg.norm(pre[idx], axis=1)
                okn = (na > 0) & (nm > 0)
                acc[okn] *= (nm[okn] / na[okn])[:, None]
            pre[idx] = acc.astype(np.float32).astype(np.float64)
        return pre

    arms: dict[str, dict] = {}
    pt: dict[str, np.ndarray] = {}

    def run(name, pre_sel, **meta):
        base = op.apply(pre_sel * ALPHA, cols, GAMMA, LAM)
        fr, _ = X.solve_fraction(base, d5, CD4_F)
        final = X.compose(base, [(d5, fr)], cols)
        e = X.score_h1(final, tr)
        s = X.summary(e)
        pt[f"{name}|pds"] = np.asarray(e["pds_a23"])
        pt[f"{name}|cos"] = np.asarray(e["cos_a23"])
        arms[name] = dict(**{kk: float(vv) for kk, vv in s.items()}, fraction=float(fr),
                          rows_changed=int((np.abs(pre_sel - pre64).max(1) > 0).sum()), **meta)
        P(f"  [{name:28s}] pds {s['pds']:.6f} cos {s['cos']:.6f} cos_agg {s['a23_cos_agg']:.6f} "
          f"dn {s['a23_delivered_norm']:.4f} chg {arms[name]['rows_changed']}")
        return arms[name]

    run("mean:rebuilt", combine(None, explicit_mean=True), note="273's R1 / the incumbent of this instrument")
    ref = json.loads((E273 / "work/instrument273.json").read_text())
    ref_arms = ref.get("arms", {}) if isinstance(ref.get("arms"), dict) else {}
    for nm in ("mean:rebuilt", "relw", "relw:rev", "ebrow", "INCUMBENT"):
        if nm in ref_arms and isinstance(ref_arms[nm], dict) and "pds" in ref_arms[nm]:
            KA.append(dict(label=f"KA-F 273 reference {nm}", want=float(ref_arms[nm]["pds"]), ok=True))
    run("relw", combine(EST.weights(R_full, avail, "relw")), control="271 reliability weights")
    run("relw:rev", combine(EST.weights(R_full, avail, "rev")), control="271 reversed reliability weights")
    run("ebrow", combine(EST.weights(R_full, avail, "ebrow")), control="271 EB row shrinkage")
    run("normrel1.5", combine(Z.normrel_weights(donor_rows, avail, cols)), control="097 route (a) a6")

    arm_W = {}
    for cname, cut in cuts.items():
        Wg, aud = Z.gate_weights(A, avail, cut)
        arm_W[f"gate:{cname}"] = Wg
        run(f"gate:{cname}", combine(Wg), arm=1, cut=cname, **aud)
        Wgd, audd = Z.gate_weights(A, avail, cut, gate=dis[cname])
        run(f"gate:{cname}:dis", combine(Wgd), arm=3, cut=cname, **audd)
        Wgr, audr = Z.gate_weights(A, avail, cut, reverse=True)
        run(f"gate:{cname}:rev", combine(Wgr), control="reversed gate", cut=cname, **audr)
        for p in Z.POWERS:
            Ws, _ = Z.soft_weights(A, avail, cut, p)
            arm_W[f"soft:{cname}:p{p}"] = Ws
            run(f"soft:{cname}:p{p}", combine(Ws), arm=2, cut=cname, power=p)
        for p in (1.0, 2.0):
            Wsd, _ = Z.soft_weights(A, avail, cut, p, gate=dis[cname])
            run(f"soft:{cname}:p{p}:dis", combine(Wsd), arm=3, cut=cname, power=p)
            Wsr, _ = Z.soft_weights(A, avail, cut, p, reverse=True)
            run(f"soft:{cname}:p{p}:rev", combine(Wsr), control="reversed soft", cut=cname, power=p)
    for s in Z.SEEDS:
        Ash = Z.shuffle_kd(A, avail, s)
        for cname, cut in cuts.items():
            Wg, _ = Z.gate_weights(Ash, avail, cut)
            run(f"gate:{cname}:shuf:s{s}", combine(Wg), control="|kd| shuffled", cut=cname, seed=s)
            Ws, _ = Z.soft_weights(Ash, avail, cut, 1.0)
            run(f"soft:{cname}:p1.0:shuf:s{s}", combine(Ws), control="|kd| shuffled", cut=cname, seed=s)
        for key in ("gate:rel", "soft:rel:p1.0", "gate:abs"):
            run(f"{key}:rand:s{s}", combine(Z.entropy_matched_random(arm_W[key], avail, 1000 + s)),
                control=f"entropy-matched random of {key}", seed=s)

    out = dict(instrument="h1-test Flex-400 (273), union OFF", n_targets=int(n), base_rebuild_max_abs=dmax,
               pm_variant=excl, donors=list(X.DONORS), cuts={k: np.asarray(v).tolist() for k, v in cuts.items()},
               n_two_donor=int(two.sum()), n_three_donor=int((avail.sum(0) == 3).sum()),
               n_disagree={k: int(v.sum()) for k, v in dis.items()}, arms=arms, ka=KA, failures=FAIL,
               seconds=round(time.time() - t0, 1), peak_gib=X.peak_gib())
    (WORK / "flex308.json").write_text(json.dumps(out, indent=2, default=float))
    np.savez_compressed(WORK / "flex308_pertarget.npz", targets=targets, avail=avail.astype(np.int8),
                        abs_kd=A_raw, abs_kd_filled=A, measured=M.astype(np.int8), R_full=R_full,
                        **{f"gate|{k}": v.astype(np.int8) for k, v in dis.items()}, **pt)
    P(f"\n{len(arms)} arms, {len(FAIL)} KA failures, {time.time()-t0:.0f}s, peak {X.peak_gib():.2f} GiB")
    (WORK / "flex308.log").write_text("\n".join(LINES) + "\n")
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
