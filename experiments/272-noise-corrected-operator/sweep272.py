"""272 -- construction A on ONE PLOLO cancer fold: the donors' noise spectrum in the fold operator's eigenbasis
(KA-0 first: 271's in-sample anti-correlation reproduced), the Wiener-damped operator at kappa in {0.25, 0.5, 1, 2, 4},
its shuffled-s and mirrored-s controls, the hard truncation reading, through cd4svd5's construction (263's staging:
CD4 deflated k = 5 re-solved at f = 0.34, the H1-pair union re-solved by the carrying-row rule at g = 0.28, one add,
one preserve; 188's exact nmae identity; per-target pds / cosine arrays).  MEASUREMENT ONLY.

    PYTHONDONTWRITEBYTECODE=1 uv run python sweep272.py orion_hct116
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/272-noise-corrected-operator"
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E257 = ROOT / "experiments/257-soft-deflation-price"
E267 = ROOT / "experiments/267-replicate-validated-oracle"
sys.path.insert(0, str(E263))
sys.path.insert(0, str(HERE))
import lib263 as Y  # noqa: E402   (imports lib257 -> lib246 -> plolo / resid218 / compose218)
import op272 as O  # noqa: E402

Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK
plolo, L, X = Y.plolo, Y.L, Y.X
E062 = plolo.E062
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main(held: str) -> None:
    t0 = time.time()
    op_sha = O.sha256_self()
    P(f"[{held}] op272.py sha256 {op_sha}")
    st = Y.stage(held)
    ident = X.Identity(st)
    n = st["n"]
    targets = st["targets"]
    wc, sc = st["wcols"], st["sc"]
    avail = st["avail"]
    donors = list(st["sources"])
    op = st["op"]
    op_meta = json.loads(Path(L.op_path_for(held)).with_suffix(".json").read_text())
    P(f"[{held}] staged {n} targets, donors {donors} ({st['donors']}), PM rows {st['pm_n']}, operator {op.name} {op.sha256[:12]} "
      f"n_cols {op.n_cols} basis rows {op_meta['basis_rows_by_source']}, KA-0 pre max|d| {st['ka0']['pre_float32_max_abs_diff']:.1e}, "
      f"{time.time() - t0:.0f}s, peak {Y.peak_gib():.2f} GiB")
    failures: list[str] = []
    ka_log: list[dict] = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka_log.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
        if not ok:
            failures.append(f"{label}: got {got!r} want {want!r} |d| {d:.3e}")
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")
        return ok

    check("KA-0a build_baseline(per-donor rows) float32 == cached pre (263)", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)

    # ------------------------------------------------------------------ the half rows and KA-0 (271's anti-correlation)
    diag271 = json.loads((ROOT / "experiments/271-donor-reliability-estimation/work/diagwht271.json").read_text())[held]
    PANEL = list(plolo.PANEL)
    H = {line: O.load_half_rows(line, PANEL) for line in donors}          # all 300 panel targets (271's target set)
    Hf = {line: O.load_half_rows(line, targets) for line in donors}       # the fold's targets
    for j, line in enumerate(donors):
        check(f"KA-0b availability of {line}'s half rows == the fold's donor availability", int((Hf[line]["ok"] != avail[j]).sum()), 0, 0)
    hf = {}
    anti = {}
    for line in donors:
        h = H[line]
        ok = h["ok"]
        WA = op.apply(h["A"] * ALPHA, wc, GAMMA, LAM)[:, wc]
        WB = op.apply(h["B"] * ALPHA, wc, GAMMA, LAM)[:, wc]
        c = O.row_cos(WA[ok], WB[ok])
        pi = O.derangement(int(ok.sum()), 1)
        null = O.row_cos(WA[ok], WB[ok][pi])
        craw = O.row_cos(h["A"][ok][:, wc], h["B"][ok][:, wc])
        anti[line] = dict(n=int(ok.sum()), wht_median=float(np.nanmedian(c)), wht_share_le0=float(np.nanmean(c <= 0)), wht_null_median=float(np.nanmedian(null)),
                          wht_null_sd=float(np.nanstd(null)), raw_median=float(np.nanmedian(craw)), raw_share_le0=float(np.nanmean(craw <= 0)))
        ref = diag271[line]
        check(f"KA-0 {line}: in-sample whitened split-half cosine median == 271 diagwht ({ref['wht']['median']:+.4f})", anti[line]["wht_median"], ref["wht"]["median"], 1e-3)
        check(f"KA-0 {line}: raw split-half median == 271 ({ref['raw']['median']:+.4f})", anti[line]["raw_median"], ref["raw"]["median"], 1e-3)
        hf[line] = O.half_factor_check(h, wc)
        P(f"  half factor {line}: mean||d||^2 {hf[line]['mean_sq_d']:.4f} vs implied noise energy (half - real) {hf[line]['implied_noise_energy']:.4f} -> ratio {hf[line]['ratio_d_over_implied']:.4f}; "
          f"||half||/||real|| {hf[line]['median_norm_ratio_half_over_real']:.4f}, cos((A+B)/2, real) {hf[line]['median_cos_halfmean_real']:.4f}, "
          f"signal share of the real row's energy {hf[line]['signal_share_of_real_energy']:.4f}")
    # the held line's halves through the same operator (out-of-sample), for the record
    Hh = O.load_half_rows(held, PANEL)
    WA = op.apply(Hh["A"] * ALPHA, wc, GAMMA, LAM)[:, wc]
    WB = op.apply(Hh["B"] * ALPHA, wc, GAMMA, LAM)[:, wc]
    okh = Hh["ok"]
    c = O.row_cos(WA[okh], WB[okh])
    anti[held] = dict(role="held (out-of-sample)", n=int(okh.sum()), wht_median=float(np.nanmedian(c)), wht_share_le0=float(np.nanmean(c <= 0)))
    check(f"KA-0 {held} (held, out-of-sample): whitened split-half median == 271 ({diag271[held]['wht']['median']:+.4f})", anti[held]["wht_median"], diag271[held]["wht"]["median"], 1e-3)
    del WA, WB, Hh

    # ------------------------------------------------------------------ the noise spectrum in the fitted eigenbasis
    raw_medians = {line: H[line]["raw_median_real"] for line in donors}
    unit, renorm = O.basis_unit_factor(raw_medians)
    weights = [float(op_meta["basis_rows_by_source"][line]) for line in donors]
    D_list = [O.half_differences(H[line])[:, wc] * unit for line in donors]      # every carried panel row of each donor
    t1 = time.time()
    n_i, s, per_donor_n = O.noise_spectrum(op.evals, op.evecs, D_list, weights)
    # units KA: the donors' FULL panel rows in basis units projected on the eigenbasis vs lambda (in-sample rows of the fit)
    Pfull = np.concatenate([H[line]["real"][H[line]["ok"]][:, wc] * unit for line in donors])
    e_raw_full, _ = O.whitened_energy(Pfull, op.evecs, op.evals, GAMMA)
    lam_hat = e_raw_full / Pfull.shape[0]
    order = np.argsort(op.evals)[::-1]
    ratio_dec = [float(np.median((lam_hat / op.evals)[order][k * op.n_cols // 10:(k + 1) * op.n_cols // 10])) for k in range(10)]
    P(f"  units check: median (panel-row variance / lambda) by rank decile {[round(x, 3) for x in ratio_dec]} (in-sample rows; O(1) expected; "
      f"unit factor ref/EFFECT_SCALE {unit:.6f}, renorm {renorm})")
    e_raw, e_wht = O.whitened_energy(np.asarray(st["pre"], dtype=np.float64)[:, wc] * ALPHA, op.evecs, op.evals, GAMMA)
    rep = O.spectrum_report(op.evals, s, n_i, e_raw, e_wht, GAMMA)
    rep.update(half_factor=O.HALF_FACTOR, half_factor_check=hf, unit_factor=unit, renorm=renorm, basis_weights=dict(zip(donors, weights)),
               lambda_hat_over_lambda_by_decile=ratio_dec, seconds=round(time.time() - t1, 1))
    P(f"  spectrum: p {rep['p']}, floored {rep['n_lambda_floored']}; noise-dominated (s < 0.5) {rep['n_noise_dominated']} ({100*rep['share_noise_dominated']:.1f} %), "
      f"s = 0 on {rep['n_s_zero']} ({100*rep['share_s_zero']:.1f} %); whitened energy of the prediction rows in noise-dominated directions "
      f"{100*rep['whitened_energy_share_noise_dominated']:.1f} % (s = 0: {100*rep['whitened_energy_share_s_zero']:.1f} %; raw energy there {100*rep['raw_energy_share_noise_dominated']:.1f} %); "
      f"rank of the noise-dominated set: median {rep['rank_of_noise_dominated']['median']:.0f} (p10 {rep['rank_of_noise_dominated']['p10']:.0f}, p90 {rep['rank_of_noise_dominated']['p90']:.0f}, first {rep['rank_of_noise_dominated']['first']}); "
      f"mean s {rep['mean_s']:.3f}, energy-weighted {rep['energy_weighted_mean_s']:.3f}; sum n / sum lambda {rep['total_noise_over_total_lambda']:.3f}")
    for r in rep["deciles"]:
        P(f"    decile {r['decile']}: lambda [{r['lambda_min']:.2e}, {r['lambda_max']:.2e}] amp median {r['amplification_median']:.1f}x; mean s {r['mean_s']:.3f} "
          f"(share < 0.5 {r['share_s_lt_half']:.2f}, = 0 {r['share_s_zero']:.2f}); n/lambda median {r['median_noise_over_lambda']:.3f}; "
          f"raw energy {100*r['raw_energy_share']:.1f} %, whitened energy {100*r['whitened_energy_share']:.1f} %")
    s_shuf = O.shuffled_s(s)
    s_mirror = O.mirrored_s(s, op.evals)
    s_hard = O.hard_s(s)
    n_shuf = O.shuffled_s(n_i)
    n_mirror = O.mirrored_s(n_i, op.evals)
    rep["n_zero"] = int((n_i == 0).sum())
    good = (n_i > 0) & (op.evals > op.evals.max() * 1e-11)
    rep["spearman_log_n_log_lambda"] = float(np.corrcoef(np.argsort(np.argsort(np.log(n_i[good]))), np.argsort(np.argsort(np.log(op.evals[good]))))[0, 1])
    P(f"  n_i = 0 on {rep['n_zero']} directions; rank correlation of log n_i with log lambda_i {rep['spearman_log_n_log_lambda']:+.3f}")

    # ------------------------------------------------------------------ the scoring path
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    ref263 = {r["arm"]: r for r in json.loads((E263 / f"work/sweep263-{held}.json").read_text())["rows"]}
    pt263 = np.load(E263 / f"work/pertarget263-{held}.npz")
    ref257 = {r["arm"]: r for r in json.loads((E257 / f"work/sweep257-{held}.json").read_text())["rows"]}

    def run(name, dop: O.DampedOperator | None = None, *, base=None, parts=None, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        if base is None:
            base = st["whitened"] if dop is None else dop.apply(st["pre"] * ALPHA, wc, GAMMA, LAM)
        if parts is None:
            parts, aud = Y.arm_parts(st, base)
        else:
            aud = dict(cd4_fraction=float(parts[0][1]), union_fraction=float(parts[1][1]), note="STACK: un-deflated CD4 at the staged fraction")
        st_arm = dict(st, whitened=base)
        final = L.compose(st_arm, parts)
        e = L.score(st_arm, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        aud["base_median_norm_wcols"] = float(np.median(np.linalg.norm(base[:, wc], axis=1)))
        aud["base_cos_to_incumbent_base_median"] = float(np.nanmedian(O.row_cos(base[:, wc], st["whitened"][:, wc])))
        if dop is not None and not dop.is_incumbent:
            wt0 = dop.evals ** (-GAMMA / 2.0)
            wt1 = dop.weights(GAMMA, LAM)
            ratio = wt1 / wt0
            aud["weight_energy_kept"] = float((wt1 ** 2).sum() / (wt0 ** 2).sum())
            aud["mean_damping"] = float(ratio.mean())
            aud["whitened_energy_kept"] = float((e_wht * ratio ** 2).sum() / e_wht.sum())
            aud["damping_by_decile"] = [float(np.median(ratio[order][k * op.n_cols // 10:(k + 1) * op.n_cols // 10])) for k in range(10)]
        row = dict(fold=held, arm=name, n_targets=n, kappa=(None if dop is None else dop.kappa), floor_c=(None if dop is None else dop.floor_c),
                   label=(None if dop is None else dop.label),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**aud, **(extra or {})), seconds=round(time.time() - t1, 1))
        msg = (f"  {name:22s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f} "
               f"dn {m['delivered_norm']:.4f} base~inc {aud['base_cos_to_incumbent_base_median']:.4f} ({row['seconds']}s)")
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            ok = d < TOL
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(ok))
            ka_log.append(dict(label=f"{name} == {ka} on {KEYS}", got=float(e["pds"]), want=float(ka_ref["pds"]), abs_diff=float(d), tol=TOL, ok=bool(ok)))
            msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({d:.1e})"
            if not ok:
                failures.append(f"{name} vs {ka}: {d:.3e}")
        rows.append(row)
        P(msg)
        return row

    # KA-1: kappa = 0 through DampedOperator == the staged base bit for bit == 263's cd4svd5 row
    d0 = O.DampedOperator.from_operator(op, s, 0.0, "incumbent (kappa 0)")
    base0 = d0.apply(st["pre"] * ALPHA, wc, GAMMA, LAM)
    check("KA-1a DampedOperator(kappa 0).apply == plolo.Operator.apply (bit-identical base)", float(np.abs(base0 - st["whitened"]).max()), 0.0, 0.0)
    r_base = run("cd4svd5", d0, ka="263 cd4svd5", ka_ref=ref263["cd4svd5"], extra=dict(control="the incumbent: the fold operator, kappa 0"))
    check("KA-1b cd4svd5 per-target pds bit-identical to 263", float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
    run("STACK", None, parts=[(st["cd4"], st["frac_cd4"]), (st["union"], st["frac_x"])], ka="257 STACK (h1union)", ka_ref=ref257["STACK"],
        extra=dict(control="h1union's stack: the retrodiction row for the pricer"))
    check("STACK nmae_wald89 == 257", rows[-1]["nmae_wald89"], ref257["STACK"]["nmae_wald89"], TOL)

    # the arms
    for kappa in O.KAPPAS + O.KAPPAS_EXTRA:
        run(f"wiener:k{kappa:g}", O.DampedOperator.from_operator(op, s, kappa, f"wiener kappa {kappa:g}"), extra=dict(arm="A: diagonal Wiener damping in the fitted eigenbasis"))
    for kappa in O.KAPPAS:
        run(f"wiener:k{kappa:g}:shuf", O.DampedOperator.from_operator(op, s_shuf, kappa, f"shuffled s, kappa {kappa:g}"), extra=dict(control="s_i permuted across directions"))
        run(f"wiener:k{kappa:g}:mirror", O.DampedOperator.from_operator(op, s_mirror, kappa, f"mirrored s, kappa {kappa:g}"), extra=dict(control="s profile reversed along the eigenvalue rank"))
    run("hard:0.5", O.DampedOperator.from_operator(op, s_hard, 1.0, "hard truncation s >= 0.5"), extra=dict(comparison="directions with s < 0.5 removed (truncation reading)"))
    run("hard:0.5:mirror", O.DampedOperator.from_operator(op, O.mirrored_s(s_hard, op.evals), 1.0, "mirrored hard truncation"), extra=dict(control="the same count of directions removed at the top"))
    # A': the eigenvalue floor at the measured noise level (soft correction of the same finding)
    for c in O.FLOOR_C:
        run(f"floor:c{c:g}", O.DampedOperator.from_operator(op, label=f"noise floor c {c:g}", n=n_i, floor_c=c), extra=dict(arm="A': lambda' = max(lambda, c n): whiten by the out-of-sample variance"))
    run("floor:c1:shuf", O.DampedOperator.from_operator(op, label="noise floor, n shuffled", n=n_shuf, floor_c=1.0), extra=dict(control="n_i permuted across directions"))
    run("floor:c1:mirror", O.DampedOperator.from_operator(op, label="noise floor, n mirrored", n=n_mirror, floor_c=1.0), extra=dict(control="n profile reversed along the eigenvalue rank"))

    # ------------------------------------------------------------------ save
    meta = dict(fold=held, n_targets=n, targets=targets, donors=donors, donor_names=st["donors"], op272_sha256=op_sha,
                sweep_sha256=O.sha256_file(Path(__file__)), operator=dict(name=op.name, sha256=op.sha256, n_cols=op.n_cols, basis_rows=op_meta["basis_rows_by_source"]),
                half_rows={line: dict(sha256=H[line]["sha256"], raw_median_real=H[line]["raw_median_real"]) for line in donors},
                anti_correlation=anti, spectrum=rep, frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"],
                ka=ka_log, ka_pass=sum(1 for k in ka_log if k["ok"]), ka_total=len(ka_log), failures=failures,
                n_rows=len(rows), seconds=round(time.time() - t0, 1), peak_rss_gib=Y.peak_gib(), fold_cache=st["fold"]["audit"].get("cache"),
                reserve_note="the orion_hek293t half-B TRUTHS (267's reserve) were not read; every arm is scored on the fold's cache truth")
    (WORK / f"sweep272-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=str))
    ref_rel = np.asarray(np.load(E267 / f"work/half/truth-{held}.npz", allow_pickle=True)["rel|batch|cos_resid"], dtype=np.float64)
    np.savez_compressed(WORK / f"pertarget272-{held}.npz", targets=np.asarray(targets), avail=avail.astype(np.int8), ref_rel=ref_rel,
                        base262_cos_truth=np.asarray(np.load(E263 / f"work/cov263-{held}.npz")["base262_cos_truth"], dtype=np.float64), **per_target)
    np.savez_compressed(WORK / f"identity272-{held}.npz", **identity)
    np.savez_compressed(WORK / f"spectrum272-{held}.npz", evals=op.evals, n=n_i, s=s, s_shuf=s_shuf, s_mirror=s_mirror, n_shuf=n_shuf, n_mirror=n_mirror, e_raw=e_raw, e_wht=e_wht,
                        **{f"n|{line}": per_donor_n[j] for j, line in enumerate(donors)})
    (WORK / f"sweep272-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{held}] done: {len(rows)} rows, KA {meta['ka_pass']}/{meta['ka_total']}, failures {failures}, {(time.time() - t0) / 60:.1f} min, peak {Y.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1])
