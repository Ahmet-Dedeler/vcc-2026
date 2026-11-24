"""272 -- construction C, DIAGNOSTIC ONLY: the noise spectrum of the PRODUCTION operator (data/donors/
whitener_union_pruned.npz, read through plolo.Operator's registry check) in its own eigenbasis, from the THREE
production donors' batch-parity half rows (271 work/half/rows-<line>.npz), and the whitened energy of the production
prediction rows (the kolf2 fold's cached `pre`: the three-donor mean with the PerturbMatch blend, 273 panel targets)
in the noise-dominated directions.  Also the in-sample whitened split-half cosine of each production donor through
the production operator (271's diagnostic, which was run on the fold operators only).  Nothing is written outside
work/; no production file is touched.

    PYTHONDONTWRITEBYTECODE=1 uv run python prod272.py
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
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))
import plolo  # noqa: E402
import op272 as O  # noqa: E402

plolo.WORK = WORK / "plolo"
E062, E064 = plolo.E062, plolo.E064
OFFICIAL, PANEL, PANEL_SET = plolo.OFFICIAL, list(plolo.PANEL), plolo.PANEL_SET
ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main() -> None:
    t0 = time.time()
    P(f"[production] op272.py sha256 {O.sha256_self()}")
    banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in DONORS]
    wc = E064.operator_columns(banks, OFFICIAL, PANEL_SET)
    axis = plolo.axis_sha256(OFFICIAL[i] for i in wc)
    op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz", name="union_pruned")
    reg = plolo.OPERATOR_REGISTRY["union_pruned"]
    assert axis == reg["axis_sha256"], (axis, reg["axis_sha256"])
    meta = json.loads((plolo.PROD_OPS / "whitener_union_pruned.json").read_text())
    P(f"[production] operator {op.name} {op.sha256[:12]} n_cols {op.n_cols} axis {axis[:12]} (registry OK); basis rows kept {meta['rows_kept_per_donor']} "
      f"(pruned {meta['rows_dropped_per_donor']}), {time.time() - t0:.0f}s")

    H = {line: O.load_half_rows(line, PANEL) for line in DONORS}
    anti = {}
    hf = {}
    for line in DONORS:
        h = H[line]
        ok = h["ok"]
        WA = op.apply(h["A"] * ALPHA, wc, GAMMA, LAM)[:, wc]
        WB = op.apply(h["B"] * ALPHA, wc, GAMMA, LAM)[:, wc]
        c = O.row_cos(WA[ok], WB[ok])
        pi = O.derangement(int(ok.sum()), 1)
        null = O.row_cos(WA[ok], WB[ok][pi])
        craw = O.row_cos(h["A"][ok][:, wc], h["B"][ok][:, wc])
        anti[line] = dict(n=int(ok.sum()), wht_median=float(np.nanmedian(c)), wht_share_le0=float(np.nanmean(c <= 0)), wht_p10=float(np.nanpercentile(c, 10)),
                          wht_p90=float(np.nanpercentile(c, 90)), wht_null_median=float(np.nanmedian(null)), wht_null_sd=float(np.nanstd(null)),
                          raw_median=float(np.nanmedian(craw)), raw_share_le0=float(np.nanmean(craw <= 0)))
        hf[line] = O.half_factor_check(h, wc)
        P(f"  {line:14s} n {anti[line]['n']}: in-sample whitened split-half median {anti[line]['wht_median']:+.4f} (share <= 0 {anti[line]['wht_share_le0']:.3f}, "
          f"null {anti[line]['wht_null_median']:+.4f} sd {anti[line]['wht_null_sd']:.4f}); raw {anti[line]['raw_median']:+.4f}; half-factor ratio {hf[line]['ratio_d_over_implied']:.4f}; "
          f"signal share {hf[line]['signal_share_of_real_energy']:.4f}")
        del WA, WB

    raw_medians = {line: H[line]["raw_median_real"] for line in DONORS}
    unit, renorm = O.basis_unit_factor(raw_medians)
    weights = [float(meta["rows_kept_per_donor"][line]) for line in DONORS]
    D_list = [O.half_differences(H[line])[:, wc] * unit for line in DONORS]
    n_i, s, per = O.noise_spectrum(op.evals, op.evecs, D_list, weights)
    # the production prediction rows: the kolf2 fold's cached pre (three-donor mean, PM blend, 273 targets)
    z = np.load(WORK / "plolo/fold-kolf2.npz", allow_pickle=True)
    assert np.array_equal(np.asarray(z["wcols"], dtype=np.int64), wc), "kolf2 fold wcols != production wcols"
    pre = np.asarray(z["pre"], dtype=np.float64)
    e_raw, e_wht = O.whitened_energy(pre[:, wc] * ALPHA, op.evecs, op.evals, GAMMA)
    Pfull = np.concatenate([H[line]["real"][H[line]["ok"]][:, wc] * unit for line in DONORS])
    e_raw_full, _ = O.whitened_energy(Pfull, op.evecs, op.evals, GAMMA)
    lam_hat = e_raw_full / Pfull.shape[0]
    order = np.argsort(op.evals)[::-1]
    ratio_dec = [float(np.median((lam_hat / op.evals)[order][k * op.n_cols // 10:(k + 1) * op.n_cols // 10])) for k in range(10)]
    rep = O.spectrum_report(op.evals, s, n_i, e_raw, e_wht, GAMMA)
    rep.update(half_factor=O.HALF_FACTOR, half_factor_check=hf, unit_factor=unit, renorm=renorm, basis_weights=dict(zip(DONORS, weights)),
               lambda_hat_over_lambda_by_decile=ratio_dec, prediction_rows="kolf2 fold cache pre (273 targets, three-donor mean, PM blend)",
               n_zero=int((n_i == 0).sum()))
    good = (n_i > 0) & (op.evals > op.evals.max() * 1e-11)
    rep["spearman_log_n_log_lambda"] = float(np.corrcoef(np.argsort(np.argsort(np.log(n_i[good]))), np.argsort(np.argsort(np.log(op.evals[good]))))[0, 1])
    P(f"  units check: median (panel-row variance / lambda) by rank decile {[round(x, 3) for x in ratio_dec]} (unit factor {unit:.6f}, renorm {renorm})")
    P(f"  spectrum: p {rep['p']}, floored {rep['n_lambda_floored']}, n_i = 0 on {rep['n_zero']}; noise-dominated (s < 0.5) {rep['n_noise_dominated']} ({100*rep['share_noise_dominated']:.1f} %), "
      f"s = 0 on {rep['n_s_zero']} ({100*rep['share_s_zero']:.1f} %); whitened energy of the production rows in noise-dominated directions "
      f"{100*rep['whitened_energy_share_noise_dominated']:.1f} % (s = 0: {100*rep['whitened_energy_share_s_zero']:.1f} %; raw energy there {100*rep['raw_energy_share_noise_dominated']:.1f} %); "
      f"rank of the noise-dominated set median {rep['rank_of_noise_dominated']['median']:.0f} (p10 {rep['rank_of_noise_dominated']['p10']:.0f}, p90 {rep['rank_of_noise_dominated']['p90']:.0f}); "
      f"mean s {rep['mean_s']:.3f}, energy-weighted {rep['energy_weighted_mean_s']:.3f}; sum n / sum lambda {rep['total_noise_over_total_lambda']:.3f}; rank corr(log n, log lambda) {rep['spearman_log_n_log_lambda']:+.3f}")
    for r in rep["deciles"]:
        P(f"    decile {r['decile']}: lambda [{r['lambda_min']:.2e}, {r['lambda_max']:.2e}] amp median {r['amplification_median']:.1f}x; mean s {r['mean_s']:.3f} "
          f"(share < 0.5 {r['share_s_lt_half']:.2f}, = 0 {r['share_s_zero']:.2f}); n/lambda median {r['median_noise_over_lambda']:.3f}; "
          f"raw energy {100*r['raw_energy_share']:.1f} %, whitened energy {100*r['whitened_energy_share']:.1f} %")
    # what the two constructions would do to the production operator's weights (no file written)
    w0 = op.evals ** (-GAMMA / 2.0)
    fam = {}
    for kappa in O.KAPPAS:
        w1 = w0 * s ** kappa
        fam[f"wiener:k{kappa:g}"] = dict(weight_energy_kept=float((w1 ** 2).sum() / (w0 ** 2).sum()), whitened_energy_kept=float((e_wht * (w1 / w0) ** 2).sum() / e_wht.sum()))
    for c in O.FLOOR_C:
        w1 = np.maximum(op.evals, c * n_i) ** (-GAMMA / 2.0)
        fam[f"floor:c{c:g}"] = dict(weight_energy_kept=float((w1 ** 2).sum() / (w0 ** 2).sum()), whitened_energy_kept=float((e_wht * (w1 / w0) ** 2).sum() / e_wht.sum()),
                                    damping_by_decile=[float(np.median((w1 / w0)[order][k * op.n_cols // 10:(k + 1) * op.n_cols // 10])) for k in range(10)])
    for k, v in fam.items():
        P(f"  {k:12s} whitened energy kept {100*v['whitened_energy_kept']:.1f} %" + (f"; damping by decile {[round(x, 3) for x in v['damping_by_decile']]}" if "damping_by_decile" in v else ""))
    out = dict(operator=dict(name=op.name, sha256=op.sha256, path=str((plolo.PROD_OPS / "whitener_union_pruned.npz").relative_to(ROOT)), axis_sha256=axis, n_cols=op.n_cols, metadata=meta),
               donors=list(DONORS), half_rows={line: dict(sha256=H[line]["sha256"], raw_median_real=H[line]["raw_median_real"]) for line in DONORS},
               anti_correlation=anti, spectrum=rep, families=fam, op272_sha256=O.sha256_self(), seconds=round(time.time() - t0, 1), peak_rss_gib=O.peak_gib())
    (WORK / "prod272.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / "spectrum272-production.npz", evals=op.evals, n=n_i, s=s, e_raw=e_raw, e_wht=e_wht, **{f"n|{line}": per[j] for j, line in enumerate(DONORS)})
    (WORK / "prod272.log").write_text("\n".join(LOG) + "\n")
    P(f"[production] done {(time.time() - t0) / 60:.1f} min, peak {O.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main()
