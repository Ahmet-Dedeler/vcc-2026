"""273 items 1, 3, 4 -- THE INSTRUMENT (union OFF) and the three re-reads, on h1-test scored with the scorer's own statistic.

Base rows: 170's staged `cand_h1` (the builder's alpha-scaled pre-whitening rows for the 85 targets: three donors,
PerturbMatch 0.5 on K562's shared columns; gated on a15's shipped sha256) whitened by the PRODUCTION operator
`union_pruned` through `plolo.Operator.apply` (gamma 0.65, lam 0, on-axis preserve) -- exactly what 221 whitened and
emitted for the cloud scorer.  Residual: CD4 only (the H1 union is never touched).  Two CD4 objects:
  * FROZEN -- the shipped banks (panel rows only): 4 complete + 1 partial of the 85 targets; 221's cloud-scored
    candidates carried this object at the CHALLENGE-frame fraction (KA-2a reproduces 170's h1 effect at f = 0.10).
  * H1 -- the same source (the pinned 16.8 GB GWCD4i.DE_stats.h5ad) re-extracted for the 85 targets with 047's
    `extract_cd4.py` (the October D/E/F path), prepared by `E062.prepare_cd4_mean` (gene_norm) + the complete-mean3
    policy, re-solved on THIS frame at f = 0.34 by production's all-rows rule.  THE INCUMBENT of this instrument is
    base + H1-CD4 deflated k = 5 (the builder's `deflate_cd4_residual` == lib246.deflate svd) -- `cd4svd5` with the
    union off.
Arms (each CD4 residual re-solved against its own base):
  (a) 272's A' eigenvalue floor max(lambda, c n) at c in {0.25, 1, 4} on `union_pruned`, n = 272's production noise
      spectrum (spectrum272-production.npz), + its shuffled / mirrored controls at c = 1; 272's Wiener s^kappa at kappa
      0.5 / 1; each also without CD4.
  (b) 271's EB row shrinkage `ebrow` (and relw / reversed / shuffled controls) from the donors' half rows (halfrows273.py).
  (c) 249's deflation: H1-CD4 off vs svd5 (fitted on the day's rows) vs the pinned V5_abc basis (October's rule) vs the
      energy-matched random subspace (246) vs the target-shuffled residual (218).
KNOWN-ANSWER ROWS (exit 1 on failure): KA-2a (170's h1 effect), KA-3 (DampedOperator c = 0 == Operator.apply), KA-0
(the per-donor rebuild == cand_h1), KA-cd4 (the re-extracted bank rows == the frozen bank rows on the shared targets).

    PYTHONDONTWRITEBYTECODE=1 uv run python instrument273.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib273 as X  # noqa: E402

sys.path.insert(0, str(X.E271))
import est271 as EST  # noqa: E402

plolo, L, O, E062, E064 = X.plolo, X.L, X.O, X.E062, X.E064
OFFICIAL, PANEL_SET = X.OFFICIAL, X.PANEL_SET
GAMMA, LAM, ALPHA, CD4_F = X.GAMMA, X.LAM, X.ALPHA, X.CD4_F
LOG: list[str] = []
FAIL: list[str] = []
KA: list[dict] = []
CD4_H1 = X.WORK / "cd4/cd4_compact_h1"
CD4_FROZEN_ROOT = X.ROOT / "data/external/cd4t/cloud-compact/cd4_compact"


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def ka(label: str, value: float, tol: float, **extra) -> bool:
    ok = bool(value <= tol)
    KA.append(dict(label=label, value=float(value), tol=tol, ok=ok, **extra))
    P(f"  {'PASS' if ok else 'FAIL'} {label}: {value:.3e} (tol {tol:.0e})")
    if not ok:
        FAIL.append(label)
    return ok


def read_results_csv(p: Path, metric: str = "pds_cosine") -> dict[str, float]:
    out = {}
    with open(p) as fh:
        for row in csv.DictReader(fh):
            if row["metric"] == metric:
                out[row["perturbation"]] = float(row["value"])
    return out


def rand_energy_matched(res, target_frac: float, seed: int, kmax: int):
    """sweep246b.rand_energy_matched, verbatim (kmax = the carried-row count - 1)."""
    lo, hi = 1, kmax
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        d, aud = L.deflate(res, "rand", mid, seed=seed)
        if aud["energy_removed_frac"] >= target_frac:
            best = (d, aud, mid)
            hi = mid - 1
        else:
            lo = mid + 1
    if best is None:
        d, aud = L.deflate(res, "rand", kmax, seed=seed)
        best = (d, aud, kmax)
    return best


def main() -> None:
    t0 = time.time()
    st = X.staged170()
    tr = X.load_truth()
    targets = st["h1_targets"]
    assert np.array_equal(targets, tr["targets"])
    n = targets.size
    cols, cand = st["cols"], st["cand_h1"]
    tset = set(targets.tolist())
    P(f"273 instrument273.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}; {n} targets, {cols.size} operator columns, "
      f"170 pre.npz {st['sha']['pre'][:12]} cd4.npz {st['sha']['cd4'][:12]}; truth273 h1 {tr['sha'][:12]}")

    # ------------------------------------------------------------------------------------------ the operator + base
    op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz", name="union_pruned")
    base = op.apply(cand, cols, GAMMA, LAM)
    P(f"  base: median row norm on cols before {np.median(np.linalg.norm(cand[:, cols], axis=1)):.6f} -> after {np.median(np.linalg.norm(base[:, cols], axis=1)):.6f} "
      f"(axis check {op.axis_checked}); {time.time() - t0:.0f}s, peak {X.peak_gib():.2f} GiB")

    # ------------------------------------------------------------------------------------------ KA-2a: 170's h1 effect
    e170 = np.load(X.E170 / "eff/h1_g0.65.npz")["effect"]
    f010, aud010 = E062.add_post_whiten_residual(base, st["cd4_frozen"], X.CHAL_FRACTION_F010, cols)
    diff = f010 - e170
    fro = float(np.linalg.norm(diff) / np.linalg.norm(e170))
    rn_a, rn_b = np.linalg.norm(f010, axis=1), np.linalg.norm(e170, axis=1)
    min_cos = float(((f010 * e170).sum(1) / np.maximum(rn_a * rn_b, 1e-300)).min())
    med = float(np.median(rn_b))
    P("KA-2a: my composition at f = 0.10 (cand_h1 -> Operator.apply -> add_post_whiten_residual at the challenge-frame fraction) vs 170's h1_g0.65.npz")
    ka("KA-2a Frobenius ratio ||d||/||E||", fro, 1e-12)
    ka("KA-2a 1 - min row cosine", 1.0 - min_cos, 1e-12)
    ka("KA-2a max|d| / median row norm", float(np.abs(diff).max()) / med, 1e-9, median_row_norm=med)

    arms: dict[str, dict] = {}
    audits: dict[str, dict] = {}

    def add(name: str, final: np.ndarray, **meta) -> dict:
        e = X.score_h1(final, tr)
        arms[name] = e
        audits[name] = dict(meta, median_row_norm=float(np.median(np.linalg.norm(final, axis=1))))
        P(f"  [{name:22s}] pds a23 {e['pds']:.6f}  a15 {float(e['pds_a15'].mean()):.6f}  eff {float(e['pds_eff'].mean()):.6f} | own cos a23 {e['cos']:.5f} "
          f"cos_agg {e['a23_cos_agg']:.5f} | delivered {e['a23_delivered_norm']:.4f}  ({time.time() - t0:.0f}s)")
        return e

    add("BASE", base)
    add("cd4frz:f010", f010, fraction=X.CHAL_FRACTION_F010, sizing="challenge frame (221)", **{k: float(v) for k, v in aud010.items() if k != "fraction"})
    f034, aud034 = E062.add_post_whiten_residual(base, st["cd4_frozen"], X.CHAL_FRACTION_F034, cols)
    add("cd4frz:f034", f034, fraction=X.CHAL_FRACTION_F034, sizing="challenge frame (221)", **{k: float(v) for k, v in aud034.items() if k != "fraction"})

    # the RENDERED composition against the cloud's read of the EMITTED counts (the emission gap; not a to-the-digit KA)
    cloud = {"a15": read_results_csv(X.E221 / "ladder6-20260904/a15/run/results.csv"),
             "a15_seed1": read_results_csv(X.E221 / "ladder7-20260904/a15_seed1/run/results.csv"),
             "cd4f034": read_results_csv(X.E221 / "ladder8-20260904/cd4f034/run/results.csv")}
    from scipy.stats import spearmanr
    emission = {}
    for arm, key in (("cd4frz:f010", "a15"), ("cd4frz:f010", "a15_seed1"), ("cd4frz:f034", "cd4f034")):
        c = np.asarray([cloud[key][t] for t in targets])
        for tag in ("a23", "a15", "eff"):
            p = arms[arm][f"pds_{tag}"]
            emission[f"{arm}|{tag}|vs {key}"] = dict(local_mean=float(p.mean()), cloud_mean=float(c.mean()), diff_mean=float(p.mean() - c.mean()),
                                                     max_abs_per_target=float(np.abs(p - c).max()), mean_abs_per_target=float(np.abs(p - c).mean()),
                                                     n_equal=int((np.abs(p - c) < 1e-12).sum()), spearman=float(spearmanr(p, c).statistic))
        e = emission[f"{arm}|a23|vs {key}"]
        P(f"  emission gap {arm} (rendered, a23) vs cloud {key} (emitted, scored): mean {e['local_mean']:.6f} vs {e['cloud_mean']:.6f} (d {e['diff_mean']:+.6f}); "
          f"per target max|d| {e['max_abs_per_target']:.4f}, mean|d| {e['mean_abs_per_target']:.4f}, equal on {e['n_equal']}/85, Spearman {e['spearman']:+.4f}")
    seed_floor = float(np.abs(np.asarray([cloud["a15"][t] for t in targets]) - np.asarray([cloud["a15_seed1"][t] for t in targets])).mean())
    emission["cloud a15 vs a15_seed1 per-target mean|d| (the emission-seed floor)"] = seed_floor
    P(f"  the cloud's own seed floor: a15 vs a15_seed1 per-target mean|d| {seed_floor:.4f}, mean d {np.mean([cloud['a15'][t] - cloud['a15_seed1'][t] for t in targets]):+.6f}")

    # ------------------------------------------------------------------------------------------ frozen CD4 deflations
    V, vinfo = X.load_v5_abc()
    cdf = st["cd4_frozen"]
    d5f, aud5f = L.deflate(cdf, "svd", 5)
    audits["cd4frz:svd5"] = dict(aud5f, note="fitted k = 5 on the 5 carried rows (4 complete + 1 zero partial): rank <= 4 after centring, the projector removes everything")
    P(f"  frozen CD4 deflated k = 5 on {aud5f['n_carried']} carried rows: energy removed {aud5f['energy_removed_frac']:.6f} -> "
      f"{'DEGENERATE (nothing left to size)' if aud5f['energy_removed_frac'] > 1 - 1e-9 else 'sizable'}")
    dVf, audVf = X.deflate_pinned(cdf, V)
    try:
        frVf, puVf = X.solve_fraction(base, dVf, CD4_F)
        sizing_frz = "all_rows"
    except plolo.SilentZeroError:
        frVf, puVf = X.solve_fraction(base, dVf, CD4_F, sizing="carrying")
        sizing_frz = "carrying (all-rows median is 0 at 5/85 coverage; October's rule)"
    add("cd4frz:V5", X.compose(base, [(dVf, frVf)], cols), fraction=frVf, per_unit=puVf, sizing=sizing_frz, **audVf)
    try:
        frf, puf = X.solve_fraction(base, cdf, CD4_F, sizing="carrying")
        add("cd4frz:off:carry", X.compose(base, [(cdf, frf)], cols), fraction=frf, per_unit=puf, sizing="carrying")
    except plolo.SilentZeroError as exc:
        P(f"  cd4frz:off carrying-row sizing refused: {exc}")

    # ------------------------------------------------------------------------------------------ H1 CD4 (re-extracted)
    banks_h1 = [E062.Bank.load(nm, CD4_H1 / nm, require_cd4_scale=True) for nm in X.CD4_NAMES]
    banks_frz = [E062.Bank.load(nm, CD4_FROZEN_ROOT / nm, require_cd4_scale=True) for nm in X.CD4_NAMES]
    # KA-cd4: on the targets both extractions carry, the delta rows (ln fold changes) are identical on the common genes
    worst = 0.0
    n_shared = 0
    for bh, bf in zip(banks_h1, banks_frz):
        gh = {g: i for i, g in enumerate(bh.genes)}
        common = [g for g in bf.genes if g in gh]
        ih = np.asarray([gh[g] for g in common]); i_f = np.asarray([i for i, g in enumerate(bf.genes) if g in gh])
        dh, df = np.load(bh.delta_path, mmap_mode="r"), np.load(bf.delta_path, mmap_mode="r")
        sh = [t for t in bh.perts if t in bf.pert_pos]
        n_shared += len(sh)
        for t in sh:
            worst = max(worst, float(np.abs(np.asarray(dh[bh.pert_pos[t], ih]) - np.asarray(df[bf.pert_pos[t], i_f])).max()))
        P(f"  CD4 bank {bh.name}: h1 extraction {bh.shape[0]} perts x {bh.shape[1]} genes; frozen {bf.shape[0]} x {bf.shape[1]}; common genes {len(common)}; shared targets {len(sh)}")
    ka(f"KA-cd4 re-extracted CD4 delta rows == frozen bank rows on the {n_shared} shared (condition, target) rows, common genes", worst, 0.0)
    prep = E062.prepare_cd4_mean(banks_h1, list(targets), OFFICIAL, tset, "gene_norm", chunk_rows=256, max_rss_gib=40.0)
    cd4h1 = X.complete_mean3(prep)
    cnt = np.asarray(cd4h1.target_condition_count)
    P(f"  H1 CD4: condition count histogram {np.bincount(cnt, minlength=4).tolist()} (complete = {int((cnt == 3).sum())}/{n}); details "
      f"{[dict((k, v) for k, v in d.items() if k in ('condition', 'covered_targets', 'background_perturbations')) for d in prep.condition_details]}")
    fr_off, pu_off = X.solve_fraction(base, cd4h1, CD4_F)
    d5, aud5 = L.deflate(cd4h1, "svd", 5)
    fr5, pu5 = X.solve_fraction(base, d5, CD4_F)
    dV, audV = X.deflate_pinned(cd4h1, V)
    frV, puV = X.solve_fraction(base, dV, CD4_F)
    drand, audr, kr = rand_energy_matched(cd4h1, aud5["energy_removed_frac"], X.RAND_SEED, int((cnt > 0).sum()) - 1)
    frr, pur = X.solve_fraction(base, drand, CD4_F)
    dshuf = L.shuffled(cd4h1, seed=L.SHUF_SEED)
    frs, pus = X.solve_fraction(base, dshuf, CD4_F)
    P(f"  H1 CD4 sizing at f = 0.34 (all-rows rule): off {fr_off:.5f} (per unit {pu_off:.5f}); svd5 {fr5:.5f} (energy removed {aud5['energy_removed_frac']:.4f}, top-5 share of centred "
      f"{aud5['top_k_energy_frac_of_centred']:.4f}); V5_abc {frV:.5f} (energy removed {audV['energy_removed_frac']:.4f}, captured {audV['captured_frac_of_centred']:.4f}); "
      f"random {kr} dims {frr:.5f} (energy removed {audr['energy_removed_frac']:.4f}); shuffled {frs:.5f}")
    add("cd4h1:off", X.compose(base, [(cd4h1, fr_off)], cols), fraction=fr_off, per_unit=pu_off)
    add("INCUMBENT", X.compose(base, [(d5, fr5)], cols), fraction=fr5, per_unit=pu5, **aud5, note="base + H1 CD4 deflated k=5 re-solved f=0.34 (cd4svd5, union OFF)")
    add("cd4h1:V5", X.compose(base, [(dV, frV)], cols), fraction=frV, per_unit=puV, **audV)
    add("cd4h1:randE5", X.compose(base, [(drand, frr)], cols), fraction=frr, per_unit=pur, **audr)
    add("cd4h1:shuf", X.compose(base, [(dshuf, frs)], cols), fraction=frs, per_unit=pus)
    # the frozen object with the H1-frame sizing rule, for the record (what "re-solved" means at 5/85 coverage)
    P(f"  frozen CD4 at f = 0.34: challenge-frame fraction {X.CHAL_FRACTION_F034:.5f} vs carrying-row re-solve on this frame "
      f"{X.solve_fraction(base, cdf, CD4_F, sizing='carrying')[0]:.5f}")

    # ------------------------------------------------------------------------------------------ (a) the operator
    spz = np.load(X.E272 / "work/spectrum272-production.npz")
    assert np.array_equal(spz["evals"], op.evals), "272's production spectrum is not this operator's"
    nvec, svec = np.asarray(spz["n"], dtype=np.float64), np.asarray(spz["s"], dtype=np.float64)
    d0 = O.DampedOperator.from_operator(op)
    ka("KA-3 DampedOperator (c = 0, kappa = 0) == plolo.Operator.apply, max|d|", float(np.abs(d0.apply(cand, cols, GAMMA, LAM) - base).max()), 0.0)

    def op_arm(name: str, dop, with_cd4: bool = True):
        b = dop.apply(cand, cols, GAMMA, LAM)
        bc = X.row_cos(b[:, cols], base[:, cols])
        if with_cd4:
            fr, pu = X.solve_fraction(b, d5, CD4_F)
            add(name, X.compose(b, [(d5, fr)], cols), fraction=fr, per_unit=pu, base_cos_to_incumbent_base=float(np.median(bc)),
                whitened_energy_kept=float(((dop.weights(GAMMA, LAM) / op.evals ** (-GAMMA / 2)) ** 2 * spz["e_wht"]).sum() / spz["e_wht"].sum()))
        else:
            add(name, b, base_cos_to_incumbent_base=float(np.median(bc)))

    for c in (0.25, 1.0, 4.0):
        op_arm(f"floor:c{c:g}", O.DampedOperator.from_operator(op, n=nvec, floor_c=c, label=f"floor c={c}"))
    op_arm("floor:c1:noCD4", O.DampedOperator.from_operator(op, n=nvec, floor_c=1.0), with_cd4=False)
    op_arm("floor:c1:shuf", O.DampedOperator.from_operator(op, n=O.shuffled_s(nvec, O.SHUF_SEED), floor_c=1.0))
    op_arm("floor:c1:mirror", O.DampedOperator.from_operator(op, n=O.mirrored_s(nvec, op.evals), floor_c=1.0))
    for k in (0.5, 1.0):
        op_arm(f"wiener:k{k:g}", O.DampedOperator.from_operator(op, s=svec, kappa=k, label=f"wiener kappa={k}"))
    op_arm("wiener:k1:noCD4", O.DampedOperator.from_operator(op, s=svec, kappa=1.0), with_cd4=False)
    op_arm("wiener:k1:shuf", O.DampedOperator.from_operator(op, s=O.shuffled_s(svec, O.SHUF_SEED), kappa=1.0))

    # ------------------------------------------------------------------------------------------ (b) the donors
    rows = {ln: np.load(X.WORK / "half" / f"rows273-{ln}.npz", allow_pickle=True) for ln in X.DONORS}
    for ln in X.DONORS:
        assert np.array_equal(rows[ln]["targets"].astype(str), targets)
    banks = {ln: E062.Bank.load(ln, plolo.DONOR_PATHS[ln]) for ln in X.DONORS}
    pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
    no_panel = np.zeros(0, dtype=np.int64)
    pm_prep = E062.prepare_main_source(pm_bank, list(targets), OFFICIAL, no_panel, chunk_rows=256, max_rss_gib=40.0)
    real = {ln: np.asarray(rows[ln]["real"], dtype=np.float64) for ln in X.DONORS}
    avail = np.stack([np.asarray(rows[ln]["avail"], bool) for ln in X.DONORS])          # (J, n)
    # the PerturbMatch 0.5 blend on K562's shared columns (E064.shared_k562_pm_columns excludes the A/B/C panel from the
    # blend; the builder's `apply_perturbmatch_k562_replacement` may not -- KA-0 tries both and reports which reproduces)
    both = avail[0] & np.asarray(pm_prep.target_available, bool)
    variants = {}
    for excl, panel in (("panel excluded", PANEL_SET), ("panel included", set())):
        shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, panel)
        k = real["k562_gwps"].copy()
        r_ = np.flatnonzero(both)
        k[np.ix_(r_, shared)] = (1 - plolo.PM_BLEND) * real["k562_gwps"][np.ix_(r_, shared)] + plolo.PM_BLEND * np.asarray(pm_prep.effects)[np.ix_(r_, shared)]
        src = [E062.PreparedSource(name="k562_gwps_pm50", effects=k, target_available=avail[0].copy(), mean_on_axis=np.asarray(rows["k562_gwps"]["mean_on_axis"]),
                                   raw_median_norm=float(rows["k562_gwps"]["raw_median_real"]), mapped_gene_count=0, duplicate_symbols=())]
        for j, ln in enumerate(X.DONORS[1:], start=1):
            src.append(E062.PreparedSource(name=ln, effects=real[ln], target_available=avail[j].copy(), mean_on_axis=np.asarray(rows[ln]["mean_on_axis"]),
                                           raw_median_norm=float(rows[ln]["raw_median_real"]), mapped_gene_count=0, duplicate_symbols=()))
        pre64, votes = E062.build_baseline(src)
        recon = pre64 * ALPHA
        d = np.abs(recon - cand)
        variants[excl] = dict(max_abs=float(d.max()), rel_to_max=float(d.max() / np.abs(cand).max()),
                              max_abs_off_panel=float(d[:, np.setdiff1d(np.arange(len(OFFICIAL)), X.PANEL_POS)].max()), shared_cols=int(shared.size), k_rows=k, pre64=pre64)
        P(f"  KA-0 rebuild of cand_h1 from the per-donor rows (PM blend {excl}, {shared.size} shared cols, {int(both.sum())} blended rows): max|d| {d.max():.3e} "
          f"(rel {variants[excl]['rel_to_max']:.1e}); off the A/B/C panel columns {variants[excl]['max_abs_off_panel']:.3e}")
    best = min(variants, key=lambda k_: variants[k_]["max_abs"])
    # H-0 (MEASURED, not a regression check): does the fold instrument's row construction (E062 prepare_main_source +
    # E064 PerturbMatch blend = plolo.build_fold's path, the rows 271 built its estimators on) reproduce the production
    # builder's rows?  Recorded with the per-row geometry; the (b) arms are read against their own mean row (R1).
    recon = variants[best]["pre64"] * ALPHA
    rc = X.row_cos(recon, cand)
    nr = np.linalg.norm(cand, axis=1) / np.linalg.norm(recon, axis=1)
    pm_genes = set(pm_bank.genes)
    pm_cols = np.asarray([i for i, g in enumerate(OFFICIAL) if g in pm_genes], dtype=np.int64)
    non_pm = np.setdiff1d(np.arange(len(OFFICIAL)), pm_cols)
    dd = np.abs(recon - cand)
    h0 = dict(label="H-0 (measured) alpha x build_baseline(E062-path donor rows, E064 PerturbMatch blend) vs the builder's cand_h1", variant=best,
              rel_to_max=variants[best]["rel_to_max"], row_cos_min=float(rc.min()), row_cos_median=float(np.median(rc)),
              rows_identical=int((dd.max(1) < 1e-9).sum()), rows_identical_no_k562=int((dd.max(1)[~avail[0]] < 1e-9).sum()), n_no_k562=int((~avail[0]).sum()),
              max_abs_on_non_pm_columns=float(dd[:, non_pm].max()), max_abs_on_pm_columns=float(dd[:, pm_cols].max()), norm_ratio_median=float(np.median(nr)),
              reproduces=bool(variants[best]["rel_to_max"] < 1e-9))
    KA.append(dict(label=h0["label"], value=h0["rel_to_max"], tol=1e-9, ok=h0["reproduces"], fatal=False,
                   **{k_: v for k_, v in h0.items() if k_ not in ("label", "rel_to_max", "reproduces")}))
    P(f"  {h0['label']}: {'REPRODUCES' if h0['reproduces'] else 'DIFFERS'} -- rows identical {h0['rows_identical']}/85 ({h0['rows_identical_no_k562']}/{h0['n_no_k562']} of the rows "
      f"without a K562 slot); row cosine min {h0['row_cos_min']:.4f} median {h0['row_cos_median']:.4f}; norm ratio median {h0['norm_ratio_median']:.4f}; "
      f"max|d| on non-PerturbMatch columns {h0['max_abs_on_non_pm_columns']:.2e} vs on PerturbMatch columns {h0['max_abs_on_pm_columns']:.3f}")
    k_rows = variants[best]["k_rows"]
    pre64 = variants[best]["pre64"]
    donor_rows = [k_rows, real["orion_hct116"], real["orion_hek293t"]]

    # reliabilities R_j on the mapped columns minus the 85 target genes (est271's definitions; the K562 halves un-blended, as 271)
    R_half = np.full((3, n), np.nan)
    ncells = np.zeros((3, n), dtype=np.int64)
    for j, ln in enumerate(X.DONORS):
        z = rows[ln]
        ok = np.asarray(z["avail_A"], bool) & np.asarray(z["avail_B"], bool)
        R_half[j] = np.where(ok, np.asarray(z["r_half_raw"], dtype=np.float64), np.nan)
        ncells[j] = np.asarray(z["ncells_full"])
    R_full = EST.spearman_brown(R_half)
    for j, ln in enumerate(X.DONORS):
        v = R_half[j][avail[j]]
        P(f"  donor {ln:14s}: {int(avail[j].sum())}/85 targets; half-row reliability median {np.nanmedian(v):+.4f} (share <= 0 {np.nanmean(v <= 0):.3f}); "
          f"SB median {np.nanmedian(R_full[j][avail[j]]):+.4f}; cells/target median {int(np.median(ncells[j][avail[j]]))}")
    two_plus = avail.sum(0) >= 2
    P(f"  targets carried by >= 2 donors: {int(two_plus.sum())}/85 (three donors: {int((avail.sum(0) == 3).sum())})")

    def combine(W: np.ndarray, *, explicit_mean: bool = False) -> np.ndarray:
        """lib263.combine: the mean row where W is NaN, else sum_j W[t, j] row_j[t] NORM-MATCHED to the mean row, float32-cast."""
        pre = pre64.copy()
        if explicit_mean:
            acc = np.zeros_like(pre); c_ = np.zeros(n, dtype=np.int64)
            for j in range(3):
                m = avail[j]; acc[m] += donor_rows[j][m]; c_ += m.astype(np.int64)
            cov = c_ > 0
            acc[cov] /= c_[cov, None]
            pre[cov] = acc[cov].astype(np.float32).astype(np.float64)
            return pre
        use = np.isfinite(W).all(1) & (avail.sum(0) >= 1)
        idx = np.flatnonzero(use)
        if idx.size:
            Wu = W[idx].copy(); Wu[~avail[:, idx].T] = 0.0
            s_ = Wu.sum(1); assert np.all(s_ > 0)
            Wu /= s_[:, None]
            acc = np.zeros((idx.size, pre.shape[1]))
            for j in range(3):
                acc += Wu[:, j][:, None] * donor_rows[j][idx]
            na, nm = np.linalg.norm(acc, axis=1), np.linalg.norm(pre[idx], axis=1)
            ok = (na > 0) & (nm > 0)
            acc[ok] *= (nm[ok] / na[ok])[:, None]
            pre[idx] = acc.astype(np.float32).astype(np.float64)
        return pre

    def donor_arm(name: str, pre_sel: np.ndarray, **meta):
        b = op.apply(pre_sel * ALPHA, cols, GAMMA, LAM)
        fr, pu = X.solve_fraction(b, d5, CD4_F)
        add(name, X.compose(b, [(d5, fr)], cols), fraction=fr, per_unit=pu, rows_changed=int((np.abs(pre_sel - pre64).max(1) > 0).sum()),
            base_cos_to_incumbent_base=float(np.median(X.row_cos(b[:, cols], base[:, cols]))), **meta)

    donor_arm("mean:rebuilt", combine(None, explicit_mean=True), note="R1: the explicit equal-weight mean of the rebuilt rows (== INCUMBENT if KA-0 holds)")
    W_eb = EST.weights(R_full, avail, "ebrow")
    donor_arm("ebrow", combine(W_eb), max_weight_median=float(np.nanmedian(np.nanmax(W_eb, axis=1))))
    W_ebs = EST.weights(EST.spearman_brown(EST.shuffle_within_donor(R_half, avail, EST.SHUF_SEED)), avail, "ebrow")
    donor_arm("ebrow:shuf", combine(W_ebs))
    W_rw = EST.weights(R_full, avail, "relw")
    donor_arm("relw", combine(W_rw), max_weight_median=float(np.nanmedian(np.nanmax(W_rw, axis=1))))
    donor_arm("relw:rev", combine(EST.weights(R_full, avail, "rev")))
    donor_arm("relw:shuf", combine(EST.weights(EST.spearman_brown(EST.shuffle_within_donor(R_half, avail, EST.SHUF_SEED)), avail, "relw")))

    # ------------------------------------------------------------------------------------------ save
    arrays = {}
    for name, e in arms.items():
        for k_, v in e.items():
            if isinstance(v, np.ndarray):
                arrays[f"{name}|{k_}"] = v
    arrays["targets"] = targets
    arrays["cd4h1_count"] = cnt
    arrays["cd4frz_count"] = np.asarray(cdf.target_condition_count)
    arrays["donor_avail"] = avail
    arrays["donor_R_half"] = R_half
    arrays["donor_ncells"] = ncells
    arrays["base_cos_eff"] = arms["BASE"]["cos_eff"]
    np.savez_compressed(X.WORK / "instrument273.npz", **arrays)
    summ = {name: X.summary(e) for name, e in arms.items()}
    out = dict(ka=KA, failures=FAIL, arms=summ, audits=audits, emission_gap=emission, v5_abc=vinfo, operator=dict(name=op.name, sha256=op.sha256, n_cols=op.n_cols),
               cd4_h1=dict(path=str(CD4_H1.relative_to(X.ROOT)), condition_details=prep.condition_details, count_hist=np.bincount(cnt, minlength=4).tolist(),
                           fractions=dict(off=fr_off, svd5=fr5, V5=frV, randE5=frr, shuf=frs), deflation=aud5, pinned=audV, random=dict(audr, k=kr)),
               cd4_frozen=dict(count_hist=np.bincount(np.asarray(cdf.target_condition_count), minlength=4).tolist(), fraction_f010=X.CHAL_FRACTION_F010,
                               fraction_f034=X.CHAL_FRACTION_F034, svd5=aud5f, V5=audVf, V5_fraction=frVf, V5_sizing=sizing_frz),
               donors=dict(avail=avail.sum(1).tolist(), two_plus=int(two_plus.sum()), reliability_median={ln: float(np.nanmedian(R_half[j][avail[j]])) for j, ln in enumerate(X.DONORS)},
                           ka0_variants={k_: {kk: vv for kk, vv in v.items() if kk not in ("k_rows", "pre64")} for k_, v in variants.items()}, pm_blend_variant=best),
               seconds=round(time.time() - t0, 1), peak_rss_gib=X.peak_gib())
    (X.WORK / "instrument273.json").write_text(json.dumps(out, indent=1, default=float))
    (X.WORK / "instrument273.log").write_text("\n".join(LOG) + "\n")
    P(f"\ndone: {len(FAIL)} KA failures, {len(arms)} arms, {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB")
    for f in FAIL:
        P(f"  *** {f}")
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
