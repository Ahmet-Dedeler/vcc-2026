"""271 -- the measurement on ONE PLOLO cancer fold: donor-row reliabilities (item 1) and the noise-aware estimators
(item 2, est271.py) through 263's construction on the cd4svd5 base, scored on the fold's CACHE truth (the deep
reference; the reserved orion_hek293t half-B truths are never read).  MEASUREMENT ONLY.

    PYTHONDONTWRITEBYTECODE=1 uv run python sweep271.py orion_hct116

Rows (every row rides cd4svd5's downstream: whitening through plolo.Operator.apply, CD4 deflated k = 5 re-solved at
f = 0.34, the H1-pair union re-solved by the carrying-row rule at g = 0.28, one add, one preserve; est271's docstring):
  cd4svd5, mean_explicit, eqR:nm  known answers (263's R1; equal reliabilities through the estimator path)
  STACK                            h1union's stack (un-deflated CD4 at f = 0.34), the retrodiction row for the pricer
  relw / relw_eb / relw_glob / relw_wht / ebrow / ebrow_eb / wiener_own / wiener_gene   the estimators (all :nm)
  a13:nm / a13:lit                 the baseline (d)
  relw:shuf / relw:rev / ebrow:shuf / wiener_own:shuf                                 the controls
  random:all:s0:nm, oracle_wht:all:nm                                                 263's reference rows (KA)
Per row: PLOLO endpoints (resid218.score_effects), 188's exact nmae identity (lib257.Identity), top-89 precision,
delivered norm, the fractions; per-target pds / own-cosine arrays, the identity arrays, the weights and the
reliability arrays to npz.
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/271-donor-reliability-estimation"
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E257 = ROOT / "experiments/257-soft-deflation-price"
E267 = ROOT / "experiments/267-replicate-validated-oracle"
sys.path.insert(0, str(E263))
sys.path.insert(0, str(HERE))
import lib263 as Y  # noqa: E402   (imports lib257 -> lib246 -> plolo / resid218 / compose218)
import est271 as E  # noqa: E402

Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK
plolo, L, X = Y.plolo, Y.L, Y.X
E062, E064 = plolo.E062, plolo.E064
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main(held: str) -> None:
    t0 = time.time()
    est_sha = E.sha256_self()
    P(f"[{held}] est271.py sha256 {est_sha}")
    st = Y.stage(held)
    ident = X.Identity(st)
    n, J = st["n"], len(st["rows"])
    targets = st["targets"]
    wc, sc = st["wcols"], st["sc"]
    avail = st["avail"]
    donors = list(st["sources"])
    truth = st["truth"]
    P(f"[{held}] staged {n} targets, donors {donors} ({st['donors']}), PM rows {st['pm_n']}, KA-0 pre max|d| "
      f"{st['ka0']['pre_float32_max_abs_diff']:.1e}, {time.time() - t0:.0f}s, peak {Y.peak_gib():.2f} GiB")
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

    check("KA-0 build_baseline(per-donor rows) float32 == cached pre", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)

    # ------------------------------------------------------------------------------ the half rows per donor
    pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
    HA, HB, NCF, NCA, NCB = [], [], [], [], []
    for j, line in enumerate(donors):
        z = np.load(WORK / f"half/rows-{line}.npz", allow_pickle=True)
        tz = [str(t) for t in z["targets"]]
        idx = np.asarray([tz.index(t) for t in targets], dtype=np.int64)
        A = z["A"][idx].astype(np.float64)
        B = z["B"][idx].astype(np.float64)
        R = z["real"][idx].astype(np.float64)
        avh = z["avail_A"][idx] & z["avail_B"][idx]
        check(f"KA-C availability of {line}'s half rows == the fold's donor availability", int((avh != avail[j]).sum()), 0, 0)
        # the production row == the real half-builder row: bit-identical (float32) off the PerturbMatch columns
        prod = st["rows"][j]
        if line == "k562_gwps":
            shared = E064.shared_k562_pm_columns(E062.Bank.load(line, plolo.DONOR_PATHS[line]), pm_bank, plolo.OFFICIAL, plolo.PANEL_SET)
            mask = np.ones(prod.shape[1], bool)
            mask[np.asarray(shared, dtype=np.int64)] = False
            d = float(np.abs(prod[avail[j]][:, mask].astype(np.float32) - R[avail[j]][:, mask].astype(np.float32)).max())
            dpm = float(np.abs(prod[avail[j]][:, ~mask] - R[avail[j]][:, ~mask]).max())
            check(f"KA-C {line} production row == real half-builder row OFF the {int((~mask).sum())} PerturbMatch columns (float32)", d, 0.0, 0.0)
            P(f"     ({line}: on the PM-blended columns the production row differs by max {dpm:.3f}, as built)")
        else:
            d = float(np.abs(prod[avail[j]].astype(np.float32) - R[avail[j]].astype(np.float32)).max())
            check(f"KA-C {line} production row == real half-builder row (float32)", d, 0.0, 0.0)
        HA.append(A)
        HB.append(B)
        NCF.append(np.asarray(z["ncells_full"][idx], dtype=np.int64))
        NCA.append(np.asarray(z["ncells_A"][idx], dtype=np.int64))
        NCB.append(np.asarray(z["ncells_B"][idx], dtype=np.int64))
    del pm_bank

    # ------------------------------------------------------------------------------ item 1: reliabilities
    r_raw = np.full((J, n), np.nan)
    r_wht = np.full((J, n), np.nan)
    r_eb = np.full((J, n), np.nan)
    eb_audit = {}
    WA, WB = [], []
    for j, line in enumerate(donors):
        ok = avail[j]
        r = E.row_cos(HA[j], HB[j])
        r_raw[j] = np.where(ok, r, np.nan)
        # KA-B: reproduces 267's rel|batch|cos_resid (the same quantity, the donor line's own fold order)
        tr = np.load(E267 / f"work/half/truth-{line}.npz", allow_pickle=True)
        t267 = {str(t): i for i, t in enumerate(tr["targets"])}
        r267 = np.asarray(tr["rel|batch|cos_resid"], dtype=np.float64)
        both = [t for t in targets if t in t267]
        d = float(np.nanmax(np.abs(np.array([r_raw[j][targets.index(t)] for t in both]) - np.array([r267[t267[t]] for t in both]))))
        check(f"KA-B {line} split-half reliability == 267 rel|batch|cos_resid on {len(both)} targets", d, 0.0, 1e-5)
        # whitened halves (each half whitened alone through the fold's operator; cosine on the operator columns)
        wa = st["op"].apply(HA[j] * ALPHA, wc, GAMMA, LAM)[:, wc]
        wb = st["op"].apply(HB[j] * ALPHA, wc, GAMMA, LAM)[:, wc]
        r_wht[j] = np.where(ok, E.row_cos(wa, wb), np.nan)
        WA.append(wa.astype(np.float32))
        WB.append(wb.astype(np.float32))
        eb = E.eb_shrink_half_reliability(r_raw[j], HA[j], HB[j], ok)
        r_eb[j] = eb.pop("r_shrunk")
        eb_audit[line] = eb
        rr = r_raw[j][ok]
        P(f"  reliability {line:14s} n {int(ok.sum())}: raw half-cos median {np.median(rr):+.4f} (p10 {np.percentile(rr, 10):+.4f}, p90 {np.percentile(rr, 90):+.4f}, "
          f"share<=0 {np.mean(rr <= 0):.3f}) -> SB full-depth median {np.median(E.spearman_brown(rr)):+.4f}; whitened median {np.nanmedian(r_wht[j][ok]):+.4f} "
          f"(share<=0 {np.mean(r_wht[j][ok] <= 0):.3f}); EB: rbar {eb['rbar']:+.4f} s2_null {eb['s2_null']:.2e} (null mean {eb['null_mean']:+.4f}) "
          f"var_r {eb['var_r']:.2e} kappa {eb['kappa']:.3f}")
    R_raw, R_wht, R_eb = E.spearman_brown(r_raw), E.spearman_brown(r_wht), E.spearman_brown(r_eb)
    r_shuf = E.shuffle_within_donor(r_raw, avail)
    R_shuf = E.spearman_brown(r_shuf)
    # the reference's own reliability on this fold (267, batch split) for the strata
    trh = np.load(E267 / f"work/half/truth-{held}.npz", allow_pickle=True)
    assert [str(t) for t in trh["targets"]] == list(targets)
    ref_rel = np.asarray(trh["rel|batch|cos_resid"], dtype=np.float64)

    # ------------------------------------------------------------------------------ the scoring path
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    weights: dict[str, np.ndarray] = {}
    ref263 = {r["arm"]: r for r in json.loads((E263 / f"work/sweep263-{held}.json").read_text())["rows"]}
    pt263 = np.load(E263 / f"work/pertarget263-{held}.npz")
    pk263 = np.load(E263 / f"work/picks263-{held}.npz")
    ref257 = {r["arm"]: r for r in json.loads((E257 / f"work/sweep257-{held}.json").read_text())["rows"]}

    def run(name, *, W=None, pre_sel=None, explicit_mean=False, post=None, parts=None, norm_match=True, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        if pre_sel is None:
            if explicit_mean:
                pre_sel = Y.combine(st, None, explicit_mean=True)
            elif W is not None:
                pre_sel = Y.combine(st, W, norm_match=norm_match)
        if pre_sel is None:
            pre_sel, base = st["pre"], st["whitened"]
        else:
            base = Y.base_from(st, pre_sel)
        if parts is None:
            parts, aud = Y.arm_parts(st, base)
        else:
            aud = dict(cd4_fraction=float(parts[0][1]), union_fraction=float(parts[1][1]), note="STACK: un-deflated CD4 at the staged fraction")
        st_arm = dict(st, whitened=base)
        final = L.compose(st_arm, parts)
        post_aud = {}
        if post is not None:
            final, post_aud = post(final)
        e = L.score(st_arm, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        if W is not None:
            weights[name] = W.astype(np.float32)
        n_changed = int((np.abs(pre_sel - st["pre"]).max(1) > 0).sum())
        aud["base_median_norm_wcols"] = float(np.median(np.linalg.norm(base[:, wc], axis=1)))
        row = dict(fold=held, arm=name, n_targets=n, n_rows_changed=n_changed, norm_match=bool(norm_match),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**aud, **post_aud, **(extra or {})), seconds=round(time.time() - t1, 1))
        msg = (f"  {name:18s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f} "
               f"dn {m['delivered_norm']:.4f} chg {n_changed:3d} ({row['seconds']}s)")
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

    # ------------------------------------------------------------------ known-answer rows
    r_base = run("cd4svd5", ka="263 cd4svd5", ka_ref=ref263["cd4svd5"], extra=dict(control="the incumbent: equal weights"))
    check("cd4svd5 per-target pds bit-identical to 263", float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
    run("mean_explicit", explicit_mean=True, ka="cd4svd5 (263's R1)", ka_ref=r_base)
    check("R1 mean_explicit per-target pds bit-identical to cd4svd5", float(np.abs(per_target["mean_explicit|pds"] - per_target["cd4svd5|pds"]).max()), 0.0, 0.0)
    W_eq = E.weights(np.full((J, n), 0.5), avail, "relw")
    r_eq = run("eqR:nm", W=W_eq, extra=dict(control="equal reliabilities through the estimator path (norm-matched mean)"))
    d_eq = max(abs(float(r_eq[k]) - float(r_base[k])) for k in KEYS)
    check("equal reliabilities through weights()+combine(norm_match) == cd4svd5 (1e-6: float32 re-cast of the norm-matched mean)", d_eq, 0.0, 1e-6)
    run("STACK", parts=[(st["cd4"], st["frac_cd4"]), (st["union"], st["frac_x"])], ka="257 STACK (h1union)", ka_ref=ref257["STACK"],
        extra=dict(control="h1union's stack: the retrodiction row for the pricer"))
    check("STACK nmae_wald89 == 257", rows[-1]["nmae_wald89"], ref257["STACK"]["nmae_wald89"], TOL)
    # 263's reference rows
    S_wht = np.stack([Y.row_cos(st["wht_rows"][j][:, sc], truth) for j in range(J)], 1)
    for j in range(J):
        S_wht[~avail[j], j] = np.nan
    gate_all = avail.sum(0) >= 2
    p_orc = Y.pick_from_scores(S_wht, avail, gate_all)
    check("oracle picks == 263 oracle_wht:all:nm", int(np.abs(p_orc - pk263["oracle_wht:all:nm"]).max()), 0, 0)
    run("oracle_wht:all:nm", W=Y.one_hot(st, p_orc), ka="263 oracle_wht:all:nm", ka_ref=ref263["oracle_wht:all:nm"], extra=dict(control="263's whitened oracle (the ceiling; uses the truth)"))
    p_rnd = np.asarray(pk263["random:all:s0:nm"], dtype=np.int64)
    run("random:all:s0:nm", W=Y.one_hot(st, p_rnd), ka="263 random:all:s0:nm", ka_ref=ref263["random:all:s0:nm"], extra=dict(control="263's random pick, seed 0"))

    # ------------------------------------------------------------------ the estimators (a), (b)
    run("relw", W=E.weights(R_raw, avail, "relw"), extra=dict(estimator="(a) reliability-weighted mean, per-target SB reliabilities"))
    run("relw_eb", W=E.weights(R_eb, avail, "relw"), extra=dict(estimator="(a) with EB-shrunk reliabilities", kappa={d: eb_audit[d]["kappa"] for d in donors}))
    run("relw_glob", W=E.weights(R_raw, avail, "glob"), extra=dict(estimator="(a) with donor-level mean reliabilities (scalar donor weights)",
                                                                   Rbar={d: float(np.nanmean(R_raw[j][avail[j]])) for j, d in enumerate(donors)}))
    run("relw_wht", W=E.weights(R_wht, avail, "relw"), extra=dict(estimator="(a) with reliabilities measured after the operator"))
    run("ebrow", W=E.weights(R_raw, avail, "ebrow"), extra=dict(estimator="(b) EB row shrinkage toward the consensus, then the mean"))
    run("ebrow_eb", W=E.weights(R_eb, avail, "ebrow"), extra=dict(estimator="(b) with EB-shrunk reliabilities"))
    # controls for (a), (b)
    run("relw:shuf", W=E.weights(R_shuf, avail, "relw"), extra=dict(control="(a) with reliabilities permuted across targets within donor"))
    run("relw:rev", W=E.weights(R_raw, avail, "rev"), extra=dict(control="(a) with reversed weights (1 - R) / R"))
    run("ebrow:shuf", W=E.weights(R_shuf, avail, "ebrow"), extra=dict(control="(b) with permuted reliabilities"))

    # ------------------------------------------------------------------ (c) the diagonal filters
    dAB = [HA[j] - HB[j] for j in range(J)]
    rows_f = [E.wiener_own(st["rows"][j], dAB[j]) for j in range(J)]
    kept = {donors[j]: float(np.mean((rows_f[j][avail[j]] != 0) & (st["rows"][j][avail[j]] != 0)) / max(np.mean(st["rows"][j][avail[j]] != 0), 1e-12)) for j in range(J)}
    en = {donors[j]: float((rows_f[j][avail[j]] ** 2).sum() / (st["rows"][j][avail[j]] ** 2).sum()) for j in range(J)}
    run("wiener_own", pre_sel=E.combine_filtered(st, rows_f), extra=dict(estimator="(c) per-gene Wiener with the row's own half difference", share_coords_kept=kept, energy_kept=en))
    rows_g, aud_g = [], {}
    for j in range(J):
        rf, a = E.wiener_gene(st["rows"][j], dAB[j], NCF[j], avail[j])
        rows_g.append(rf)
        aud_g[donors[j]] = a
    run("wiener_gene", pre_sel=E.combine_filtered(st, rows_g), extra=dict(estimator="(c) per-gene Wiener, gene-pooled noise scaled by cell count", per_donor=aud_g))
    rng = np.random.default_rng(E.SHUF_SEED)
    rows_s = []
    for j in range(J):
        d_s = dAB[j].copy()
        idx = np.flatnonzero(avail[j])
        d_s[idx] = dAB[j][idx][rng.permutation(idx.size)]
        rows_s.append(E.wiener_own(st["rows"][j], d_s))
    run("wiener_own:shuf", pre_sel=E.combine_filtered(st, rows_s), extra=dict(control="(c) with the half differences permuted across targets within donor"))

    # ------------------------------------------------------------------ (d) a13's damping on the final effect
    alpha, aud_a = E.a13_alpha(st["rows"], avail)
    P(f"  a13 key: {aud_a['n_active']} active rows, median graded fraction {aud_a['median_graded_fraction']:.4f}; alpha median {np.median(alpha):.4f}, "
      f"alpha on the zero-agreement block {float(np.median(alpha[alpha < 0.999])):.4f}")
    run("a13:nm", post=lambda F: (E.apply_alpha(F, alpha, norm_match=True), dict(a13=aud_a)), extra=dict(baseline="(d) a13 per-gene agreement damping, gamma 0.25, rows re-normalised"))
    run("a13:lit", post=lambda F: (E.apply_alpha(F, alpha, norm_match=False), dict(a13=aud_a)), extra=dict(baseline="(d) a13 literal (amplitude-changing)"))

    # ------------------------------------------------------------------ save
    meta = dict(fold=held, n_targets=n, targets=targets, donors=donors, donor_names=st["donors"], J=J, est271_sha256=est_sha,
                sweep_sha256=sha256_file(Path(__file__)), frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"],
                eb=eb_audit, a13=aud_a, ka=ka_log, ka_pass=sum(1 for k in ka_log if k["ok"]), ka_total=len(ka_log), failures=failures,
                n_rows=len(rows), seconds=round(time.time() - t0, 1), peak_rss_gib=Y.peak_gib(), fold_cache=st["fold"]["audit"].get("cache"),
                reserve_note="the orion_hek293t half-B TRUTHS (267's reserve) were not read; HEK293T's half-B bank rows enter only as a DONOR's own replicate")
    (WORK / f"sweep271-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=str))
    np.savez_compressed(WORK / f"pertarget271-{held}.npz", targets=np.asarray(targets), avail=avail.astype(np.int8),
                        r_raw=r_raw, r_wht=r_wht, r_eb=r_eb, r_shuf=r_shuf, R_raw=R_raw, R_wht=R_wht, R_eb=R_eb,
                        ncells_full=np.stack(NCF), ncells_A=np.stack(NCA), ncells_B=np.stack(NCB), ref_rel=ref_rel,
                        base262_cos_truth=np.asarray(np.load(E263 / f"work/cov263-{held}.npz")["base262_cos_truth"], dtype=np.float64),
                        cos_truth_wht=S_wht, alpha_a13_median_row=np.median(alpha, axis=1),
                        **per_target, **{f"W|{k}": v for k, v in weights.items()})
    np.savez_compressed(WORK / f"identity271-{held}.npz", **identity)
    np.savez_compressed(WORK / f"halfwht271-{held}.npz", **{f"WA|{donors[j]}": WA[j] for j in range(J)}, **{f"WB|{donors[j]}": WB[j] for j in range(J)})
    (WORK / f"sweep271-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{held}] done: {len(rows)} rows, KA {meta['ka_pass']}/{meta['ka_total']}, failures {failures}, {(time.time() - t0) / 60:.1f} min, peak {Y.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1])
