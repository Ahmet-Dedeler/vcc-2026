"""263 sweep -- every arbiter x gating cell of the per-target donor selection on ONE PLOLO fold, on the `cd4svd5` base.

    uv run python sweep263.py --fold <k562_gwps | orion_hct116 | orion_hek293t | kolf2> [--quick]

Rows (every cell rides `cd4svd5`'s downstream: CD4 deflated k = 5 re-solved at f = 0.34, the H1-pair union re-solved
by the carrying-row rule at g = 0.28, one add, one preserve; see lib263):
  cd4svd5                 the incumbent, the cached mean base (KA vs 262 diag262-<fold>.json / .npz at 1e-9 / bit-identical)
  mean_explicit           R1: the equal-weight mean rebuilt from the per-donor rows through the selection path (0.0e+00)
  oracle_wht:<g>          R2, the CEILING: per target the donor whose WHITENED row has the higher cosine to the held truth
  oracle_raw:<g>          the same with the un-whitened rows
  oracle3_wht:<g>         the abstaining oracle: the best of {donor rows, the mean row} by whitened cosine to the truth
  <arb>:<space>:<g>       THE ARMS: per target the donor the arbiter scores higher; <space> raw = the arbiter row against
                          the donors' un-whitened rows on the operator columns, wht = against their whitened rows
  <arb>:<space>:<g>:rev   the REVERSE control (the donor the arbiter scores LOWER)
  random:<g>:s<k>         the random-pick control, three seeds
  <arb>:<space>:all:sm<b> continuous weighting, w ∝ softmax(beta * standardised score), beta in {1, 3}
Gating <g>: tau0 = only targets whose donors' un-whitened rows have cosine < 0; tau0.1 = < 0.1; all = every target carried
by >= 2 donors.  Arbiters: cd4 (the un-deflated CD4 row), cd4d (the deflated CD4 row), union (the H1-pair union row,
carried targets only), kolf2 (KOLF2's own row, cancer folds only; on the kolf2 fold `majority`: keep the two donors that
agree, drop the odd one out), norm (the donor with the larger row norm), selfcons (kNN self-consistency, k = 5), basal
(the donor whose control profile is closer to the held line's, a fold-constant pick; raw only).
Per row: PLOLO endpoints, 188's exact nmae identity (wald89, frame A), top-89 precision, delivered norm, the fractions,
the picks; per-target pds / cosine arrays and the identity arrays to npz; the covariates and arbiter scores to npz.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib263 as Y  # noqa: E402

L, X, plolo = Y.L, Y.X, Y.plolo
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
GATES = ("tau0", "tau0.1", "all")
BETAS = (1.0, 3.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    held, quick = args.fold, args.quick
    t0 = time.time()
    ref262 = json.loads((Y.E262 / f"work/diag262-{held}.json").read_text())
    pt262 = np.load(Y.E262 / f"work/diag262-{held}.npz")

    st = Y.stage(held)
    ident = X.Identity(st)
    n, J = st["n"], len(st["rows"])
    targets = st["targets"]
    wc, sc = st["wcols"], st["sc"]
    avail = st["avail"]
    print(f"[{held}] staged {n} targets, donors {st['donors']} (votes hist "
          f"{ {int(k): int(v) for k, v in zip(*np.unique(avail.sum(0), return_counts=True))} }), PM rows {st['pm_n']}, "
          f"frac_cd4 {st['frac_cd4']:.6f} frac_x {st['frac_x']:.6f}; KA-0 pre float32 max|d| {st['ka0']['pre_float32_max_abs_diff']:.1e}, "
          f"votes {st['ka0']['votes_max_abs_diff']} ({time.time()-t0:.0f}s, peak {Y.peak_gib():.2f} GiB)", flush=True)

    failures: list[str] = []
    ka_log: list[dict] = []
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    picks: dict[str, np.ndarray] = {}

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka_log.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        if not ok:
            failures.append(f"{label}: got {got!r} want {want!r} |d| {d:.3e}")
        return ok

    check("KA-0 build_baseline(per-donor rows) float32 == cached pre", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)
    check("KA-0 votes == per-donor availability", st["ka0"]["votes_max_abs_diff"], 0.0, 0.0)

    # ------------------------------------------------------------------ covariates and arbiter scores
    R = [r[:, wc] for r in st["rows"]]                       # un-whitened donor rows on the operator columns
    Wh = [w[:, wc] for w in st["wht_rows"]]                  # the same rows whitened (each donor alone)
    truth = st["truth"]
    cos_truth_raw = np.stack([Y.row_cos(st["rows"][j][:, sc], truth) for j in range(J)], 1)       # (n, J)
    cos_truth_wht = np.stack([Y.row_cos(st["wht_rows"][j][:, sc], truth) for j in range(J)], 1)
    for j in range(J):
        cos_truth_raw[~avail[j], j] = np.nan
        cos_truth_wht[~avail[j], j] = np.nan
    mean_cos_truth_raw = Y.row_cos(st["pre"][:, sc], truth)
    mean_cos_truth_wht = Y.row_cos(st["whitened"][:, sc], truth)               # 262's base·truth
    pairs = [(a, b) for a in range(J) for b in range(a + 1, J)]
    dd_raw_pairs = np.stack([Y.row_cos(R[a], R[b]) for a, b in pairs], 1)      # (n, P)
    dd_wht_pairs = np.stack([Y.row_cos(Wh[a], Wh[b]) for a, b in pairs], 1)
    for p, (a, b) in enumerate(pairs):
        m = ~(avail[a] & avail[b])
        dd_raw_pairs[m, p] = np.nan
        dd_wht_pairs[m, p] = np.nan
    with np.errstate(invalid="ignore"):
        dd_raw = np.nanmin(dd_raw_pairs, axis=1) if J > 2 else dd_raw_pairs[:, 0]
        dd_wht = np.nanmin(dd_wht_pairs, axis=1) if J > 2 else dd_wht_pairs[:, 0]
    two = avail.sum(0) >= 2
    gates = {"tau0": two & (dd_raw < 0.0), "tau0.1": two & (dd_raw < 0.1), "all": two.copy()}
    print(f"[{held}] gates: " + ", ".join(f"{g} {int(m.sum())}" for g, m in gates.items()) +
          f"; donor-donor raw cos median {np.nanmedian(dd_raw):+.4f} (whitened {np.nanmedian(dd_wht):+.4f})", flush=True)

    def ext_scores(A):                     # an external (un-whitened) row against each donor's rows, both spaces
        Aw = A[:, wc]
        live = np.linalg.norm(Aw, axis=1) > 0
        s_raw = np.stack([Y.row_cos(Aw, R[j]) for j in range(J)], 1)
        s_wht = np.stack([Y.row_cos(Aw, Wh[j]) for j in range(J)], 1)
        s_raw[~live] = np.nan
        s_wht[~live] = np.nan
        for j in range(J):
            s_raw[~avail[j], j] = np.nan
            s_wht[~avail[j], j] = np.nan
        return dict(raw=s_raw, wht=s_wht)

    cd4_row = np.where(np.asarray(st["cd4"].gene_votes) > 0, np.asarray(st["cd4"].mean_effect), 0.0)
    d5_row = np.where(np.asarray(st["d5"].gene_votes) > 0, np.asarray(st["d5"].mean_effect), 0.0)
    u_row = np.where(np.asarray(st["union"].gene_votes) > 0, np.asarray(st["union"].mean_effect), 0.0)
    S: dict[str, dict[str, np.ndarray]] = {}
    S["cd4"] = ext_scores(cd4_row)
    S["cd4d"] = ext_scores(d5_row)
    S["union"] = ext_scores(u_row)
    if st["kolf_rows"] is not None:
        k_row = np.where(st["kolf_avail"][:, None], st["kolf_rows"], 0.0)
        S["kolf2"] = ext_scores(k_row)
    norm_raw = np.stack([np.linalg.norm(R[j], axis=1) for j in range(J)], 1)
    norm_wht = np.stack([np.linalg.norm(Wh[j], axis=1) for j in range(J)], 1)
    for j in range(J):
        norm_raw[~avail[j], j] = np.nan
        norm_wht[~avail[j], j] = np.nan
    S["norm"] = dict(raw=norm_raw, wht=norm_wht)

    def selfcons(space_rows):
        out = np.full((n, J), np.nan)
        for j in range(J):
            others = [k for k in range(J) if k != j]
            if len(others) == 1:
                other, oav = space_rows[others[0]], avail[others[0]]
            else:                          # the mean of the other donors' rows over those carrying the target
                acc = np.zeros_like(space_rows[0])
                cnt = np.zeros(n)
                for k in others:
                    acc[avail[k]] += space_rows[k][avail[k]]
                    cnt += avail[k]
                oav = cnt > 0
                other = np.where(oav[:, None], acc / np.maximum(cnt, 1)[:, None], 0.0)
            out[:, j] = Y.knn_self_consistency(space_rows[j], avail[j], other, oav, k=Y.KNN)
        return out
    S["selfcons"] = dict(raw=selfcons(R), wht=selfcons(Wh))
    basal = np.array([np.corrcoef(st["b_held"][wc], st["b_donor"][j][wc])[0, 1] for j in range(J)])
    S["basal"] = dict(raw=np.where(avail.T, basal[None, :], np.nan))
    print(f"[{held}] basal Pearson(held, donor) on the operator columns: " +
          ", ".join(f"{st['donors'][j]} {basal[j]:.4f}" for j in range(J)), flush=True)

    # agreement of each arbiter's pick with the whitened oracle's pick, all targets with both defined
    orc_pick_all = Y.pick_from_scores(cos_truth_wht, avail, gates["all"])
    agree = {}
    for arb, sp in S.items():
        for space, s in sp.items():
            p = Y.pick_from_scores(s, avail, gates["all"])
            ok = (p >= 0) & (orc_pick_all >= 0)
            agree[f"{arb}:{space}"] = dict(n=int(ok.sum()), agree=float((p[ok] == orc_pick_all[ok]).mean()) if ok.any() else float("nan"))
    print(f"[{held}] arbiter == whitened-oracle pick: " + ", ".join(f"{k} {v['agree']*100:.0f}% (n {v['n']})" for k, v in agree.items()), flush=True)

    cov = dict(targets=np.asarray(targets), avail=avail.astype(np.int8), cos_truth_raw=cos_truth_raw, cos_truth_wht=cos_truth_wht,
               mean_cos_truth_raw=mean_cos_truth_raw, mean_cos_truth_wht=mean_cos_truth_wht, dd_raw=dd_raw, dd_wht=dd_wht,
               dd_raw_pairs=dd_raw_pairs, dd_wht_pairs=dd_wht_pairs, norm_raw=norm_raw, norm_wht=norm_wht, basal=basal,
               base262_cos_truth=pt262["cov|base_cos_truth"], cd4svd5_262_pds=pt262["cd4svd5|pds"], cd4svd5_262_own=pt262["cd4svd5|own"],
               **{f"gate|{g}": m.astype(np.int8) for g, m in gates.items()},
               **{f"S|{arb}|{space}": s for arb, sp in S.items() for space, s in sp.items()})

    # ------------------------------------------------------------------ the runner
    def run(name, *, W=None, explicit_mean=False, pick=None, extra=None, ka=None, ka_ref=None, norm_match=False):
        t1 = time.time()
        if explicit_mean or W is not None:
            pre_sel = Y.combine(st, W, explicit_mean=explicit_mean, norm_match=norm_match)
            wh = Y.base_from(st, pre_sel)
        else:
            pre_sel, wh = st["pre"], st["whitened"]
        e, m = Y.evaluate(st, wh, ident)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        if pick is None and W is not None:
            pick = np.where(np.isfinite(W).all(1), np.argmax(np.nan_to_num(W, nan=-1.0), axis=1), -1)
        if pick is not None:
            picks[name] = np.asarray(pick, dtype=np.int8)
        n_sel = int((pick >= 0).sum()) if pick is not None else 0
        n_changed = int((np.abs(pre_sel - st["pre"]).max(1) > 0).sum())
        row = dict(fold=held, arm=name, n_targets=n, n_selected=n_sel, n_rows_changed=n_changed, norm_match=bool(norm_match),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**m["audit"], **(extra or {})), seconds=round(time.time() - t1, 1))
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
               f"nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} sel {n_sel} chg {n_changed} ({row['seconds']}s)")
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            ok = d < TOL
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(ok))
            msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({d:.1e})"
            if not ok:
                failures.append(f"{name} vs {ka}: {d:.3e}")
        rows.append(row)
        print(msg, flush=True)
        return row

    # ------------------------------------------------------------------ known-answer rows
    r_base = run("cd4svd5", ka="262 cd4svd5 (== 257 k5_s1 == 261 R1)", ka_ref=ref262["arms"]["cd4svd5"])
    check("cd4svd5 per-target pds bit-identical to 262", float(np.abs(per_target["cd4svd5|pds"] - pt262["cd4svd5|pds"]).max()), 0.0, 0.0)
    check("cd4svd5 per-target cos bit-identical to 262", float(np.abs(per_target["cd4svd5|cos"] - pt262["cd4svd5|own"]).max()), 0.0, 0.0)
    check("cd4svd5 CD4 fraction == 262 frac_cd4_resolved", r_base["audit"]["cd4_fraction"], ref262["meta"]["frac_cd4_resolved"], 0.0)
    check("cd4svd5 union fraction == 262 frac_x", r_base["audit"]["union_fraction"], ref262["meta"]["frac_x"], 0.0)
    pre_expl = Y.combine(st, None, explicit_mean=True)
    check("R1 explicit mean (sum / count over the per-donor rows, float32) == cached pre, bit for bit",
          float(np.abs(pre_expl - st["pre"]).max()), 0.0, 0.0)
    r_expl = run("mean_explicit", explicit_mean=True, ka="cd4svd5 (R1)", ka_ref=r_base)
    check("R1 mean_explicit per-target pds bit-identical to cd4svd5", float(np.abs(per_target["mean_explicit|pds"] - per_target["cd4svd5|pds"]).max()), 0.0, 0.0)
    for key in ("nmae_wald89", "top89_precision", "delivered_norm"):
        check(f"R1 mean_explicit.{key} == cd4svd5", r_expl[key], r_base[key], 0.0)
    # the whitened base of the trivial selection is the staged base bit for bit
    check("R1 base_from(explicit mean) == staged whitened, bit for bit",
          0.0 if np.array_equal(Y.base_from(st, pre_expl), st["whitened"]) else 1.0, 0.0, 0.0)

    def run2(name, **kw):
        r = run(name, **kw)
        run(name + ":nm", norm_match=True, **kw)
        return r

    # ------------------------------------------------------------------ oracle rows (the ceiling), per gating
    for g in GATES:
        p = Y.pick_from_scores(cos_truth_wht, avail, gates[g])
        run2(f"oracle_wht:{g}", W=Y.one_hot(st, p), pick=p, extra=dict(oracle="whitened row cosine to the held truth", gate=g))
        if quick:
            continue
        p = Y.pick_from_scores(cos_truth_raw, avail, gates[g])
        run2(f"oracle_raw:{g}", W=Y.one_hot(st, p), pick=p, extra=dict(oracle="un-whitened row cosine to the held truth", gate=g))
        # abstaining oracle: the best of {each donor row, the mean row}
        p3 = np.full(n, -1, dtype=np.int64)
        for t in np.flatnonzero(gates[g]):
            s = np.where(avail[:, t], cos_truth_wht[t], -np.inf)
            j = int(np.argmax(s))
            if np.isfinite(s[j]) and s[j] > mean_cos_truth_wht[t]:
                p3[t] = j
        run2(f"oracle3_wht:{g}", W=Y.one_hot(st, p3), pick=p3, extra=dict(oracle="best of {donor rows, mean row} by whitened cosine", gate=g))

    # ------------------------------------------------------------------ the arms and their controls
    arb_list = list(S.keys())
    if quick:
        arb_list = ["cd4"]
    for arb in arb_list:
        for space, s in S[arb].items():
            if quick and space != "raw":
                continue
            for g in GATES:
                p = Y.pick_from_scores(s, avail, gates[g])
                run2(f"{arb}:{space}:{g}", W=Y.one_hot(st, p), pick=p, extra=dict(arbiter=arb, space=space, gate=g))
                pr = Y.pick_from_scores(s, avail, gates[g], reverse=True)
                run2(f"{arb}:{space}:{g}:rev", W=Y.one_hot(st, pr), pick=pr, extra=dict(arbiter=arb, space=space, gate=g, control="reverse"))
            if not quick:
                for beta in BETAS:
                    Wsm = Y.softmax_weights(st, s, beta, gates["all"])
                    run2(f"{arb}:{space}:all:sm{beta:g}", W=Wsm, extra=dict(arbiter=arb, space=space, gate="all", beta=beta, weighting="softmax"))
    if J == 3:                                          # kolf2: the majority rule (drop the odd one out)
        for space, pairs_cos in (("raw", dd_raw_pairs), ("wht", dd_wht_pairs)):
            for g in GATES:
                for rev in (False, True):
                    W = np.full((n, J), np.nan)
                    p = np.full(n, -1, dtype=np.int64)
                    for t in np.flatnonzero(gates[g]):
                        if not avail[:, t].all() or not np.isfinite(pairs_cos[t]).all():
                            continue
                        c = {pairs[k]: pairs_cos[t, k] for k in range(len(pairs))}
                        a = np.array([np.mean([c[pr] for pr in pairs if j in pr]) for j in range(J)])
                        odd = int(np.argmax(a) if rev else np.argmin(a))
                        W[t] = 1.0
                        W[t, odd] = 0.0
                        p[t] = odd
                    run2(f"majority:{space}:{g}{':rev' if rev else ''}", W=W, pick=p,
                        extra=dict(arbiter="majority", space=space, gate=g, control="reverse" if rev else None,
                                   note="pick = the DROPPED donor"))
    for g in GATES:
        for k in range(3):
            if quick and (g != "all" or k > 0):
                continue
            p = Y.random_pick(avail, gates[g], Y.SHUF_SEED + 101 * k + len(g))
            run2(f"random:{g}:s{k}", W=Y.one_hot(st, p), pick=p, extra=dict(control="random", gate=g, seed=k))

    n_ok = sum(1 for r in ka_log if r["ok"])
    meta = dict(fold=held, n_targets=n, targets=targets, donors=st["donors"], sources=st["sources"], J=J, pm_rows=st["pm_n"],
                frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"], deflation=st["aud5"], gates={g: int(m.sum()) for g, m in gates.items()},
                dd_raw_median=float(np.nanmedian(dd_raw)), dd_wht_median=float(np.nanmedian(dd_wht)), basal=basal.tolist(),
                agree_with_oracle=agree, n_rows=len(rows), ka_rows=ka_log, ka_pass=n_ok, ka_total=len(ka_log), failures=failures,
                quick=quick, seconds=round(time.time() - t0, 1), peak_rss_gib=Y.peak_gib(), fold_cache=st["fold"]["audit"].get("cache"))
    tag = "quick263" if quick else "sweep263"
    (Y.WORK / f"{tag}-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=float))
    np.savez_compressed(Y.WORK / f"pertarget263-{held}{'-quick' if quick else ''}.npz", **per_target)
    np.savez_compressed(Y.WORK / f"identity263-{held}{'-quick' if quick else ''}.npz", **identity)
    np.savez_compressed(Y.WORK / f"picks263-{held}{'-quick' if quick else ''}.npz", **picks)
    np.savez_compressed(Y.WORK / f"cov263-{held}{'-quick' if quick else ''}.npz", **cov)
    print(f"[{held}] {len(rows)} rows, known-answer {n_ok}/{len(ka_log)} (+{sum(1 for r in rows if 'ka' in r and r['ka']['ok'])}/"
          f"{sum(1 for r in rows if 'ka' in r)} row KAs), {len(failures)} failures, {meta['seconds']}s, peak RSS {meta['peak_rss_gib']:.2f} GiB",
          flush=True)
    for f in failures:
        print(f"  *** {f}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
