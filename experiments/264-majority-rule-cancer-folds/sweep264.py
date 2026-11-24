"""264 sweep -- ONE PLOLO fold per process: the three-donor majority rule with kolf2 as the third donor (cancer folds)
or with the fold's own three donors (kolf2, the known-answer fold), every row on cd4svd5's downstream.

    PYTHONDONTWRITEBYTECODE=1 uv run python sweep264.py --fold <k562_gwps | orion_hct116 | orion_hek293t | kolf2> [--quick] [--kolf-weight equal|matched]
`--kolf-weight matched` (sensitivity, cancer folds; outputs tagged `sweep264m-*`): kolf2's rows scaled by ONE fold-constant scalar
so their un-preserved whitened median norm equals the mean of the two production donors' (the operator, fitted without kolf2,
amplifies its rows ~50x against ~6x; the raw-space picks are scale-invariant, so they are identical to the equal-mean run).

Rows (all through lib263's `evaluate`: CD4 deflated k = 5 re-solved at f = 0.34, the H1-pair union re-solved by the
carrying-row rule at g = 0.28, one add, one preserve; 188's exact nmae identity; the scorer's per-target rank statistic):
  cd4svd5                     (c) the incumbent: the cached two-donor base (R1: == 263's row at 0.0e+00)
  mean3                       (a) the THREE-donor equal mean, `E062.build_baseline` over [donor 1, donor 2, kolf2]
  mean3:nm                    (a) with every changed row rescaled to the cached two-donor row's norm (shape-only form of adding kolf2)
  mean3:unrestricted          (a) with kolf2's columns outside the operator axis kept (sensitivity; cancer folds only)
  majority:raw:<g>[:nm|:nmc]  (b) drop the odd donor out by un-whitened cosine on the operator columns; :nm = norm-matched to
                              the three-row mean (THE construction); :nmc = norm-matched to the cached two-donor row (shape-only vs c)
  majority:wht:all:nm         the same rule on the rows whitened alone (263: loses on kolf2)
  reverse:raw:<g>[:nm]        CONTROL: drop the donor that agrees MOST (keep the two least agreeing)
  random:<g>:s<k>[:nm]        CONTROL: drop a random donor, three seeds
  drop:<donor>:all[:nm]       CONTROL: always drop that donor on every three-donor target (a global weight)
  only:kolf2:all:nm           kolf2 ALONE on the three-donor targets, norm-matched (what the third donor is worth by itself)
  oracle_drop:<space>:all:nm  CEILING reference: drop the donor whose row has the lowest cosine to the held truth
Gating <g>: all = targets carried by all three donors; tau0 = additionally the minimum pairwise un-whitened cosine < 0.
Per row: PLOLO endpoints, nmae89, delivered norm, the fractions, the dropped donor; per-target pds / cosine and the
identity arrays to npz.  Known-answer rows: R1 (cd4svd5 == 263), R2 (kolf2: majority:raw:all:nm == 263 at 0.0e+00, plus
the as-is / tau0 / wht / reverse rows), R3 (kolf2's prepared donor rows == the kolf2 fold's truth rows, hashed).
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
import lib264 as Z  # noqa: E402

Y, X, L, plolo = Z.Y, Z.X, Z.L, Z.plolo
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
GATES = ("all", "tau0")
# 263's names for the rows my rows must reproduce on the kolf2 fold
KA263_KOLF2 = {"majority:raw:all:nm": "majority:raw:all:nm", "majority:raw:all": "majority:raw:all",
               "majority:raw:tau0:nm": "majority:raw:tau0:nm", "majority:raw:tau0": "majority:raw:tau0",
               "majority:wht:all:nm": "majority:wht:all:nm",
               "reverse:raw:all:nm": "majority:raw:all:rev:nm", "reverse:raw:all": "majority:raw:all:rev",
               "reverse:raw:tau0:nm": "majority:raw:tau0:rev:nm"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--kolf-weight", default="equal", choices=("equal", "matched"),
                    help="equal = the mandate's three-donor equal mean; matched = kolf2 scaled by a fold-constant scalar so its "
                         "whitened amplitude matches the production donors' (sensitivity; picks identical)")
    args = ap.parse_args()
    held, quick, kw = args.fold, args.quick, args.kolf_weight
    t0 = time.time()
    sw263 = json.loads((Z.E263 / f"work/sweep263-{held}.json").read_text())
    ref263 = {r["arm"]: r for r in sw263["rows"]}
    pt263 = np.load(Z.E263 / f"work/pertarget263-{held}.npz")

    st = Z.stage3(held, kolf_weight=kw)
    ident = X.Identity(st)
    n, J = st["n"], len(st["rows3"])
    targets = st["targets"]
    wc, sc = st["wcols"], st["sc"]
    avail3 = st["avail3"]
    three = avail3.all(0)
    print(f"[{held}] staged {n} targets; donors {st['donors3']}; votes3 {st['aud3']['votes3_hist']}; three-donor targets {int(three.sum())}; "
          f"rows changed vs cached pre {st['aud3']['n_changed_vs_cached_pre']}; filler rows {st['aud3']['n_filler_rows']}; "
          f"pre_a == cached pre: {st['aud3']['pre_a_equals_cached_pre']} ({time.time()-t0:.0f}s, peak {Z.peak_gib():.2f} GiB)", flush=True)
    if held != "kolf2":
        a = st["aud3"]
        print(f"[{held}] R3 kolf2 rows: {a['n_common']} targets in common with the kolf2 fold, rows bit-identical to its cached truth "
              f"{a['rows_bit_identical']}/{a['n_common']} (max|d| {a['max_abs_diff_float32']:.1e}); sha mine {a['sha256_mine'][:12]} theirs "
              f"{a['sha256_kolf2_fold_truth'][:12]}; columns outside the operator axis {a['n_cols_outside_operator']}/{a['n_carried_cols']} carrying "
              f"{100*a['energy_frac_outside_operator']:.3f} % of kolf2's energy; kolf2 weight {a['kolf2_weight']:.6f} ({kw})", flush=True)
    for nm, v in st["whitening_audit"].items():
        print(f"[{held}] whitening audit {nm:40s} n {v['n']:3d} raw med {v['raw_median']:.3f} amp x{v['amplification']:.2f} energy on gate A raw {100*v['energy_on_gateA_raw']:.1f} % "
              f"/ wht {100*v['energy_on_gateA_wht']:.1f} %, bottom eigen-decile {100*v['energy_bottom_eigen_decile_wht']:.1f} %, delivered/effect {v['delivered_over_effect']:.3f}", flush=True)

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

    check("KA-0 (263) build_baseline(per-donor rows) float32 == cached pre", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)
    if held == "kolf2":
        check("KA-0 kolf2: the three-donor base (a) rebuilt here == cached pre, bit for bit", 0.0 if st["aud3"]["pre_a_equals_cached_pre"] else 1.0, 0.0, 0.0)
    else:
        check("R3 kolf2 donor rows == 263's kolf2 fold truth rows (float32, all common targets)", st["aud3"]["rows_bit_identical"], st["aud3"]["n_common"], 0)

    # ------------------------------------------------------------------ covariates: pair cosines, gates, truth cosines
    R = [r[:, wc] for r in st["rows3"]]
    Wh = [w[:, wc] for w in st["wht_rows3"]]
    truth = st["truth"]
    dd_raw_pairs = Z.pair_cosines(R, avail3)
    dd_wht_pairs = Z.pair_cosines(Wh, avail3)
    with np.errstate(invalid="ignore"):
        dd_raw = np.nanmin(dd_raw_pairs, axis=1)
    gates = {"all": three.copy(), "tau0": three & (dd_raw < 0.0)}
    cos_truth_raw = np.stack([Z.row_cos(st["rows3"][j][:, sc], truth) for j in range(J)], 1)
    cos_truth_wht = np.stack([Z.row_cos(st["wht_rows3"][j][:, sc], truth) for j in range(J)], 1)
    for j in range(J):
        cos_truth_raw[~avail3[j], j] = np.nan
        cos_truth_wht[~avail3[j], j] = np.nan
    base_a_cos_truth = Z.row_cos(st["whitened_a"][:, sc], truth)
    base_c_cos_truth = Z.row_cos(st["whitened"][:, sc], truth)
    print(f"[{held}] gates: all {int(gates['all'].sum())}, tau0 {int(gates['tau0'].sum())}; min pairwise raw cosine median {np.nanmedian(dd_raw):+.4f}; "
          f"pair medians " + ", ".join(f"{st['donors3'][a]}~{st['donors3'][b]} {np.nanmedian(dd_raw_pairs[:, k]):+.4f}" for k, (a, b) in enumerate(Z.PAIRS)), flush=True)

    # ------------------------------------------------------------------ the runner
    def run(name, *, W=None, pre=None, pick=None, extra=None, ka=None, ka_ref=None, norm_match=False, norm_ref=None):
        t1 = time.time()
        if pre is not None:
            pre_sel = pre
            wh = Y.base_from(st, pre_sel)
        elif W is not None:
            pre_sel = Z.combine3(st, W, norm_match=norm_match, norm_ref=norm_ref)
            wh = Y.base_from(st, pre_sel)
        else:
            pre_sel, wh = st["pre"], st["whitened"]
        e, m = Y.evaluate(st, wh, ident)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        if pick is not None:
            picks[name] = np.asarray(pick, dtype=np.int8)
        n_sel = int((pick >= 0).sum()) if pick is not None else 0
        row = dict(fold=held, arm=name, n_targets=n, n_selected=n_sel,
                   n_rows_changed_vs_a=int((np.abs(pre_sel - st["pre_a"]).max(1) > 0).sum()),
                   n_rows_changed_vs_c=int((np.abs(pre_sel - st["pre"]).max(1) > 0).sum()), norm_match=bool(norm_match),
                   base_median_norm_pre_wcols=float(np.median(np.linalg.norm(pre_sel[:, wc], axis=1))),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**m["audit"], **(extra or {})), seconds=round(time.time() - t1, 1))
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
               f"nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} sel {n_sel} chgA {row['n_rows_changed_vs_a']} chgC {row['n_rows_changed_vs_c']} ({row['seconds']}s)")
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            d2 = max(abs(float(ka_ref[k]) - float(row[k])) for k in ("nmae_wald89", "delivered_norm") if k in ka_ref)
            ok = d < TOL and d2 < TOL
            row["ka"] = dict(ref=ka, max_abs_diff=float(max(d, d2)), ok=bool(ok))
            msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({max(d, d2):.1e})"
            if not ok:
                failures.append(f"{name} vs {ka}: {max(d, d2):.3e}")
        rows.append(row)
        print(msg, flush=True)
        return row

    # ------------------------------------------------------------------ R1: the incumbent through this path
    r_c = run("cd4svd5", ka="263 cd4svd5 (== 262 == 257 k5_s1 == 261 R1)", ka_ref=ref263["cd4svd5"])
    check("R1 cd4svd5 per-target pds bit-identical to 263", float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
    check("R1 cd4svd5 per-target cos bit-identical to 263", float(np.abs(per_target["cd4svd5|cos"] - pt263["cd4svd5|cos"]).max()), 0.0, 0.0)

    # ------------------------------------------------------------------ (a): the three-donor mean
    r_a = run("mean3", pre=st["pre_a"], extra=dict(construction="E062.build_baseline over the three donors, float32-cast"),
              **(dict(ka="cd4svd5 (kolf2: the cached pre IS the three-donor mean)", ka_ref=r_c) if held == "kolf2" else {}))
    if held == "kolf2":
        check("KA kolf2 mean3 per-target pds == cd4svd5", float(np.abs(per_target["mean3|pds"] - per_target["cd4svd5|pds"]).max()), 0.0, 0.0)
    else:
        Wall = np.full((n, J), np.nan)
        chg = np.abs(st["pre_a"] - st["pre"]).max(1) > 0
        Wall[chg & (avail3.sum(0) >= 1)] = 1.0
        run("mean3:nm", W=Wall, norm_match=True, norm_ref=st["pre"], extra=dict(construction="three-donor mean rescaled to the cached two-donor row norm on every changed covered row"))
        if not quick and st["kolf_unrestricted"] is not None:
            src_u = list(st["sources3"][:-1]) + [st["kolf_unrestricted"]]
            pre_u = Y.E062.build_baseline(src_u)[0].astype(np.float32).astype(np.float64)
            run("mean3:unrestricted", pre=pre_u, extra=dict(construction="as mean3 with kolf2's columns outside the operator axis kept (pass through Operator.apply un-whitened)"))

    # ------------------------------------------------------------------ (b) the majority rule and its controls
    def run_nm(name, drop, extra=None, nmc=False, asis=True):
        W = Z.drop_weights(n, J, drop)
        out = {}
        if asis:
            out["asis"] = run(name, W=W, pick=drop, extra=extra, **({"ka": f"263 {KA263_KOLF2[name]}", "ka_ref": ref263[KA263_KOLF2[name]]} if held == "kolf2" and name in KA263_KOLF2 and KA263_KOLF2[name] in ref263 else {}))
        nm = name + ":nm"
        out["nm"] = run(nm, W=W, pick=drop, extra=extra, norm_match=True,
                        **({"ka": f"263 {KA263_KOLF2[nm]}", "ka_ref": ref263[KA263_KOLF2[nm]]} if held == "kolf2" and nm in KA263_KOLF2 and KA263_KOLF2[nm] in ref263 else {}))
        if nmc and held != "kolf2":
            out["nmc"] = run(name + ":nmc", W=W, pick=drop, extra=dict(**(extra or {}), norm_ref="cached two-donor row"), norm_match=True, norm_ref=st["pre"])
        return out

    for g in GATES:
        d = Z.majority_drop(dd_raw_pairs, avail3, gates[g])
        run_nm(f"majority:raw:{g}", d, extra=dict(rule="majority", space="raw", gate=g, note="pick = the DROPPED donor"), nmc=(g == "all"))
        if held == "kolf2":
            check(f"R2 majority:raw:{g}:nm per-target pds bit-identical to 263", float(np.abs(per_target[f"majority:raw:{g}:nm|pds"] - pt263[f"majority:raw:{g}:nm|pds"]).max()), 0.0, 0.0)
        if quick:
            continue
        dr = Z.majority_drop(dd_raw_pairs, avail3, gates[g], reverse=True)
        run_nm(f"reverse:raw:{g}", dr, extra=dict(rule="reverse", space="raw", gate=g, control="reverse"), asis=(g == "all"))
        for k in range(3):
            drnd = Z.random_drop(avail3, gates[g], Z.SHUF_SEED + 101 * k + len(g))
            run_nm(f"random:{g}:s{k}", drnd, extra=dict(control="random", gate=g, seed=k), asis=(g == "all" and k == 0))
    if not quick:
        dw = Z.majority_drop(dd_wht_pairs, avail3, gates["all"])
        run_nm("majority:wht:all", dw, extra=dict(rule="majority", space="wht", gate="all"), asis=False)
        for j in range(J):
            run_nm(f"drop:{st['donors3'][j]}:all", Z.fixed_drop(avail3, gates["all"], j), extra=dict(control="global_drop", dropped=st["donors3"][j], gate="all"))
        if held != "kolf2":                                    # kolf2 ALONE on the three-donor targets (what the third donor is worth by itself)
            Wk = np.full((n, J), np.nan)
            Wk[gates["all"]] = 0.0
            Wk[gates["all"], J - 1] = 1.0
            run("only:kolf2:all:nm", W=Wk, pick=np.where(gates["all"], J - 1, -1), norm_match=True, extra=dict(control="kolf2 alone", gate="all"))
        for space, ct in (("raw", cos_truth_raw), ("wht", cos_truth_wht)):
            run_nm(f"oracle_drop:{space}:all", Z.oracle_drop(ct, avail3, gates["all"]), extra=dict(oracle=f"drop the donor with the lowest {space} cosine to the held truth", gate="all"), asis=False)

    cov = dict(targets=np.asarray(targets), avail3=avail3.astype(np.int8), dd_raw_pairs=dd_raw_pairs, dd_wht_pairs=dd_wht_pairs, dd_raw=dd_raw,
               cos_truth_raw=cos_truth_raw, cos_truth_wht=cos_truth_wht, base_a_cos_truth=base_a_cos_truth, base_c_cos_truth=base_c_cos_truth,
               votes3=st["votes3"], **{f"gate|{g}": m.astype(np.int8) for g, m in gates.items()})
    n_ok = sum(1 for r in ka_log if r["ok"])
    aud3 = {k: v for k, v in st["aud3"].items() if k not in ("kolf2_targets_not_in_fold", "fold_targets_kolf2_lacks")}
    meta = dict(fold=held, n_targets=n, targets=targets, donors=st["donors3"], sources2=st["sources"], J=J, pm_rows=st["pm_n"],
                frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"], deflation=st["aud5"], gates={g: int(m.sum()) for g, m in gates.items()},
                dd_raw_median=float(np.nanmedian(dd_raw)), pair_medians={f"{st['donors3'][a]}~{st['donors3'][b]}": float(np.nanmedian(dd_raw_pairs[:, k])) for k, (a, b) in enumerate(Z.PAIRS)},
                aud3=aud3, kolf2_targets_not_in_fold=st["aud3"].get("kolf2_targets_not_in_fold"), fold_targets_kolf2_lacks=st["aud3"].get("fold_targets_kolf2_lacks"),
                n_rows=len(rows), ka_rows=ka_log, ka_pass=n_ok, ka_total=len(ka_log), failures=failures, quick=quick,
                seconds=round(time.time() - t0, 1), peak_rss_gib=Z.peak_gib(), fold_cache=st["fold"]["audit"].get("cache"),
                operator=dict(name=st["op"].name, sha256=st["op"].sha256, n_cols=st["op"].n_cols), kolf_weight_mode=kw,
                whitening_audit=st["whitening_audit"])
    tag = ("quick264" if quick else "sweep264") + ("m" if kw == "matched" else "")
    (Z.WORK / f"{tag}-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=float))
    sfx = ("m" if kw == "matched" else "") + "-" + held + ("-quick" if quick else "")
    np.savez_compressed(Z.WORK / f"pertarget264{sfx}.npz", **per_target)
    np.savez_compressed(Z.WORK / f"identity264{sfx}.npz", **identity)
    np.savez_compressed(Z.WORK / f"picks264{sfx}.npz", **picks)
    np.savez_compressed(Z.WORK / f"cov264{sfx}.npz", **cov)
    print(f"[{held}] {len(rows)} rows, known-answer {n_ok}/{len(ka_log)} (+{sum(1 for r in rows if 'ka' in r and r['ka']['ok'])}/"
          f"{sum(1 for r in rows if 'ka' in r)} row KAs), {len(failures)} failures, {meta['seconds']}s, peak RSS {meta['peak_rss_gib']:.2f} GiB", flush=True)
    for f in failures:
        print(f"  *** {f}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
