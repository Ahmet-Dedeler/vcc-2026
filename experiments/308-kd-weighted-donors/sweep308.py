"""308 -- one PLOLO cancer fold: every |kd|-weighted-donor arm x every control on the cd4svd5 base.  MEASUREMENT.

    PYTHONDONTWRITEBYTECODE=1 uv run python sweep308.py <k562_gwps | orion_hct116 | orion_hek293t>

Every row rides cd4svd5's downstream verbatim (lib263: whitening through plolo.Operator.apply, CD4 deflated k = 5
re-solved at f = 0.34, the H1-pair union re-solved by the carrying-row rule at g = 0.28, one add, one preserve;
scored by lib246.score, 188's exact nmae identity through lib257.Identity).  Nothing is written outside
experiments/308-kd-weighted-donors/.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
import lib308 as Z  # noqa: E402

Y, E = Z.Y, Z.E
L, X, plolo = Y.L, Y.X, Y.plolo
E062, E064 = plolo.E062, plolo.E064
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
E263, E271 = Z.E263, Z.E271
E262 = ROOT / "experiments/262-where-pds-is-lost"
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main(held: str) -> None:
    t0 = time.time()
    P(f"[{held}] lib308 sha256 {Z.sha256_self()}")
    st = Y.stage(held)
    ident = X.Identity(st)
    n, J = st["n"], len(st["rows"])
    targets = list(st["targets"])
    donors = list(st["sources"])
    wc, sc = st["wcols"], st["sc"]
    avail = st["avail"]
    P(f"[{held}] staged {n} targets, donors {donors} ({st['donors']}), PM rows {st['pm_n']}, "
      f"KA-0 pre max|d| {st['ka0']['pre_float32_max_abs_diff']:.1e}, {time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB")

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
    check("KA-0 votes == per-donor availability", st["ka0"]["votes_max_abs_diff"], 0.0, 0.0)

    # ---------------------------------------------------------------- the feature on this fold
    feat = Z.load_feature()
    ff = Z.fold_feature(feat, targets, donors)
    A_raw, M = ff["A"], ff["M"]
    A = Z._fill_missing(A_raw, M, avail)
    two = avail.sum(0) >= 2
    P(f"[{held}] feature coverage on the fold: " + ", ".join(
        f"{d} measured {int((M[j] & avail[j]).sum())}/{int(avail[j].sum())} carried" for j, d in enumerate(donors)))
    cuts = dict(abs=ff["cut_abs"], rel=ff["cut_rel"])
    dis = {}
    for cname, cut in cuts.items():
        works = (A >= cut[:, None]) & avail
        fails = (A < cut[:, None]) & avail
        d_ = two & (works.sum(0) >= 1) & (fails.sum(0) >= 1)
        allfail = two & (works.sum(0) == 0)
        dis[cname] = d_
        P(f"[{held}/{cname}] cut {np.round(cut, 4).tolist()}: per-donor pass "
          + ", ".join(f"{d} {int(works[j].sum())}/{int(avail[j].sum())}" for j, d in enumerate(donors))
          + f"; two-donor targets {int(two.sum())}; DISAGREEING {int(d_.sum())}; all-fail (gate falls back) {int(allfail.sum())}")

    # ---------------------------------------------------------------- 271's reliabilities (the decisive control)
    pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
    HA, HB = [], []
    for j, line in enumerate(donors):
        z = np.load(E271 / f"work/half/rows-{line}.npz", allow_pickle=True)
        tz = [str(t) for t in z["targets"]]
        idx = np.asarray([tz.index(t) for t in targets], dtype=np.int64)
        HA.append(z["A"][idx].astype(np.float64))
        HB.append(z["B"][idx].astype(np.float64))
        R_ = z["real"][idx].astype(np.float64)
        prod = st["rows"][j]
        if line == "k562_gwps":
            shared = E064.shared_k562_pm_columns(E062.Bank.load(line, plolo.DONOR_PATHS[line]), pm_bank,
                                                 plolo.OFFICIAL, plolo.PANEL_SET)
            mask = np.ones(prod.shape[1], bool)
            mask[np.asarray(shared, dtype=np.int64)] = False
            d = float(np.abs(prod[avail[j]][:, mask].astype(np.float32) - R_[avail[j]][:, mask].astype(np.float32)).max())
        else:
            d = float(np.abs(prod[avail[j]].astype(np.float32) - R_[avail[j]].astype(np.float32)).max())
        check(f"KA-C {line} production row == 271's half-builder real row (float32)", d, 0.0, 0.0)
    del pm_bank
    r_raw = np.full((J, n), np.nan)
    for j in range(J):
        r_raw[j] = np.where(avail[j], E.row_cos(HA[j], HB[j]), np.nan)
    R_raw = E.spearman_brown(r_raw)
    W_relw = E.weights(R_raw, avail, "relw")
    W_relg = E.weights(R_raw, avail, "glob")
    del HA, HB

    # ---------------------------------------------------------------- the runner
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    wsave: dict[str, np.ndarray] = {}
    ref262 = json.loads((E262 / f"work/diag262-{held}.json").read_text())
    pt262 = np.load(E262 / f"work/diag262-{held}.npz")
    ref271 = {r["arm"]: r for r in json.loads((E271 / f"work/sweep271-{held}.json").read_text())["rows"]}

    def run(name, *, W=None, explicit_mean=False, norm_match=True, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        if explicit_mean or W is not None:
            pre_sel = Y.combine(st, W, explicit_mean=explicit_mean, norm_match=norm_match)
            base = Y.base_from(st, pre_sel)
        else:
            pre_sel, base = st["pre"], st["whitened"]
        e, m = Y.evaluate(st, base, ident)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        if W is not None:
            wsave[name] = W.astype(np.float32)
        n_changed = int((np.abs(pre_sel - st["pre"]).max(1) > 0).sum())
        ent = Z.entropy_of(W) if W is not None else np.full(n, np.nan)
        row = dict(fold=held, arm=name, n_targets=n, n_rows_changed=n_changed, norm_match=bool(norm_match),
                   n_weighted=int(np.isfinite(W).all(1).sum()) if W is not None else 0,
                   mean_entropy=float(np.nanmean(ent)) if np.isfinite(ent).any() else None,
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                   delivered_norm=m["delivered_norm"], audit=dict(**m["audit"], **(extra or {})),
                   seconds=round(time.time() - t1, 1))
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
               f"nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} w {row['n_weighted']} chg {n_changed} ({row['seconds']}s)")
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            ok = d < TOL
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(ok))
            msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({d:.1e})"
            if not ok:
                failures.append(f"{name} vs {ka}: {d:.3e}")
        rows.append(row)
        P(msg)
        return row

    # ---------------------------------------------------------------- known-answer rows
    r_base = run("cd4svd5", ka="262 cd4svd5 (== 257 k5_s1 == 261 R1 == 271)", ka_ref=ref262["arms"]["cd4svd5"])
    check("cd4svd5 per-target pds bit-identical to 262",
          float(np.abs(per_target["cd4svd5|pds"] - pt262["cd4svd5|pds"]).max()), 0.0, 0.0)
    check("cd4svd5 pds == 271's cd4svd5", r_base["pds"], ref271["cd4svd5"]["pds"], 0.0)
    check("cd4svd5 nmae_wald89 == 271's cd4svd5", r_base["nmae_wald89"], ref271["cd4svd5"]["nmae_wald89"], 0.0)
    run("mean_explicit", explicit_mean=True, ka="cd4svd5 (263 R1)", ka_ref=r_base)
    # P = 0 in the soft family is the equal-weight mean through the weighting path
    W0, _ = Z.soft_weights(A, avail, cuts["rel"], 0.0)
    run("soft:rel:p0", W=W0, ka="cd4svd5 (equal weights through the weight path)", ka_ref=r_base, extra=dict(power=0.0))

    # ---------------------------------------------------------------- the controls that must be beaten
    r_relw = run("relw", W=W_relw, extra=dict(control="271 reliability-weighted mean, w ∝ R/(1-R)"),
                 ka="271 relw", ka_ref=ref271["relw"])
    run("relw_glob", W=W_relg, extra=dict(control="271 donor-level scalar reliability (019/142 class)"),
        ka="271 relw_glob", ka_ref=ref271["relw_glob"])
    W_nr = Z.normrel_weights(st["rows"], avail, wc)
    run("normrel1.5", W=W_nr, extra=dict(control="097 route (a) a6, normrel^-1.5"))

    # ---------------------------------------------------------------- the arms
    arm_W: dict[str, np.ndarray] = {}
    for cname, cut in cuts.items():
        Wg, aud = Z.gate_weights(A, avail, cut)
        arm_W[f"gate:{cname}"] = Wg
        run(f"gate:{cname}", W=Wg, extra=dict(arm=1, cut=cname, **aud))
        Wg_d, aud_d = Z.gate_weights(A, avail, cut, gate=dis[cname])
        arm_W[f"gate:{cname}:dis"] = Wg_d
        run(f"gate:{cname}:dis", W=Wg_d, extra=dict(arm=3, cut=cname, **aud_d))
        Wg_r, aud_r = Z.gate_weights(A, avail, cut, reverse=True)
        run(f"gate:{cname}:rev", W=Wg_r, extra=dict(control="reversed gate: keep the FAILING donors", cut=cname, **aud_r))
        for p in Z.POWERS:
            Ws, auds = Z.soft_weights(A, avail, cut, p)
            arm_W[f"soft:{cname}:p{p}"] = Ws
            run(f"soft:{cname}:p{p}", W=Ws, extra=dict(arm=2, cut=cname, power=p, **auds))
        # arm 3 and the reverse for the best-behaved power only (p = 1 and p = 2), to keep the cell count honest
        for p in (1.0, 2.0):
            Ws_d, _ = Z.soft_weights(A, avail, cut, p, gate=dis[cname])
            arm_W[f"soft:{cname}:p{p}:dis"] = Ws_d
            run(f"soft:{cname}:p{p}:dis", W=Ws_d, extra=dict(arm=3, cut=cname, power=p))
            Ws_r, _ = Z.soft_weights(A, avail, cut, p, reverse=True)
            run(f"soft:{cname}:p{p}:rev", W=Ws_r, extra=dict(control="reversed soft weight", cut=cname, power=p))

    # ---------------------------------------------------------------- shuffled and entropy-matched controls
    for s in Z.SEEDS:
        Ash = Z.shuffle_kd(A, avail, s)                       # permute within donor over the CARRIED targets
        for cname, cut in cuts.items():
            Wg, aud = Z.gate_weights(Ash, avail, cut)
            run(f"gate:{cname}:shuf:s{s}", W=Wg, extra=dict(control="|kd| shuffled across targets within donor",
                                                            cut=cname, seed=s, **aud))
            Ws, _ = Z.soft_weights(Ash, avail, cut, 1.0)
            run(f"soft:{cname}:p1.0:shuf:s{s}", W=Ws, extra=dict(control="|kd| shuffled", cut=cname, seed=s, power=1.0))
    for s in Z.SEEDS:
        for key in ("gate:rel", "soft:rel:p1.0", "gate:abs"):
            Wr = Z.entropy_matched_random(arm_W[key], avail, 1000 + s)
            run(f"{key}:rand:s{s}", W=Wr, extra=dict(control=f"random weights, entropy matched to {key}", seed=s))

    # ---------------------------------------------------------------- the literal (amplitude-changing) reading
    for key in ("gate:rel", "soft:rel:p1.0", "gate:abs"):
        run(f"{key}:lit", W=arm_W[key], norm_match=False, extra=dict(note="literal, not norm-matched"))

    # ---------------------------------------------------------------- save
    meta = dict(fold=held, donors=donors, n_targets=n, lib308_sha256=Z.sha256_self(),
                cuts={k: np.asarray(v).tolist() for k, v in cuts.items()},
                n_two_donor=int(two.sum()),
                n_disagree={k: int(v.sum()) for k, v in dis.items()},
                feature_measured={d: int((M[j] & avail[j]).sum()) for j, d in enumerate(donors)},
                peak_gib=peak_gib(), seconds=round(time.time() - t0, 1), failures=failures, ka=ka_log)
    (WORK / f"sweep308-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=2))
    np.savez_compressed(WORK / f"pertarget308-{held}.npz", targets=np.asarray(targets), avail=avail.astype(np.int8),
                        abs_kd=A_raw, abs_kd_filled=A, measured=M.astype(np.int8),
                        cells=ff["cells"], R_raw=R_raw, **{f"gate|{k}": v.astype(np.int8) for k, v in dis.items()},
                        **per_target)
    np.savez_compressed(WORK / f"identity308-{held}.npz", **identity)
    np.savez_compressed(WORK / f"weights308-{held}.npz", **wsave)
    P(f"\n[{held}] {len(rows)} rows, {len(failures)} KA failures, {time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB")
    if failures:
        P("FAILURES:\n  " + "\n  ".join(failures))
    (WORK / f"sweep308-{held}.log").write_text("\n".join(LINES) + "\n")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main(sys.argv[1])
