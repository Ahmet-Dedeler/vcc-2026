"""286 STAGE B -- ONE PLOLO FOLD PER PROCESS: the cross-product estimators through 271's plumbing on the cd4svd5 base,
scored on the fold's CACHE truth.  MEASUREMENT ONLY.

    PYTHONDONTWRITEBYTECODE=1 uv run python sweepB286.py orion_hct116

Every row arm replaces the donors' un-whitened production rows by their filtered versions (fit286.py, per line) and
then rides cd4svd5 verbatim: whitening through `plolo.Operator.apply`, CD4 deflated k = 5 re-solved at f = 0.34, the
H1-pair union re-solved by the carrying-row rule at g = 0.28, one add, one preserve (lib263 / lib246), 188's exact
`nmae` identity through `lib257.Identity`.  E3 replaces the OPERATOR instead (pooled signal covariance + floor, in the
operator's place); it is evaluated last, after the fitted operator has been released.

Known-answer rows: cd4svd5 == 263; mean_explicit == cd4svd5; STACK == 257 (the BR249 retrodiction's input);
relw / ebrow == 271; the per-line filtered E0 rows == the fold's staged donor rows bit for bit (KA-C).
"""
from __future__ import annotations

import gc
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/286-crossproduct-estimator"
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E257 = ROOT / "experiments/257-soft-deflation-price"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E272 = ROOT / "experiments/272-noise-corrected-operator"
for p in (E263, E271, E272, HERE):
    sys.path.insert(0, str(p))
import lib263 as Y  # noqa: E402
import est271 as E271M  # noqa: E402
import op272 as O  # noqa: E402
import est286 as E  # noqa: E402

Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK
plolo, L, X = Y.plolo, Y.L, Y.X
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
KEYS = ("cos_agg", "a23_cos", "pds")
TOL = 1e-9
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main(held: str) -> None:
    t0 = time.time()
    P(f"[{held}] est286.py sha256 {E.sha256_self()}")
    st = Y.stage(held)
    ident = X.Identity(st)
    n, J = st["n"], len(st["rows"])
    targets = list(st["targets"])
    wc, sc, avail = st["wcols"], st["sc"], st["avail"]
    donors = list(st["sources"])
    P(f"[{held}] staged {n} targets, donors {donors} ({st['donors']}), PM rows {st['pm_n']}, wcols {wc.size}, "
      f"KA-0 pre max|d| {st['ka0']['pre_float32_max_abs_diff']:.1e}, {time.time()-t0:.0f}s, peak {rss():.2f} GiB")
    failures: list[str] = []
    ka_log: list[dict] = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka_log.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
        if not ok:
            failures.append(f"{label}: |d| {d:.3e}")
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")
        return ok

    check("KA-0 build_baseline(per-donor rows) float32 == cached pre", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)

    # ------------------------------------------------------------------ the per-line filtered rows
    filt = {}
    for j, line in enumerate(donors):
        z = np.load(WORK / "filt" / f"{line}.npz", allow_pickle=True)
        tz = [str(t) for t in z["targets"]]
        idx = np.asarray([tz.index(t) for t in targets], dtype=np.int64)
        filt[line] = {k: np.asarray(z[k])[idx] for k in z.files
                      if k not in ("topvec", "lam", "nu", "gains", "cols", "targets")}
        filt[line]["_meta"] = json.loads((WORK / "filt" / f"{line}.json").read_text())
        d = float(np.abs(filt[line]["E0"].astype(np.float64) - st["rows"][j]).max())
        check(f"KA-C {line} fit286 E0 rows == the fold's staged production rows", d, 0.0, 1e-6)
        dav = int((np.asarray(filt[line]["avail"], bool) != avail[j]).sum())
        check(f"KA-C {line} availability == the fold's donor availability", dav, 0, 0)
    arm_names = sorted(k for k in filt[donors[0]] if k not in ("E0", "avail", "_meta"))
    P(f"  filtered arms available: {arm_names}")

    # 271's half rows (for the E4 baseline and 272's noise spectrum)
    HA, HB = [], []
    for line in donors:
        z = np.load(E271 / f"work/half/rows-{line}.npz", allow_pickle=True)
        tz = [str(t) for t in z["targets"]]
        i2 = np.asarray([tz.index(t) for t in targets], dtype=np.int64)
        HA.append(z["A"][i2].astype(np.float64))
        HB.append(z["B"][i2].astype(np.float64))

    # ------------------------------------------------------------------ the scoring path (sweep271.run, verbatim)
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    ref263 = {r["arm"]: r for r in json.loads((E263 / f"work/sweep263-{held}.json").read_text())["rows"]}
    pt263 = np.load(E263 / f"work/pertarget263-{held}.npz")
    ref257 = {r["arm"]: r for r in json.loads((E257 / f"work/sweep257-{held}.json").read_text())["rows"]}
    ref271 = {r["arm"]: r for r in json.loads((E271 / f"work/sweep271-{held}.json").read_text())["rows"]}

    def run(name, *, W=None, pre_sel=None, explicit_mean=False, base=None, parts=None, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        if base is None:
            if pre_sel is None:
                if explicit_mean:
                    pre_sel = Y.combine(st, None, explicit_mean=True)
                elif W is not None:
                    pre_sel = Y.combine(st, W, norm_match=True)
            base = st["whitened"] if pre_sel is None else Y.base_from(st, pre_sel)
        if parts is None:
            parts, aud = Y.arm_parts(st, base)
        else:
            aud = dict(cd4_fraction=float(parts[0][1]), union_fraction=float(parts[1][1]), note="STACK")
        st_arm = dict(st, whitened=base)
        final = L.compose(st_arm, parts)
        e = L.score(st_arm, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        chg = 0 if pre_sel is None else int((np.abs(pre_sel - st["pre"]).max(1) > 0).sum())
        aud["base_median_norm_wcols"] = float(np.median(np.linalg.norm(base[:, wc], axis=1)))
        row = dict(fold=held, arm=name, n_targets=n, n_rows_changed=chg,
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**aud, **(extra or {})), seconds=round(time.time() - t1, 1))
        msg = (f"  {name:26s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
               f"nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} chg {chg:3d} ({row['seconds']}s)")
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            ok = d < TOL
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(ok))
            ka_log.append(dict(label=f"{name} == {ka} on {KEYS}", got=float(e["pds"]), want=float(ka_ref["pds"]),
                               abs_diff=float(d), tol=TOL, ok=bool(ok)))
            msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({d:.1e})"
            if not ok:
                failures.append(f"{name} vs {ka}: {d:.3e}")
        rows.append(row)
        P(msg)
        return row

    # ---------------- known answers
    r_base = run("cd4svd5", ka="263 cd4svd5", ka_ref=ref263["cd4svd5"], extra=dict(control="the incumbent"))
    check("cd4svd5 per-target pds bit-identical to 263",
          float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
    run("mean_explicit", explicit_mean=True, ka="cd4svd5", ka_ref=r_base)
    run("STACK", parts=[(st["cd4"], st["frac_cd4"]), (st["union"], st["frac_x"])], ka="257 STACK", ka_ref=ref257["STACK"],
        extra=dict(control="h1union's stack: the BR249 retrodiction row"))
    # 271's two estimator rows through THIS plumbing
    r_raw = np.full((J, n), np.nan)
    for j in range(J):
        r_raw[j] = np.where(avail[j], E271M.row_cos(HA[j], HB[j]), np.nan)
    R_raw = E271M.spearman_brown(r_raw)
    run("relw", W=E271M.weights(R_raw, avail, "relw"), ka="271 relw", ka_ref=ref271["relw"])
    run("ebrow", W=E271M.weights(R_raw, avail, "ebrow"), ka="271 ebrow", ka_ref=ref271["ebrow"])
    # E4: 271's diagonal per-gene Wiener (the already-failed baseline)
    rows_f = [E271M.wiener_own(st["rows"][j], HA[j] - HB[j]) for j in range(J)]
    run("E4:wiener_own", pre_sel=E271M.combine_filtered(st, rows_f), ka="271 wiener_own", ka_ref=ref271["wiener_own"],
        extra=dict(baseline="271 (c) diagonal per-gene Wiener"))
    del rows_f
    gc.collect()

    # ---------------- the cross-product arms (both `:nm` and `:lit`)
    def combine_filtered(key: str, *, norm_match: bool) -> np.ndarray:
        pre = st["pre"].copy()
        acc = np.zeros_like(pre)
        cnt = np.zeros(n, dtype=np.int64)
        for j, line in enumerate(donors):
            mk = avail[j]
            acc[mk] += filt[line][key][mk].astype(np.float64)
            cnt += mk.astype(np.int64)
        cov = cnt > 0
        acc[cov] /= cnt[cov, None]
        if norm_match:
            acc[cov] = E.norm_match(acc[cov], pre[cov])
        pre[cov] = acc[cov].astype(np.float32).astype(np.float64)
        return pre

    for key in arm_names:
        for tag, nm in ((":nm", True), (":lit", False)):
            if tag == ":lit" and (":" in key.split("phi")[-1] or key.startswith("e2")):
                continue                                    # the literal reading only for the primary E1 cells
            run(key + tag, pre_sel=combine_filtered(key, norm_match=nm),
                extra=dict(estimator=key, norm_matched=nm))

    # ---------------- mechanism: principal angles between S-hat's leading directions and the fitted operator's
    mech: dict = {}
    ev_ord = np.argsort(st["op"].evals)[::-1]
    Vop = np.asarray(st["op"].evecs[:, ev_ord[:100]], dtype=np.float64)
    pos_w = {int(c): i for i, c in enumerate(wc.tolist())}
    for line in donors:
        zf = np.load(WORK / "filt" / f"{line}.npz", allow_pickle=True)
        cj, TV = zf["cols"], np.asarray(zf["topvec"], dtype=np.float64)
        keep = np.asarray([i for i, c in enumerate(cj) if int(c) in pos_w], dtype=np.int64)
        dest = np.asarray([pos_w[int(cj[i])] for i in keep], dtype=np.int64)
        A = np.zeros((wc.size, TV.shape[1]))
        A[dest] = TV[keep]
        A, _ = np.linalg.qr(A)
        pa = E.principal_angles(A, Vop, 100)
        lam_j = np.asarray(zf["lam"], dtype=np.float64)
        mech[line] = dict(principal_cos_top100=[float(x) for x in pa[:10]], mean_principal_cos=float(pa.mean()),
                          median_principal_cos=float(np.median(pa)),
                          top100_share_of_positive_energy=float(np.sort(np.clip(lam_j, 0, None))[::-1][:100].sum()
                                                                / np.clip(lam_j, 0, None).sum()),
                          share_negative=float((lam_j < 0).mean()),
                          clipped_energy=float(np.abs(lam_j[lam_j < 0]).sum() / np.abs(lam_j).sum()))
        P(f"  mechanism {line}: principal cosines of S-hat's top-100 vs the operator's top-100 -- median "
          f"{mech[line]['median_principal_cos']:.4f}, mean {mech[line]['mean_principal_cos']:.4f}, top "
          f"{pa[0]:.4f}; S-hat top-100 carry {100*mech[line]['top100_share_of_positive_energy']:.1f} % of its "
          f"positive energy; {100*mech[line]['share_negative']:.1f} % of its eigenvalues are negative")
    del Vop
    gc.collect()

    # ---------------- E5: 272's eigenvalue floor in the FITTED operator's basis (c = 1)
    try:
        rmeds = {ln: json.loads((E271 / f"work/half/rows-{ln}.json").read_text())["raw_median_real"] for ln in donors}
        bf, _ = O.basis_unit_factor(rmeds)
        Ds = [O.HALF_FACTOR * 0.5 * (HA[j] - HB[j])[avail[j]][:, wc] * bf for j in range(J)]
        wts = [float(np.sum(avail[j])) for j in range(J)]
        nspec, s_i, _ = O.noise_spectrum(st["op"].evals, st["op"].evecs, Ds, wts)
        dop = O.DampedOperator.from_operator(st["op"], None, 0.0, "e5-floor-c1", n=nspec, floor_c=1.0)
        base5 = dop.apply(st["pre"] * ALPHA, wc, GAMMA, LAM)
        run("E5:opfloor_c1", base=base5, extra=dict(baseline="272 A' eigenvalue floor at the noise level, c = 1",
                                                    share_noise_dominated=float((s_i < 0.5).mean())))
        del Ds, base5, dop, nspec, s_i
        gc.collect()
    except Exception as exc:                                                             # noqa: BLE001
        P(f"  E5 SKIPPED: {exc!r}")
        failures.append(f"E5 skipped: {exc!r}")

    # ---------------- E3: the signal-covariance whitener in the OPERATOR's place (operator released first)
    e3_aud = {}
    try:
        del st["op"], st["wht_rows"]
        gc.collect()
        P(f"  operator released; rss {rss():.2f} GiB -- pooling S-hat / N-hat on the {wc.size} operator columns")
        pos_w = {int(c): i for i, c in enumerate(wc.tolist())}
        Sp = np.zeros((wc.size, wc.size), dtype=np.float32)
        Np = np.zeros((wc.size, wc.size), dtype=np.float32)
        wsum = 0.0
        for line in donors:
            cj = np.load(WORK / "moments" / f"{line}_cols.npy")
            keep = np.asarray([i for i, c in enumerate(cj) if int(c) in pos_w], dtype=np.int64)
            dest = np.asarray([pos_w[int(cj[i])] for i in keep], dtype=np.int64)
            wj = float(json.loads((WORK / "moments" / f"{line}.json").read_text())["n_rows"])
            for buf, nm2 in ((Sp, "S"), (Np, "N")):
                A = np.load(WORK / "moments" / f"{line}_{nm2}.npy")
                buf[np.ix_(dest, dest)] += np.float32(wj) * A[np.ix_(keep, keep)]
                del A
                gc.collect()
            wsum += wj
            e3_aud[line] = dict(cols_used=int(keep.size), of=int(cj.size), weight=wj)
        Sp /= np.float32(wsum)
        Np /= np.float32(wsum)
        lam3, V3, a3 = E.psd_eigh(Sp, log=P)
        del Sp
        gc.collect()
        nu3 = E.direction_noise(Np, V3)
        del Np
        gc.collect()
        e3_aud["spectrum"] = a3
        e3_aud["share_lam_below_nu"] = float((np.clip(lam3, 0, None) < nu3).mean())
        # the mechanism reading: where the INCUMBENT's whitened prediction energy sits relative to S-hat's spectrum
        Zw = np.asarray(st["whitened"][:, wc], dtype=np.float32) @ V3
        ew = np.asarray((Zw * Zw).sum(0), dtype=np.float64)
        del Zw
        gc.collect()
        below = np.clip(lam3, 0, None) < nu3
        ordl = np.argsort(lam3)[::-1]
        e3_aud["whitened_energy_share_where_lam_below_nu"] = float(ew[below].sum() / ew.sum())
        e3_aud["whitened_energy_share_top100_of_S"] = float(ew[ordl[:100]].sum() / ew.sum())
        e3_aud["whitened_energy_share_negative_lam"] = float(ew[lam3 < 0].sum() / ew.sum())
        P(f"  mechanism (pooled S-hat on the operator columns): {100*e3_aud['share_lam_below_nu']:.1f} % of the "
          f"directions have lam+ < nu, and they carry {100*e3_aud['whitened_energy_share_where_lam_below_nu']:.1f} % "
          f"of the incumbent's whitened prediction energy; S-hat's top-100 carry "
          f"{100*e3_aud['whitened_energy_share_top100_of_S']:.1f} %; the negative-eigenvalue directions "
          f"{100*e3_aud['whitened_energy_share_negative_lam']:.1f} %")
        for c in E.E3_CS:
            w3 = E.e3_weights(lam3, nu3, float(c))
            T = E.apply_e3(st["pre"][:, wc] * ALPHA, V3, w3)
            base3 = O.preserve(st["pre"] * ALPHA, T, wc)
            run(f"E3:c{c:g}", base=base3, extra=dict(estimator="signal-covariance whitener in the operator's place",
                                                     c=float(c), share_lam_below_nu=e3_aud["share_lam_below_nu"]))
            del T, base3
            gc.collect()
        np.savez_compressed(WORK / f"e3spec-{held}.npz", lam=lam3, nu=nu3, topvec=np.ascontiguousarray(V3[:, ::-1][:, :100]))
        del V3, lam3, nu3
        gc.collect()
    except Exception as exc:                                                             # noqa: BLE001
        P(f"  E3 SKIPPED: {exc!r}")
        failures.append(f"E3 skipped: {exc!r}")

    meta = dict(fold=held, n_targets=n, targets=targets, donors=donors, donor_names=st["donors"], J=J,
                est286_sha256=E.sha256_self(), sweep_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                wcols=int(wc.size), frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"], e3=e3_aud,
                filt_meta={ln: filt[ln]["_meta"] for ln in donors}, mechanism=mech, ka=ka_log,
                ka_pass=sum(1 for k in ka_log if k["ok"]), ka_total=len(ka_log), failures=failures,
                n_rows=len(rows), seconds=round(time.time() - t0, 1), peak_rss_gib=rss(),
                reserve_note="267's reserved orion_hek293t half-B truths were not read; the Broad-16 reserve lines "
                             "were not read")
    (WORK / f"sweepB286-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=str))
    np.savez_compressed(WORK / f"pertargetB286-{held}.npz", targets=np.asarray(targets), avail=avail.astype(np.int8),
                        base262_cos_truth=np.asarray(np.load(E263 / f"work/cov263-{held}.npz")["base262_cos_truth"]),
                        ref_rel=np.asarray(np.load(ROOT / f"experiments/267-replicate-validated-oracle/work/half/truth-{held}.npz",
                                                   allow_pickle=True)["rel|batch|cos_resid"]),
                        **per_target)
    np.savez_compressed(WORK / f"identityB286-{held}.npz", **identity)
    (WORK / f"sweepB286-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{held}] done: {len(rows)} rows, KA {meta['ka_pass']}/{meta['ka_total']}, failures {failures}, "
      f"{(time.time()-t0)/60:.1f} min, peak {rss():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1])
