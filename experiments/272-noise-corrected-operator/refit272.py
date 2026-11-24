"""272 -- construction B on ONE fold (orion_hct116): the cross-fitted operators.  Two modes, two processes:

  fit      (HEAVY, gated, serial; ONE operator per process -- E062's guard compares the process's sticky PEAK RSS
           with its limit, and the eigendecomposition's peak of ~14.5 GiB would trip every later fit in the same
           process, which is what happened on the first attempt)
           `fit-full <fold>`     every row of both FULL banks  -> KA-2: must reproduce the cached held operator
                                 (eigenvalues to 1e-6 relative, axis_sha256 equal, and the apply-equivalence on the fold's pre),
           `fit-half <fold> A|B` rebuilds the two donors' batch-parity half banks for that side from 267's sums
                                 (271's write_half_bank, verbatim in op272; KA-H against 271's saved rows) and fits
                                 the operator on every row of both half banks; `cleanup <fold>` removes the delta files.
           Writes work/refit/{full,halfA,halfB}-<fold>.{npz,json} and fit-{full,halfA,halfB}-<fold>.json.
  measure  (light) the anti-correlation test and the arms through cd4svd5's construction (263's staging), one
           operator resident at a time:
             cos(W_A A_t, W_A B_t) and cos(W_B A_t, W_B B_t) for each donor (one half in-sample, one out) against
             the fully in-sample cos(W A_t, W B_t) (KA-0) and the held line's halves through each operator;
             arms: cd4svd5 (KA), refit_full (KA-2 in score terms), op_A, op_B, xfit_avg (the production pre through
             the average of the two half operators, calibrated on pre), half:inc / half:xfit / half:infit /
             half:avgop / half:full (the half-row constructions: each half through the operator that never saw it,
             the operator that saw it, the average, the cached operator, the full refit; matched baseline half:inc).

    PYTHONDONTWRITEBYTECODE=1 uv run python refit272.py fit-full orion_hct116
    PYTHONDONTWRITEBYTECODE=1 uv run python refit272.py fit-half orion_hct116 A
    PYTHONDONTWRITEBYTECODE=1 uv run python refit272.py fit-half orion_hct116 B
    PYTHONDONTWRITEBYTECODE=1 uv run python refit272.py measure orion_hct116
    PYTHONDONTWRITEBYTECODE=1 uv run python refit272.py cleanup orion_hct116
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/272-noise-corrected-operator"
WORK = HERE / "work"
REFIT = WORK / "refit"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
sys.path.insert(0, str(E263))
sys.path.insert(0, str(HERE))
import lib263 as Y  # noqa: E402
import op272 as O  # noqa: E402

Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK
plolo, L, X = Y.plolo, Y.L, Y.X
E062 = plolo.E062
OFFICIAL, PANEL = plolo.OFFICIAL, list(plolo.PANEL)
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def panel_positions():
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    return np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)


# ============================================================================================ fit (one operator per process)
MAX_RSS_FIT = 24.0     # E062's guard on the PEAK RSS of this process; the machine gate (>= 20 GiB free+inactive) is checked before launch


def _fold_axis(held: str):
    spec = plolo.FOLDS[held]
    donors = list(spec["sources"])
    reg = plolo.OPERATOR_REGISTRY[f"held_{held}"]
    banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in donors]
    wc = plolo.E064.operator_columns(banks, OFFICIAL, plolo.PANEL_SET)
    axis = plolo.axis_sha256(OFFICIAL[i] for i in wc)
    assert axis == reg["axis_sha256"], (axis, reg["axis_sha256"])
    return donors, reg, banks, wc


def _checker(ka: list, failures: list):
    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
        if not ok:
            failures.append(label)
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")
    return check


def _raw_full(banks, check, pp) -> dict:
    """062's source_scale = prepare_main_source's raw_median_norm on the FULL banks (KA against 271's rows json)."""
    raw_full = {}
    for b in banks:
        prep = E062.prepare_main_source(b, PANEL, OFFICIAL, pp, chunk_rows=256, max_rss_gib=40.0)
        raw_full[b.name] = float(prep.raw_median_norm)
        j271 = json.loads((E271 / f"work/half/rows-{b.name}.json").read_text())
        check(f"raw median of {b.name} == 271 rows json", raw_full[b.name], j271["raw_median_real"], 1e-9)
        del prep
    return raw_full


def fit_full(held: str) -> None:
    t0 = time.time()
    REFIT.mkdir(parents=True, exist_ok=True)
    donors, reg, banks, wc = _fold_axis(held)
    cached_npz = ROOT / reg["path"]
    P(f"[fit-full {held}] op272.py sha256 {O.sha256_self()}; donors {donors}; cached operator {reg['sha256'][:12]}")
    ka: list[dict] = []
    failures: list[str] = []
    check = _checker(ka, failures)
    raw_full = _raw_full(banks, check, panel_positions())
    P(f"  raw medians {raw_full}, {time.time() - t0:.0f}s, peak {O.peak_gib():.2f} GiB")
    out_full = REFIT / f"full-{held}.npz"
    meta_full = O.fit_operator(banks, raw_full, wc, OFFICIAL, out_full, E062=E062, log=P, max_rss_gib=MAX_RSS_FIT)
    with np.load(cached_npz) as z:
        ev_c = np.sort(np.asarray(z["evals"], dtype=np.float64))
    with np.load(out_full) as z:
        ev_r = np.sort(np.asarray(z["evals"], dtype=np.float64))
    rel = np.abs(ev_r - ev_c) / np.maximum(np.abs(ev_c), 1e-300)
    big = ev_c > ev_c.max() * 1e-8
    check("KA-2a eigenvalues of the full-row refit == cached held operator (max rel diff, eigenvalues > 1e-8 max)", float(rel[big].max()), 0.0, 1e-6)
    check("KA-2a' the same on the whole spectrum incl. the floored tail (max rel diff)", float(rel.max()), 0.0, 1e-3)
    check("KA-2b axis_sha256 of the refit == registry", 0.0 if meta_full["axis_sha256"] == reg["axis_sha256"] else 1.0, 0.0, 0.0)
    check("KA-2c wcols_sha256 of the refit == 062's fit-time record", 0.0 if meta_full["wcols_sha256"] == reg["wcols_sha256"] else 1.0, 0.0, 0.0)
    check("KA-2d max eigenvalue == registry", meta_full["max_eigenvalue"], reg["max_eigenvalue"], 1e-6)
    P(f"  KA-2 eigenvalue comparison: max rel {rel[big].max():.2e} on {int(big.sum())} directions > 1e-8 max; median rel {np.median(rel):.2e}; "
      f"cached min {ev_c.min():.3e} refit min {ev_r.min():.3e}; {time.time() - t0:.0f}s")
    rec = dict(fold=held, donors=donors, cached=dict(path=reg["path"], sha256=reg["sha256"], axis_sha256=reg["axis_sha256"]), raw_medians_full=raw_full,
               operator=meta_full, eigen_comparison=dict(max_rel_big=float(rel[big].max()), n_big=int(big.sum()), max_rel_all=float(rel.max()), median_rel=float(np.median(rel))),
               ka=ka, ka_pass=sum(1 for k in ka if k["ok"]), ka_total=len(ka), failures=failures, op272_sha256=O.sha256_self(),
               seconds=round(time.time() - t0, 1), peak_rss_gib=O.peak_gib())
    (REFIT / f"fit-full-{held}.json").write_text(json.dumps(rec, indent=1, default=float))
    (REFIT / f"fit-full-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[fit-full {held}] done: KA {rec['ka_pass']}/{rec['ka_total']}, failures {failures}, {(time.time() - t0) / 60:.1f} min, peak {O.peak_gib():.2f} GiB")


def fit_half(held: str, side: str) -> None:
    t0 = time.time()
    REFIT.mkdir(parents=True, exist_ok=True)
    donors, reg, banks, wc = _fold_axis(held)
    P(f"[fit-half {held} {side}] op272.py sha256 {O.sha256_self()}; donors {donors}")
    ka: list[dict] = []
    failures: list[str] = []
    check = _checker(ka, failures)
    pp = panel_positions()
    raw_full = _raw_full(banks, check, pp)
    del banks
    raw_half = {}
    half_banks = []
    metas_bank = {}
    for n in donors:
        z271 = np.load(E271 / f"work/half/rows-{n}.npz", allow_pickle=True)
        assert [str(t) for t in z271["targets"]] == PANEL
        bd = WORK / "half" / n / f"bank-batch{side}"
        if not (bd / "delta.npy").exists():
            sm = O.Sums(n)
            metas_bank[n] = O.write_half_bank(sm, side, bd)
            del sm
        else:
            metas_bank[n] = json.loads((bd / "meta.json").read_text())
        hb = E062.Bank.load(f"{n}-batch{side}", bd)
        prep = E062.prepare_main_source(hb, PANEL, OFFICIAL, pp, chunk_rows=256, max_rss_gib=40.0)
        rows = (prep.effects * (prep.raw_median_norm / raw_full[n])).astype(np.float32)
        check(f"KA-H {n} half {side}: rebuilt half bank through prepare_main_source == 271's saved rows (float32)", float(np.abs(rows - z271[side]).max()), 0.0, 0.0)
        check(f"KA-H {n} half {side}: raw median == 271", prep.raw_median_norm, float(z271[f"raw_median_{side}"]), 1e-9)
        raw_half[hb.name] = float(prep.raw_median_norm)
        half_banks.append(hb)
        m = metas_bank[n]
        P(f"  half bank {n} {side}: {m['n_perts']} perts x {m['n_genes']} genes, cells/pert median {m['cells_per_pert']['median']:.0f}, raw median {prep.raw_median_norm:.6f}, "
          f"{time.time() - t0:.0f}s, peak {O.peak_gib():.2f} GiB")
        del prep, rows, z271
    wc_half = plolo.E064.operator_columns(half_banks, OFFICIAL, plolo.PANEL_SET)
    check("half banks' operator axis == the fold's axis (index digest)", 0.0 if np.array_equal(wc_half, wc) else 1.0, 0.0, 0.0)
    out = REFIT / f"half{side}-{held}.npz"
    meta = O.fit_operator(half_banks, raw_half, wc, OFFICIAL, out, E062=E062, log=P, max_rss_gib=MAX_RSS_FIT)
    rec = dict(fold=held, side=side, donors=donors, raw_medians_full=raw_full, raw_medians_half=raw_half, half_banks=metas_bank, operator=meta,
               ka=ka, ka_pass=sum(1 for k in ka if k["ok"]), ka_total=len(ka), failures=failures, op272_sha256=O.sha256_self(),
               seconds=round(time.time() - t0, 1), peak_rss_gib=O.peak_gib())
    (REFIT / f"fit-half{side}-{held}.json").write_text(json.dumps(rec, indent=1, default=float))
    (REFIT / f"fit-half{side}-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[fit-half {held} {side}] done: KA {rec['ka_pass']}/{rec['ka_total']}, failures {failures}, {(time.time() - t0) / 60:.1f} min, peak {O.peak_gib():.2f} GiB")


def cleanup(held: str) -> None:
    for n in plolo.FOLDS[held]["sources"]:
        for side in ("A", "B"):
            p = WORK / "half" / n / f"bank-batch{side}" / "delta.npy"
            if p.exists():
                p.unlink()
                print(f"removed {p.relative_to(ROOT)}", flush=True)


def fit_record(held: str) -> dict:
    """The three per-operator fit records merged into the shape `measure` reads."""
    full = json.loads((REFIT / f"fit-full-{held}.json").read_text())
    hA = json.loads((REFIT / f"fit-halfA-{held}.json").read_text())
    hB = json.loads((REFIT / f"fit-halfB-{held}.json").read_text())
    ka = full["ka"] + hA["ka"] + hB["ka"]
    return dict(fold=held, donors=full["donors"], cached=full["cached"], raw_medians_full=full["raw_medians_full"],
                raw_medians_half={"A": hA["raw_medians_half"], "B": hB["raw_medians_half"]},
                operators={"full": full["operator"], "halfA": hA["operator"], "halfB": hB["operator"]}, eigen_comparison=full["eigen_comparison"],
                ka=ka, ka_pass=sum(1 for k in ka if k["ok"]), ka_total=len(ka), failures=full["failures"] + hA["failures"] + hB["failures"],
                seconds=full["seconds"] + hA["seconds"] + hB["seconds"], peak_rss_gib=max(full["peak_rss_gib"], hA["peak_rss_gib"], hB["peak_rss_gib"]))


# ============================================================================================ measure
def measure(held: str) -> None:
    t0 = time.time()
    fitrec = fit_record(held)
    (REFIT / f"fit-{held}.json").write_text(json.dumps(fitrec, indent=1, default=float))
    st = Y.stage(held)
    ident = X.Identity(st)
    n = st["n"]
    targets = st["targets"]
    wc, sc = st["wcols"], st["sc"]
    avail = st["avail"]
    donors = list(st["sources"])
    pre = np.asarray(st["pre"], dtype=np.float64)
    P(f"[measure {held}] op272.py sha256 {O.sha256_self()}; staged {n} targets, donors {donors}; fit KA {fitrec['ka_pass']}/{fitrec['ka_total']}; peak {Y.peak_gib():.2f} GiB")
    failures: list[str] = []
    ka_log: list[dict] = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka_log.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
        if not ok:
            failures.append(f"{label}: |d| {d:.3e}")
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")

    check("KA-0a build_baseline(per-donor rows) float32 == cached pre", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)
    H = {line: O.load_half_rows(line, targets) for line in donors}
    for j, line in enumerate(donors):
        check(f"availability of {line}'s half rows == the fold's donor availability", int((H[line]["ok"] != avail[j]).sum()), 0, 0)
    Hh = O.load_half_rows(held, targets)
    # the half-row prediction: mean over the carrying donors of (A + B)/2; targets no donor carries keep the cached pre row
    pre_half = pre.copy()
    cnt = np.zeros(n, dtype=np.int64)
    acc = np.zeros_like(pre)
    for line in donors:
        ok = H[line]["ok"]
        acc[ok] += 0.5 * (H[line]["A"][ok] + H[line]["B"][ok])
        cnt += ok.astype(np.int64)
    cov = cnt > 0
    pre_half[cov] = acc[cov] / cnt[cov, None]
    c_half_pre = O.row_cos(pre_half[:, wc], pre[:, wc])
    P(f"  pre_half vs cached pre: {int(cov.sum())}/{n} targets covered by a half row; row cosine median {np.nanmedian(c_half_pre):.4f} (min {np.nanmin(c_half_pre):.4f}); "
      f"the K562 rows are un-blended here (PerturbMatch 0.5 on {st['pm_n']} rows in the cache)")

    # ---- transforms through each operator, one resident at a time ----
    names = {"inc": None, "full": REFIT / f"full-{held}.npz", "A": REFIT / f"halfA-{held}.npz", "B": REFIT / f"halfB-{held}.npz"}
    T: dict[str, dict[str, np.ndarray]] = {}
    scale: dict[str, dict[str, float]] = {}
    inputs = {"pre": pre[:, wc] * ALPHA, "pre_half": pre_half[:, wc] * ALPHA, "held_A": Hh["A"][:, wc] * ALPHA, "held_B": Hh["B"][:, wc] * ALPHA}
    for line in donors:
        inputs[f"A|{line}"] = H[line]["A"][:, wc] * ALPHA
        inputs[f"B|{line}"] = H[line]["B"][:, wc] * ALPHA
    for key, path in names.items():
        t1 = time.time()
        if key == "inc":
            dop = O.DampedOperator.from_operator(st["op"], label="incumbent (cached held operator)")
        else:
            dop = O.load_fitted(path, key)
        T[key] = {k: dop.transform(v, GAMMA, LAM) for k, v in inputs.items()}
        # calibration scalars: production's preserve for the reference rows (pre, pre_half)
        scale[key] = {ref: float(np.median(np.linalg.norm(inputs[ref], axis=1)) / np.median(np.linalg.norm(T[key][ref], axis=1))) for ref in ("pre", "pre_half")}
        if key == "inc":
            base_inc = O.preserve(pre * ALPHA, T["inc"]["pre"], wc)
            check("KA-1 transform + preserve through the cached operator == the staged base (bit-identical)", float(np.abs(base_inc - st["whitened"]).max()), 0.0, 0.0)
            del st["op"], st["wht_rows"]
        P(f"  transforms through {key}: {len(inputs)} matrices, scale on pre {scale[key]['pre']:.6f}, {time.time() - t1:.0f}s, peak {Y.peak_gib():.2f} GiB")
        del dop

    # ---- KA-2 in apply terms: the full refit's whitened base vs the cached operator's ----
    base_full = O.preserve(pre * ALPHA, T["full"]["pre"], wc)
    rel = np.linalg.norm((base_full - st["whitened"])[:, wc], axis=1) / np.linalg.norm(st["whitened"][:, wc], axis=1)
    check("KA-2e full-row refit applied to the fold's pre == the cached operator's base (max row rel diff)", float(rel.max()), 0.0, 1e-4)
    P(f"  KA-2e: median row rel diff {np.median(rel):.2e}, row cosine min {np.nanmin(O.row_cos(base_full[:, wc], st['whitened'][:, wc])):.9f}")

    # ---- the anti-correlation test ----
    anti = {}
    diag271 = json.loads((E271 / "work/diagwht271.json").read_text())[held]
    for line in donors:
        ok = H[line]["ok"]
        rec = {}
        for key in ("inc", "full", "A", "B"):
            WA, WB = T[key][f"A|{line}"][ok], T[key][f"B|{line}"][ok]
            c = O.row_cos(WA, WB)
            pi = O.derangement(int(ok.sum()), 1)
            null = O.row_cos(WA, WB[pi])
            rec[key] = dict(n=int(ok.sum()), median=float(np.nanmedian(c)), mean=float(np.nanmean(c)), share_le0=float(np.nanmean(c <= 0)),
                            p10=float(np.nanpercentile(c, 10)), p90=float(np.nanpercentile(c, 90)), null_median=float(np.nanmedian(null)), null_sd=float(np.nanstd(null)))
        craw = O.row_cos(H[line]["A"][ok][:, wc], H[line]["B"][ok][:, wc])
        rec["raw"] = dict(median=float(np.nanmedian(craw)), share_le0=float(np.nanmean(craw <= 0)))
        anti[line] = rec
        P(f"  anti-correlation {line:14s} (n {rec['inc']['n']}): raw {rec['raw']['median']:+.4f}; cached op (both halves in-sample) {rec['inc']['median']:+.4f} "
          f"(share<=0 {rec['inc']['share_le0']:.3f}); full refit {rec['full']['median']:+.4f}; W_A (A in, B out) {rec['A']['median']:+.4f} (share<=0 {rec['A']['share_le0']:.3f}, "
          f"null sd {rec['A']['null_sd']:.3f}); W_B (B in, A out) {rec['B']['median']:+.4f} (share<=0 {rec['B']['share_le0']:.3f})")
    okh = Hh["ok"]
    rec = {}
    for key in ("inc", "full", "A", "B"):
        c = O.row_cos(T[key]["held_A"][okh], T[key]["held_B"][okh])
        rec[key] = dict(n=int(okh.sum()), median=float(np.nanmedian(c)), share_le0=float(np.nanmean(c <= 0)), p10=float(np.nanpercentile(c, 10)), p90=float(np.nanpercentile(c, 90)))
    anti[held] = dict(role="held (out-of-sample for every operator)", **rec)
    P(f"  held line {held}: whitened split-half median through cached {rec['inc']['median']:+.4f} / full refit {rec['full']['median']:+.4f} / W_A {rec['A']['median']:+.4f} / W_B {rec['B']['median']:+.4f}")
    # KA-0 on the fold's target set is not 271's exact set (271 used all 300 panel rows); the sweep carries the exact KA-0.  Here: sign and size.
    for line in donors:
        check(f"KA-0 (fold targets) {line}: cached-operator anti-correlation within 0.01 of 271's", anti[line]["inc"]["median"], diag271[line]["wht"]["median"], 1e-2)

    # ---- the arms ----
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    ref263 = {r["arm"]: r for r in json.loads((E263 / f"work/sweep263-{held}.json").read_text())["rows"]}
    pt263 = np.load(E263 / f"work/pertarget263-{held}.npz")

    def run(name, base, *, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        parts, aud = Y.arm_parts(st, base)
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
        row = dict(fold=held, arm=name, n_targets=n, **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**aud, **(extra or {})), seconds=round(time.time() - t1, 1))
        msg = (f"  {name:14s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f} "
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

    def combine_halves(op_for_A: str, op_for_B: str) -> np.ndarray:
        """mean over the carrying donors of (A_j through op_for_A + B_j through op_for_B)/2, each transform calibrated on
        pre_half; targets no donor carries take the cached-operator row of pre_half; ONE preserve to pre_half's norm."""
        acc = np.zeros((n, wc.size))
        cnt = np.zeros(n, dtype=np.int64)
        for line in donors:
            ok = H[line]["ok"]
            acc[ok] += 0.5 * (T[op_for_A][f"A|{line}"][ok] * scale[op_for_A]["pre_half"] + T[op_for_B][f"B|{line}"][ok] * scale[op_for_B]["pre_half"])
            cnt += ok.astype(np.int64)
        cov = cnt > 0
        acc[cov] /= cnt[cov, None]
        acc[~cov] = T["inc"]["pre_half"][~cov] * scale["inc"]["pre_half"]
        return O.preserve(pre_half * ALPHA, acc, wc)

    run("cd4svd5", st["whitened"], ka="263 cd4svd5", ka_ref=ref263["cd4svd5"], extra=dict(control="the incumbent"))
    check("KA-1b cd4svd5 per-target pds bit-identical to 263", float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
    r_full = run("refit_full", base_full, extra=dict(ka="KA-2 in score terms: the full-row refit through cd4svd5's construction"))
    d_full = max(abs(float(r_full[k]) - float(ref263["cd4svd5"][k])) for k in KEYS)
    check("KA-2f refit_full == cd4svd5 on (cos_agg, macro, pds) within 1e-5 (eigensolver / platform precision)", d_full, 0.0, 1e-5)
    run("op_A", O.preserve(pre * ALPHA, T["A"]["pre"], wc), extra=dict(arm="the half-A operator alone on the production pre"))
    run("op_B", O.preserve(pre * ALPHA, T["B"]["pre"], wc), extra=dict(arm="the half-B operator alone on the production pre"))
    run("xfit_avg", O.preserve(pre * ALPHA, 0.5 * (T["A"]["pre"] * scale["A"]["pre"] + T["B"]["pre"] * scale["B"]["pre"]), wc),
        extra=dict(arm="B: the production pre through (W_A + W_B)/2, each calibrated on pre, one preserve (the mandate's cross-fitted operator)"))
    run("half:inc", O.preserve(pre_half * ALPHA, T["inc"]["pre_half"], wc), extra=dict(baseline="the half-row prediction (A+B)/2 through the cached operator: the matched baseline"))
    run("half:full", O.preserve(pre_half * ALPHA, T["full"]["pre_half"], wc), extra=dict(baseline="the half-row prediction through the full-row refit"))
    run("half:avgop", O.preserve(pre_half * ALPHA, 0.5 * (T["A"]["pre_half"] * scale["A"]["pre_half"] + T["B"]["pre_half"] * scale["B"]["pre_half"]), wc),
        extra=dict(arm="the half-row prediction through (W_A + W_B)/2"))
    run("half:xfit", combine_halves("B", "A"), extra=dict(arm="B: the PROPER cross-fit: A_t through W_B, B_t through W_A, averaged -- no half sees the operator fitted on it"))
    run("half:infit", combine_halves("A", "B"), extra=dict(control="the in-sample counterpart: A_t through W_A, B_t through W_B"))

    meta = dict(fold=held, n_targets=n, targets=targets, donors=donors, op272_sha256=O.sha256_self(), refit_sha256=O.sha256_file(Path(__file__)),
                fit=dict(ka_pass=fitrec["ka_pass"], ka_total=fitrec["ka_total"], failures=fitrec["failures"], operators={k: dict(sha256=v["sha256"], max_eigenvalue=v["max_eigenvalue"], basis_rows=v["basis_rows_by_source"], seconds=v["seconds"], peak_rss_gib=v["peak_rss_gib"]) for k, v in fitrec["operators"].items()},
                         eigen_comparison=fitrec["eigen_comparison"], seconds=fitrec["seconds"], peak_rss_gib=fitrec["peak_rss_gib"]),
                scale=scale, anti_correlation=anti, pre_half=dict(covered=int(cov.sum()), row_cos_to_pre_median=float(np.nanmedian(c_half_pre)), row_cos_to_pre_min=float(np.nanmin(c_half_pre))),
                ka=ka_log, ka_pass=sum(1 for k in ka_log if k["ok"]), ka_total=len(ka_log), failures=failures, n_rows=len(rows),
                seconds=round(time.time() - t0, 1), peak_rss_gib=Y.peak_gib())
    (WORK / f"refit272-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=str))
    np.savez_compressed(WORK / f"pertarget272-refit-{held}.npz", targets=np.asarray(targets), avail=avail.astype(np.int8), **per_target)
    np.savez_compressed(WORK / f"identity272-refit-{held}.npz", **identity)
    (WORK / f"refit272-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[measure {held}] done: {len(rows)} rows, KA {meta['ka_pass']}/{meta['ka_total']}, failures {failures}, {(time.time() - t0) / 60:.1f} min, peak {Y.peak_gib():.2f} GiB")


if __name__ == "__main__":
    mode, held = sys.argv[1], sys.argv[2]
    if mode == "fit-full":
        fit_full(held)
    elif mode == "fit-half":
        fit_half(held, sys.argv[3])
    elif mode == "measure":
        measure(held)
    elif mode == "cleanup":
        cleanup(held)
    else:
        raise SystemExit(f"unknown mode {mode}")
