"""296 -- CD4 GWCD4i (CRISPRi, primary human T cells): the knockdown-stratified construct term.

280's bank has no cells, so every side is a UMI budget T and the cell-sampling reference is a binomial
COUNT split (Poisson-only, therefore optimistic -- 285 sec 6 measured the count-minus-cell offset at
+0.0438 on the two sources that have both; it is applied nowhere here, only reported).

Per (condition, donor, target) with two guide slots, at budget T:
  (a)  count split of guide slot 0        (a2) count split of guide slot 1
  (b)  slot 0 vs slot 1                   (x)  slot 0 of t vs slot 1 of the NEXT target
Per-target values are averaged over the four donors inside a condition; the knockdown of each guide is its
own gene's column at FULL depth, pooled over the condition's donors, against that condition's NTC pool.

    PYTHONDONTWRITEBYTECODE=1 uv run python cd4_296.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path("~/Code/virtual-cell-challenge-2026/experiments/"
                            "280-gwcd4i-donor-replicates/work")))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402
import lib280 as B280  # noqa: E402

LIB = "CD4-GWCD4i"
E280 = "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work"
P = L.Logger(X.WORK / "cd4_296.log")
BUDGETS = (1e5, 2e5)
NREP = 8
CONDS = ("Rest", "Stim8hr", "Stim48hr")
MIN_COHORT = 15


def main() -> None:
    t0 = time.time()
    P(f"296 cd4_296.py [{LIB}] -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    bank = B280.Bank(E280)
    obs = bank.obs
    ax = np.load(f"{E280}/axis280.npz", allow_pickle=True)
    sym = ax["official"][ax["cols_off"]].astype(str)          # symbol of each mapped column
    col_of = {s: int(c) for s, c in zip(sym, ax["cols_cd4"])}  # symbol -> CD4-axis column
    scored_cols = np.asarray(bank.cols_cd4[bank.scored], np.int64)
    Xb = np.asarray(bank.X, dtype=np.float32)                 # (n_rows, 18129) raw counts
    tot_row = Xb.sum(1, dtype=np.float64)
    P(f"  bank {Xb.shape}; {scored_cols.size} scored columns; {len(bank.panel)} panel targets; "
      f"peak {L.peak_gib():.2f} GiB")

    tgt = obs.pert.values.astype(str)
    gid = obs.guide.values.astype(str)
    cond = obs.cond.values.astype(str)
    dv = obs.D.values.astype(str)
    isT = (obs.gtype.values.astype(str) == "targeting")
    ncell = obs.n_cells.values.astype(np.float64)

    res: dict = {"library": LIB, "budgets": {}, "kd": {}, "conds": list(CONDS)}
    store: dict[str, np.ndarray] = {}

    # ---------------------------------------------------------------- full-depth knockdown per (cond, guide)
    kd: dict[tuple[str, str], dict] = {}
    for c in CONDS:
        ctrl = sum(bank.ctrl[(c, d)] for d in B280.DORDER)
        pc = L.prof(ctrl)
        cp = ctrl / ctrl.sum()
        for g in sorted(set(gid[isT & (cond == c)])):
            rows = np.flatnonzero(isT & (cond == c) & (gid == g))
            t = tgt[rows[0]]
            j = col_of.get(t, -1)
            S = Xb[rows].sum(0, dtype=np.float64)
            tt = S.sum()
            if j < 0 or tt <= 0:
                kd[(c, g)] = dict(target=t, lfc=np.nan, kfrac=np.nan, umis=float(tt),
                                  cells=float(ncell[rows].sum()))
                continue
            p = S[j] / tt
            kd[(c, g)] = dict(target=t, lfc=float(np.log1p(L.TS * p) - pc[j]),
                              kfrac=float(1.0 - p / cp[j]) if cp[j] > 0 else np.nan,
                              umis=float(tt), cells=float(ncell[rows].sum()))
    res["kd"] = {f"{c}|{g}": v for (c, g), v in kd.items()}
    fin = [v["lfc"] for v in kd.values() if np.isfinite(v["lfc"])]
    P(f"  full-depth on-target lfc measurable on {len(fin)}/{len(kd)} (condition, guide) rows; "
      f"median {np.median(fin):+.3f}")

    # cross-estimator check (KA-5b): my |lfc| against the authors' own t-statistic
    kdt = pd.read_csv(X.ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24/"
                              "guide_kd_efficiency.suppl_table.csv", index_col=0)
    mine, theirs = [], []
    for (c, g), v in kd.items():
        r = kdt[(kdt.index.astype(str) == g) & (kdt.culture_condition == c)]
        if len(r) == 1 and np.isfinite(v["lfc"]):
            mine.append(abs(v["lfc"])); theirs.append(-float(r.t_statistic.iloc[0]))
    rho_x, p_x = X.spearman(mine, theirs)
    res["ka5b"] = dict(n=len(mine), spearman=rho_x, p=p_x, rule="> +0.3")
    P(f"KA-5b my |on-target lfc| vs the authors' -t_statistic on {len(mine)} (cond, guide) rows: "
      f"Spearman {rho_x:+.3f} (p {p_x:.2g}) -> {'PASS' if rho_x > 0.3 else 'FAIL'}")

    # ---------------------------------------------------------------- the ladder
    for T in BUDGETS:
        key = f"{T:.0e}"
        res["budgets"][key] = {}
        for c in CONDS:
            keys = ("a", "a2", "b", "x")
            per: dict[str, dict[str, list]] = {k: {} for k in keys}
            perp: dict[str, dict[str, list]] = {k: {} for k in keys}
            cells_side: list[float] = []
            cmean_acc: dict[str, list] = {}
            cross_acc: list[tuple[list, np.ndarray, np.ndarray]] = []
            for d in B280.DORDER:
                ctrl = bank.ctrl[(c, d)]
                sub = np.flatnonzero(isT & (cond == c) & (dv == d))
                by_t: dict[str, list[int]] = {}
                for i in sub:
                    by_t.setdefault(tgt[i], []).append(int(i))
                pairs = []
                for t, rr in sorted(by_t.items()):
                    if len(rr) < 2:
                        continue
                    rr = sorted(rr, key=lambda i: gid[i])[:2]
                    if tot_row[rr[0]] >= 2 * T and tot_row[rr[1]] >= 2 * T:
                        pairs.append((t, rr[0], rr[1]))
                if len(pairs) < MIN_COHORT:
                    continue
                tl = [t for t, _, _ in pairs]
                cells_side.append(float(np.median([ncell[r0] * T / tot_row[r0] for _, r0, _ in pairs])))
                cosr = {k: np.full((NREP, len(pairs)), np.nan) for k in keys}
                pdsr = {k: np.full((NREP, len(pairs)), np.nan) for k in keys}
                cmean = np.full((NREP, 2 * len(pairs)), np.nan)
                cr, cu = [], []
                for rep in range(NREP):
                    r = np.random.default_rng(L.SEED + 7919 * rep + int(T) + len(c) + len(d))
                    half = {}
                    for slot, gi in ((0, 1), (1, 2)):
                        h2 = np.stack([L.thin(Xb[p[gi]].astype(np.float64), 2 * T / tot_row[p[gi]], r)
                                       for p in pairs])
                        h1 = r.binomial(np.rint(h2).astype(np.int64), 0.5).astype(np.float64)
                        half[slot] = (h1, h2 - h1)
                    V = {slot: np.stack([L.thin(Xb[p[gi]].astype(np.float64), T / tot_row[p[gi]], r)
                                         for p in pairs]) for slot, gi in ((0, 1), (1, 2))}
                    sides = {"a": half[0], "a2": half[1], "b": (V[0], V[1]),
                             "x": (V[0], np.roll(V[1], 1, axis=0))}
                    for k in keys:
                        v = L.views(sides[k][0], sides[k][1], ctrl, ctrl, scored_cols)
                        cosr[k][rep] = v["cos_resid"]
                        pdsr[k][rep] = v["pds"]
                    CM = np.vstack([V[0], V[1]])
                    DC = L.delta(CM, ctrl)[:, scored_cols]
                    mu = DC.mean(0)
                    nn = np.linalg.norm(DC, axis=1)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        cmean[rep] = (DC @ mu) / (nn * np.linalg.norm(mu))
                    RC = DC - mu[None, :]
                    for Mx, acc in ((RC, cr), (DC, cu)):
                        q = np.linalg.norm(Mx, axis=1)
                        with np.errstate(invalid="ignore", divide="ignore"):
                            acc.append((Mx @ Mx.T) / np.outer(q, q))
                    del CM, DC, RC
                for k in keys:
                    m = np.nanmean(cosr[k], 0); mp = np.nanmean(pdsr[k], 0)
                    for i, t in enumerate(tl):
                        per[k].setdefault(t, []).append(float(m[i]))
                        perp[k].setdefault(t, []).append(float(mp[i]))
                cm = np.nanmean(cmean, 0)
                labs = [f"{t}|s0" for t in tl] + [f"{t}|s1" for t in tl]
                for i, lb in enumerate(labs):
                    cmean_acc.setdefault(lb, []).append(float(cm[i]))
                cross_acc.append((labs, np.nanmean(np.stack(cr), 0).astype(np.float32),
                                  np.nanmean(np.stack(cu), 0).astype(np.float32)))
                del cr, cu
            if not per["b"]:
                continue
            targets = sorted(set(per["a"]) & set(per["b"]) & set(per["x"]))
            if len(targets) < MIN_COHORT:
                continue
            arr = {k: np.asarray([np.mean(per[k][t]) if t in per[k] else np.nan for t in targets])
                   for k in keys}
            arrp = {k: np.asarray([np.mean(perp[k][t]) if t in perp[k] else np.nan for t in targets])
                    for k in keys}
            kdw = np.asarray([min(abs(kd[(c, g)]["lfc"]) for g in sorted(
                {gid[i] for i in np.flatnonzero(isT & (cond == c) & (tgt == t))})[:2]) if t in per["b"]
                else np.nan for t in targets])
            gg = {t: sorted({gid[i] for i in np.flatnonzero(isT & (cond == c) & (tgt == t))})[:2]
                  for t in targets}
            row = dict(cond=c, n_targets=len(targets), targets=targets,
                       cells_per_side=float(np.mean(cells_side)) if cells_side else None,
                       guides={t: gg[t] for t in targets},
                       kd={t: [kd[(c, g)]["lfc"] for g in gg[t]] for t in targets},
                       kf={t: [kd[(c, g)]["kfrac"] for g in gg[t]] for t in targets},
                       **{k: dict(**L.stats(arr[k]), boot=L.boot_median(arr[k])) for k in keys},
                       rho_guide=X.rho_boot(arr["b"], arr["a"]),
                       b_minus_x=L.boot_median_diff(arr["b"], arr["x"]),
                       pds={k: float(np.nanmean(arrp[k])) for k in keys})
            res["budgets"][key][c] = row
            for k in keys:
                store[f"{key}|{c}|{k}"] = arr[k]
                store[f"{key}|{c}|pds_{k}"] = arrp[k]
            store[f"{key}|{c}|targets"] = np.asarray(targets)
            store[f"{key}|{c}|kd_weak"] = kdw
            store[f"{key}|{c}|kd_pair"] = np.asarray([[kd[(c, g)]["lfc"] for g in gg[t]] for t in targets])
            store[f"{key}|{c}|kf_pair"] = np.asarray([[kd[(c, g)]["kfrac"] for g in gg[t]] for t in targets])
            store[f"{key}|{c}|cos_to_mean_lab"] = np.asarray(sorted(cmean_acc))
            store[f"{key}|{c}|cos_to_mean"] = np.asarray([np.mean(cmean_acc[k]) for k in sorted(cmean_acc)])
            np.savez_compressed(X.WORK / f"cd4_296_cross_{key}_{c}.npz",
                                **{f"lab{i}": np.asarray(a[0]) for i, a in enumerate(cross_acc)},
                                **{f"r{i}": a[1] for i, a in enumerate(cross_acc)},
                                **{f"u{i}": a[2] for i, a in enumerate(cross_acc)})
            P(f"  [T={key}] {c:9s} {len(targets):3d} targets, ~{row['cells_per_side'] or 0:.0f} cells/side: "
              f"(a) {L.fmt(row['a']['boot'])}  (a2) {L.fmt(row['a2']['boot'])}  (b) {L.fmt(row['b']['boot'])}"
              f"  (x) {L.fmt(row['x']['boot'])}  rho_guide {X.fmt_rho(row['rho_guide'])}")
            P(f"           pds (a) {row['pds']['a']:.4f} / (b) {row['pds']['b']:.4f} / "
              f"(x) {row['pds']['x']:.4f}")

    res["seconds"] = round(time.time() - t0, 1)
    res["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(X.WORK / "cd4_296.npz", **store)
    (X.WORK / "cd4_296.json").write_text(json.dumps(res, indent=1, default=float))
    P(f"\ndone {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
