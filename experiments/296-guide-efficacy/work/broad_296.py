"""296 -- Broad-16 (Cas9 KNOCKOUT, 12 working lines): the knockdown-stratified construct term.

268's `halves/<line>.npz` carries per-target count SUMS only, so cells cannot be re-sampled; every side is
thinned to a common TOTAL-UMI budget T (285's currency for this library).  Contrasts, per (line, target):

  (a)   random cell halves A vs B        -- 285's matched reference; it MIXES the two constructs
  (a1)  binomial count split of v1       -- v1's own (Poisson-only) reliability
  (a2)  binomial count split of v2       -- v2's own
  (b)   v1 vs v2                         -- the construct term
  (x)   v1 of target t vs v2 of the NEXT target (a different gene) -- the cross-target null

The on-target knockdown of each construct is its own gene's column at FULL depth against the line's
control profile (the 10 olfactory-receptor knockouts summed, 268's `CONTROLS`).  For Cas9 a "failure" is a
construct whose target mRNA does not drop -- editing without nonsense-mediated decay reads as no knockdown,
which is a real limitation of this readout and is stated in RESULTS.

    PYTHONDONTWRITEBYTECODE=1 uv run python broad_296.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402

LIB = "Broad-16"
E268 = X.ROOT / "experiments/268-broad16-context-residuals/work"
HV = E268 / "halves"
P = L.Logger(X.WORK / "broad_296.log")
BUDGETS = (2e5, 5e5)
NREP = 8


def main() -> None:
    t0 = time.time()
    P(f"296 broad_296.py [{LIB}] -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    lines = sorted(p.stem for p in HV.glob("*.npz"))
    inv = json.load(open(E268 / "inventory.json"))
    CONTROLS = inv["controls"]
    data = {}
    for ln in lines:
        z = np.load(HV / f"{ln}.npz", allow_pickle=True)
        d = {k: np.asarray(z[k]) for k in ("targets", "genes", "ncell", "full", "A", "B", "v1", "v2")}
        d["targets"] = d["targets"].astype(str); d["genes"] = d["genes"].astype(str)
        data[ln] = d
    genes = data[lines[0]]["genes"]
    assert all(list(data[ln]["genes"]) == list(genes) for ln in lines)
    tg = sorted(set.intersection(*[set(data[ln]["targets"].tolist()) for ln in lines]))
    tg = [t for t in tg if t not in CONTROLS]
    for ln in lines:
        pos = {t: i for i, t in enumerate(data[ln]["targets"])}
        ci = [pos[c] for c in CONTROLS]
        data[ln]["ctrl"] = np.nan_to_num(data[ln]["full"][ci], nan=0.0).sum(0)
        sel = np.asarray([pos[t] for t in tg], np.int64)
        for k in ("targets", "ncell", "full", "A", "B", "v1", "v2"):
            data[ln][k] = data[ln][k][sel]
    gset = {g: i for i, g in enumerate(genes)}
    scored = np.asarray(sorted(set(range(len(genes))) - {gset[t] for t in inv["targets"] if t in gset}),
                        np.int64)
    P(f"  {len(lines)} lines, {len(tg)} shared non-control targets, {len(genes)} genes "
      f"({scored.size} scored), {len(CONTROLS)} OR controls")

    # ------------------------------------------------------------------ full-depth on-target knockdown
    kd: dict[tuple[str, str], dict] = {}
    for ln in lines:
        d = data[ln]
        pc = L.prof(d["ctrl"])
        cp = d["ctrl"] / d["ctrl"].sum()
        for i, t in enumerate(tg):
            j = gset.get(t, -1)
            for k in ("v1", "v2", "full"):
                S = np.nan_to_num(d[k][i], nan=0.0)
                tot = S.sum()
                if j < 0 or tot <= 0:
                    kd[(ln, t, k)] = dict(lfc=np.nan, kfrac=np.nan, umis=float(tot))
                    continue
                p = S[j] / tot
                kd[(ln, t, k)] = dict(lfc=float(np.log1p(L.TS * p) - pc[j]),
                                      kfrac=float(1.0 - p / cp[j]) if cp[j] > 0 else np.nan,
                                      umis=float(tot))
    nfin = sum(np.isfinite(v["lfc"]) for k, v in kd.items() if k[2] == "full")
    P(f"  on-target measurable on {nfin}/{len(lines) * len(tg)} (line, target) full rows")

    res: dict = {"library": LIB, "lines": lines, "n_targets": len(tg), "budgets": {},
                 "kd": {f"{ln}|{t}|{k}": v for (ln, t, k), v in kd.items()}}
    store: dict[str, np.ndarray] = {}

    for T in BUDGETS:
        key = f"{T:.0e}"
        rows: list[dict] = []
        cross_r_all, cross_u_all, cmean_all, clab_all = [], [], [], []
        for ln in lines:
            d = data[ln]
            S = {k: np.nan_to_num(d[k], nan=0.0) for k in ("A", "B", "v1", "v2")}
            tot = {k: S[k].sum(1) for k in S}
            keep = np.flatnonzero((tot["A"] >= T) & (tot["B"] >= T)
                                  & (tot["v1"] >= 2 * T) & (tot["v2"] >= 2 * T))
            if keep.size < 15:
                P(f"  [T={key}] {ln}: only {keep.size} targets clear the budget, skipped")
                continue
            K = keep.size
            keys = ("a", "a1", "a2", "b", "x")
            cosr = {k: np.full((NREP, K), np.nan) for k in keys}
            pdsr = {k: np.full((NREP, K), np.nan) for k in keys}
            cmean = np.full((NREP, 2 * K), np.nan)
            cr, cu = [], []
            pc = L.prof(d["ctrl"])
            for rep in range(NREP):
                r = np.random.default_rng(L.SEED + 101 * rep + int(T))
                M = {k: np.stack([L.thin(S[k][i], T / tot[k][i], r) for i in keep]) for k in ("A", "B")}
                # v1 / v2 at 2T, split into two independent binomial halves of T each
                H = {}
                for k in ("v1", "v2"):
                    h2 = np.stack([L.thin(S[k][i], 2 * T / tot[k][i], r) for i in keep])
                    h1 = r.binomial(np.rint(h2).astype(np.int64), 0.5).astype(np.float64)
                    H[k] = (h1, h2 - h1)
                # (b) and (x) must sit at the SAME per-side budget T as (a), (a1), (a2):
                # v1 and v2 are each thinned to T independently (NOT the 2T rows summed back).
                V = {k: np.stack([L.thin(S[k][i], T / tot[k][i], r) for i in keep]) for k in ("v1", "v2")}
                sides = {"a": (M["A"], M["B"]),
                         "a1": (H["v1"][0], H["v1"][1]),
                         "a2": (H["v2"][0], H["v2"][1]),
                         "b": (V["v1"], V["v2"])}
                sides["x"] = (V["v1"], np.roll(V["v2"], 1, axis=0))
                for k in keys:
                    v = L.views(sides[k][0], sides[k][1], d["ctrl"], d["ctrl"], scored)
                    cosr[k][rep] = v["cos_resid"]
                    pdsr[k][rep] = v["pds"]
                CM = np.vstack([V["v1"], V["v2"]])                   # v1 rows then v2 rows, each at T
                DC = L.delta(CM, d["ctrl"])[:, scored]
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
            m = {k: np.nanmean(cosr[k], 0) for k in keys}
            mp = {k: np.nanmean(pdsr[k], 0) for k in keys}
            for i, ix in enumerate(keep):
                t = tg[ix]
                rows.append(dict(line=ln, target=t, **{k: float(m[k][i]) for k in keys},
                                 **{f"pds_{k}": float(mp[k][i]) for k in keys},
                                 kd_v1=kd[(ln, t, "v1")]["lfc"], kd_v2=kd[(ln, t, "v2")]["lfc"],
                                 kf_v1=kd[(ln, t, "v1")]["kfrac"], kf_v2=kd[(ln, t, "v2")]["kfrac"],
                                 kd_full=kd[(ln, t, "full")]["lfc"],
                                 ncell=int(d["ncell"][ix, 0]), umi_v1=kd[(ln, t, "v1")]["umis"],
                                 umi_v2=kd[(ln, t, "v2")]["umis"]))
            cross_r_all.append(np.nanmean(np.stack(cr), 0).astype(np.float32))
            cross_u_all.append(np.nanmean(np.stack(cu), 0).astype(np.float32))
            cmean_all.append(np.nanmean(cmean, 0))
            clab_all.append(np.asarray([f"{ln}|{tg[ix]}|v1" for ix in keep]
                                       + [f"{ln}|{tg[ix]}|v2" for ix in keep]))
            P(f"  [T={key}] {ln:9s} n_t {K:3d}  (a) {np.nanmedian(m['a']):+.4f}  (b) {np.nanmedian(m['b']):+.4f}"
              f"  (x) {np.nanmedian(m['x']):+.4f}  rho_guide {np.nanmedian(m['b']) / np.nanmedian(m['a']):+.3f}")
            del cr, cu
        if not rows:
            continue
        arr = {k: np.asarray([r[k] for r in rows], float) for k in
               ("a", "a1", "a2", "b", "x", "pds_a", "pds_b", "pds_x", "kd_v1", "kd_v2", "kf_v1", "kf_v2")}
        res["budgets"][key] = dict(
            n_rows=len(rows), n_lines=len(set(r["line"] for r in rows)),
            **{k: dict(**L.stats(arr[k]), boot=L.boot_median(arr[k])) for k in ("a", "a1", "a2", "b", "x")},
            rho_guide=X.rho_boot(arr["b"], arr["a"]),
            pds={k: float(np.nanmean(arr[f"pds_{k}"])) for k in ("a", "b", "x")},
            rows=rows)
        for k in arr:
            store[f"{key}|{k}"] = arr[k]
        store[f"{key}|line"] = np.asarray([r["line"] for r in rows])
        store[f"{key}|target"] = np.asarray([r["target"] for r in rows])
        store[f"{key}|cross_labels"] = np.concatenate(clab_all)
        store[f"{key}|cos_to_mean"] = np.concatenate(cmean_all)
        # block-diagonal stack: cross-construct cosines are only defined WITHIN a line
        np.savez_compressed(X.WORK / f"broad_296_cross_{key}.npz",
                            **{f"{i}": c for i, c in enumerate(cross_r_all)},
                            **{f"u{i}": c for i, c in enumerate(cross_u_all)},
                            lines=np.asarray([r for r in sorted(set(x["line"] for x in rows))]),
                            labels=np.concatenate(clab_all))
        v = res["budgets"][key]
        P(f"  [T={key}] POOLED {len(rows)} (line,target) rows over {v['n_lines']} lines: "
          f"(a) {L.fmt(v['a']['boot'])}  (a1) {L.fmt(v['a1']['boot'])}  (a2) {L.fmt(v['a2']['boot'])}  "
          f"(b) {L.fmt(v['b']['boot'])}  (x) {L.fmt(v['x']['boot'])}  rho_guide {X.fmt_rho(v['rho_guide'])}")
        P(f"           pds (a) {v['pds']['a']:.4f} / (b) {v['pds']['b']:.4f} / (x) {v['pds']['x']:.4f}")

    res["seconds"] = round(time.time() - t0, 1)
    res["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(X.WORK / "broad_296.npz", **store)
    (X.WORK / "broad_296.json").write_text(json.dumps(res, indent=1, default=float))
    P(f"\ndone {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
