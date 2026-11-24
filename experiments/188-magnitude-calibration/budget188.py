"""188 — the decisive held-out screen: which magnitude SHAPE minimises L1 at a fixed L2 norm?

Every candidate below preserves each row's L2 norm in the builder's space exactly, so they
differ ONLY in where the magnitude sits.  They are then read on the held line that never
contributed to the map, on two endpoints that are the two things the six members split into:

  `nmae*`  per row, `sum|x - r| / sum|r|` over the row's top-2000 |r| genes -- the exact
           functional `de_wilcoxon_lfc_nmae` computes (169 §5 identity (I)), evaluated in
           log2 fold-change space, on the reference the fold actually holds.  It is a FOURTH
           `nmae` surface and it is NOT board-validated; 162 classifies `nmae` as not locally
           priceable (3/6 signs) and both existing harnesses missed a13 by 0.075 raw IN THE
           WRONG DIRECTION.  It is reported as a third opinion, never as the answer.
  `pds`    the exact rank statistic, computed on the DELIVERED effect (156: an effect-level
           score measures an object the board never sees).

CANDIDATES
----------
  identity            a15 + a23 as shipped
  power q             169's arm, `sign(h)|h|^q`, row norm restored
  median-map          the unconstrained L1 optimum (M), row norm restored
  budget-map          the L2-budget-constrained L1 optimum, which hits the norm itself
  abs-map             the mandate's literal "conditional median of the true MAGNITUDE"
  de-map              (M) fitted only on each row's top-2000 |r| genes

KNOWN-ANSWER ROWS
-----------------
  * `identity` must reproduce the cached PLOLO a23 rows bit-for-bit (`plolo169.json`);
  * `power q=1.00` must be bit-identical to `identity`;
  * every candidate's delivered row norm ratio `rA` must be 1.000 by construction.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "experiments/169-sharpening-six-member"))
sys.path.insert(0, os.path.join(REPO, "experiments/162-retrodicting-validator"))
from mapfit import (Map, apply_map, apply_map_lfc, fit_map, fit_map_budget,  # noqa: E402
                    to_lfc)
from sharp import sharpen, sharpen_lfc                                        # noqa: E402

FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
N_DE = 2000
TS = 50_000.0
PDS_SPAN = 0.95228 - 0.50020
POWERS = (1.05, 1.10, 1.15, 1.25, 1.40)


# ---------------------------------------------------------------- renderer (162 verbatim)
def mass_consistent_effect(gated, b, on, mode="minnorm", iters=64, rtol=1e-12):
    h = np.asarray(gated, dtype=np.float64).copy()
    for _ in range(iters):
        live = on & (b + h > 0.0)
        residual = float(np.clip(np.expm1(b + h), 0.0, None).sum() - TS)
        if abs(residual) <= rtol * TS:
            break
        w = np.where(live, np.exp(b + h), 0.0)
        denom = float(w @ w) if mode == "minnorm" else float(w.sum())
        if not np.isfinite(denom) or denom <= 0.0:
            break
        h = h + (-residual / denom) * w
    return h


def deliver(effects, b, gate):
    out = np.empty_like(effects)
    for i in range(effects.shape[0]):
        g = mass_consistent_effect(np.where(gate, effects[i], 0.0), b, gate, "minnorm")
        p = np.clip(np.expm1(b + g), 0.0, None) / TS
        p = p / float(p.sum())
        out[i] = np.log1p(TS * p) - b
    return out


def pds_of(P, T):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


def nmae_star(Lp, Lr, de):
    """`mean_p( sum_{g in DE_p} |lfc_hat - lfc_ref| / sum |lfc_ref| )` — identity (I)'s object."""
    out = []
    for i in range(Lp.shape[0]):
        s = np.flatnonzero(de[i])
        den = np.abs(Lr[i, s]).sum()
        if den <= 0 or s.size < 10:
            continue
        out.append(float(np.abs(Lp[i, s] - Lr[i, s]).sum() / den))
    return float(np.mean(out))


def load(fold):
    z = np.load(os.path.join(WORK, f"pair-{fold}.npz"), allow_pickle=True)
    return dict(pred=z["pred"].astype(np.float64), truth=z["truth"].astype(np.float64),
                b=z["b_ctrl"].astype(np.float64), scored=z["scored"])


def main() -> None:
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    cg = z["gA"] & z["gB"] & z["gC"]
    bA, gA = z["bA"], z["gA"]
    F = {f: load(f) for f in FOLDS}
    cols = {}
    for f in FOLDS:
        m = np.zeros(cg.size, bool)
        m[F[f]["scored"]] = True
        cols[f] = np.flatnonzero(m & cg)

    # ---- fitting pairs, both spaces --------------------------------------------------
    D = {}
    for space in ("h", "lfc"):
        for f in FOLDS:
            c = cols[f]
            P, R = F[f]["pred"][:, c], F[f]["truth"][:, c]
            if space == "lfc":
                one = np.ones(c.size, bool)
                P, R = to_lfc(P, F[f]["b"][c], one), to_lfc(R, F[f]["b"][c], one)
            live = (P != 0) & (R != 0)
            de = np.zeros(R.shape, bool)
            k = min(N_DE, R.shape[1])
            np.put_along_axis(de, np.argpartition(-np.abs(R), k - 1, axis=1)[:, :k], True, 1)
            D[(space, f)] = (np.abs(P)[live], (R * np.sign(P))[live], np.abs(R)[live],
                             de[live])

    maps, mus = {}, {}
    for space in ("h", "lfc"):
        for held in FOLDS:
            x = np.concatenate([D[(space, f)][0] for f in FOLDS if f != held])
            y = np.concatenate([D[(space, f)][1] for f in FOLDS if f != held])
            a = np.concatenate([D[(space, f)][2] for f in FOLDS if f != held])
            d = np.concatenate([D[(space, f)][3] for f in FOLDS if f != held])
            maps[(space, "median", held)] = fit_map(x, y)
            maps[(space, "abs", held)] = fit_map(x, a)
            maps[(space, "de", held)] = fit_map(x[d], y[d])
            maps[(space, "budget", held)], mus[(space, held)] = fit_map_budget(x, y)
            maps[(space, "debudget", held)], mus[(space, "de", held)] = \
                fit_map_budget(x[d], y[d])
    print("budget-map multipliers mu (the quantile level rises as 1/2 + mu*|p|):")
    for space in ("h", "lfc"):
        print(f"  {space:<4}" + "  ".join(f"{h}={mus[(space, h)]:.3f}" for h in FOLDS)
              + "   |  DE-restricted: "
              + "  ".join(f"{h}={mus[(space, 'de', h)]:.3f}" for h in FOLDS))

    # ---- sign accuracy and |r| information, per |p| bin ------------------------------
    print(f"\n{'=' * 104}\nWHAT OUR PREDICTED MAGNITUDE ACTUALLY CARRIES, per decile of |p| "
          f"(pooled over all four folds, h space)")
    x = np.concatenate([D[("h", f)][0] for f in FOLDS])
    y = np.concatenate([D[("h", f)][1] for f in FOLDS])
    a = np.concatenate([D[("h", f)][2] for f in FOLDS])
    d = np.concatenate([D[("h", f)][3] for f in FOLDS])
    q = np.quantile(x, np.linspace(0, 1, 11))
    print(f"  {'decile':>7}{'|p| mid':>10}{'P(sign ok)':>12}{'median |r|':>12}"
          f"{'median r*s':>12}{'|p|/|r| med':>13}{'P(sign ok|DE)':>15}{'median|r| DE':>14}")
    signrows = []
    for i in range(10):
        s = (x >= q[i]) & (x < q[i + 1] if i < 9 else x <= q[10])
        sd = s & d
        row = (float(np.median(x[s])), float((y[s] > 0).mean()), float(np.median(a[s])),
               float(np.median(y[s])), float(np.median(x[s] / np.maximum(a[s], 1e-12))),
               float((y[sd] > 0).mean()), float(np.median(a[sd])))
        signrows.append(row)
        print(f"  {i+1:>7}{row[0]:>10.5f}{row[1]:>12.4f}{row[2]:>12.5f}{row[3]:>12.5f}"
              f"{row[4]:>13.4f}{row[5]:>15.4f}{row[6]:>14.5f}")

    # ---- held-out evaluation ---------------------------------------------------------
    print(f"\n{'=' * 104}\nHELD-OUT SCREEN — map fitted on three lines, read on the fourth. "
          f"Challenge frame A, delivered.\nEvery arm preserves the builder-space row L2 norm "
          f"exactly, so `mse` is neutral by construction.")
    ref = {}
    p169 = json.loads(open(os.path.join(REPO,
        "experiments/169-sharpening-six-member/work/plolo169.json")).read())
    for r in p169:
        ref[(r["fold"], r["arm"], round(r["q"], 2), r["space"])] = r["pds"]

    ARMS = [("identity", None, None)]
    ARMS += [(f"power{p:.2f}", "h", ("power", p)) for p in POWERS]
    ARMS += [(f"powerL{p:.2f}", "lfc", ("power", p)) for p in POWERS]
    for space in ("h", "lfc"):
        for kind in ("median", "abs", "de", "budget", "debudget"):
            ARMS.append((f"{kind}{'L' if space == 'lfc' else ''}", space, ("map", kind)))

    res = {}
    for held in FOLDS:
        P0 = F[held]["pred"]
        c = cols[held]
        Rh = F[held]["truth"]
        one = np.ones(cg.size, bool)
        # truth in log2 fold change, in the frame the truth was measured in (its own line)
        LR = to_lfc(Rh, F[held]["b"], gA)
        de = np.zeros((Rh.shape[0], c.size), bool)
        k = min(N_DE, c.size)
        np.put_along_axis(de, np.argpartition(-np.abs(LR[:, c]), k - 1, axis=1)[:, :k],
                          True, 1)
        Tsc = Rh[:, F[held]["scored"]]
        for name, space, spec in ARMS:
            if spec is None:
                Q = np.where(gA[None, :], P0, 0.0)
            elif spec[0] == "power":
                Q = (sharpen(P0, spec[1], gA, mode="row") if space == "h"
                     else sharpen_lfc(P0, spec[1], bA, gA, mode="row"))
            else:
                mm = maps[(space, spec[1], held)]
                Q = (apply_map(P0, mm, gA, mode="row") if space == "h"
                     else apply_map_lfc(P0, mm, bA, gA, restore="h"))
            rA = float(np.median(np.linalg.norm(Q, axis=1)) /
                       np.median(np.linalg.norm(np.where(gA[None, :], P0, 0.0), axis=1)))
            Dl = deliver(Q, bA, gA)
            LP = to_lfc(Dl, bA, gA)
            res[(held, name)] = dict(
                pds=pds_of(Dl[:, F[held]["scored"]], Tsc),
                nmae=nmae_star(LP[:, c], LR[:, c], de),
                l1=float(np.abs(LP[:, c] - LR[:, c]).sum() / np.abs(LR[:, c]).sum()),
                rA=rA,
                over=float(np.abs(LP[:, c])[de][np.abs(LP[:, c])[de] >
                                              np.abs(LR[:, c])[de]].sum() /
                           max(np.abs(LP[:, c])[de].sum(), 1e-12)),
            )
        print(f"  [{held}] done", flush=True)

    # ---- known-answer gate ------------------------------------------------------------
    print(f"\n--- KNOWN-ANSWER GATE: `identity` and `power1.00` against 169's cached PLOLO "
          f"chalA rows ---")
    worst = 0.0
    for held in FOLDS:
        want = ref.get((held, "a23", 1.00, "chalA/delivered"))
        got = res[(held, "identity")]["pds"]
        worst = max(worst, abs(got - want))
        print(f"  {held:<14} 169 {want:.6f}   here {got:.6f}   |d| {abs(got-want):.2e}")
    print(f"  worst |d pds| = {worst:.2e}")
    if worst > 1e-6:
        sys.exit("ABORT: identity does not reproduce 169's PLOLO chalA rows")

    # ---- report ------------------------------------------------------------------------
    print(f"\n{'=' * 104}\nHELD-OUT RESULTS, 4-fold unweighted mean (162's pre-committed "
          f"predictor), challenge frame A")
    print(f"  {'arm':<14}{'d pds raw':>11}{'d pds scaled':>14}{'nmae*':>9}{'d nmae*':>10}"
          f"{'L1 all':>9}{'d L1':>9}{'rA':>8}{'DE L1 over k=1':>16}")
    base = {f: res[(f, "identity")] for f in FOLDS}
    rows = {}
    for name, _, _ in ARMS:
        dp = float(np.mean([res[(f, name)]["pds"] - base[f]["pds"] for f in FOLDS]))
        nm = float(np.mean([res[(f, name)]["nmae"] for f in FOLDS]))
        dn = float(np.mean([res[(f, name)]["nmae"] - base[f]["nmae"] for f in FOLDS]))
        l1 = float(np.mean([res[(f, name)]["l1"] for f in FOLDS]))
        dl = float(np.mean([res[(f, name)]["l1"] - base[f]["l1"] for f in FOLDS]))
        ra = float(np.mean([res[(f, name)]["rA"] for f in FOLDS]))
        ov = float(np.mean([res[(f, name)]["over"] for f in FOLDS]))
        rows[name] = dict(dpds=dp, dpds_scaled=dp / PDS_SPAN, nmae=nm, dnmae=dn,
                          l1=l1, dl1=dl, rA=ra, over=ov,
                          per_fold_pds={f: res[(f, name)]["pds"] - base[f]["pds"]
                                        for f in FOLDS},
                          per_fold_nmae={f: res[(f, name)]["nmae"] - base[f]["nmae"]
                                         for f in FOLDS})
        print(f"  {name:<14}{dp:>+11.6f}{dp/PDS_SPAN:>+14.6f}{nm:>9.5f}{dn:>+10.5f}"
              f"{l1:>9.5f}{dl:>+9.5f}{ra:>8.5f}{ov:>16.4f}")

    print(f"\nper-fold `pds` and `nmae*` deltas for the shortlist")
    for name in ("power1.15", "powerL1.15", "median", "medianL", "budget", "budgetL",
                 "debudget", "debudgetL"):
        if name not in rows:
            continue
        print(f"  {name:<12} pds " + " ".join(
            f"{f[:6]}={rows[name]['per_fold_pds'][f]:+.5f}" for f in FOLDS)
            + "   nmae* " + " ".join(
            f"{f[:6]}={rows[name]['per_fold_nmae'][f]:+.5f}" for f in FOLDS))

    with open(os.path.join(WORK, "budget188.json"), "w") as fh:
        json.dump(dict(rows=rows, signrows=signrows,
                       mu={f"{k}": v for k, v in mus.items()},
                       maps={f"{s}|{k}|{h}": maps[(s, k, h)].to_dict()
                             for (s, k, h) in maps}), fh, indent=1, default=float)
    print(f"\nwrote {WORK}/budget188.json")


if __name__ == "__main__":
    sys.exit(main())
