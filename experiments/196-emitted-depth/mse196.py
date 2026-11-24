"""196 — `mse` against emitted depth, both halves, through 096's route.

`expr_mse_unbiased_capped_norm`'s numerator is

    N_p = (1/G_p) ( ||b_hat_p - b_p||^2  -  rho * min(C_hat_p, C_p)  -  C_p )

and DEPTH moves both terms, in opposite directions:

  PROFILE half   integer rounding adds `eps = log1p(TS q) - log1p(TS p)` to the delivered
                 pseudobulk, so `||b_hat - b||^2` carries `||eps||^2` plus a cross term with
                 the (large) prediction error.  `||eps||^2` is computed EXACTLY on the
                 challenge from 196's reconstruction; the cross term needs a truth and is
                 therefore measured on PLOLO, where one exists, as the ratio between the
                 realised change in `||b_hat - truth||^2` and `||eps||^2`.
  DISPERSION half `C_hat_p` is the delete-one jackknife of OUR emitted noise.  Deeper cells
                 mean more counts per cell, so the Poisson component of that noise falls as
                 1/m -- but the gamma component does NOT (Var = mu + phi mu^2 with phi = 4.0,
                 and `mu` grows with `m`, so the second term is depth-INVARIANT after the
                 (TS/S)^2 normalisation).  Which one dominates decides the sign, and 191
                 measured a15 at `C_hat/C_p` = 1.094 with only 5.6% of margin above 185's
                 knee, so this has to be measured rather than argued.

185's route: `claim = mean_p min(C_hat_p, K C_p)`, `raw mse = 1 + (s^2 A - 2 s X + Kc -
min(s^2 B, claim)) / sumD` with 096's board-solved constants.  a15 sits at claim 0.24581
against a knee of `s^2 B = 0.23268`, so `mse` is FLAT until `claim` drops 5.6%.

    MAXT=6 DEPTHS=1,2,4,8,16 uv run python mse196.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
E188 = ROOT / "experiments/188-magnitude-calibration"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(E188))
import recon196                                                # noqa: E402
import vcc2026.rc_emit as RC                                   # noqa: E402
from vcc2026.donor_transfer import emit_counts                 # noqa: E402

TS = 50_000.0
n, PHI, DISPERSION = 400, 4.0, 0.1
DEPTHS = tuple(float(x) for x in os.environ.get("DEPTHS", "1,2,4,8,16").split(","))
MAXT = int(os.environ.get("MAXT", "6"))
ARMS = tuple(os.environ.get("ARMS", "quasi").split(","))
CB = {}
NL_1X = 400 * 20_109

bs = np.load(ROOT / "experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X_, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
S15 = 0.982 / 1.0796
KNEE = S15 * S15 * B


def raw_mse(s, claim):
    return 1 + (s * s * A - 2 * s * X_ + Kc - min(s * s * B, claim)) / sD


def score(raw):
    return min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)


def to_block(v, col, ip, G):
    Bm = np.zeros((n, G), np.int64)
    for r in range(n):
        s = slice(ip[r], ip[r + 1])
        Bm[r, col[s]] = v[s]
    return Bm


def plolo_cross():
    """The cross term, measured where a truth exists.  Returns U(m)-U(inf) and ||eps||^2."""
    from budget188 import FOLDS, deliver, load, mass_consistent_effect
    z = np.load(E188 / "work/chalframes.npz")
    bA, gA = z["bA"], z["gA"]
    out = {}
    for f in FOLDS:
        F = load(f)
        P0 = np.where(gA[None, :], F["pred"], 0.0)
        inf = deliver(P0, bA, gA)
        T = F["truth"]
        c = F["scored"]
        u_inf = float(np.mean(((inf - T)[:, c] ** 2).sum(1) / c.size))
        for m in DEPTHS:
            D = np.empty_like(P0)
            for i in range(P0.shape[0]):
                g = mass_consistent_effect(np.where(gA, P0[i], 0.0), bA, gA, "minnorm")
                p = np.clip(np.expm1(bA + g), 0.0, None) / TS
                p = p / float(p.sum())
                Tq = np.floor(NL_1X * m * p)
                D[i] = np.log1p(TS * (Tq / Tq.sum())) - bA
            eps = D - inf
            u = float(np.mean(((D - T)[:, c] ** 2).sum(1) / c.size))
            e2 = float(np.mean((eps[:, c] ** 2).sum(1) / c.size))
            out[(f, m)] = dict(dU=u - u_inf, eps2=e2, ratio=(u - u_inf) / e2 if e2 else np.nan)
        print(f"  [{f}] " + "  ".join(
            f"m{m:g}: dU {out[(f,m)]['dU']:.3e} eps2 {out[(f,m)]['eps2']:.3e} "
            f"ratio {out[(f,m)]['ratio']:.3f}" for m in DEPTHS), flush=True)
    return out


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction
    t0 = time.time()
    print(f"knee s^2 B = {KNEE:.5f}   cap K = {K:.5f}   a15 raw mse "
          f"{raw_mse(S15, 0.24581):.5f}  score {score(raw_mse(S15, 0.24581)):.5f}   "
          f"(185's known-answer row: 0.88228 / 0.11112)")

    print("\n--- PROFILE HALF, cross term on PLOLO (dU / ||eps||^2; 1.0 = pure eps energy) ---")
    cross = plolo_cross()

    print("\n--- PROFILE HALF, rounding energy on the CHALLENGE (exact, per row, /G_p) ---")
    prof = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        G = D["genes"].size
        sel = np.linspace(0, D["sums"].shape[0] - 1, 24).astype(int)
        for m in DEPTHS:
            acc = []
            for i in sel:
                R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
                p = R.x / R.x.sum()
                Tq = R.totals(m).astype(np.float64)
                eps = np.log1p(TS * Tq / Tq.sum()) - np.log1p(TS * p)
                self_i = np.flatnonzero(D["genes"] == D["targets"][i])
                if self_i.size:
                    eps[self_i[0]] = 0.0
                acc.append(float((eps ** 2).sum() / (G - 1)))
            prof[(ctx, m)] = float(np.mean(acc))
        print(f"  ctx {ctx}: " + "  ".join(f"m{m:g} {prof[(ctx,m)]:.3e}" for m in DEPTHS),
              flush=True)

    print("\n--- DISPERSION HALF, C_hat on the CHALLENGE at each depth ---")
    chat = {}
    cp = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        Ad = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        M = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        ids = np.asarray(Ad.obs["ntc_id"].astype(str).values)
        cps = [float(jackknife_correction(sp.csr_matrix(M[np.flatnonzero(ids == nm)]),
                                          np.zeros(int((ids == nm).sum()), np.intp), 1, TS,
                                          chunk=64)[0])
               for nm in dict.fromkeys(ids.tolist())]
        cp[ctx] = float(np.mean(cps))
        if "rc_hybrid" in ARMS:
            CB[ctx] = RC.control_block_from_counts(M, min_cpm=5.0)
        del Ad, M
        G = D["genes"].size
        sel = np.linspace(0, D["sums"].shape[0] - 1, MAXT).astype(int)
        for i in sel:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
            sv = np.full(G, 1.0)
            sv[order[:min(8000, int(D["gate"].sum()))]] = 0.25
            for m in DEPTHS:
                for arm in ARMS:
                    rng = np.random.default_rng(90_000 + int(i))
                    if arm == "quasi":
                        v, c, ip = emit_counts(R.profile(m), np.full(n, R.emit_L(m)),
                                               DISPERSION, rng, mode="quasi", spread=sv,
                                               overdispersion=PHI)
                        blk = sp.csr_matrix((v.astype(np.float64), c.astype(np.int64),
                                             np.asarray(ip)), shape=(n, G))
                    else:
                        prof = R.profile(m)
                        keyv = np.where(D["gate"], np.log1p(TS * prof) - D["b"], 0.0)
                        blkd, _ = RC.rc_hybrid_block(
                            prof, np.full(n, R.emit_L(m)), DISPERSION, rng,
                            pass2_rng=np.random.default_rng(90_000 + int(i)), spread=sv,
                            overdispersion=PHI, control=CB[ctx], key=keyv,
                            jitter_seed=1000 + int(i), emit_counts=emit_counts)
                        blk = sp.csr_matrix(blkd.astype(np.float64))
                        del blkd
                    chat.setdefault((arm, ctx, m), []).append(float(jackknife_correction(
                        blk, np.zeros(n, np.intp), 1, TS, chunk=64)[0]))
                    del blk
        for arm in ARMS:
            print(f"  ctx {ctx} [{arm}]: C_p {cp[ctx]:.4f}  " + "  ".join(
                f"m{m:g} Chat {np.mean(chat[(arm,ctx,m)]):.3f} "
                f"(/C_p {np.mean(chat[(arm,ctx,m)])/cp[ctx]:.3f})"
                for m in DEPTHS) + f"  [{time.time()-t0:.0f}s]", flush=True)
        CB.pop(ctx, None)

    print("\n=== THE SIX-MEMBER PRICE OF DEPTH ON `mse` ===")
    print(f"{'depth':>7}{'Chat/C_p':>10}{'rows<cap':>10}{'claim':>10}{'branch':>11}"
          f"{'d claim mse':>13}{'round energy':>14}{'d prof mse':>12}{'raw mse':>10}"
          f"{'score':>9}{'d overall':>11}")
    rows = []
    base_raw = raw_mse(S15, 0.24581)
    for arm in ARMS:
      print(f"--- emitter `{arm}` ---")
      for m in DEPTHS:
        ratios, claims, below = [], [], 0
        for ctx in "ABC":
            ch = np.asarray(chat[(arm, ctx, m)])
            cap = PRED_TRACE_CAP_K * cp[ctx]
            ratios.append(float(ch.mean() / cp[ctx]))
            claims.append(float(np.mean(np.minimum(ch, cap)) / cp[ctx]))
            below += int((ch < cap).sum())
        ratio = float(np.mean(ratios))
        claim = float(np.mean(claims)) * K            # in 096's normalised units
        branch = "SATURATED" if claim >= KNEE else "linear"
        d_claim_mse = (min(KNEE, claim) - min(KNEE, 0.24581)) / sD * -1.0
        energy = float(np.mean([prof[(ctx, m)] for ctx in "ABC"]))
        energy1 = float(np.mean([prof[(ctx, 1.0)] for ctx in "ABC"]))
        # the cross term, as a multiple of the pure rounding energy, from PLOLO
        xr = float(np.mean([v["ratio"] for k, v in cross.items() if k[1] == m]))
        d_prof_mse = (energy - energy1) * xr / sD
        raw = base_raw + d_claim_mse + d_prof_mse
        rows.append(dict(arm=arm, depth=m, ratio=ratio, claim=claim, branch=branch,
                         d_claim_mse=d_claim_mse, energy=energy, cross_ratio=xr,
                         d_prof_mse=d_prof_mse, raw_mse=raw, score=score(raw),
                         d_overall=(score(raw) - score(base_raw)) / 6, rows_below_cap=below))
        print(f"{m:>7g}{ratio:>10.3f}{below:>10d}{claim:>10.5f}{branch:>11}"
              f"{d_claim_mse:>+13.5f}{energy:>14.3e}{d_prof_mse:>+12.5f}{raw:>10.5f}"
              f"{score(raw):>9.5f}{(score(raw)-score(base_raw))/6:>+11.5f}")
    json.dump(dict(rows=rows, C_p=cp,
                   chat={f"{k[0]}|{k[1]}|{k[2]:g}": v for k, v in chat.items()},
                   prof={f"{k[0]}|{k[1]:g}": v for k, v in prof.items()},
                   cross={f"{k[0]}|{k[1]:g}": v for k, v in cross.items()}),
              open(WORK / "mse196.json", "w"), indent=1, default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
