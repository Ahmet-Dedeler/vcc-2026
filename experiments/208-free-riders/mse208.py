"""208 — `mse` for the 2x2 rider grid, through 096's `claim` route, on the real challenge panel.

`mse` is the binding check.  a15 has only **5.6%** of margin above the knee
(`claim` 0.24581 against `s^2 B` = 0.23268) and 196 showed the whole 0.111 `mse` score
vanishes the moment `C_hat` leaves saturation.  Neither rider is exempt: the rounding rule
moves the delivered per-gene totals themselves, and the spread policy moves the emitted
cell-to-cell scatter, which IS what `C_hat` reads.

EIGHT configurations, because two of them are known-answer rows and six are the live grid:

  emitter `quasi`      a15's own law.  `floor/tier` must reproduce 196's `m = 1` row
                       (`C_hat/C_p` 1.1258, `claim` 0.24653, raw `mse` 0.88228) and
                       `floor/flat` must reproduce 200's flat-1.00 row (`claim` 0.24653,
                       raw `mse` 0.88228 bit-equal to a15's).
  emitter `rc_callset` the INCUMBENT.  Same 2x2, and these are the numbers that price the
                       riders marginally to what is on the board.

The profile half (the rounding energy, and PLOLO's cross term that converts it into raw
`mse`) depends on the rule only -- `emit_counts` computes `totals` before `spread` is read --
and is measured once per rule.

    MAXT=8 uv run python mse208.py
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
WORK.mkdir(exist_ok=True)
E188 = ROOT / "experiments/188-magnitude-calibration"
E196 = ROOT / "experiments/196-emitted-depth"
sys.path.insert(0, str(E196))
sys.path.insert(0, str(E188))
sys.path.insert(0, str(ROOT / "src"))
import recon196  # noqa: E402
import vcc2026.rc_emit as RC  # noqa: E402
from ham196 import hamilton  # noqa: E402
from mse196 import DISPERSION, KNEE, PHI, S15, TS, n, raw_mse, score, sD  # noqa: E402
from mse196 import K as CAP_K  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

MAXT = int(os.environ.get("MAXT", "8"))
NL_1X = 400 * 20_109
RULES = ("floor", "ham")
SPREADS = ("tier", "flat")
EMITTERS = tuple(os.environ.get("EMITTERS", "quasi,rc_callset").split(","))
EFFECTS = ROOT / "experiments/178-composed-riders/work/eff.a15.npy"
BUNDLE = ROOT / "data/challenge"
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90
KEEP_TOP, KEEP_SPREAD = 8000, 0.25
# 196's own m = 1 quasi row, and 200's flat-1.00 row.
KA = {"quasi|floor|tier": dict(ratio=1.125770, claim=0.246529, raw=0.882283),
      "quasi|floor|flat": dict(claim=0.246529, raw=0.882283)}


def plolo_cross():
    """The cross term: how much of the rounding error's own energy lands in squared error."""
    from budget188 import FOLDS, deliver, load, mass_consistent_effect
    z = np.load(E188 / "work/chalframes.npz")
    bA, gA = z["bA"], z["gA"]
    out = {}
    for f in FOLDS:
        F = load(f)
        P0 = np.where(gA[None, :], F["pred"], 0.0)
        inf = deliver(P0, bA, gA)
        T, c = F["truth"], F["scored"]
        u_inf = float(np.mean(((inf - T)[:, c] ** 2).sum(1) / c.size))
        for rule in RULES:
            D = np.empty_like(P0)
            for i in range(P0.shape[0]):
                g = mass_consistent_effect(np.where(gA, P0[i], 0.0), bA, gA, "minnorm")
                p = np.clip(np.expm1(bA + g), 0.0, None) / TS
                p = p / float(p.sum())
                x = NL_1X * p
                Tq = np.floor(x) if rule == "floor" else hamilton(x).astype(np.float64)
                D[i] = np.log1p(TS * (Tq / Tq.sum())) - bA
            eps = D - inf
            u = float(np.mean(((D - T)[:, c] ** 2).sum(1) / c.size))
            e2 = float(np.mean((eps[:, c] ** 2).sum(1) / c.size))
            out[(f, rule)] = dict(dU=u - u_inf, eps2=e2,
                                  ratio=(u - u_inf) / e2 if e2 else np.nan)
        print("  [%s] " % f + "  ".join(
            f"{r}: dU {out[(f,r)]['dU']:.3e} eps2 {out[(f,r)]['eps2']:.3e} "
            f"ratio {out[(f,r)]['ratio']:.3f}" for r in RULES), flush=True)
    return out


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction
    t0 = time.time()
    base_raw = raw_mse(S15, 0.24581)
    print(f"knee s^2 B = {KNEE:.5f}   cap K = {CAP_K:.5f}   PRED_TRACE_CAP_K "
          f"{PRED_TRACE_CAP_K}   a15 raw mse {base_raw:.5f}  score {score(base_raw):.5f}"
          f"   (185's known-answer row: 0.88228 / 0.11112)")

    print("\n--- PROFILE HALF, cross term on PLOLO ---")
    cross = plolo_cross()

    print("\n--- PROFILE HALF, rounding energy on the CHALLENGE (exact, per row, /G_p) ---")
    prof = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        G = D["genes"].size
        sel = np.linspace(0, D["sums"].shape[0] - 1, 24).astype(int)
        for rule in RULES:
            acc = []
            for i in sel:
                R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
                p = R.x / R.x.sum()
                Tq = (R.totals(1.0) if rule == "floor" else hamilton(R.x)).astype(np.float64)
                eps = np.log1p(TS * Tq / Tq.sum()) - np.log1p(TS * p)
                si = np.flatnonzero(D["genes"] == D["targets"][i])
                if si.size:
                    eps[si[0]] = 0.0
                acc.append(float((eps ** 2).sum() / (G - 1)))
            prof[(ctx, rule)] = float(np.mean(acc))
        print(f"  ctx {ctx}: " + "  ".join(f"{r} {prof[(ctx,r)]:.3e}" for r in RULES),
              flush=True)

    # ---- dispersion half ---------------------------------------------------------------
    effects = np.load(EFFECTS)
    perts300 = np.asarray([ln.split(",")[0] for ln in
                           (BUNDLE / "pert_counts.csv").read_text().splitlines()[1:]])
    pos300 = {str(p): i for i, p in enumerate(perts300)}
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)

    print("\n--- DISPERSION HALF, C_hat on the CHALLENGE, paired per block ---")
    chat, cp, nnz, maxcell = {}, {}, {}, {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        G = D["genes"].size
        assert effects.shape[1] == G
        Ad = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        M = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        ids = np.asarray(Ad.obs["ntc_id"].astype(str).values)
        cps = [float(jackknife_correction(sp.csr_matrix(M[np.flatnonzero(ids == nm)]),
                                          np.zeros(int((ids == nm).sum()), np.intp), 1, TS,
                                          chunk=64)[0])
               for nm in dict.fromkeys(ids.tolist())]
        cp[ctx] = float(np.mean(cps))
        control = (RC.control_block_from_counts(M, min_cpm=5.0)
                   if "rc_callset" in EMITTERS else None)
        del Ad, M
        print(f"  ctx {ctx}: C_p {cp[ctx]:.4f}"
              + (f", scorer gate {int(control.keep.sum())}/{G}" if control else "")
              + f"  [{time.time()-t0:.0f}s]", flush=True)

        sel = np.linspace(0, D["sums"].shape[0] - 1, MAXT).astype(int)
        for i in sel:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            Lp = R.emit_L(1.0)
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
            sv = {"flat": np.full(G, 1.0), "tier": np.full(G, 1.0)}
            sv["tier"][order[:min(KEEP_TOP, int(D["gate"].sum()))]] = KEEP_SPREAD
            key = np.where(D["gate"], effects[pos300[str(D["targets"][i])]], 0.0)
            for rule in RULES:
                Tw = R.totals(1.0) if rule == "floor" else hamilton(R.x)
                prof_p = (Tw + 0.5) / (n * Lp)
                for spd in SPREADS:
                    for em in EMITTERS:
                        rng = np.random.default_rng(90_000 + int(i))
                        if em == "quasi":
                            v, c, ip = emit_counts(prof_p, np.full(n, Lp), DISPERSION, rng,
                                                   mode="quasi", spread=sv[spd],
                                                   overdispersion=PHI)
                            blk = sp.csr_matrix((v.astype(np.float64), c.astype(np.int64),
                                                 np.asarray(ip)), shape=(n, G))
                        else:
                            B, _dg = RC.rc_hybrid_block(
                                prof_p, np.full(n, Lp), DISPERSION, rng, pass2_rng=rng,
                                spread=sv[spd], overdispersion=PHI, control=control,
                                key=key, jitter_seed=1000 + int(i),
                                emit_counts=emit_counts, call_set=call_set)
                            blk = sp.csr_matrix(B.astype(np.float64))
                            del B
                        got = np.asarray(blk.sum(0)).ravel().astype(np.int64)
                        assert np.array_equal(got, Tw), (em, rule, spd, ctx, i)
                        k = (ctx, em, rule, spd)
                        chat.setdefault(k, []).append(float(jackknife_correction(
                            blk, np.zeros(n, np.intp), 1, TS, chunk=64)[0]))
                        nnz.setdefault(k, []).append(int(blk.nnz))
                        maxcell.setdefault(k, []).append(
                            int(np.asarray(blk.sum(1)).ravel().max()))
                        del blk
        for em in EMITTERS:
            print(f"    {em:<11} " + "  ".join(
                f"{r}/{s} Chat {np.mean(chat[(ctx,em,r,s)]):.3f} "
                f"(/C_p {np.mean(chat[(ctx,em,r,s)])/cp[ctx]:.3f})"
                for r in RULES for s in SPREADS) + f"  [{time.time()-t0:.0f}s]", flush=True)
        del control

    # ---- assemble ----------------------------------------------------------------------
    print("\n=== THE PRICE OF THE TWO RIDERS ON `mse` ===")
    print(f"{'emitter':<11}{'rule':>6}{'spread':>7}{'Chat/C_p':>10}{'rows<cap':>10}"
          f"{'claim':>10}{'branch':>11}{'d claim':>10}{'round E':>11}{'d prof':>10}"
          f"{'raw mse':>10}{'score':>9}{'d overall':>11}{'nnz/blk':>12}{'max cell':>10}")
    rows = {}
    e_floor = float(np.mean([prof[(c, "floor")] for c in "ABC"]))
    for em in EMITTERS:
        for rule in RULES:
            for spd in SPREADS:
                ratios, claims, below = [], [], 0
                for ctx in "ABC":
                    ch = np.asarray(chat[(ctx, em, rule, spd)])
                    cap = PRED_TRACE_CAP_K * cp[ctx]
                    ratios.append(float(ch.mean() / cp[ctx]))
                    claims.append(float(np.mean(np.minimum(ch, cap)) / cp[ctx]))
                    below += int((ch < cap).sum())
                ratio = float(np.mean(ratios))
                claim = float(np.mean(claims)) * CAP_K
                branch = "SATURATED" if claim >= KNEE else "linear"
                d_claim = (min(KNEE, claim) - min(KNEE, 0.24581)) / sD * -1.0
                energy = float(np.mean([prof[(c, rule)] for c in "ABC"]))
                xr = float(np.mean([v["ratio"] for k, v in cross.items() if k[1] == rule]))
                d_prof = (energy - e_floor) * xr / sD
                raw = base_raw + d_claim + d_prof
                mn = float(np.mean([np.mean(nnz[(c, em, rule, spd)]) for c in "ABC"]))
                mx = float(np.max([np.max(maxcell[(c, em, rule, spd)]) for c in "ABC"]))
                r = dict(emitter=em, rule=rule, spread=spd, ratio=ratio, claim=claim,
                         branch=branch, d_claim_mse=d_claim, energy=energy,
                         cross_ratio=xr, d_prof_mse=d_prof, raw_mse=raw, score=score(raw),
                         d_overall=(score(raw) - score(base_raw)) / 6,
                         rows_below_cap=below, nnz=mn, max_cell=mx)
                rows[f"{em}|{rule}|{spd}"] = r
                print(f"{em:<11}{rule:>6}{spd:>7}{ratio:>10.3f}{below:>10d}{claim:>10.5f}"
                      f"{branch:>11}{d_claim:>+10.5f}{energy:>11.3e}{d_prof:>+10.5f}"
                      f"{raw:>10.5f}{score(raw):>9.5f}"
                      f"{(score(raw)-score(base_raw))/6:>+11.5f}{mn:>12,.0f}{mx:>10,.0f}")

    print("\nKNOWN-ANSWER ROWS (rule 10)")
    bad = []
    for k, want in KA.items():
        if k not in rows:
            continue
        g = rows[k]
        for fld, tol in (("ratio", 2e-2), ("claim", 5e-3), ("raw", 5e-3)):
            if fld not in want:
                continue
            got = g["raw_mse"] if fld == "raw" else g[fld]
            ok = abs(got - want[fld]) <= tol
            print(f"  {k:<22} {fld:<7} {got:>10.6f}  expect {want[fld]:<10.6f} "
                  f"{'PASS' if ok else 'FAIL'}")
            if not ok:
                bad.append(f"{k}.{fld}")
    if bad and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {bad}")
    print("  known-answer rows PASS")

    json.dump(dict(rows=rows, C_p=cp, MAXT=MAXT, knee=KNEE, cap_k=CAP_K,
                   base_raw=base_raw, base_score=score(base_raw),
                   chat={"|".join(k): v for k, v in chat.items()},
                   nnz={"|".join(k): v for k, v in nnz.items()},
                   max_cell={"|".join(k): v for k, v in maxcell.items()},
                   prof={f"{k[0]}|{k[1]}": v for k, v in prof.items()},
                   cross={f"{k[0]}|{k[1]}": v for k, v in cross.items()}),
              open(WORK / "mse208.json", "w"), indent=1, default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
