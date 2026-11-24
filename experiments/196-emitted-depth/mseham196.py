"""196 — `mse` for the ROUNDING RULE at production's depth, both halves, through 096's route.

`mse196.py` swept the depth and found the cliff: at `m >= 2` the emitted per-cell noise falls
far enough that `C_hat` drops below `185`'s knee, `claim` leaves saturation, and raw `mse`
goes 0.88228 -> 1.00780, i.e. the member's whole 0.111 score.  That cliff is a property of
DEPTH: it is the counts per cell that set the Poisson component of `C_hat`.

The rounding RULE moves each gene's panel total by at most one count at the SAME depth, so
the emitted cell-to-cell spread — and therefore `C_hat`, and therefore `claim` — should not
move at all.  This file measures that instead of assuming it, which is Phase 3 falsifier 3.

Paired by construction: the same blocks, the same seeds, the same spread vector, the same
`L`.  The `floor` row is the known-answer row — it must reproduce `mse196`'s own `m = 1`
numbers (`C_hat/C_p` 1.1258, `claim` 0.24653, raw `mse` 0.88228, score 0.11112).

    MAXT=8 uv run python mseham196.py
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
from ham196 import hamilton                                    # noqa: E402
from mse196 import KNEE, S15, TS, DISPERSION, PHI, n, raw_mse, score, sD   # noqa: E402
from mse196 import K as CAP_K                                  # noqa: E402
from vcc2026.donor_transfer import emit_counts                 # noqa: E402

MAXT = int(os.environ.get("MAXT", "8"))
NL_1X = 400 * 20_109
RULES = ("floor", "ham")
KA = dict(ratio=1.125770, claim=0.246529, raw=0.882283)         # mse196's own m = 1 row


def rule_totals(R, rule):
    return R.totals(1.0) if rule == "floor" else hamilton(R.x)


def plolo_cross():
    """The cross term, measured where a truth exists.  U(rule) - U(inf) and ||eps||^2."""
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
    print(f"knee s^2 B = {KNEE:.5f}   cap K = {CAP_K:.5f}   a15 raw mse {base_raw:.5f}  "
          f"score {score(base_raw):.5f}   (185's known-answer row: 0.88228 / 0.11112)")

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
                Tq = rule_totals(R, rule).astype(np.float64)
                eps = np.log1p(TS * Tq / Tq.sum()) - np.log1p(TS * p)
                si = np.flatnonzero(D["genes"] == D["targets"][i])
                if si.size:
                    eps[si[0]] = 0.0
                acc.append(float((eps ** 2).sum() / (G - 1)))
            prof[(ctx, rule)] = float(np.mean(acc))
        print(f"  ctx {ctx}: " + "  ".join(f"{r} {prof[(ctx,r)]:.3e}" for r in RULES),
              flush=True)

    print("\n--- DISPERSION HALF, C_hat on the CHALLENGE, paired per block ---")
    chat, cp, nnz = {}, {}, {}
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
        del Ad, M
        G = D["genes"].size
        sel = np.linspace(0, D["sums"].shape[0] - 1, MAXT).astype(int)
        for i in sel:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            Lp = R.emit_L(1.0)
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
            sv = np.full(G, 1.0)
            sv[order[:min(8000, int(D["gate"].sum()))]] = 0.25
            for rule in RULES:
                Tw = rule_totals(R, rule)
                rng = np.random.default_rng(90_000 + int(i))
                v, c, ip = emit_counts((Tw + 0.5) / (n * Lp), np.full(n, Lp), DISPERSION,
                                       rng, mode="quasi", spread=sv, overdispersion=PHI)
                blk = sp.csr_matrix((v.astype(np.float64), c.astype(np.int64),
                                     np.asarray(ip)), shape=(n, G))
                assert np.array_equal(np.asarray(blk.sum(0)).ravel().astype(np.int64), Tw)
                chat.setdefault((ctx, rule), []).append(float(jackknife_correction(
                    blk, np.zeros(n, np.intp), 1, TS, chunk=64)[0]))
                nnz.setdefault((ctx, rule), []).append(int(blk.nnz))
                del blk
        print(f"  ctx {ctx}: C_p {cp[ctx]:.4f}  " + "  ".join(
            f"{r} Chat {np.mean(chat[(ctx,r)]):.3f} (/C_p "
            f"{np.mean(chat[(ctx,r)])/cp[ctx]:.3f}) nnz {np.mean(nnz[(ctx,r)]):,.0f}"
            for r in RULES) + f"  [{time.time()-t0:.0f}s]", flush=True)

    print("\n=== THE PRICE OF THE ROUNDING RULE ON `mse` ===")
    print(f"{'rule':>7}{'Chat/C_p':>10}{'rows<cap':>10}{'claim':>10}{'branch':>11}"
          f"{'d claim mse':>13}{'round energy':>14}{'d prof mse':>12}{'raw mse':>10}"
          f"{'score':>9}{'d overall':>11}{'nnz/block':>12}")
    rows = []
    for rule in RULES:
        ratios, claims, below = [], [], 0
        for ctx in "ABC":
            ch = np.asarray(chat[(ctx, rule)])
            cap = PRED_TRACE_CAP_K * cp[ctx]
            ratios.append(float(ch.mean() / cp[ctx]))
            claims.append(float(np.mean(np.minimum(ch, cap)) / cp[ctx]))
            below += int((ch < cap).sum())
        ratio = float(np.mean(ratios))
        claim = float(np.mean(claims)) * CAP_K
        branch = "SATURATED" if claim >= KNEE else "linear"
        d_claim = (min(KNEE, claim) - min(KNEE, 0.24581)) / sD * -1.0
        energy = float(np.mean([prof[(ctx, rule)] for ctx in "ABC"]))
        e_floor = float(np.mean([prof[(ctx, "floor")] for ctx in "ABC"]))
        xr = float(np.mean([v["ratio"] for k, v in cross.items() if k[1] == rule]))
        d_prof = (energy - e_floor) * xr / sD
        raw = base_raw + d_claim + d_prof
        mnnz = float(np.mean([np.mean(nnz[(ctx, rule)]) for ctx in "ABC"]))
        rows.append(dict(rule=rule, ratio=ratio, claim=claim, branch=branch,
                         d_claim_mse=d_claim, energy=energy, cross_ratio=xr,
                         d_prof_mse=d_prof, raw_mse=raw, score=score(raw),
                         d_overall=(score(raw) - score(base_raw)) / 6,
                         rows_below_cap=below, nnz=mnnz))
        print(f"{rule:>7}{ratio:>10.3f}{below:>10d}{claim:>10.5f}{branch:>11}"
              f"{d_claim:>+13.5f}{energy:>14.3e}{d_prof:>+12.5f}{raw:>10.5f}"
              f"{score(raw):>9.5f}{(score(raw)-score(base_raw))/6:>+11.5f}{mnnz:>12,.0f}")

    f = rows[0]
    print("\nKNOWN-ANSWER ROW (rule 10): the `floor` row must reproduce mse196's m = 1 row.")
    bad = []
    for k, w, g in (("Chat/C_p", KA["ratio"], f["ratio"]), ("claim", KA["claim"], f["claim"]),
                    ("raw mse", KA["raw"], f["raw_mse"])):
        tol = 2e-2 if k == "Chat/C_p" else 5e-3
        ok = abs(g - w) <= tol
        print(f"  {k:<10} {g:>10.6f}  expect {w:<10.6f} {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(k)
    if bad and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {bad}")
    print("  known-answer rows PASS")
    json.dump(dict(rows=rows, C_p=cp, MAXT=MAXT,
                   chat={f"{k[0]}|{k[1]}": v for k, v in chat.items()},
                   nnz={f"{k[0]}|{k[1]}": v for k, v in nnz.items()},
                   prof={f"{k[0]}|{k[1]}": v for k, v in prof.items()},
                   cross={f"{k[0]}|{k[1]}": v for k, v in cross.items()}),
              open(WORK / "mseham196.json", "w"), indent=1, default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
