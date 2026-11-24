"""211 Phase 0/1 — the emitted jackknife `C_hat` per CONTEXT, for the three board emitters, and
the per-context `claim` each of them buys.

208's Phase-4 by-product ("rc_callset sits at C_hat/C_p 2.940 against a15's 1.125, so the mse
cliff is 2.6x further away") is the premise of this experiment.  It is reproduced here
independently, on the same instrument, and then decomposed PER CONTEXT — because 187's cliff is
per context and the pooled ratio hides which context binds.

The identity the whole thing turns on (096 (a), 187's per-context correction):

    MSE_c - 1 = ( s^2 A_c + Kc_c - min(claim_c, s^2 B_c) - 2 s X ) / sumD
    claim_c   = sum_p w_p min( C_hat_p , PRED_TRACE_CAP_K * C_p )      PRED_TRACE_CAP_K = 1.0

`C_hat_p` is OURS (the delete-one jackknife of our own emitted 400 cells).  `C_p` is the
REFERENCE's, and `PRED_TRACE_CAP_K = 1.0` makes `K_c = sum_p w_p C_p` a HARD CAP on `claim_c`
that no emission can exceed.  So the amplitude cliff `s*_c = alpha * sqrt(claim_c / B_c)` has an
emission-independent ceiling at `alpha * sqrt(K_c / B_c)`.  This script measures how much of that
ceiling each emitter already collects.

Three emitters, one rule (`floor` — production's), one spread policy (`tier` — production's):

  quasi        a15 / a23's emission law.  Known-answer: 208's `quasi|floor|tier` row, C_hat/C_p
               1.125457, and per context 1.2316 / 1.0211 / 1.1236.
  rc_hybrid    the 2026-09-03 11:20 board point (`rc_hybrid_block`, call_set=None).
  rc_callset   the CURRENT INCUMBENT (`rc_hybrid_block` + the prescribed CallSet).  Known-answer:
               208's 2.939562, per context 2.8027 / 2.6999 / 3.3161.

Two row sets, both from the same per-index deterministic emissions (seed = f(index) only, so the
8-row subset reproduces 208 exactly whatever else is measured alongside it):
    MAXT = 8   -> reproduces 208's published numbers to the digit
    MAXT = 24  -> the row count 180 used for its per-context claim, for the cliff table

FALSIFIER built in: for every (ctx, target) the three emitters MUST deliver bit-identical
per-gene panel totals.  That is what pins `A`, `B`, `X` and `Kc` across the family and is why
the board reproduced `mse` to sixteen significant figures on both rc arms.  If it fails, 187's
solve does not transfer and this experiment stops.

    uv run python chat211.py            # ~8 min, peak RSS ~6 GB

Writes work/chat211.{json,log}.  Nothing outside experiments/211-amplitude-reopened/.
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
E196 = ROOT / "experiments/196-emitted-depth"
sys.path.insert(0, str(E196))
sys.path.insert(0, str(ROOT / "src"))
import recon196  # noqa: E402
import vcc2026.rc_emit as RC  # noqa: E402
from mse196 import DISPERSION, KNEE, PHI, S15, TS, n, raw_mse, score, sD  # noqa: E402
from mse196 import K as CAP_K  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

EFFECTS = ROOT / "experiments/178-composed-riders/work/eff.a15.npy"
BUNDLE = ROOT / "data/challenge"
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90
KEEP_TOP, KEEP_SPREAD = 8000, 0.25
EMITTERS = ("quasi", "rc_hybrid", "rc_callset")
MAXT_SMALL, MAXT_BIG = 8, 24

# 208's published rows (work/mse208.json), the known-answer gate.
KA208 = {"quasi": dict(pooled=1.125457, per_ctx={"A": 1.2316, "B": 1.0211, "C": 1.1236}),
         "rc_callset": dict(pooled=2.939562, per_ctx={"A": 2.8027, "B": 2.6999, "C": 3.3161})}
# 191/196 measured rc_hybrid at C_hat/C_p 2.786 on 196's own instrument (12 rows) — a soft check.
KA_RCH = 2.786


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction
    t0 = time.time()
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))  # noqa: E731

    P(f"PRED_TRACE_CAP_K = {PRED_TRACE_CAP_K}  -> claim_c <= K_c = sum_p w_p C_p, a HARD cap")
    P(f"knee s^2 B = {KNEE:.5f}   pooled cap K = {CAP_K:.5f}   sumD = {sD:.6f}   "
      f"s15 = {S15:.5f}   a15 raw mse {raw_mse(S15, 0.24581):.5f} score {score(raw_mse(S15, 0.24581)):.5f}")

    effects = np.load(EFFECTS)
    perts300 = np.asarray([ln.split(",")[0] for ln in
                           (BUNDLE / "pert_counts.csv").read_text().splitlines()[1:]])
    pos300 = {str(p): i for i, p in enumerate(perts300)}
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)

    chat, cp, totals_match = {}, {}, {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        G = D["genes"].size
        assert effects.shape[1] == G
        Ad = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad")
        M = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        ids = np.asarray(Ad.obs["ntc_id"].astype(str).values)
        names = list(dict.fromkeys(ids.tolist()))
        cps = [float(jackknife_correction(sp.csr_matrix(M[np.flatnonzero(ids == nm)]),
                                          np.zeros(int((ids == nm).sum()), np.intp), 1, TS,
                                          chunk=64)[0]) for nm in names]
        cp[ctx] = float(np.mean(cps))
        control = RC.control_block_from_counts(M, min_cpm=5.0)
        del Ad, M
        P(f"\nctx {ctx}: C_p {cp[ctx]:.4f} over {len(cps)} NTC constructs "
          f"(sd {np.std(cps):.4f}), scorer gate {int(control.keep.sum())}/{G}  [{time.time()-t0:.0f}s]")

        nrows = D["sums"].shape[0]
        sel8 = np.linspace(0, nrows - 1, MAXT_SMALL).astype(int)
        sel24 = np.linspace(0, nrows - 1, MAXT_BIG).astype(int)
        union = sorted(set(sel8.tolist()) | set(sel24.tolist()))
        for i in union:
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            Lp = R.emit_L(1.0)
            M1 = D["sums"][i].astype(np.int64)
            pcpm = 1e6 * M1 / max(int(M1.sum()), 1)
            eff = np.abs(np.log2((pcpm + 1e-9) / (D["cpm"] + 1e-9)))
            order = np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
            sv = np.full(G, 1.0)
            sv[order[:min(KEEP_TOP, int(D["gate"].sum()))]] = KEEP_SPREAD
            key = np.where(D["gate"], effects[pos300[str(D["targets"][i])]], 0.0)
            Tw = R.totals(1.0)                                    # `floor`, production's rule
            prof_p = (Tw + 0.5) / (n * Lp)
            for em in EMITTERS:
                rng = np.random.default_rng(90_000 + int(i))
                if em == "quasi":
                    v, c, ip = emit_counts(prof_p, np.full(n, Lp), DISPERSION, rng,
                                           mode="quasi", spread=sv, overdispersion=PHI)
                    blk = sp.csr_matrix((v.astype(np.float64), c.astype(np.int64),
                                         np.asarray(ip)), shape=(n, G))
                else:
                    B, _dg = RC.rc_hybrid_block(
                        prof_p, np.full(n, Lp), DISPERSION, rng, pass2_rng=rng,
                        spread=sv, overdispersion=PHI, control=control, key=key,
                        jitter_seed=1000 + int(i), emit_counts=emit_counts,
                        call_set=(call_set if em == "rc_callset" else None))
                    blk = sp.csr_matrix(B.astype(np.float64))
                    del B
                got = np.asarray(blk.sum(0)).ravel().astype(np.int64)
                # FALSIFIER: every emitter must deliver the SAME per-gene panel totals.
                totals_match[(ctx, em, int(i))] = bool(np.array_equal(got, Tw))
                chat.setdefault((ctx, em), {})[int(i)] = float(jackknife_correction(
                    blk, np.zeros(n, np.intp), 1, TS, chunk=64)[0])
                del blk
        for em in EMITTERS:
            v8 = np.array([chat[(ctx, em)][int(i)] for i in sel8])
            v24 = np.array([chat[(ctx, em)][int(i)] for i in sel24])
            P(f"  {em:<11} MAXT8  Chat {v8.mean():8.3f}  /C_p {v8.mean()/cp[ctx]:6.4f}  "
              f"min/C_p {(v8/cp[ctx]).min():6.4f}  rows<cap {(v8 < cp[ctx]).sum()}/{v8.size}   "
              f"|  MAXT24 Chat {v24.mean():8.3f}  /C_p {v24.mean()/cp[ctx]:6.4f}  "
              f"min/C_p {(v24/cp[ctx]).min():6.4f}  rows<cap {(v24 < cp[ctx]).sum()}/{v24.size}"
              f"  [{time.time()-t0:.0f}s]")
        del control, D

    bad_totals = [k for k, v in totals_match.items() if not v]
    P(f"\n=== FALSIFIER: per-gene panel totals identical across the three emitters ===")
    P(f"  {len(totals_match) - len(bad_totals)}/{len(totals_match)} (ctx, emitter, target) blocks "
      f"deliver EXACTLY the `floor` totals.  {'PASS' if not bad_totals else 'FAIL ' + str(bad_totals[:5])}")
    P("  -> A_c, B_c, X and Kc_c are functions of the delivered pseudobulk alone and are therefore")
    P("     IDENTICAL across a15 / rc_hybrid / rc_callset.  Only `claim_c` can move.")
    if bad_totals:
        raise SystemExit("ABORT: the pseudobulk is not shared; 187's solve does not transfer.")

    # ---- per-context ratios, claim fractions, and the KA gates ---------------------------
    out = {"C_p": cp, "PRED_TRACE_CAP_K": float(PRED_TRACE_CAP_K), "cap_K_pooled": CAP_K,
           "knee": KNEE, "sumD": sD, "chat": {}, "rows": {}}
    for (ctx, em), d in chat.items():
        out["chat"][f"{ctx}|{em}"] = {str(k): v for k, v in sorted(d.items())}

    P(f"\n=== per-context `claim` fraction of the cap, q_c = mean_p min(C_hat_p, C_p) / C_p ===")
    P(f"  {'emitter':<11}{'rows':>6}" + "".join(f"{'q_' + c:>10}" for c in "ABC")
      + f"{'q pooled':>11}{'Chat/C_p pooled':>17}")
    for maxt, tag in ((MAXT_SMALL, "MAXT8"), (MAXT_BIG, "MAXT24")):
        P(f"  --- {tag} ---")
        for em in EMITTERS:
            qs, rs = {}, {}
            for ctx in "ABC":
                nrows = len(chat[(ctx, em)])
                sel = np.linspace(0, max(chat[(ctx, em)].keys()), maxt).astype(int)
                v = np.array([chat[(ctx, em)][int(i)] for i in sel])
                qs[ctx] = float(np.mean(np.minimum(v, cp[ctx])) / cp[ctx])
                rs[ctx] = float(v.mean() / cp[ctx])
            qp = float(np.mean(list(qs.values())))
            rp = float(np.mean(list(rs.values())))
            out["rows"][f"{em}|{tag}"] = dict(q=qs, ratio=rs, q_pooled=qp, ratio_pooled=rp,
                                              claim_pooled=qp * CAP_K)
            P(f"  {em:<11}{maxt:>6}" + "".join(f"{qs[c]:>10.6f}" for c in "ABC")
              + f"{qp:>11.6f}{rp:>17.4f}")

    P(f"\n=== KNOWN-ANSWER GATES against 208's published rows (rule 10) ===")
    ok = True
    for em, want in KA208.items():
        got = out["rows"][f"{em}|MAXT8"]
        g = abs(got["ratio_pooled"] - want["pooled"]) < 1e-4
        ok &= g
        P(f"  {em:<11} pooled C_hat/C_p {got['ratio_pooled']:.6f}  expect {want['pooled']:.6f}"
          f"  {'PASS' if g else 'FAIL'}")
        for c in "ABC":
            gc = abs(got["ratio"][c] - want["per_ctx"][c]) < 5e-4
            ok &= gc
            P(f"      ctx {c}  {got['ratio'][c]:.4f}  expect {want['per_ctx'][c]:.4f}  "
              f"{'PASS' if gc else 'FAIL'}")
    rch = out["rows"]["rc_hybrid|MAXT8"]["ratio_pooled"]
    P(f"  rc_hybrid   pooled C_hat/C_p {rch:.4f}  (196/191 recorded {KA_RCH} on a 12-row "
      f"instrument; soft check, {'consistent' if abs(rch - KA_RCH) < 0.5 else 'DIVERGENT'})")
    out["ka_pass"] = bool(ok)
    P(f"  -> {'ALL 208 GATES PASS' if ok else 'GATE FAILURE — the premise is unverified'}")

    import resource
    out["peak_rss_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    P(f"  peak RSS (getrusage, this process) {out['peak_rss_gb']:.2f} GB")

    (WORK / "chat211.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "chat211.log").write_text("\n".join(lines))
    P(f"\ndone [{time.time()-t0:.0f}s]")
    if not ok and not os.environ.get("NOABORT"):
        raise SystemExit("ABORT (rule 10)")


if __name__ == "__main__":
    main()
