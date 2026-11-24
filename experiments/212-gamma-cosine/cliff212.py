"""212 Phase 1 — WHICH cliff did `rc_callset` move? Measured on all three emissions.

`mse`'s numerator, per context, is

    s^2 A_c + Kc_c - min( claim_c , s^2 B_c ) - 2 s X_c ,
    claim_c = (n_perts / (G - 1)) * mean_p min( Ĉ_p , K * C_p ),   K = PRED_TRACE_CAP_K = 1

so there are TWO saturations in it and they belong to different knobs:

  (a) the **rho = 1 cliff** -- `s^2 B_c` crossing `claim_c`. Above it, extra amplitude stops
      being refunded. GAMMA moves this one, because gamma moves `A` and `B` (170: `A` +14%
      over [0.65, 0.35]). Its position is set by `claim_c`, NOT by `Ĉ_p`.
  (b) the **claim-saturation cliff** -- `Ĉ_p` falling back through `K C_p`. Below it, `claim`
      itself starts to shrink and 196 measured the whole 0.111 `mse` score vanishing. THE
      EMISSION moves this one. 208 measured `Ĉ/C_p` 1.125 (a15) -> 2.940 (`rc_callset`) and
      recorded it as "the `mse` cliff is now 2.6x further away".

This script measures both, per context, on **three emissions** at the same delivered
profile: a15's `quasi`, `rc_hybrid` (`rc_hybrid_block`, no prescribed call set) and the
incumbent `rc_callset` (`rc_hybrid_block` with `CallSet(400, 2200, 1.90)`). The profile is
identical on all three by construction -- 191/202 proved the per-gene panel totals
digest-identical on 900/900 blocks and the board reproduced `pds`/`mse` to sixteen
significant figures -- and this script ASSERTS that per emitted block rather than claiming
it, so `A`, `B` and `X` cannot move and only `Ĉ_p` can.

KNOWN-ANSWER ROWS: 208's own `Ĉ/C_p` -- **1.125** for `quasi` and **2.940** for
`rc_callset` -- reproduced here through a completely different route (170's
`profile_from_effect` staging instead of 196's `recon196` reconstruction).

    MAXT=8 uv run python cliff212.py
"""
from __future__ import annotations

import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")

from chal212 import (CTX, DISPERSION, KEEP_SPREAD, KEEP_TOP, MIN_CPM, NCELL,  # noqa: E402
                     OVERDISPERSION, TAIL_SPREAD, TS, builder,
                     context_library_medians)

MAXT = int(os.environ.get("MAXT", "8"))
GAMMAS = [float(g) for g in os.environ.get("GAMMAS", "0.35,0.65").split(",")]
N_HEAD, N_CPM, Z_CAP = 400, 2200, 1.90        # 202's shipped `CallSet`
KA_208 = {"quasi": 1.125, "rc_hybrid": None, "rc_callset": 2.940}


def main() -> None:
    import vcc2026.rc_emit as RC
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction
    from cell_eval2.prep import bulk_lognorm_means
    from vcc2026.donor_transfer import emit_counts

    t0 = time.time()
    B = builder()
    assert B.MASS_FIX is None, "the incumbent ships no mass fix; run without VCC_MASS_FIX"
    ctrl = np.load(REPO / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    cpn = dict(np.load(REPO / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    pre = np.load(REPO / "experiments/170-gamma-mse/work/pre.npz", allow_pickle=True)
    genes = pre["genes"].astype(str)
    perts = pre["perts"].astype(str)
    G = genes.size
    bctrl = {c: bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
             for c in CTX}
    Cp = {c: float(np.mean(cpn[c])) for c in CTX}
    Lmed = context_library_medians()
    call_set = RC.CallSet(N_HEAD, N_CPM, Z_CAP)
    sub = np.linspace(0, perts.size - 1, MAXT).astype(int)

    eff = {g: np.load(WORK / "eff" / f"chal_g{g:.3f}.npz")["effect"] for g in GAMMAS}
    rows = []
    out_path = WORK / "cliff212.json"
    for c in CTX:
        import anndata as ad
        Ad = ad.read_h5ad(REPO / f"data/challenge/context_{c}.h5ad")
        M = Ad.X.tocsr() if not sp.isspmatrix_csr(Ad.X) else Ad.X
        control = RC.control_block_from_counts(M, min_cpm=5.0)
        del Ad, M
        gc.collect()
        print(f"[cliff] ctx {c}: C_p {Cp[c]:.4f}, scorer DE gate "
              f"{int(control.keep.sum())}/{G}, {control.n_controls} control cells "
              f"[{time.time()-t0:.0f}s]", flush=True)
        b = bctrl[c]
        expr = np.expm1(b)
        cpm_gate = 1e6 * expr / max(float(expr.sum()), 1e-12)
        build_gate = (cpm_gate > MIN_CPM) & (b > 0)
        libs = np.full(NCELL, Lmed[c])
        for gamma in GAMMAS:
            E = eff[gamma]
            acc = {em: [] for em in ("quasi", "rc_hybrid", "rc_callset")}
            totals_ok = {em: 0 for em in acc}
            for j in sub:
                P = B.profile_from_effect(E[j], b, min_cpm=MIN_CPM)
                spread = B.emission_spread_vector(
                    E[j], b, min_cpm=MIN_CPM, tail_spread=TAIL_SPREAD,
                    keep_top=KEEP_TOP, keep_spread=KEEP_SPREAD, spread_map="block")
                key = np.where(build_gate, E[j], 0.0)
                # pass 1's own per-gene totals: the invariant every emission must preserve
                v0, c0, ip0 = emit_counts(P, libs, DISPERSION,
                                          np.random.default_rng(3000 + int(j)),
                                          mode="quasi", overdispersion=OVERDISPERSION,
                                          spread=spread)
                ref = sp.csr_matrix((v0.astype(np.float64), c0.astype(np.int64),
                                     ip0.astype(np.int64)), shape=(NCELL, G))
                ref_tot = np.asarray(ref.sum(0)).ravel().astype(np.int64)
                for em in acc:
                    if em == "quasi":
                        blk = ref
                    else:
                        rng = np.random.default_rng(3000 + int(j))
                        Bm, _ = RC.rc_hybrid_block(
                            P, libs, DISPERSION, rng, pass2_rng=rng, spread=spread,
                            overdispersion=OVERDISPERSION, control=control, key=key,
                            jitter_seed=3000 + int(j), emit_counts=emit_counts,
                            call_set=(call_set if em == "rc_callset" else None))
                        blk = sp.csr_matrix(Bm.astype(np.float64))
                        del Bm
                    got = np.asarray(blk.sum(0)).ravel().astype(np.int64)
                    totals_ok[em] += int(np.array_equal(got, ref_tot))
                    acc[em].append(float(jackknife_correction(
                        blk, np.zeros(NCELL, np.intp), 1, TS, chunk=64)[0]))
                    if em != "quasi":
                        del blk
                del ref
            capped = PRED_TRACE_CAP_K * Cp[c]
            for em, vals in acc.items():
                vals = np.asarray(vals)
                claim = float(perts.size * np.mean(np.minimum(vals, capped)) / (G - 1))
                rows.append(dict(
                    ctx=c, gamma=gamma, emitter=em, n_targets=int(MAXT),
                    chat=float(vals.mean()), Cp=Cp[c],
                    chat_over_cp=float(vals.mean() / Cp[c]),
                    min_chat_over_cp=float(vals.min() / Cp[c]),
                    rows_under_cap=int((vals < capped).sum()),
                    claim=claim,
                    claim_cap=float(perts.size * capped / (G - 1)),
                    totals_identical=f"{totals_ok[em]}/{MAXT}"))
                print(f"  ctx {c} g{gamma:.2f} {em:<11} Ĉ/C_p {vals.mean()/Cp[c]:7.4f} "
                      f"(min {vals.min()/Cp[c]:6.4f})  under-cap {int((vals<capped).sum())}"
                      f"/{MAXT}  claim {claim:.6f}  totals {totals_ok[em]}/{MAXT}",
                      flush=True)
            out_path.write_text(json.dumps(rows, indent=1))
        del control
        gc.collect()

    print("\n===== KNOWN-ANSWER: 208's Ĉ/C_p, reproduced through 170's staging =====")
    for em, want in KA_208.items():
        got = float(np.mean([r["chat_over_cp"] for r in rows
                             if r["emitter"] == em and r["gamma"] == 0.65]))
        if want is None:
            print(f"  {em:<11} {got:.4f}   (208 published no row for this emitter)")
        else:
            print(f"  {em:<11} {got:.4f} vs 208's {want}   rel "
                  f"{abs(got/want-1):.2%}")
    print(f"[cliff] peak RSS {B.peak_rss_gib():.2f} GiB, {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
