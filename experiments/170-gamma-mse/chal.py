"""170 — translate the gamma sweep onto the CHALLENGE panel, through 096's `claim` route.

`A`, `B` and `claim` are properties of OUR OWN submission and the challenge's own control
cells, so all three are computable at every gamma with no access to the answers (096 §(a)).
`X` (alignment) and `sumD` (the reference's corrected effect) are not, so:

  * `sumD` is the board-solved constant -- a property of the data, identical for everyone.
  * `X` is solved from a15's OWN board point at gamma 0.65 (one equation, one unknown), and
    then carried across gamma by the RATIO measured on h1-test. That is rule 1 done properly:
    transfer the mechanism, not the magnitude.

The identity, applied PER CONTEXT (the 2026-08-29 correction to MSE-MECHANISM: `min(claim, B)`
is a per-context quantity and averaging it is wrong once any context sits on the rho cliff --
164 measured a15 at rho = 0.9466, i.e. BELOW it):

    MSE_c - 1 = ( A_c - 2 X_c + Kc_c - min(claim_c, B_c) ) / sumD_c

KNOWN-ANSWER GATE: run at gamma 0.65 with MASS_FIX unset it must reproduce a15's board raw
`mse`, and with MASS_FIX=minnorm it must reproduce a23's. Those are two independent board
points that differ in exactly one field.

    VCC_MASS_FIX=minnorm uv run python chal.py            # a23 geometry (the arm's base)
    uv run python chal.py                                 # a15 geometry
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments/096-mse-floor"))

from h1kit import (  # noqa: E402
    DISPERSION, KEEP_SPREAD, KEEP_TOP, MIN_CPM, NCELL, OVERDISPERSION, TAIL_SPREAD, TS,
    builder)

CTX = ("A", "B", "C")
NSUB = 24                      # targets per context used for `claim` (096 used 16)
GRID = [0.35, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]

# Board-anchored constants.
SUMD = 0.229141                # MSE-MECHANISM's 2026-08-29 corrected solve
MSE_MAP = (0.987862879, 0.950162566)     # score_mse = clip((b - raw)/span, 0, 1), 096 §(d)
A15_SCORE_MSE, A23_SCORE_MSE = 0.1200566151966586, 0.121581


def raw_from_score(s):
    return MSE_MAP[0] - s * MSE_MAP[1]


def context_library_medians():
    """Median control library per context, straight off the challenge bundle's own cells."""
    import h5py

    cache = WORK / "chal_libs.json"
    if cache.exists():
        return {k: float(v) for k, v in json.loads(cache.read_text()).items()}
    out = {}
    for c in CTX:
        with h5py.File(REPO / f"data/challenge/context_{c}.h5ad", "r") as h:
            ip = h["X/indptr"][:]
            data = h["X/data"]
            tot = np.empty(ip.size - 1, dtype=np.float64)
            step = 2000
            for lo in range(0, ip.size - 1, step):
                hi = min(lo + step, ip.size - 1)
                blk = np.asarray(data[int(ip[lo]):int(ip[hi])], dtype=np.float64)
                tot[lo:hi] = np.add.reduceat(blk, (ip[lo:hi] - ip[lo]).astype(np.int64))
        out[c] = float(np.median(tot))
    cache.write_text(json.dumps(out, indent=1))
    return out


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget, _row_weights
    from cell_eval2.moments import jackknife_correction
    from cell_eval2.prep import bulk_lognorm_means
    from vcc2026.donor_transfer import emit_counts

    B = builder()
    massfix = B.MASS_FIX
    tag = "minnorm" if massfix else "none"
    print(f"[chal] MASS_FIX={tag}  ALPHA={B.ALPHA}  CD4_TARGET_F={B.CD4_TARGET_F}", flush=True)

    ctrl = np.load(REPO / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    cpn = dict(np.load(REPO / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    pre = np.load(WORK / "pre.npz", allow_pickle=True)
    genes = pre["genes"].astype(str)
    perts = pre["perts"].astype(str)
    assert np.array_equal(genes, ctrl["genes"].astype(str))

    bctrl = {c: bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
             for c in CTX}
    # the identity's Kc is `sum_p w_p C_ctrl` = n_perts * C_ctrl / (G - 1), NOT
    # C_ctrl itself. 096's table: C_ctrl 0.3141/0.3229/0.3656 -> Kc_c
    # 0.005084/0.005228/0.005919.
    Kc = {c: 300 * float(ctrl[f"jk_{c}"]) / 18532 for c in CTX}
    Cp = {c: float(np.mean(cpn[c])) for c in CTX}
    Lmed = context_library_medians()
    print(f"[chal] Kc {[round(Kc[c], 6) for c in CTX]}  C_p proxy "
          f"{[round(Cp[c], 3) for c in CTX]}  median control library "
          f"{[int(Lmed[c]) for c in CTX]}", flush=True)

    # the scorer's own row weights and own-target columns, on the 300-target panel
    gpos = {g: i for i, g in enumerate(genes)}
    cols = {i: gpos[p] for i, p in enumerate(perts) if p in gpos}
    w = _row_weights(perts.size, genes.size, cols)
    sub = np.linspace(0, perts.size - 1, NSUB).astype(int)

    out = []
    for gamma in GRID:
        t0 = time.time()
        eff = np.load(WORK / "eff" / f"chal_g{gamma:.2f}.npz")["effect"]
        per = {}
        for c in CTX:
            b = bctrl[c]
            P = np.empty_like(eff)
            for i in range(eff.shape[0]):
                P[i] = B.profile_from_effect(eff[i], b, min_cpm=MIN_CPM)
            pm = np.log1p(TS * P)                      # the predicted bulk the scorer reads
            dh = pm - b[None, :]
            a = np.einsum("ij,ij->i", dh, dh)
            for k, col in cols.items():
                a[k] -= dh[k, col] ** 2
            A = float((a * w).sum())
            bud = float(_across_pert_budget(pm, cols, w))
            # `claim` on a subsample: min(C_hat, C_p) saturates for almost every row at a15's
            # policy (164 measured claim at 99.7% of the cap), so a 24-target draw resolves
            # it. `emit_counts(mode="quasi")` uses only `median(lib_sizes)`, so the context's
            # own median control library is the right and arm-independent choice.
            libs = np.full(NCELL, Lmed[c])
            vals = []
            for j in sub:
                spread = B.emission_spread_vector(
                    eff[j], b, min_cpm=MIN_CPM, tail_spread=TAIL_SPREAD,
                    keep_top=KEEP_TOP, keep_spread=KEEP_SPREAD, spread_map="block")
                v, cl, ip = emit_counts(P[j], libs, DISPERSION,
                                        np.random.default_rng(3000 + int(j)), mode="quasi",
                                        overdispersion=OVERDISPERSION, spread=spread)
                Xm = csr_matrix((v.astype(np.float64), cl.astype(np.int64),
                                 ip.astype(np.int64)), shape=(NCELL, genes.size))
                vals.append(jackknife_correction(Xm, np.zeros(NCELL, np.intp), 1, TS,
                                                 chunk=64)[0])
            chat = float(np.mean(vals))
            claim = float(perts.size * np.mean(
                np.minimum(vals, PRED_TRACE_CAP_K * Cp[c])) / (genes.size - 1))
            per[c] = dict(A=A, B=bud, claim=claim, chat=chat, Kc=Kc[c],
                          rho=min(1.0, bud / claim) if claim else 1.0)
            del P, pm, dh
        row = dict(gamma=gamma, per={c: per[c] for c in CTX},
                   A=float(np.mean([per[c]["A"] for c in CTX])),
                   B=float(np.mean([per[c]["B"] for c in CTX])),
                   claim=float(np.mean([per[c]["claim"] for c in CTX])),
                   chat=float(np.mean([per[c]["chat"] for c in CTX])),
                   secs=round(time.time() - t0, 1))
        out.append(row)
        print(f"[chal] gamma {gamma:.2f}  A {row['A']:.6f}  B {row['B']:.6f}  "
              f"claim {row['claim']:.6f}  C_hat {row['chat']:.2f}  "
              f"rho {[round(per[c]['rho'], 4) for c in CTX]}  [{row['secs']:.0f}s]",
              flush=True)
        with open(WORK / f"chal_{tag}.json", "w") as fh:
            json.dump(out, fh, indent=1)

    # --- known-answer check on the claim route itself: 164 measured a15's challenge-panel
    #     `claim` at 0.24581 against the 0.24653 cap, with rho = 0.9466.
    base = [r for r in out if r["gamma"] == 0.65][0]
    cap = float(np.mean([perts.size * Cp[c] / (genes.size - 1) for c in CTX]))
    print(f"\n[chal] CLAIM KNOWN-ANSWER: gamma 0.65 claim {base['claim']:.5f} vs 164's "
          f"0.24581 (cap {cap:.5f} vs 0.24653); rho "
          f"{[round(base['per'][c]['rho'], 4) for c in CTX]} vs 164's 0.9466")

    # --- solve X from the board point at gamma 0.65, then report the gate
    target_raw = raw_from_score(A23_SCORE_MSE if massfix else A15_SCORE_MSE)
    num = float(np.mean([base["per"][c]["A"] + base["per"][c]["Kc"]
                         - min(base["per"][c]["claim"], base["per"][c]["B"]) for c in CTX]))
    X = (num - SUMD * (target_raw - 1.0)) / 2.0
    print(f"\n[chal] board anchor: {'a23' if massfix else 'a15'} score_mse "
          f"{A23_SCORE_MSE if massfix else A15_SCORE_MSE} -> raw {target_raw:.6f}")
    print(f"[chal] solved X at gamma 0.65 = {X:.6f}   (096's v11/v12 solve gave 0.019240; "
          f"a15 differs from v11 in alpha, CD4 and the pruned operator)")
    with open(WORK / f"chal_{tag}_solve.json", "w") as fh:
        json.dump({"X_at_065": X, "sumD": SUMD, "target_raw": target_raw,
                   "anchor": "a23" if massfix else "a15"}, fh, indent=1)


if __name__ == "__main__":
    main()
