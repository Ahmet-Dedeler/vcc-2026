"""211 Phase 1 — `A_c`, `B_c` and `X` measured DIRECTLY on the two shipped artifacts.

187's cliff table was computed on **a23's** profile (board raw `mse` 0.872353).  The current
incumbent chain — a15 -> `rc_hybrid` -> `rc_callset` — carries **a15's** profile (board raw
0.8738045, reproduced to sixteen significant figures on both rc arms).  Rather than transport
187's a23 terms, this script measures `A_c` and `B_c` on the artifacts themselves, with
`cell_eval2`'s own `_across_pert_budget` and `_row_weights`, so the cliff is located on the
object that is actually on the board.

    A_c  = sum_p w_p ( ||dh_p||^2 - dh[p, target_p]^2 )          w_p = 1/(G-1)
    B_c  = _across_pert_budget(mu_hat, cols, w)                  the SUBMISSION's own budget
    Kc_c = the reference's control term (096 section (e), unchanged by us)
    X    = solved from the arm's own board raw `mse` (one equation, per 187/096)

    cliff_c = alpha * sqrt( claim_c / B_c )      and      claim_c <= claim_cap_c = sum_p w_p C_p

Two artifacts, one code path, so the RATIO is the measurement that matters:

  submissions/202-rc_callset.h5ad   the incumbent, board 0.1918638367, raw `mse` 0.8738045
  submissions/210-cd4f034.h5ad      built, not yet scored; 210 measured its delivered
                                    `A` at 0.97970x on `amp210.py`'s (unweighted, all-gene)
                                    definition.  Whether `B` tracks `A` is the open question
                                    that decides whether cd4f034 opens amplitude headroom, and
                                    210 did not measure it.

KNOWN-ANSWER ROWS (rule 10):
  1. amp210's own unweighted all-gene `A`: inc 14.708266, arm 14.409672, ratio 0.9796989.
  2. the per-context `A_c` must land within ~2% of 187/180's a23 values
     (0.234830 / 0.224778 / 0.234941) — a15 and a23 differ by the mass fix, board raw `mse`
     0.873805 vs 0.872353, so a small offset is expected and a large one is a bug.
  3. `B_c/A_c` must sit in the 0.993-0.995 band every measurement of this object has produced.

    uv run python bx211.py

Reads two artifacts out of submissions/ and writes NOTHING there.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
E210 = ROOT / "experiments/210-cd4-build"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(E210))
from profile210 import control_frames, read_obs_blocks, N_GENES  # noqa: E402
from vcc2026.donor_transfer import TS  # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402

ARMS = {"rc_callset": ROOT / "submissions/202-rc_callset.h5ad",
        "cd4f034": ROOT / "submissions/210-cd4f034.h5ad"}
BOARD_RAW = {"rc_callset": 0.8738045054832962}          # 197 section 5 / 185 BOARD-RESULT
SUMD = 0.229141                                          # 187 [KA-i], the C_p-free LS solve
KC = {"A": 0.005083962657824914, "B": 0.005227704077951962, "C": 0.005919139136949927}  # 096 (e)
CP46 = {"A": 14.439751080335085, "B": 14.805734135005, "C": 16.439650143279888}          # 180
A23_187 = {"A": 0.2348299009419585, "B": 0.22477784721536484, "C": 0.23494104987562134}
B23_187 = {"A": 0.23359310250430515, "B": 0.22362398630982533, "C": 0.2337275546353769}
KA_AMP210 = dict(inc=14.708265772725685, arm=14.409672249221215, ratio=0.9796989306476792)


def totals_of(path: Path, blocks, indptr) -> np.ndarray:
    out = np.zeros((len(blocks), N_GENES), dtype=np.float64)
    with h5py.File(path, "r") as handle:
        data, indices = handle["X"]["data"], handle["X"]["indices"]
        for k, (_key, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            out[k] = np.bincount(indices[lo:hi].astype(np.int64),
                                 weights=data[lo:hi].astype(np.float64), minlength=N_GENES)
    return out


def main() -> None:
    t0 = time.time()
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))  # noqa: E731

    frames = control_frames()
    with h5py.File(ROOT / "submissions/202-rc_callset.h5ad", "r") as h:
        gene_names = np.asarray(h["var"]["_index"][:]).astype(str)
    gpos = {g: i for i, g in enumerate(gene_names)}
    P(f"gene axis {gene_names.size}  [{time.time()-t0:.0f}s]")

    out, eff_store = {}, {}
    for name, path in ARMS.items():
        blocks, indptr, shape = read_obs_blocks(path)
        keys = [k for k, _, _ in blocks]
        tot = totals_of(path, blocks, indptr)
        P(f"\n{name}: {len(blocks)} blocks, shape {shape}  [{time.time()-t0:.0f}s]")
        eff = np.empty_like(tot)
        for k, (ctx, _t) in enumerate(keys):
            eff[k] = np.log1p(TS * tot[k] / tot[k].sum()) - frames[ctx]["b_ctrl"]
        eff_store[name] = (keys, eff)
        per = {}
        for ctx in "ABC":
            rows = [i for i, k in enumerate(keys) if k[0] == ctx]
            mu = np.log1p(TS * tot[rows] / tot[rows].sum(1, keepdims=True))
            dh = mu - frames[ctx]["b_ctrl"][None, :]
            tgts = [keys[i][1] for i in rows]
            cols = {j: gpos[t] for j, t in enumerate(tgts) if t in gpos}
            w = _row_weights(len(rows), N_GENES, cols)
            a = np.einsum("ij,ij->i", dh, dh)
            for j, col in cols.items():
                a[j] -= dh[j, col] ** 2
            A = float((a * w).sum())
            B = float(_across_pert_budget(mu, cols, w))
            per[ctx] = dict(A=A, B=B, n_rows=len(rows), n_resolved=len(cols),
                            weight_sum=float(w.sum()), Kc=KC[ctx],
                            claim_cap=float(w.sum()) * CP46[ctx])
            P(f"  ctx {ctx}: rows {len(rows)} resolved {len(cols)}  A {A:.6f}  B {B:.6f}  "
              f"B/A {B/A:.5f}  weight_sum {w.sum():.8f}  claim_cap {per[ctx]['claim_cap']:.6f}")
        out[name] = per
        del tot

    # ---- KA 1: amp210's unweighted, all-gene A on the SAME two artifacts -------------------
    P(f"\n=== KA-1: amp210's unweighted all-gene `A = mean ||row||^2` ===")
    amp = {}
    for name in ARMS:
        _keys, eff = eff_store[name]
        amp[name] = float(np.mean(np.linalg.norm(eff, axis=1) ** 2))
    ratio = amp["cd4f034"] / amp["rc_callset"]
    ok1 = (abs(amp["rc_callset"] - KA_AMP210["inc"]) < 1e-4
           and abs(amp["cd4f034"] - KA_AMP210["arm"]) < 1e-4
           and abs(ratio - KA_AMP210["ratio"]) < 1e-6)
    P(f"  rc_callset {amp['rc_callset']:.6f} (210: {KA_AMP210['inc']:.6f})   "
      f"cd4f034 {amp['cd4f034']:.6f} (210: {KA_AMP210['arm']:.6f})   "
      f"ratio {ratio:.7f} (210: {KA_AMP210['ratio']:.7f})   {'PASS' if ok1 else 'FAIL'}")

    # ---- KA 2/3 and the scored ratios -------------------------------------------------------
    P(f"\n=== KA-2: per-context `A_c` against 187/180's a23 values (a15 vs a23 = the mass fix) ===")
    P(f"  {'ctx':>4}{'A(rc_callset)':>15}{'A(a23, 187)':>13}{'rel':>9}"
      f"{'B(rc_callset)':>15}{'B(a23, 187)':>13}{'rel':>9}")
    ok2 = True
    for ctx in "ABC":
        r = out["rc_callset"][ctx]
        dA, dB = r["A"] / A23_187[ctx] - 1, r["B"] / B23_187[ctx] - 1
        ok2 &= abs(dA) < 0.03 and abs(dB) < 0.03
        P(f"  {ctx:>4}{r['A']:>15.6f}{A23_187[ctx]:>13.6f}{dA:>+9.2%}"
          f"{r['B']:>15.6f}{B23_187[ctx]:>13.6f}{dB:>+9.2%}")
    P(f"  {'PASS' if ok2 else 'FAIL'} (tolerance 3%: a15 and a23 are different profiles)")

    P(f"\n=== KA-3: B/A in the 0.993-0.995 band, and the SCORED per-context A and B ratios ===")
    P(f"  {'ctx':>4}{'A ratio':>10}{'B ratio':>10}{'B/A rc':>10}{'B/A cd4':>10}")
    ok3 = True
    ratios = {}
    for ctx in "ABC":
        r, a = out["rc_callset"][ctx], out["cd4f034"][ctx]
        ratios[ctx] = dict(A=a["A"] / r["A"], B=a["B"] / r["B"])
        ok3 &= 0.990 < r["B"] / r["A"] < 0.999 and 0.990 < a["B"] / a["A"] < 0.999
        P(f"  {ctx:>4}{ratios[ctx]['A']:>10.5f}{ratios[ctx]['B']:>10.5f}"
          f"{r['B']/r['A']:>10.5f}{a['B']/a['A']:>10.5f}")
    mA = float(np.mean([ratios[c]["A"] for c in "ABC"]))
    mB = float(np.mean([ratios[c]["B"] for c in "ABC"]))
    P(f"  mean over contexts: A ratio {mA:.5f}   B ratio {mB:.5f}   "
      f"B/A of the ratios {mB/mA:.5f}   {'PASS' if ok3 else 'FAIL'}")
    P(f"  -> cd4f034's `B` {'TRACKS' if abs(mB - mA) < 0.002 else 'does NOT track'} its `A`; "
      f"210 measured the `A` side only.")

    # ---- X on the incumbent, from its own board raw `mse` -----------------------------------
    P(f"\n=== `X` solved from rc_callset's board raw `mse` {BOARD_RAW['rc_callset']:.10f} ===")
    P("  (claim_c held at its cap: chat211 measures 0 rows under the cap on all three contexts)")
    num = float(np.mean([out["rc_callset"][c]["A"] + KC[c]
                         - min(out["rc_callset"][c]["claim_cap"], out["rc_callset"][c]["B"])
                         for c in "ABC"]))
    X = (num - (BOARD_RAW["rc_callset"] - 1.0) * SUMD) / 2.0
    P(f"  X = {X:.6f}   (187 solved 0.017930 on a23's profile; 170 solved 0.017766 on a15's)")

    import resource
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    P(f"\npeak RSS (getrusage, this process) {rss:.2f} GB")

    res = dict(per_context=out, amp_unweighted=amp, amp_ratio=ratio, ratios=ratios,
               A_ratio_mean=mA, B_ratio_mean=mB, X_rc_callset=X, sumD=SUMD, Kc=KC, cp46=CP46,
               peak_rss_gb=rss, ka=dict(ka1=bool(ok1), ka2=bool(ok2), ka3=bool(ok3)))
    (WORK / "bx211.json").write_text(json.dumps(res, indent=1, default=float))
    (WORK / "bx211.log").write_text("\n".join(lines))
    P(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
