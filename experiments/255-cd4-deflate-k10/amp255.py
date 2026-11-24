"""255 falsifier 5 -- delivered amplitude off the PACKAGE (vs cd4svd5, the incumbent, and vs h1union
for the record), and the predicted board raw `mse` from 197's identity at the arm's OWN terms
(243 section 5.2: never the floor).  249's `amp249.py` with the incumbent changed.

Streams the arm's per-block panel totals (214's `block_totals`, imported from 243's
`totals243.py`), computes 211's per-context (A_c, B_c, Kc_c, cap_c) with 243's `terms`.

KNOWN-ANSWER ROWS: (1) the shipped `submissions/222-h1union-g028.h5ad` streamed the same way must
equal 243's `invert243.json` (A_c, B_c, cap_c) to 1e-9 on all three contexts and 197's identity at
those terms must reproduce its board raw `mse` 0.8172630; (2) `submissions/249-cd4svd5-f034.h5ad`
streamed the same way must equal 249's `amp249.json` per-context terms to 1e-9.

Predicted raw `mse` for the arm: raw = 1 + (mean_c[A_c + Kc_c - min(cap_c, B_c)] - 2 X)/sumD with
X = cos * sqrt(sumD) * mean_c sqrt(A_c) and cos = cos_h1union * (1 + k_cos * rel_agg), rel_agg the
three-cancer-fold PLOLO aggregate read vs h1union (repro255), at k_cos in {0, 0.94, 1.00, 1.36}
(243's classes); the same for cd4svd5 (its own terms, its own rel_agg) so the k 5 -> 10 step reads
as a difference of two predictions.  Also the count caps for falsifier 6 (nnz, max counts per cell).

    uv run python amp255.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E243 = ROOT / "experiments/243-transfer-factor"
E249 = ROOT / "experiments/249-cd4-deflate"
sys.path.insert(0, str(E243))
from totals243 import block_totals  # noqa: E402
from invert243 import bctrl, terms  # noqa: E402
from lib243 import CTX, SUMD, rss_gb  # noqa: E402

ARM = ROOT / "submissions/255-cd4svd10-f034.h5ad"
INC5 = ROOT / "submissions/249-cd4svd5-f034.h5ad"
INC = ROOT / "submissions/222-h1union-g028.h5ad"
INV = json.loads((E243 / "work/invert243.json").read_text())
REPRO = json.loads((WORK / "repro255-summary.json").read_text())
AMP249 = json.loads((E249 / "work/amp249.json").read_text())
K_COS = (0.0, 0.94, 1.00, 1.36)


def stream(tag, path):
    out = WORK / f"totals-{tag}.npz"
    if out.exists():
        z = np.load(out, allow_pickle=True)
        return {k: z[k] for k in z.files}
    t0 = time.time()
    r = block_totals(path)
    np.savez_compressed(out, **r)
    print(f"  streamed {path.name}: nnz {r['nnz']} max_cell {r['max_cell']} [{time.time()-t0:.0f}s] "
          f"rss {rss_gb():.2f} GB", flush=True)
    return r


def raw_for(per, cos):
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["cap"], per[c]["B"]) for c in CTX]))
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
    X = cos * np.sqrt(SUMD) * msa
    return 1.0 + (num - 2.0 * X) / SUMD


def main() -> None:
    genes, b = bctrl()
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    # ---- known-answer row 1: the shipped incumbent's terms reproduce 243's ------------------
    inc = stream("h1union", INC)
    assert np.array_equal(np.asarray(inc["genes"]).astype(str), genes)
    per_inc = terms(np.asarray(inc["totals"]), np.asarray(inc["blocks"]).astype(str), genes, b)
    ref = INV["arms"]["h1union"]
    worst = max(abs(per_inc[c][k] - ref["per"][c][k]) for c in CTX for k in ("A", "B", "cap"))
    P(f"KA1: h1union (A, B, cap) streamed here vs 243 invert243.json: max|d| {worst:.2e} "
      f"{'PASS' if worst < 1e-9 else 'FAIL'}; nnz {int(inc['nnz'])}")
    assert worst < 1e-9
    m222 = json.loads((ROOT / "submissions/222-h1union-g028.h5ad.meta.json").read_text())
    assert int(inc["nnz"]) == m222["output"]["nnz"]

    # ---- known-answer row 2: cd4svd5's terms reproduce 249's amp249.json -------------------
    inc5 = stream("cd4svd5", INC5)
    assert np.array_equal(np.asarray(inc5["genes"]).astype(str), genes)
    per5 = terms(np.asarray(inc5["totals"]), np.asarray(inc5["blocks"]).astype(str), genes, b)
    worst5 = max(abs(per5[c][k] - AMP249["per_arm"][c][k]) for c in CTX for k in ("A", "B", "cap", "Kc"))
    P(f"KA2: cd4svd5 (A, B, cap, Kc) streamed here vs 249 amp249.json: max|d| {worst5:.2e} "
      f"{'PASS' if worst5 < 1e-9 else 'FAIL'}; nnz {int(inc5['nnz'])} (249 meta {AMP249['nnz']})")
    assert worst5 < 1e-9 and int(inc5["nnz"]) == AMP249["nnz"]

    # ---- the arm ---------------------------------------------------------------------------
    arm = stream("cd4svd10", ARM)
    assert np.array_equal(np.asarray(arm["genes"]).astype(str), genes)
    meta = json.loads((ROOT / "submissions/255-cd4svd10-f034.h5ad.meta.json").read_text())
    assert int(arm["nnz"]) == meta["output"]["nnz"], "streamed nnz != meta nnz"
    per = terms(np.asarray(arm["totals"]), np.asarray(arm["blocks"]).astype(str), genes, b)
    P(f"\ncd4svd10 vs cd4svd5 (and h1union), delivered terms (211's convention):")
    for c in CTX:
        P(f"  {c}: A {per[c]['A']:.6f} (cd4svd5 {per5[c]['A']:.6f}, x{per[c]['A']/per5[c]['A']:.5f}; h1union {per_inc[c]['A']:.6f}, x{per[c]['A']/per_inc[c]['A']:.5f})  "
          f"B {per[c]['B']:.6f} (cd4svd5 {per5[c]['B']:.6f})  cap {per[c]['cap']:.6f}  "
          f"leg {'claim' if per[c]['cap'] < per[c]['B'] else 'B'}")
    amp5 = float(np.mean([np.sqrt(per[c]["A"] / per5[c]["A"]) for c in CTX]))
    ampB5 = float(np.mean([np.sqrt(per[c]["B"] / per5[c]["B"]) for c in CTX]))
    amp = float(np.mean([np.sqrt(per[c]["A"] / per_inc[c]["A"]) for c in CTX]))
    ampB = float(np.mean([np.sqrt(per[c]["B"] / per_inc[c]["B"]) for c in CTX]))
    P(f"  delivered amplitude sqrt(A_arm/A_cd4svd5) - 1, mean over contexts: {100*(amp5-1):+.4f}%   "
      f"(B-term {100*(ampB5-1):+.4f}%)   falsifier 5 bound +-0.5% vs cd4svd5: {'PASS' if abs(amp5-1) <= 0.005 else 'FAIL'}")
    P(f"  vs h1union: sqrt(A_arm/A_h1union) - 1 = {100*(amp-1):+.4f}%   (B-term {100*(ampB-1):+.4f}%); "
      f"cd4svd5 vs h1union was {AMP249['amplitude_pct']:+.4f}%")

    # ---- the predicted board raw mse at each arm's own terms ---------------------------------
    cos_inc = ref["inv"]["cos"]
    raw_inc = ref["board"]["raw"]["mse"]
    rel = REPRO["arm"]["rel_agg"]
    rel5 = REPRO["arm5"]["rel_agg"]
    chk = raw_for(per_inc, cos_inc)
    P(f"\n197's identity at h1union's own terms reproduces its board raw mse: {chk:.7f} vs {raw_inc:.7f} "
      f"(|d| {abs(chk-raw_inc):.1e})")
    P(f"PLOLO aggregate cosine read vs h1union (3 cancer folds, repro255): k=10 {100*rel:+.4f}%, k=5 {100*rel5:+.4f}% "
      f"(249: {100*AMP249['rel_agg_plolo']:+.4f}%); h1union measured cos {cos_inc:.6f}")
    P(f"  {'k_cos':>6s}{'arm cos':>10s}{'d cos':>10s}{'pred raw mse k=10':>19s}{'pred k=5':>10s}{'k10 - k5':>10s}{'d vs h1union':>14s}{'mse leg vs cd4svd5':>20s}")
    preds, preds5 = {}, {}
    for k in K_COS:
        cos = cos_inc * (1 + k * rel)
        cos5 = cos_inc * (1 + k * rel5)
        raw = raw_for(per, cos)
        raw5 = raw_for(per5, cos5)
        preds[str(k)] = raw
        preds5[str(k)] = raw5
        P(f"  {k:>6.2f}{cos:>10.6f}{cos-cos_inc:>+10.6f}{raw:>19.5f}{raw5:>10.5f}{raw-raw5:>+10.5f}{raw-raw_inc:>+14.5f}{-(raw-raw5)/0.950162566/6:>+20.5f}")
    P(f"  (k_cos 0 is the amplitude/budget leg alone; 249's own table for cd4svd5 read 0.81809 / 0.81902 / 0.81908 / 0.81943; "
      f"the mse leg column converts the k 5 -> 10 raw difference to overall units at the scorer's mse span 0.950162566)")

    # ---- the count caps (falsifier 6) -------------------------------------------------------
    P(f"\ncaps: nnz {int(arm['nnz']):,} = {100*int(arm['nnz'])/4.75e9:.2f}% of 4.75e9; max counts per cell "
      f"{int(arm['max_cell']):,} vs 1,000,000 ({1e6/int(arm['max_cell']):.1f}x under); cd4svd5 nnz {int(inc5['nnz']):,} "
      f"({100*(int(arm['nnz'])/int(inc5['nnz'])-1):+.3f}%), h1union nnz {int(inc['nnz']):,} ({100*(int(arm['nnz'])/int(inc['nnz'])-1):+.3f}%)")
    out = dict(per_arm=per, per_cd4svd5=per5, per_h1union=per_inc,
               amplitude_ratio_vs_cd4svd5=amp5, amplitude_pct_vs_cd4svd5=100*(amp5-1), B_ratio_vs_cd4svd5=ampB5,
               amplitude_ratio_vs_h1union=amp, amplitude_pct_vs_h1union=100*(amp-1), B_ratio_vs_h1union=ampB,
               rel_agg_plolo_k10=rel, rel_agg_plolo_k5=rel5, cos_inc=cos_inc, raw_inc=raw_inc,
               predicted_raw_mse_k10=preds, predicted_raw_mse_k5=preds5,
               nnz=int(arm["nnz"]), max_cell=int(arm["max_cell"]), nnz_cd4svd5=int(inc5["nnz"]), nnz_h1union=int(inc["nnz"]),
               falsifier5_pass=bool(abs(amp5-1) <= 0.005), ka1_max_abs_diff=worst, ka2_max_abs_diff=worst5, peak_rss_gb=rss_gb())
    (WORK / "amp255.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "amp255.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'amp255.json'}  peak rss {rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
