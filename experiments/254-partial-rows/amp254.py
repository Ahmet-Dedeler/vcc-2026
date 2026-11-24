"""254 falsifier 5 -- delivered amplitude off the PACKAGE (243's stream), the arm against cd4svd5
(the mandate's incumbent) and against h1union, and the predicted board raw `mse` from 197's
identity at the arm's OWN terms (243 section 5.2: never the floor).

249's `amp249.py` with a second incumbent.  KNOWN-ANSWER ROWS: h1union's re-streamed (A_c, B_c, cap_c)
== 243's `invert243.json` to 1e-9 and 197's identity reproduces its board raw `mse` 0.8172630;
cd4svd5's re-streamed per-block totals == 249's `work/totals-cd4svd5.npz` bit for bit and its
amplitude vs h1union reproduces 249's -0.409%.

    uv run python amp254.py
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

ARM = ROOT / "submissions/254-cd4svd5-partial-f034.h5ad"
SVD5 = ROOT / "submissions/249-cd4svd5-f034.h5ad"
INC = ROOT / "submissions/222-h1union-g028.h5ad"
INV = json.loads((E243 / "work/invert243.json").read_text())
REPRO = json.loads((WORK / "repro254-summary.json").read_text())
K_COS = (0.0, 0.94, 1.00, 1.36)


def stream(tag, path):
    out = WORK / f"totals-{tag}.npz"
    if out.exists():
        z = np.load(out, allow_pickle=True)
        return {k: z[k] for k in z.files}
    t0 = time.time()
    r = block_totals(path)
    np.savez_compressed(out, **r)
    print(f"  streamed {path.name}: nnz {r['nnz']} max_cell {r['max_cell']} [{time.time()-t0:.0f}s] rss {rss_gb():.2f} GB", flush=True)
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

    # ---- known-answer rows -------------------------------------------------------------------
    inc = stream("h1union", INC)
    assert np.array_equal(np.asarray(inc["genes"]).astype(str), genes)
    per_inc = terms(np.asarray(inc["totals"]), np.asarray(inc["blocks"]).astype(str), genes, b)
    ref = INV["arms"]["h1union"]
    worst = max(abs(per_inc[c][k] - ref["per"][c][k]) for c in CTX for k in ("A", "B", "cap"))
    P(f"KA: h1union (A, B, cap) streamed here vs 243 invert243.json: max|d| {worst:.2e} {'PASS' if worst < 1e-9 else 'FAIL'}; nnz {int(inc['nnz'])}")
    assert worst < 1e-9
    cos_inc = ref["inv"]["cos"]
    raw_inc = ref["board"]["raw"]["mse"]
    chk = raw_for(per_inc, cos_inc)
    P(f"KA: 197's identity at h1union's own terms reproduces its board raw mse: {chk:.7f} vs {raw_inc:.7f} (|d| {abs(chk-raw_inc):.1e})")
    svd5 = stream("cd4svd5", SVD5)
    z249 = np.load(E249 / "work/totals-cd4svd5.npz", allow_pickle=True)
    same249 = np.array_equal(np.asarray(svd5["totals"]), np.asarray(z249["totals"])) and int(svd5["nnz"]) == int(z249["nnz"])
    P(f"KA: cd4svd5 per-block totals streamed here == 249's totals-cd4svd5.npz bit for bit: {same249}")
    assert same249
    per_svd5 = terms(np.asarray(svd5["totals"]), np.asarray(svd5["blocks"]).astype(str), genes, b)
    amp_svd5 = float(np.mean([np.sqrt(per_svd5[c]["A"] / per_inc[c]["A"]) for c in CTX]))
    P(f"KA: cd4svd5 delivered amplitude vs h1union {100*(amp_svd5-1):+.4f}% (249 section 9: -0.409%)")
    assert abs(100 * (amp_svd5 - 1) + 0.409) < 0.002

    # ---- the arm ------------------------------------------------------------------------------
    arm = stream("composed", ARM)
    assert np.array_equal(np.asarray(arm["genes"]).astype(str), genes)
    meta = json.loads((ROOT / "submissions/254-cd4svd5-partial-f034.h5ad.meta.json").read_text())
    assert int(arm["nnz"]) == meta["output"]["nnz"], "streamed nnz != meta nnz"
    per = terms(np.asarray(arm["totals"]), np.asarray(arm["blocks"]).astype(str), genes, b)
    P(f"\ncomposed vs cd4svd5 and h1union, delivered terms (211's convention):")
    for c in CTX:
        P(f"  {c}: A {per[c]['A']:.6f} (cd4svd5 {per_svd5[c]['A']:.6f} x{per[c]['A']/per_svd5[c]['A']:.5f}; h1union {per_inc[c]['A']:.6f} x{per[c]['A']/per_inc[c]['A']:.5f})  "
          f"B {per[c]['B']:.6f} (cd4svd5 {per_svd5[c]['B']:.6f})  cap {per[c]['cap']:.6f}  leg {'claim' if per[c]['cap'] < per[c]['B'] else 'B'}")
    amp_s = float(np.mean([np.sqrt(per[c]["A"] / per_svd5[c]["A"]) for c in CTX]))
    ampB_s = float(np.mean([np.sqrt(per[c]["B"] / per_svd5[c]["B"]) for c in CTX]))
    amp_h = float(np.mean([np.sqrt(per[c]["A"] / per_inc[c]["A"]) for c in CTX]))
    ampB_h = float(np.mean([np.sqrt(per[c]["B"] / per_inc[c]["B"]) for c in CTX]))
    P(f"  delivered amplitude sqrt(A_arm/A_cd4svd5) - 1, mean over contexts: {100*(amp_s-1):+.4f}%   (B-term {100*(ampB_s-1):+.4f}%)   "
      f"falsifier 5 bound +-0.5% vs cd4svd5: {'PASS' if abs(amp_s-1) <= 0.005 else 'FAIL'}")
    P(f"  delivered amplitude sqrt(A_arm/A_h1union) - 1: {100*(amp_h-1):+.4f}%   (B-term {100*(ampB_h-1):+.4f}%)")

    # ---- the predicted board raw mse at the arm's own terms -----------------------------------
    rel_stack = REPRO["composed_resolved"]["rel_agg"]      # composed vs h1union's stack, 3 cancer folds
    rel_svd5 = REPRO["svd5"]["rel_agg"]
    rel_incr = REPRO["increment_over_svd5"]["rel_agg"]
    P(f"\nPLOLO aggregate cosine reads (3 cancer folds): composed vs stack {100*rel_stack:+.4f}%, cd4svd5 vs stack {100*rel_svd5:+.4f}%, "
      f"increment {100*rel_incr:+.4f}%; h1union measured cos {cos_inc:.6f}")
    P(f"  {'k_cos':>6s}{'arm cos':>10s}{'pred raw mse ARM':>18s}{'d vs h1union':>14s}{'pred raw mse cd4svd5':>22s}{'ARM - cd4svd5':>15s}{'mse leg, overall':>18s}")
    preds = {}
    for k in K_COS:
        cos_arm = cos_inc * (1 + k * rel_stack)
        cos_s = cos_inc * (1 + k * rel_svd5)
        raw_arm = raw_for(per, cos_arm)
        raw_s = raw_for(per_svd5, cos_s)
        preds[str(k)] = dict(raw_arm=raw_arm, raw_svd5=raw_s, d_vs_h1union=raw_arm - raw_inc, d_vs_svd5=raw_arm - raw_s)
        P(f"  {k:>6.2f}{cos_arm:>10.6f}{raw_arm:>18.5f}{raw_arm-raw_inc:>+14.5f}{raw_s:>22.5f}{raw_arm-raw_s:>+15.5f}{-(raw_arm-raw_s)/0.950162566/6:>+18.5f}")
    P(f"  (k_cos 0 is the amplitude/budget leg alone; 249 predicted cd4svd5 at 0.8190-0.8194)")
    P(f"\ncaps: nnz {int(arm['nnz']):,} = {100*int(arm['nnz'])/4.75e9:.2f}% of 4.75e9; max counts per cell {int(arm['max_cell']):,} vs 1,000,000 "
      f"({1e6/int(arm['max_cell']):.1f}x under); cd4svd5 nnz {int(svd5['nnz']):,} ({100*(int(arm['nnz'])/int(svd5['nnz'])-1):+.3f}%)")
    out = dict(per_arm=per, per_svd5=per_svd5, per_inc=per_inc, amplitude_vs_svd5=amp_s, amplitude_pct_vs_svd5=100*(amp_s-1),
               B_ratio_vs_svd5=ampB_s, amplitude_vs_h1union=amp_h, amplitude_pct_vs_h1union=100*(amp_h-1),
               rel_agg=dict(composed_vs_stack=rel_stack, svd5_vs_stack=rel_svd5, increment=rel_incr), cos_inc=cos_inc, raw_inc=raw_inc,
               predicted_raw_mse=preds, nnz=int(arm["nnz"]), max_cell=int(arm["max_cell"]), nnz_svd5=int(svd5["nnz"]),
               falsifier5_pass=bool(abs(amp_s-1) <= 0.005), peak_rss_gb=rss_gb())
    (WORK / "amp254.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "amp254.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'amp254.json'}  peak rss {rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
