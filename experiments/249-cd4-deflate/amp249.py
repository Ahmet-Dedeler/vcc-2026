"""249 falsifier 5 -- delivered amplitude off the PACKAGE, and the predicted board raw `mse`
from 197's identity at the arm's OWN terms (243 section 5.2: never the floor).

Streams the arm's per-block panel totals (214's `block_totals`, imported from 243's
`totals243.py`), computes 211's per-context (A_c, B_c, Kc_c, cap_c) with 243's `terms`, and
reads the delivered amplitude against `h1union`'s terms in 243's `invert243.json`.

KNOWN-ANSWER ROW: the incumbent `submissions/222-h1union-g028.h5ad` is streamed the same way and
its (A_c, B_c, cap_c) must equal 243's to 1e-9 on all three contexts.

Predicted raw `mse` for the arm: raw = 1 + (mean_c[A_c + Kc_c - min(cap_c, B_c)] - 2 X)/sumD with
X = cos * sqrt(sumD) * mean_c sqrt(A_c) and cos = cos_h1union * (1 + k_cos * rel_agg), rel_agg the
three-cancer-fold PLOLO aggregate read (repro249), at k_cos in {0 (amplitude leg alone), 0.94,
1.00, 1.36} (243's classes).  Also the count caps for falsifier 6 (nnz, max counts per cell).

    uv run python amp249.py
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
sys.path.insert(0, str(E243))
from totals243 import block_totals  # noqa: E402
from invert243 import bctrl, terms  # noqa: E402
from lib243 import CTX, SUMD, rss_gb  # noqa: E402

ARM = ROOT / "submissions/249-cd4svd5-f034.h5ad"
INC = ROOT / "submissions/222-h1union-g028.h5ad"
INV = json.loads((E243 / "work/invert243.json").read_text())
REPRO = json.loads((WORK / "repro249-summary.json").read_text())
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

    # ---- known-answer row: the incumbent's terms reproduce 243's ------------------------
    inc = stream("h1union", INC)
    assert np.array_equal(np.asarray(inc["genes"]).astype(str), genes)
    per_inc = terms(np.asarray(inc["totals"]), np.asarray(inc["blocks"]).astype(str), genes, b)
    ref = INV["arms"]["h1union"]
    worst = max(abs(per_inc[c][k] - ref["per"][c][k]) for c in CTX for k in ("A", "B", "cap"))
    P(f"KA: h1union (A, B, cap) streamed here vs 243 invert243.json: max|d| {worst:.2e} "
      f"{'PASS' if worst < 1e-9 else 'FAIL'}; nnz {int(inc['nnz'])}")
    assert worst < 1e-9
    m222 = json.loads((ROOT / "submissions/222-h1union-g028.h5ad.meta.json").read_text())
    assert int(inc["nnz"]) == m222["output"]["nnz"]

    # ---- the arm ---------------------------------------------------------------------------
    arm = stream("cd4svd5", ARM)
    assert np.array_equal(np.asarray(arm["genes"]).astype(str), genes)
    meta = json.loads((ROOT / "submissions/249-cd4svd5-f034.h5ad.meta.json").read_text())
    assert int(arm["nnz"]) == meta["output"]["nnz"], "streamed nnz != meta nnz"
    per = terms(np.asarray(arm["totals"]), np.asarray(arm["blocks"]).astype(str), genes, b)
    P(f"\ncd4svd5 vs h1union, delivered terms (211's convention):")
    for c in CTX:
        P(f"  {c}: A {per[c]['A']:.6f} (h1union {per_inc[c]['A']:.6f}, x{per[c]['A']/per_inc[c]['A']:.5f})  "
          f"B {per[c]['B']:.6f} (h1union {per_inc[c]['B']:.6f})  cap {per[c]['cap']:.6f}  "
          f"leg {'claim' if per[c]['cap'] < per[c]['B'] else 'B'}")
    amp = float(np.mean([np.sqrt(per[c]["A"] / per_inc[c]["A"]) for c in CTX]))
    ampB = float(np.mean([np.sqrt(per[c]["B"] / per_inc[c]["B"]) for c in CTX]))
    P(f"  delivered amplitude sqrt(A_arm/A_inc) - 1, mean over contexts: {100*(amp-1):+.4f}%   "
      f"(B-term {100*(ampB-1):+.4f}%)   falsifier 5 bound +-0.5%: {'PASS' if abs(amp-1) <= 0.005 else 'FAIL'}")

    # ---- the predicted board raw mse at the arm's own terms ---------------------------------
    cos_inc = ref["inv"]["cos"]
    raw_inc = ref["board"]["raw"]["mse"]
    rel = REPRO["arm"]["rel_agg"]
    chk = raw_for(per_inc, cos_inc)
    P(f"\n197's identity at h1union's own terms reproduces its board raw mse: {chk:.7f} vs {raw_inc:.7f} "
      f"(|d| {abs(chk-raw_inc):.1e})")
    P(f"PLOLO aggregate cosine read (3 cancer folds, repro249): {100*rel:+.4f}%; h1union measured cos {cos_inc:.6f}")
    P(f"  {'k_cos':>6s}{'arm cos':>10s}{'d cos':>10s}{'pred raw mse':>14s}{'d raw':>10s}{'d overall (mse leg)':>22s}")
    preds = {}
    for k in K_COS:
        cos = cos_inc * (1 + k * rel)
        raw = raw_for(per, cos)
        preds[str(k)] = raw
        P(f"  {k:>6.2f}{cos:>10.6f}{cos-cos_inc:>+10.6f}{raw:>14.5f}{raw-raw_inc:>+10.5f}{-(raw-raw_inc)/0.950162566/6:>+22.5f}")
    P(f"  (k_cos 0 is the amplitude/budget leg alone; the prereg's 'raw mse ~ 0.8183' is k_cos = 1 with h1union's "
      f"terms; falsifier 2 fires above {raw_inc + 0.0030:.5f})")

    # ---- the count caps (falsifier 6) -------------------------------------------------------
    P(f"\ncaps: nnz {int(arm['nnz']):,} = {100*int(arm['nnz'])/4.75e9:.2f}% of 4.75e9; max counts per cell "
      f"{int(arm['max_cell']):,} vs 1,000,000 ({1e6/int(arm['max_cell']):.1f}x under); h1union nnz {int(inc['nnz']):,} "
      f"({100*(int(arm['nnz'])/int(inc['nnz'])-1):+.3f}%)")
    out = dict(per_arm=per, per_inc=per_inc, amplitude_ratio=amp, amplitude_pct=100*(amp-1), B_ratio=ampB,
               rel_agg_plolo=rel, cos_inc=cos_inc, raw_inc=raw_inc, predicted_raw_mse=preds,
               nnz=int(arm["nnz"]), max_cell=int(arm["max_cell"]), nnz_inc=int(inc["nnz"]),
               falsifier5_pass=bool(abs(amp-1) <= 0.005), peak_rss_gb=rss_gb())
    (WORK / "amp249.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "amp249.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'amp249.json'}  peak rss {rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
