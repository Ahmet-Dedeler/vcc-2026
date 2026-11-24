"""207 Phase 4 falsifier 2 — is what CD4 supplies LINEAGE-SPECIFIC?

The decisive test the mandate asks for (larger gain on context A/Jurkat than on B/C) CANNOT be
run: the CD4 banks on disk are panel-restricted (283/287/284 perturbations, all inside the 2026
panel) and `data/donors/nadig_jurkat` -- the only T-lineage line we hold -- covers **0 of the
300 panel targets**, because the panel excludes common essentials by construction (166/206).
There is no surface on which a T-lineage RECIPIENT and the CD4 donor share a target.

What CAN be measured, and is: how close CD4's own effect is to each of the four held recipient
lines' truths, and how close our three-donor consensus is.  If CD4 sits at a similar distance
from every lineage, the direction it supplies is lineage-GENERIC.

    uv run python lineage207.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
from plolo import FOLDS, build_fold, cd4_for  # noqa: E402
from sweep207 import shuffled_cd4  # noqa: E402

LINEAGE = {"k562_gwps": "CML / erythroleukemia (myeloid)",
           "orion_hct116": "colorectal carcinoma (epithelial)",
           "orion_hek293t": "embryonic kidney (epithelial)",
           "kolf2": "iPSC (pluripotent)"}


def agg_cos(P, T):
    x = float((P * T).sum()); a = float((P * P).sum()); d = float((T * T).sum())
    return x / np.sqrt(a * d) if a > 0 and d > 0 else 0.0


def macro_cos(P, T):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    return float(np.mean(np.sum(A * B, axis=1)))


def main():
    out = {}
    print(f"{'held line':<16}{'lineage':<34}{'cos(CD4,truth)':>16}{'cos(pre,truth)':>16}"
          f"{'cos(shufCD4,truth)':>20}{'CD4/pre':>9}")
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        sc = fold["scored"]
        T = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        pre = np.asarray(fold["pre"], dtype=np.float64)[:, sc]
        cd4 = cd4_for(fold["targets"])
        carry = np.asarray(cd4.target_condition_count) == 3
        C = np.asarray(cd4.mean_effect, dtype=np.float64)[:, sc]
        S = np.asarray(shuffled_cd4(cd4).mean_effect, dtype=np.float64)[:, sc]
        r = dict(n_carry=int(carry.sum()),
                 cd4_agg=agg_cos(C[carry], T[carry]), pre_agg=agg_cos(pre[carry], T[carry]),
                 shuf_agg=agg_cos(S[carry], T[carry]),
                 cd4_macro=macro_cos(C[carry], T[carry]),
                 pre_macro=macro_cos(pre[carry], T[carry]),
                 cd4_pre=agg_cos(C[carry], pre[carry]))
        out[held] = r
        print(f"{held:<16}{LINEAGE[held]:<34}{r['cd4_agg']:>16.6f}{r['pre_agg']:>16.6f}"
              f"{r['shuf_agg']:>20.6f}{r['cd4_agg']/r['pre_agg']:>9.3f}")
    v = [out[f]["cd4_agg"] / out[f]["pre_agg"] for f in FOLDS]
    print(f"\nCD4/consensus alignment ratio across four unrelated non-T lineages: "
          f"{min(v):.3f} - {max(v):.3f}  (cv {np.std(v)/np.mean(v):.3f})")
    (WORK / "lineage207.json").write_text(json.dumps(out, indent=2))
    print(f"wrote {WORK/'lineage207.json'}")


if __name__ == "__main__":
    main()
