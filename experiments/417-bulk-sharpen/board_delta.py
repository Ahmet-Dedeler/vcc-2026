"""Print the six-member raw and scaled deltas between two captured board entries (data/board-entries/<id>.json).
usage: uv run python experiments/417-bulk-sharpen/board_delta.py <new_id> [<ref_id> (default 410b 4WgVi6FIE4iVw9PmRUP3)]"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}
SC = {"pds": "score_pds", "mse": "score_mse", "nmae": "score_nmae", "fid": "score_fid", "reach": "score_reach", "jac": "score_jac"}


def find(x):
    if isinstance(x, dict):
        if "score_avg" in x:
            return x
        for v in x.values():
            y = find(v)
            if y:
                return y
    if isinstance(x, list):
        for v in x:
            y = find(v)
            if y:
                return y


def load(i):
    return find(json.loads((ROOT / f"data/board-entries/{i}.json").read_text()))


a = load(sys.argv[1])
b = load(sys.argv[2] if len(sys.argv) > 2 else "4WgVi6FIE4iVw9PmRUP3")
print(f"overall {a['score_avg']:.7f} vs {b['score_avg']:.7f}  d {a['score_avg'] - b['score_avg']:+.7f}")
for k in RAW:
    ra, rb = a.get(RAW[k]), b.get(RAW[k])
    sa, sb = a.get(SC[k]), b.get(SC[k])
    ds = (sa - sb) if (sa is not None and sb is not None) else float("nan")
    print(f"  {k:5} raw {ra:.6f} vs {rb:.6f} (d {ra - rb:+.6f})   scaled d {ds:+.6f}  -> overall {ds / 6:+.6f}")
