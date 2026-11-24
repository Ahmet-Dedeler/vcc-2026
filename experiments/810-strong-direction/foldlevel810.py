"""810 -- fold-level, direction-only composition read of the candidate constructions (relative aggregate cosine r, %, at the
incumbent's realised norms), gated by tier: 707's U blend (b 0.25) vs the base -> pre swap (b 0.5 / 1). Light, local.
Calibrates the composition read against 707's engine-realised r on A/B/C (decomp810 realise KA-706 record)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyse810 as A  # noqa: E402

def fold_r(c, pt, cos_new, mask):
    live = pt["rd"] > 0
    w = np.sqrt(pt["ra"] * pt["rd"])
    base = (w[live] * pt["cos_fin"][live]).sum()
    new = pt["cos_fin"].copy()
    new[mask] = np.nan_to_num(cos_new[mask])
    return 100 * ((w[live] * new[live]).sum() / base - 1)

def main():
    res, PT = A.main()
    out = {}
    for st, (c, pt) in PT.items():
        strong = c["tier"] == 2
        ag = pt.get("cos_src01", np.full(len(strong), np.nan))
        hi_ag = strong & (ag > np.nanmedian(ag[strong]))
        r = {"U0.25 strong": fold_r(c, pt, pt["blend_U_b0.25"], strong),
             "pre0.5 strong": fold_r(c, pt, pt["swapbase_pre_b0.5"], strong),
             "pre1 strong": fold_r(c, pt, pt["swapbase_pre_b1"], strong),
             "pre1 strong&agree": fold_r(c, pt, pt["swapbase_pre_b1"], hi_ag),
             "pre1 all": fold_r(c, pt, pt["swapbase_pre_b1"], np.ones_like(strong)),
             "pre0.5 all": fold_r(c, pt, pt["swapbase_pre_b0.5"], np.ones_like(strong))}
        out[st] = r
    for grp, sel in (("W", lambda s: s.startswith("L")), ("ABC cancer", lambda s: s.startswith("ABC") and "kolf2" not in s), ("kolf2", lambda s: "kolf2" in s)):
        ks = list(next(iter(out.values())).keys())
        v = {k: np.array([out[s][k] for s in out if sel(s)]) for k in ks}
        print(grp, " | ".join(f"{k} {v[k].mean():+.2f}% ({(v[k] > 0).sum()}/{v[k].size})" for k in ks))
    for s in out:
        if s.startswith("ABC"):
            print(s, {k: round(x, 3) for k, x in out[s].items()})
    (HERE / "work" / "foldlevel810.json").write_text(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
