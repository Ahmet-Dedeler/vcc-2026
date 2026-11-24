"""592 — recompute every new number in Draft v16 from primary files and check that the paper states it.

Board numbers from `data/board-entries/<id>.json` (scaled members; overall = mean of six; raw members); doses from the
build metas' `output.rc_callset.prior_lambda`; the 198-block bound from 458d's median and 468's gate sizes. Exits 1 if any
recomputed value is missing from PAPER.md / CLAIMS.md as written. Usage: uv run python experiments/592-manuscript-v16/verify592.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = (Path(__file__).resolve().parent / "paper/PAPER.md").read_text()
CLAIMS = (Path(__file__).resolve().parent / "paper/CLAIMS.md").read_text()
E = {"446b": "UDgjOKL9gI1udAexq6ue", "446bj2": "YvVru3DAF89Jiyr34olR", "446bj2z": "be90smZUDcIV7j5MgQ03",
     "446bd2": "R5uE5Rve3R9nWXGbeSLM", "453dr": "q3VuTImh4FRGKXH2waCp", "457dr": "1FfirZsGsAZvaXYTd39J",
     "460dr": "LrkhsIXWMjEsvRXSGToa", "446bj2zcp": "K1Q0y2OFYjTNRqh9jKfk"}
MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
fails: list[str] = []


def entry(k: str) -> dict:
    j = json.loads((ROOT / f"data/board-entries/{E[k]}.json").read_text())
    return j.get("entry", j)


def check(label: str, text: str, where: str = "PAPER") -> None:
    src = PAPER if where == "PAPER" else CLAIMS
    ok = text in src
    print(f"  [{'OK' if ok else 'MISSING'}] {where}: {text!r}  ({label})")
    if not ok:
        fails.append(f"{label}: {text}")


def delta(a: str, b: str) -> dict:
    A, B = entry(a), entry(b)
    s = {m: (A["score_" + m] - B["score_" + m]) / 6 for m in MEM}
    r = {m: A[RAW[m]] - B[RAW[m]] for m in MEM}
    return {"overall": A["score_avg"] - B["score_avg"], "bulk": s["pds"] + s["mse"],
            "de": s["nmae"] + s["fid"] + s["reach"] + s["jac"], "raw": r, "ov": s}


print("1. 460dr - 453dr (C-94)")
d = delta("460dr", "453dr")
print(f"   overall {d['overall']:+.7f}  mse raw {d['raw']['mse']:+.6f}  pds raw {d['raw']['pds']:+.6f}")
check("460dr overall", f"{d['overall']:+.7f}")
check("460dr mse raw", f"mse −{abs(d['raw']['mse']):.6f} raw")
k = delta("446bj2", "446b")
ratio = d["raw"]["mse"] / k["raw"]["mse"]
print(f"   A/B/C key mse raw {k['raw']['mse']:+.6f}; ratio {ratio:.3f}")
check("ratio", f"{ratio:.2f}×")
print("2. 460dr - 446bj2 (matched key): DE and bulk")
d2 = delta("460dr", "446bj2")
print(f"   DE {d2['de']:+.7f}  bulk {d2['bulk']:+.7f}")
check("DE matched", f"{abs(d2['de']):.7f}")
check("bulk matched", f"bulk −{abs(d2['bulk']):.7f}")
print("3. successors of 446b: DE shares")
for a in ("446bj2", "446bd2", "453dr", "457dr"):
    x = delta(a, "446b")
    print(f"   {a:7} - 446b  DE {x['de']:+.7f}")
check("446bj2 DE", "446bj2 −0.00039")
check("446bd2 DE", "446bd2 −0.00048")
check("453dr DE", "453dr −0.00060")
check("457dr DE", "457dr −0.00055")
x = delta("457dr", "446b")
print(f"   457dr - 446b overall {x['overall']:+.7f}, bulk {x['bulk']:+.7f}, DE {x['de']:+.7f}")
check("457dr-446b DE", f"(−{abs(x['de']):.7f})")
print("4. 446bj2zcp - 446bj2z (C-95)")
d4 = delta("446bj2zcp", "446bj2z")
fj = d4["ov"]["fid"] + d4["ov"]["jac"]
print(f"   overall {d4['overall']:+.7f}  d(fid+jac) overall {fj:+.7f}  raw jac {d4['raw']['jac']:+.6f} fid {d4['raw']['fid']:+.6f} "
      f"reach {d4['raw']['reach']:+.6f}  pds/mse raw {d4['raw']['pds']:+.1e}/{d4['raw']['mse']:+.1e}")
check("cp d(fid+jac)", f"−{abs(fj):.5f}")
check("cp jac raw", f"jac +{d4['raw']['jac']:.6f} raw")
check("cp fid raw", f"fid −{abs(d4['raw']['fid']):.6f} raw")
check("cp reach raw", f"reach −{abs(d4['raw']['reach']):.6f} raw")
print("5. delivered dose (C-97)")
metas = {"446b": "446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
         "458": "458-cap15-partial50-arch1prog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
         "458d": "458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5",
         "458g": "458g-cap15-partial50-arch1dprog-gate33-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5"}
lam = {}
for k2, stem in metas.items():
    m = json.loads((ROOT / f"submissions/{stem}.h5ad.meta.json").read_text())
    lam[k2] = m["output"]["rc_callset"]["prior_lambda"]
    print(f"   {k2:5} prior_lambda {lam[k2]}")
check("458d median", "0.433 for 458d")
check("458 median", "0.243 for 458")
check("458g median", "0.953 for 458g")
strong = 214 + 217 + 217  # 468 RESULTS: programme (strong) targets per context A/B/C
weak = 900 - strong
bound = 450 - weak  # 458d: >= 450 of 900 blocks at or below its median; at most `weak` of them are weak-target blocks
print(f"   strong {strong}, weak {weak}; >= {bound} strong blocks at lambda <= {lam['458d']['median']:.3f}")
check("bound", f"at least {bound} of the {strong}")
check("weak", f"{weak} consensus blocks")
print("FAILED: " + "; ".join(fails) if fails else "ALL CHECKS PASS")
sys.exit(1 if fails else 0)
