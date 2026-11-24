"""The 10-05 rule, FIXED 2026-09-27 (before any 09-28 entry exists). Applied by chain1005.sh at 10-04 20:00Z.
  s1 = 496a (446bj2 + the uncalled-gene sign hook, VCC_DE_MDE_SIGN=clean4:0.05; DE-only, pds/mse must be bit-identical to
       446bj2): the hook's board read. Bins in PREREGISTRATION-board.md.
  s2 = the first READY, not-yet-scored of: 497dr (only if the cap and the per-target form both held: C != FAILS and P = S_P2,
       as in rule1004), then 498sp2, 498p2g (each only if C != FAILS and its own 479 bin was at least PARTIAL), then
       469dr, 461dr, 460dr, 458e.
The P / C bins are rule1004's (imported). READY and "not yet scored" (novelty) are checked by the chain.
usage: rule1005.py <453dr id|none> <479p2g id|none> <479sp2 id|none> <495 id|none>  ->  "<s1 list>|<s2 list>|P=..,C=..,O=.."
"""
import importlib.util
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("r1004", ROOT / "experiments/498-combined-bulk/rule1004.py")
r = importlib.util.module_from_spec(spec)
argv, sys.argv = sys.argv, [sys.argv[0]] + (sys.argv[1:5] + ["none"] * 4)[:4]
buf = io.StringIO()
with redirect_stdout(buf):
    spec.loader.exec_module(r)          # rule1004 computes P, C, O and the bins at import time
sys.argv = argv
P, C, O, bp, bs = r.P, r.C, r.O, r.bp, r.bs
s2 = []
if C != "FAILS" and P == "S_P2" and O != "BAD":
    s2.append("497dr")
if C != "FAILS":
    for tag, bb in (("498sp2", bs), ("498p2g", bp)):
        if bb is not None and bb[0] in ("TRANSFERS", "PARTIAL"):
            s2.append(tag)
s2 += ["460dr"]   # 2026-09-30 06:38Z amendment (before any 10-01..10-05 entry): the Arc-programme arms 469dr/461dr/458e
                   # read HARM on 09-30 (458d); they no longer fill a slot
print(f"496a {' '.join(s2)}|{' '.join(s2)}|P={P},C={C},O={O}")
