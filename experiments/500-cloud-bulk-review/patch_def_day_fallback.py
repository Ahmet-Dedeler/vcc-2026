"""Replacement for the embedded python in scripts/def_day.sh (the not-converged fallback, ~line 148).

Today that block picks the safe build nearest rho = 0.9793, which is exp(2 * sd * Phi^{-1}(q)).
This block picks the safe build nearest rho*(e) = exp(e * sd * Phi^{-1}(q)), with e the log-log
slope of that context's builds. One build for a context keeps the exponent-2 anchor.

Drop-in: the script still receives the same "${POINTS[@]}" arguments
(`<ntclevel.json>:<D=f,E=f,F=f>`) and prints the chosen artifact path.
"""
import json
import math
import os
import sys
from collections import defaultdict

SD = 0.01
# Phi^{-1}(0.31 / (0.31 + 1.78)). Inlined so this block stays on system python3,
# matching def_day.sh's `python3` heredoc (no scipy).
Z = -1.0436427253969442


def rho_target(e: float) -> float:
    return math.exp(e * SD * Z)


by_ctx = defaultdict(list)  # context -> [(f, rho, path)]
for spec in sys.argv[1:]:
    path, fs = spec.rsplit(":", 1)
    per = json.load(open(path))["per_context"]
    h5 = path[: -len(".ntclevel.json")] + ".h5ad"
    if not os.path.exists(h5):
        continue
    if max(v["rho"] for v in per.values()) > 1.0:
        continue
    for item in fs.split(","):
        c, fv = item.split("=")
        by_ctx[c].append((float(fv), float(per[c]["rho"]), path))

# one target rho per context
target = {}
for c, rows in by_ctx.items():
    rows = sorted(set(rows))
    if len(rows) < 2:
        target[c] = rho_target(2.0)
        continue
    anchor = rho_target(2.0)
    lo = [r for r in rows if r[1] < anchor]
    hi = [r for r in rows if r[1] > anchor]
    if lo and hi:
        p0, p1 = max(lo, key=lambda r: r[1]), min(hi, key=lambda r: r[1])
    else:
        p0, p1 = sorted(rows, key=lambda r: abs(math.log(r[1] / anchor)))[:2]
    e = math.log(p1[1] / p0[1]) / math.log(p1[0] / p0[0])
    target[c] = rho_target(e)

best = None
seen = {}
for spec in sys.argv[1:]:
    path = spec.rsplit(":", 1)[0]
    h5 = path[: -len(".ntclevel.json")] + ".h5ad"
    if not os.path.exists(h5):
        continue
    per = json.load(open(path))["per_context"]
    if max(v["rho"] for v in per.values()) > 1.0:
        continue
    gap = max(abs(math.log(per[c]["rho"] / target[c])) for c in per)
    if best is None or gap < best[0]:
        best = (gap, path[: -len(".ntclevel.json")])
        seen = {c: per[c]["rho"] for c in per}
print(best[1] if best else "")
if best:
    sys.stderr.write(
        "fallback gap {g:.5f} targets {t} rhos {r}\n".format(
            g=best[0],
            t={c: round(v, 4) for c, v in target.items()},
            r={c: round(v, 4) for c, v in seen.items()},
        )
    )
