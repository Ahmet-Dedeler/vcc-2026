"""448: the day-of re-level loop, one secant step per build, per context.

The render maps the builder factor f to the across-perturbation spread nonlinearly: rho(f) ~ f^e. On A/B/C e ~ 1.86
(compressive), but 448's efr rehearsal measured e = 1.89 / 1.90 / 2.30 for D / E / F between f = 1 and the first levelled
factors, and F's LOCAL exponent near f 1.24-1.29 was 2.71 (expansive). One pass is not enough, one global exponent is not
enough, and even a per-context exponent from far away (f = 1 -> f1) misses F by 3 %. So step like a secant: take each
context's two measured builds nearest the target (a bracketing pair when there is one), read the local exponent in log-log
space, and solve for the EV target rho* = exp(2 * sd * Phi^-1(q)), q = 0.31 / (0.31 + 1.78).

usage: uv run python relevel_iter.py <build.ntclevel.json>:<D=f,E=f,F=f> [<...>:<...> ...] [--sd 0.01] [--tol 0.01]
  Pass every levelled build so far, the un-levelled one as "<lvl0.ntclevel.json>:D=1,E=1,F=1". Order does not matter.
  With the un-levelled build alone, a context's next factor is ntclevel's own EV scale (the first step, 448's lvl1).
  Prints, per context, the points used, the local exponent and the next factor, then the VCC_BULK_CTX_SCALE string for the
  next build, or CONVERGED (every context within --tol of rho* in rho, which is tol/2 in scale) with the factors to ship.
  Never ship a build with any context above its cliff (rho > 1); the loop approaches from below when it can.
"""
import argparse
import json
import math

from scipy.stats import norm

ap = argparse.ArgumentParser()
ap.add_argument("points", nargs="+", help="<ntclevel.json>:<D=f,E=f,F=f>")
ap.add_argument("--sd", type=float, default=0.01)
ap.add_argument("--tol", type=float, default=0.01, help="relative tolerance on rho (0.01 = 0.5 %% in scale)")
a = ap.parse_args()
q = 0.31 / (0.31 + 1.78)
rho_star = math.exp(2.0 * a.sd * norm.ppf(q))

pts = {}    # context -> list of (f, rho)
scale = {}  # context -> {f: ntclevel's EV scale at that build}, for the first step (one build only)
for spec in a.points:
    path, fs = spec.rsplit(":", 1)
    per = json.load(open(path))["per_context"]
    for c, v in (p.split("=") for p in fs.split(",")):
        pts.setdefault(c, []).append((float(v), float(per[c]["rho"])))
        scale.setdefault(c, {})[float(v)] = float(per[c]["scale"])

print(f"EV target rho* = {rho_star:.4f} (sd {a.sd}, q {q:.3f}); tolerance |rho/rho* - 1| <= {a.tol}")
nxt, done = {}, True
for c in sorted(pts):
    P = sorted(set(pts[c]))
    best = min(P, key=lambda p: abs(math.log(p[1] / rho_star)))
    if abs(best[1] / rho_star - 1.0) <= a.tol and best[1] <= 1.0:
        nxt[c] = best[0]
        print(f"  {c}: converged at f {best[0]:.4f} (rho {best[1]:.4f})")
        continue
    done = False
    if len(P) < 2:
        # the first step (only the un-levelled build so far): ntclevel's own EV scale, relative to that build (448's lvl1)
        f = P[0][0] * scale[c][P[0][0]]
        nxt[c] = f
        print(f"  {c}: one build ({P[0][0]:.4f}, {P[0][1]:.4f}); ntclevel EV scale {scale[c][P[0][0]]:.4f}  next f {f:.4f}"
              + ("   WARNING outside the builder bound [0.7, 1.4]" if not 0.7 <= f <= 1.4 else ""))
        continue
    lo = [p for p in P if p[1] < rho_star]
    hi = [p for p in P if p[1] > rho_star]
    if lo and hi:  # bracketing pair nearest the target on each side
        p0, p1 = max(lo, key=lambda p: p[1]), min(hi, key=lambda p: p[1])
        how = "bracket"
    else:  # the two nearest the target
        p0, p1 = sorted(P, key=lambda p: abs(math.log(p[1] / rho_star)))[:2]
        how = "extrapolate"
    e = math.log(p1[1] / p0[1]) / math.log(p1[0] / p0[0])
    f = p1[0] * (rho_star / p1[1]) ** (1.0 / e)
    nxt[c] = f
    print(f"  {c}: {how} ({p0[0]:.4f}, {p0[1]:.4f}) - ({p1[0]:.4f}, {p1[1]:.4f})  local exponent {e:.3f}  next f {f:.4f}"
          + ("   WARNING outside the builder bound [0.7, 1.4]" if not 0.7 <= f <= 1.4 else ""))
s = ",".join(f"{c}={v:.4f}" for c, v in nxt.items())
print(("CONVERGED: ship " if done else "next build: ") + "VCC_BULK_CTX_SCALE=" + s)
