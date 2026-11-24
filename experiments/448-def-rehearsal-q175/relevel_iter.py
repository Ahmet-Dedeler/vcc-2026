"""448: the day-of re-level loop, one secant step per build, per context.

The render maps the builder factor f to the across-perturbation spread nonlinearly: rho(f) ~ f^e. On A/B/C e ~ 1.86
(compressive), but 448's efr rehearsal measured e = 1.89 / 1.90 / 2.30 for D / E / F between f = 1 and the first levelled
factors, and F's LOCAL exponent near f 1.24-1.29 was 2.71 (expansive). One pass is not enough, one global exponent is not
enough, and even a per-context exponent from far away (f = 1 -> f1) misses F by 3 %. So step like a secant: take each
context's two measured builds nearest the target (a bracketing pair when there is one), read the local exponent in log-log
space, and solve for the EV target rho* = exp(2 * sd * Phi^-1(q)), q = 0.31 / (0.31 + 1.78).

usage: uv run python relevel_iter.py <build.ntclevel.json>:<D=f,E=f,F=f> [<...>:<...> ...] [--sd 0.01] [--tol 0.01] [--rho-scale 1]
  --rho-scale s (530, default 1, must be in [0.9, 1.1]) targets rho* x s. Omitted, the target is rho*
  and the arithmetic is the unscaled formula (no multiply).
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
ap.add_argument("--rho-scale", type=float, default=1.0, help="530: target rho* x this (default 1)")
a = ap.parse_args()
if not 0.9 <= a.rho_scale <= 1.1:
    raise SystemExit(f"--rho-scale must be in [0.9, 1.1] (got {a.rho_scale})")
q = 0.31 / (0.31 + 1.78)
base = math.exp(2.0 * a.sd * norm.ppf(q))
rho_star = base if a.rho_scale == 1.0 else base * a.rho_scale

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
    # 485: a flat or repeated-f pair makes e ~ 0 or a ZeroDivision. Do not print a VCC line;
    # exit 2 so def_day.sh falls through to the nearest safe build. |e| < 0.05 is far below
    # the 1.89–2.73 exponents measured on the 452 rehearsal, so the happy path is unchanged.
    # Finite factors outside [0.7, 1.4] are still printed (the shell's inbound check stops them).
    try:
        e = math.log(p1[1] / p0[1]) / math.log(p1[0] / p0[0])
        f = p1[0] * (rho_star / p1[1]) ** (1.0 / e)
    except (ZeroDivisionError, ValueError, OverflowError):
        e, f = float("nan"), float("nan")
    if not (math.isfinite(e) and abs(e) >= 0.05 and math.isfinite(f) and f > 0):
        print(f"  {c}: NO_STEP {how} ({p0[0]:.4f}, {p0[1]:.4f}) - ({p1[0]:.4f}, {p1[1]:.4f})  "
              f"local exponent {e!r} is unusable")
        print("NO_STEP")
        raise SystemExit(2)
    nxt[c] = f
    print(f"  {c}: {how} ({p0[0]:.4f}, {p0[1]:.4f}) - ({p1[0]:.4f}, {p1[1]:.4f})  local exponent {e:.3f}  next f {f:.4f}"
          + ("   WARNING outside the builder bound [0.7, 1.4]" if not 0.7 <= f <= 1.4 else ""))
s = ",".join(f"{c}={v:.4f}" for c, v in nxt.items())
print(("CONVERGED: ship " if done else "next build: ") + "VCC_BULK_CTX_SCALE=" + s)
