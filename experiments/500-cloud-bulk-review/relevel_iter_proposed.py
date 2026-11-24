"""Proposed replacement for experiments/448-def-rehearsal-q175/relevel_iter.py.

The only behaviour change: once a context has two builds, the EV target is
rho*(e) = exp(e * sd * Phi^{-1}(q)) with the local exponent e just measured.
exp(2 * sd * Phi^{-1}(q)) remains the target only while e is unknown (one build),
which is ntclevel's sqrt cliff. See RESULTS.md finding 1.

Not applied to the production script. Run it the same way as relevel_iter.py.
"""
import argparse
import json
import math

from scipy.stats import norm

ap = argparse.ArgumentParser()
ap.add_argument("points", nargs="+", help="<ntclevel.json>:<D=f,E=f,F=f>")
ap.add_argument("--sd", type=float, default=0.01)
ap.add_argument("--tol", type=float, default=0.01, help="relative tolerance on rho (0.01 = 0.5% in scale at e=2)")
a = ap.parse_args()
q = 0.31 / (0.31 + 1.78)
z = float(norm.ppf(q))
rho_e2 = math.exp(2.0 * a.sd * z)


def rho_target(e: float) -> float:
    return math.exp(float(e) * a.sd * z)


pts = {}
scale = {}
for spec in a.points:
    path, fs = spec.rsplit(":", 1)
    per = json.load(open(path))["per_context"]
    for c, v in (p.split("=") for p in fs.split(",")):
        pts.setdefault(c, []).append((float(v), float(per[c]["rho"])))
        scale.setdefault(c, {})[float(v)] = float(per[c]["scale"])

print(f"EV hedge z={z:.4f} (sd {a.sd}, q {q:.3f}); e=2 rho*={rho_e2:.4f}; tolerance |rho/rho*(e) - 1| <= {a.tol}")
nxt, done = {}, True
for c in sorted(pts):
    P = sorted(set(pts[c]))
    if len(P) < 2:
        if abs(P[0][1] / rho_e2 - 1.0) <= a.tol and P[0][1] <= 1.0:
            nxt[c] = P[0][0]
            print(f"  {c}: converged at f {P[0][0]:.4f} (rho {P[0][1]:.4f}, e unknown, rho*(2) {rho_e2:.4f})")
            continue
        done = False
        f = P[0][0] * scale[c][P[0][0]]
        nxt[c] = f
        print(f"  {c}: one build ({P[0][0]:.4f}, {P[0][1]:.4f}); ntclevel EV scale {scale[c][P[0][0]]:.4f}  next f {f:.4f}"
              + ("   WARNING outside the builder bound [0.7, 1.4]" if not 0.7 <= f <= 1.4 else ""))
        continue
    lo = [p for p in P if p[1] < rho_e2]
    hi = [p for p in P if p[1] > rho_e2]
    if lo and hi:
        p0, p1 = max(lo, key=lambda p: p[1]), min(hi, key=lambda p: p[1])
        how = "bracket"
    else:
        p0, p1 = sorted(P, key=lambda p: abs(math.log(p[1] / rho_e2)))[:2]
        how = "extrapolate"
    e = math.log(p1[1] / p0[1]) / math.log(p1[0] / p0[0])
    target = rho_target(e)
    best = min(P, key=lambda p: abs(math.log(p[1] / target)))
    if abs(best[1] / target - 1.0) <= a.tol and best[1] <= 1.0:
        nxt[c] = best[0]
        print(f"  {c}: converged at f {best[0]:.4f} (rho {best[1]:.4f}, e {e:.3f}, rho*(e) {target:.4f})")
        continue
    done = False
    f = p1[0] * (target / p1[1]) ** (1.0 / e)
    nxt[c] = f
    print(f"  {c}: {how} ({p0[0]:.4f}, {p0[1]:.4f}) - ({p1[0]:.4f}, {p1[1]:.4f})  local exponent {e:.3f}  "
          f"rho*(e) {target:.4f}  next f {f:.4f}"
          + ("   WARNING outside the builder bound [0.7, 1.4]" if not 0.7 <= f <= 1.4 else ""))
s = ",".join(f"{c}={v:.4f}" for c, v in nxt.items())
print(("CONVERGED: ship " if done else "next build: ") + "VCC_BULK_CTX_SCALE=" + s)
