"""Route (b), part 2: extend the functional form, combine keys, and PUT AN ERROR BAR ON IT.

perpert.py found `donor_cos` -- the plain cosine between the two donors' effect vectors for a
perturbation, which is NOT in 096's candidate set -- transferring at 1.206x where 096's best
per-perturbation key transferred at 1.035x.  Before believing that:

  1. gamma=3 was the EDGE of the grid, so the form is untuned.
  2. n = 85 / 132, and the multiplier is a max over a grid, so the selection alone inflates it.
     A PERMUTATION null on the score -- same pipeline, same grid, same argmax, score shuffled --
     prices exactly that.  This is the check 096's key never got.
  3. a bootstrap over perturbations gives the sampling error on the multiplier itself.
"""
from __future__ import annotations
import json, sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
D = json.load(open(f"{ROOT}/experiments/097-direction/work/perpert.json"))
RNG = np.random.default_rng(0)
GAM = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]


def W_of(z, x, a, w):
    n = float((w * z * x).sum()); d = float((w * z * z * a).sum())
    return (n * n / d) if d > 0 else 0.0


def rank01(s):
    return (np.argsort(np.argsort(s)).astype(float) + 0.5) / len(s)


def curve(sc, x, a, w, base, gams=GAM):
    rk = rank01(sc)
    return np.array([np.sqrt(W_of(rk ** g, x, a, w) / base) for g in gams])


def topk_curve(sc, x, a, w, base, fracs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9)):
    rk = rank01(sc)
    return np.array([np.sqrt(W_of((rk >= 1 - f).astype(float), x, a, w) / base) for f in fracs])


H = {}
for name in ("h1-test", "orion"):
    d = D[name]
    x, a, w = (np.array(d[k]) for k in ("x", "a", "w"))
    base, S = d["base"], {k: np.array(v) for k, v in d["scores"].items()}
    # combinations of the two donor keys that survive partialling on both surfaces
    S["cos+096"] = rank01(S["donor_cos"]) + rank01(S["donor096"])
    S["cos*096"] = rank01(S["donor_cos"]) * rank01(S["donor096"])
    S["cos+self"] = rank01(S["donor_cos"]) + rank01(S["self_norm"])
    S["cos+ndon"] = rank01(S["donor_cos"]) + rank01(S["n_donors"])
    H[name] = dict(x=x, a=a, w=w, base=base, S=S)

print("### soft power form, sqrt(W/W_base) vs gamma "
      "(h1-test row 1, orion row 2 for each score)")
hdr = "".join(f"{g:>7.2f}" for g in GAM)
print(f"  {'score':<14}{'harness':<10}{hdr}")
KEYS = ["donor_cos", "donor096", "donor_signfrac", "self_norm", "n_donors",
        "cos+096", "cos*096", "cos+self", "cos+ndon"]
for k in KEYS:
    for name in ("h1-test", "orion"):
        h = H[name]
        c = curve(h["S"][k], h["x"], h["a"], h["w"], h["base"])
        print(f"  {k if name=='h1-test' else '':<14}{name:<10}"
              + "".join(f"{v:>7.3f}" for v in c))

print("\n### hard top-fraction gate, sqrt(W/W_base)")
FR = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9)
print(f"  {'score':<14}{'harness':<10}" + "".join(f"{f:>7.2f}" for f in FR))
for k in ("donor_cos", "donor096", "cos+096"):
    for name in ("h1-test", "orion"):
        h = H[name]
        c = topk_curve(h["S"][k], h["x"], h["a"], h["w"], h["base"])
        print(f"  {k if name=='h1-test' else '':<14}{name:<10}"
              + "".join(f"{v:>7.3f}" for v in c))

# ---- permutation null on the WHOLE procedure (max over the gamma grid) -------------------
print("\n### permutation null: score shuffled, same grid, same argmax.  1000 draws")
print(f"  {'score':<14}{'harness':<10}{'observed max':>13}{'null p50':>10}{'null p95':>10}"
      f"{'null p99':>10}{'p-value':>9}")
NPERM = 1000
for k in ("donor_cos", "donor096", "donor_signfrac", "self_norm", "cos+096"):
    for name in ("h1-test", "orion"):
        h = H[name]
        x, a, w, base = h["x"], h["a"], h["w"], h["base"]
        obs = curve(h["S"][k], x, a, w, base).max()
        null = np.empty(NPERM)
        sc = h["S"][k].copy()
        rng = np.random.default_rng(12345)
        for i in range(NPERM):
            null[i] = curve(rng.permutation(sc), x, a, w, base).max()
        p = float((null >= obs).mean())
        print(f"  {k if name=='h1-test' else '':<14}{name:<10}{obs:>13.3f}"
              f"{np.percentile(null,50):>10.3f}{np.percentile(null,95):>10.3f}"
              f"{np.percentile(null,99):>10.3f}{p:>9.3f}")

# ---- bootstrap over perturbations, at a FIXED gamma chosen on the other harness -----------
print("\n### bootstrap over perturbations, gamma FIXED at 3.0 (no in-sample selection)")
print(f"  {'score':<14}{'harness':<10}{'point':>8}{'boot p5':>9}{'boot p50':>10}{'boot p95':>9}"
      f"{'P(mult>1)':>11}")
NB = 2000
for k in ("donor_cos", "donor096", "cos+096"):
    for name in ("h1-test", "orion"):
        h = H[name]
        x, a, w, sc = h["x"], h["a"], h["w"], h["S"][k]
        n = len(x)
        z = rank01(sc) ** 3.0
        pt = np.sqrt(W_of(z, x, a, w) / W_of(np.ones(n), x, a, w))
        rng = np.random.default_rng(7)
        bs = np.empty(NB)
        for i in range(NB):
            idx = rng.integers(0, n, n)
            bs[i] = np.sqrt(W_of(z[idx], x[idx], a[idx], w[idx])
                            / max(W_of(np.ones(n), x[idx], a[idx], w[idx]), 1e-30))
        print(f"  {k if name=='h1-test' else '':<14}{name:<10}{pt:>8.3f}"
              f"{np.percentile(bs,5):>9.3f}{np.percentile(bs,50):>10.3f}"
              f"{np.percentile(bs,95):>9.3f}{np.mean(bs>1):>11.3f}")
