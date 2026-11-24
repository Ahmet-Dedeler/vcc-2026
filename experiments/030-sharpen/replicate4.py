"""030e — replicate the sharpening gain on the independent 4-line LODO, and fix the order.

Whitening became credible when it replicated on a second harness with three predictors and a
mostly different gene population (exp 020). Sharpening deserves the same test before it goes
anywhere near a submission — especially since two of my explanations for it have already died.

Also settled here: whether the power transform belongs before or after whitening. They do not
commute, and there is no argument from first principles either way.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale   # noqa: E402
from vcc2026.whiten import Whitener                      # noqa: E402

LINES = {"k562_gwps": "data/donors/k562_gwps", "orion_hct116": "data/donors/orion_hct116",
         "orion_hek293t": "data/donors/orion_hek293t", "h1_2025": "data/donors/h1_2025"}

banks = {}
for n, p in LINES.items():
    d = Donor.load(p)
    _, s = d.decompose()
    gl = [str(g) for g in d.genes]
    banks[n] = dict(s=(s / _norm_scale(s)).astype(np.float32), raw=s,
                    gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                    idx={q: i for i, q in enumerate(d.perts)})
common = sorted(set.intersection(*[set(b["idx"]) for b in banks.values()]))
vcc = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
print(f"{len(common)} perturbations across 4 lines ({len(set(common) & vcc)} on the VCC "
      f"panel)\n")


def pds(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


def sharpen(M, q):
    return M if q == 1.0 else np.sign(M) * np.abs(M) ** q


res = {}
for held in banks:
    others = [n for n in banks if n != held]
    shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
    feat = [g for g in shared if g not in set(common)]
    col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
           for n in [held] + others}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][q]][col[held]]
                   for q in common]).astype(np.float64)
    sc = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
          for n in others}
    ref = float(np.median(list(sc.values())))
    parts = {n: (ref / sc[n]) * np.vstack(
        [banks[n]["s"][banks[n]["idx"][q]][col[n]] for q in common]).astype(np.float64)
        for n in others}
    P = sum(parts.values()) / len(others)
    W = Whitener.fit(np.vstack([(ref / sc[n]) * banks[n]["s"][:, col[n]].astype(np.float64)
                                for n in others]))
    print(f"fold {held:14} fitted", flush=True)
    for g in (0.0, 0.55, 0.65):
        for q in (1.0, 1.5, 2.0, 2.5):
            res.setdefault(f"after  g={g:g} q={q:g}", {})[held] = pds(sharpen(W.apply(P, g), q), T)
            if q != 1.0 and g > 0:
                res.setdefault(f"before g={g:g} q={q:g}", {})[held] = \
                    pds(W.apply(sharpen(P, q), g), T)

print(f"\n{'config':22} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:22} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/030-sharpen/replicate4.json").write_text(json.dumps(out, indent=2))
