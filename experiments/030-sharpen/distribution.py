"""Is the predicted magnitude distribution actually flatter than the truth's?

Experiment 030 found that `sign(P)*|P|^q` with q ~ 1.5-2 is worth +0.038 macro PDS, and read
the mechanism as "`gene_norm` and whitening leave the predicted magnitudes flatter than the
truth's". That is a claim about two distributions, so it should be measurable directly rather
than inferred from a tuned exponent.

If it is right, the exponent that matches the two magnitude distributions should land near the
exponent that maximises PDS -- which would turn a knob into a derived quantity, and would say
the gain is distribution *matching* rather than plain sharpening. If the matching exponent is
somewhere else entirely, the sharpening is doing something other than what I claimed.

Measured per perturbation on log|effect|, over the genes the scorer sees, at each stage of the
pipeline.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
import run as proxy_run  # noqa: E402

folds, common = proxy_run.build_folds("1:1")


def shape(M):
    """Spread of log|effect| within each perturbation, plus the top-decile mass share."""
    A = np.abs(M)
    sd, share = [], []
    for row in A:
        v = row[row > 0]
        if v.size < 100:
            continue
        sd.append(np.std(np.log(v)))
        k = max(1, v.size // 10)
        share.append(np.sort(v)[-k:].sum() / v.sum())
    return float(np.mean(sd)), float(np.mean(share))


print(f"{'stage':28} {'sd(log|e|)':>11} {'top-decile share':>17}")
acc = {}
for h, f in folds.items():
    stages = {
        "donor mean (pre-whiten)": f["P"],
        "whitened g=0.65": f["W"].apply(f["P"], 0.65),
        "whitened g=0.85": f["W"].apply(f["P"], 0.85),
        "TRUTH (held-out line)": f["T"],
    }
    for k, M in stages.items():
        s, sh = shape(M)
        acc.setdefault(k, []).append((s, sh))
for k, v in acc.items():
    s = float(np.mean([x[0] for x in v]))
    sh = float(np.mean([x[1] for x in v]))
    print(f"{k:28} {s:11.4f} {sh:17.4f}")

sd_pred = float(np.mean([x[0] for x in acc["whitened g=0.85"]]))
sd_true = float(np.mean([x[0] for x in acc["TRUTH (held-out line)"]]))
print(f"\nA power transform |e|^q multiplies sd(log|e|) by exactly q.")
print(f"Matching exponent q* = sd_true / sd_pred = {sd_true / sd_pred:.3f}")
print(f"PDS-optimal exponent measured in experiment 030: q ~ 2.0-2.5")
