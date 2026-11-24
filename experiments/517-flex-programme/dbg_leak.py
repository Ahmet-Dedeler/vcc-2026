"""517 debug: the Arc block prior's extra calls -- where are they and is the block mass balance failing?"""
import json, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
TAG = "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse"
for name, d in (("c4", "experiments/435-zs-on-coarse/work/h1"), ("g33", "experiments/467-strength-gated-programme/work/h1_gated33"),
                ("arc", "experiments/462-broad16-lab-signature/work/h1_arcfull")):
    dg = json.loads((ROOT / d / f"diag-{TAG}.json").read_text())
    st = np.load(ROOT / d / f"stats-{TAG}.npz")
    called = (st["q"] < 0.05).sum(1)
    kf = np.array([x["prior_shrink_k"] for x in dg]); uns = np.array([x.get("call_set_unsuppressible", 0) for x in dg])
    clip = np.array([x.get("clip_frac", 0) for x in dg])
    keys = [k for k in dg[0] if "prior" in k or "mde" in k]
    print(name, "called/target median", np.median(called), "mean", called.mean(), "| kf median", np.median(kf), "kf==0:", int((kf == 0).sum()),
          "kf==1:", int((kf == 1).sum()), "| unsuppressible mean", uns.mean(), "| clip_frac mean", clip.mean())
print(sorted(dg[0].keys()))
