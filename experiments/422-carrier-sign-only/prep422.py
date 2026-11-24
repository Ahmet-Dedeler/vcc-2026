"""422 prep: SIGN-ONLY carriers -- the shipped carrier magnitudes |cap15| with the signs of 416's gamma-0 design (g000),
falling back to the shipped sign where g000 is exactly 0. After 417c (board: fid/reach up, nmae +0.0089 from the magnitude
shift), this keeps the magnitude profile the board already prices and moves only the signs.
Also 's030' (signs from 416's g030). Reads 416's carriers416.npz and 323's effects.npz (cap15); writes work/carriers422.npz."""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W323 = ROOT / "experiments/323-current-family-validation/work"
with np.load(W323 / "effects.npz", allow_pickle=False) as ze:
    eff = ze["cap15"].astype(np.float64); genes = ze["genes"]; targets = ze["targets"]
with np.load(ROOT / "experiments/416-de-carrier-menu/work/carriers416.npz", allow_pickle=False) as zc:
    assert (zc["genes"] == genes).all() and (zc["targets"] == targets).all()
    g000 = zc["g000"].astype(np.float64); g030 = zc["g030"].astype(np.float64)
out, rec = {}, {}
for name, src in (("s000", g000), ("s030", g030)):
    sg = np.where(src != 0, np.sign(src), np.sign(eff))
    d = np.abs(eff) * sg
    flips = (np.sign(d) != np.sign(eff)) & (eff != 0)
    wflip = float((np.abs(eff) * flips).sum() / np.abs(eff).sum())
    assert np.allclose(np.linalg.norm(d, axis=1), np.linalg.norm(eff, axis=1))
    out[name] = d
    rec[name] = dict(frac_entries_flipped=float(flips.mean()), abs_weighted_frac_flipped=wflip,
                     row_cos_vs_cap15_q=np.quantile((d * eff).sum(1) / (np.linalg.norm(d, axis=1) * np.linalg.norm(eff, axis=1)), [0, .1, .5, .9, 1]).tolist())
    print(name, json.dumps(rec[name]), flush=True)
np.savez(HERE / "work/carriers422.npz", genes=genes, targets=targets, **out)
rec["npz_sha256"] = hashlib.sha256((HERE / "work/carriers422.npz").read_bytes()).hexdigest()
(HERE / "work/carriers422.json").write_text(json.dumps(rec, indent=1))
