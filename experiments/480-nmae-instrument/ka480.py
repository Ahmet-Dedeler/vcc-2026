"""480 step 1: known answers on h1-test with 480's own nmae code (aborts on failure).
KA-1: 426's base_s0_bps delivered lfc -> nmae 0.95237 (434's base).
KA-2: 434's analytic MDE row (c 0.1, cap 0.1, min count 0.3, our sign) -> delta -0.0098.
Also reports the emitted coarse-MDE arm (435 stats, the shipped DE construction) for later tables."""
import json
import sys

import numpy as np

from common480 import HERE, ROOT, load_h1, nmae, nmae_rows

h = load_h1()
T, cnt, genes, targets = h["T"], h["cnt"], h["genes"], h["targets"]
print("scored targets", int(np.isfinite(T).any(1).sum()), "gate cells", int(np.isfinite(T).sum()))
s426 = np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz", allow_pickle=True)
assert (s426["targets"].astype(str) == targets).all()
L = np.nan_to_num(np.asarray(s426["lfc"], float))
base = nmae(L, T)
print(f"KA-1 base_s0_bps nmae {base:.5f} (expect 0.95237)")
ok1 = abs(base - 0.95237) < 5e-5

# 434's MDE row: s = our sign; zeros take the donor-mean sign (rare), then +1
z = np.load(ROOT / "experiments/434-posterior-median-magnitude/work/allcells434.npz", allow_pickle=True)
D = z["D"].astype(float)
s = np.sign(L)
dm = np.nan_to_num(np.nanmean(D, 0))
s = np.where(s == 0, np.sign(dm), s)
s = np.where(s == 0, 1.0, s)
add = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)
mde = L + s * add[None, :]
d_mde = nmae(mde, T) - base
print(f"KA-2 analytic MDE delta {d_mde:+.5f} (expect -0.0098)")
ok2 = abs(d_mde - (-0.0098)) < 2e-4

s435 = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz",
               allow_pickle=True)
assert (s435["targets"].astype(str) == targets).all()
Lc = np.nan_to_num(np.asarray(s435["lfc"], float))
coarse = nmae(Lc, T)
print(f"emitted coarse MDE (435 prodem_mde0.1cap0.1min0.3coarse) nmae {coarse:.5f} (delta {coarse - base:+.5f})")
json.dump(dict(base=base, mde_delta=d_mde, coarse=coarse, ok1=bool(ok1), ok2=bool(ok2)),
          open(HERE / "work/ka480.json", "w"), indent=1)
if not (ok1 and ok2):
    print("KNOWN ANSWER FAILED")
    sys.exit(1)
print("known answers pass")
