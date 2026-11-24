"""Weighted vs unweighted decile means -- the last 0.003 of the pilot table."""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common189 import WORK, decile_spread, partial_rank
z = np.load(os.path.join(WORK, "pilot189.npz"), allow_pickle=True)
s = z["support"]
rHE = z["kHE"][s]/z["nHE"][s]; rHK = z["kHK"][s]/z["nHK"][s]; rEK = z["kEK"][s]/z["nEK"][s]
b = z["b_hct"][s]
mag = None
for tag, w in (("support-weighted", z["nHK"][s]), ("unweighted", None)):
    hi, lo, sp = decile_spread(rHE, rHK, w=w)
    print(f"{tag:<18} r_HE -> r_HK   top {hi:.4f} bottom {lo:.4f} spread {sp:+.4f}   (pilot 0.5515/0.5113/+0.0402)")
