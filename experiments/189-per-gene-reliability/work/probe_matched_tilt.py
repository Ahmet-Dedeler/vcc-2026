"""Matched-tilt comparison: z-score every prior so the fitting sets are compared at the SAME
perturbation strength, and read each against its own 30-seed permutation null.

Without this the `tau` differences (0.083 for HCT116+K562 against 0.397 for the pilot's pair)
confound "how good the prior is" with "how hard it tilts |H|".
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness189 import Surface189, check, paired_bootstrap
from prior189 import LINES, fit_hier, load_pairs

T = Surface189(); check(T)
base = np.log(np.maximum(np.abs(T.H), 1e-30))
ref = T.reach_vec(np.abs(T.H))
genes, names, K, N, B = load_pairs()
SETS = [("2L pilot pair HCT+HEK", ["HCT116","HEK293T"]),
        ("2L HCT+K562", ["HCT116","K562"]),
        ("3L pilot three-way", ["HCT116","HEK293T","K562"]),
        ("5L delivered", LINES)]
print(f"\n{'fitting set':<24}{'s':>5}{'reach':>9}{'d vs |H|':>10}{'SE':>8}{'P(bet)':>8}"
      f"{'null mean':>11}{'null sd':>9}{'z vs null':>10}")
for tag, sub in SETS:
    cols = [j for j,p in enumerate(names) if p.split('|')[0] in sub and p.split('|')[1] in sub]
    th,_,_ = fit_hier(K, N, cols)
    thh = T.map_gene(genes, th)
    z = (thh - thh.mean())/thh.std()                 # matched tilt strength
    for s in (0.15, 0.30):
        key = np.exp(base + s*z)
        r = T.reach(key)
        d,se,lo,hi,pb = paired_bootstrap(T.reach_vec(key), ref)
        nul = [T.reach(np.exp(base + s*np.random.default_rng(3000+q).permutation(z)))
               for q in range(30)]
        print(f"{tag:<24}{s:>5.2f}{r:>9.5f}{d:>+10.5f}{se:>8.5f}{pb:>8.3f}"
              f"{np.mean(nul):>11.5f}{np.std(nul,ddof=1):>9.5f}"
              f"{(r-np.mean(nul))/np.std(nul,ddof=1):>+10.2f}")
