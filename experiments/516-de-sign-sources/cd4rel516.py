"""516: CD4's internal replicate structure as a reliability weight. GWCD4i reports, per (target, condition), the mean
Pearson of per-gene log2FC across disjoint donor pairs (donor_correlation_all_mean) and across guides
(guide_correlation_all). Weight each CD4 condition's Stouffer term by a monotone function of that reliability, on H1
(vs production, uncalled genes) and on the folds (vs the production analog, rest genes)."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
import fold516 as FD
from common516 import price, summary_line

d = V.d; z = V.z


def weighted_stz(zz, names, rel_key, f):
    num = 0; den = 0
    for n in names:
        sc = FD.std(zz, n)
        if n.startswith("cd4_") and rel_key:
            r = np.nan_to_num(zz[f"{n}__meta_{rel_key}"].astype(float), nan=0.0)
            w = f(r)[:, None]
        else:
            w = np.ones((sc.shape[0], 1))
        num = num + w * np.nan_to_num(sc); den = den + (w ** 2) * np.isfinite(sc)
    return np.where(den > 0, num / np.sqrt(np.maximum(den, 1e-9)), np.nan)


six = ["k562_gwps", "orion_hct116", "orion_hek293t", "cd4_Rest", "cd4_Stim8hr", "cd4_Stim48hr"]
for c in ("Rest", "Stim8hr", "Stim48hr"):
    r = z[f"cd4_{c}__meta_donor_correlation_all_mean"]
    print(f"H1 CD4 {c}: donor_correlation_all_mean median {np.nanmedian(r):.3f} (IQR {np.nanpercentile(r,25):.3f}-"
          f"{np.nanpercentile(r,75):.3f}); guide corr median {np.nanmedian(z[f'cd4_{c}__meta_guide_correlation_all']):.3f}")
funcs = {"none": None, "clip(r,0)*4": lambda r: 4 * np.clip(r, 0, None), "1+4r": lambda r: np.clip(1 + 4 * r, 0, None),
         "r>0.05 only": lambda r: (r > 0.05).astype(float)}
for key in ("donor_correlation_all_mean", "guide_correlation_all"):
    for fn, f in funcs.items():
        if fn == "none" and key != "donor_correlation_all_mean":
            continue
        Z = weighted_stz(z, six, None if fn == "none" else key, f)
        C = np.where(~d["calls"] & (np.abs(np.nan_to_num(Z)) >= 0.5), np.sign(np.nan_to_num(Z)), 0.0)
        r = price(d, C, f"H1 stz6 CD4 weight {fn} by {key[:12]}", subsets=V.SUBS, verbose=False)
        for sn in V.SUBS:
            print(summary_line(r[sn]), flush=True)
        rows = []
        for name, spec in FD.FOLDS.items():
            Fo = FD.load_fold(name); donors = [x for x in spec["prod"] if x != spec["sib"]]
            A, blk = FD.analogs(Fo, spec, donors)
            Zf = weighted_stz(Fo["z"], donors + FD.CD4, None if fn == "none" else key, f)
            Cf = np.where(~blk[None, :] & (np.abs(np.nan_to_num(Zf)) >= 0.5), np.sign(np.nan_to_num(Zf)), 0.0)
            FD.price_fold(Fo, A["base_mcd4+c4q80"], Cf, f"stz CD4 weight {fn} by {key[:12]} vs mcd4+c4q80", rows)
