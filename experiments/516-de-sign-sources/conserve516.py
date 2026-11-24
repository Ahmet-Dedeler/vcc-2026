"""516: is the Stouffer-6 gain concentrated in targets whose donor rows agree with each other (a target-level
conservation key, computable on any surface without the recipient's truth)? Per-target conservation r_p = mean pairwise
Pearson (over genes) between the available standardised source rows. H1 (vs production, uncalled genes) and folds
(vs the production analog, rest genes)."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
import fold516 as FD

d = V.d; T = d["targets"]; G = d["genes"]; prod = d["prod"]; truth = d["truth"]; pool = d["pool"]; calls = d["calls"]
six = ["K", "H", "E", "R", "S8", "S48"]


def conservation(Ss):
    T_ = Ss[0].shape[0]; r = np.full(T_, np.nan)
    for k in range(T_):
        cs = []
        for i in range(len(Ss)):
            for j in range(i + 1, len(Ss)):
                a, b = Ss[i][k], Ss[j][k]; m = np.isfinite(a) & np.isfinite(b)
                if m.sum() > 200:
                    cs.append(np.corrcoef(a[m], b[m])[0, 1])
        if cs:
            r[k] = np.mean(cs)
    return r


Ss = [V.std_score(V.SRC[k], False, 3.0) for k in six]
r = conservation(Ss)
Z = V.stouffer(six, False, 3.0); S = np.sign(np.nan_to_num(Z))
H = np.where(~calls & (np.abs(np.nan_to_num(Z)) >= 0.5), S, prod)
m = pool & (prod != 0)
q = np.nanquantile(r, [1/3, 2/3])
print(f"H1 conservation r_p: median {np.nanmedian(r):.3f}, tertile cuts {q.round(3)}; corr(r_p, n_real) "
      f"{np.corrcoef(np.nan_to_num(r), np.log1p(d['n_real']))[0,1]:+.2f}")
for lo, hi, nm in ((-1, q[0], "low"), (q[0], q[1], "mid"), (q[1], 1, "high")):
    tg = (r >= lo) & (r < hi) if nm != "high" else (r >= lo)
    mm = m & tg[:, None]
    a0 = (prod[mm] == truth[mm]).mean(); a1 = (H[mm] == truth[mm]).mean()
    print(f"  H1 {nm:4s} r: {tg.sum():2d} targets, pool cells {mm.sum():6d}, prod {a0:.4f} -> stz6 {a1:.4f} ({a1-a0:+.4f}); "
          f"median n_real {np.median(d['n_real'][tg]):.0f}")
for name, spec in FD.FOLDS.items():
    Fo = FD.load_fold(name)
    donors = [x for x in spec["prod"] if x != spec["sib"]]
    srcs = donors + FD.CD4 + spec["extra"]
    A, blk = FD.analogs(Fo, spec, donors)
    base = A["base_mcd4+c4q80"]
    Ssf = [FD.std(Fo["z"], s_) for s_ in srcs]
    rf = conservation(Ssf)
    Zf = FD.stouffer(Fo["z"], srcs)
    Hf = np.where(~blk[None, :] & (np.abs(np.nan_to_num(Zf)) >= 0.5), np.sign(np.nan_to_num(Zf)), base)
    mf = Fo["pool"] & (base != 0)
    qf = np.nanquantile(rf, [1/3, 2/3])
    for lo, hi, nm in ((-1, qf[0], "low"), (qf[0], qf[1], "mid"), (qf[1], 1, "high")):
        tg = (rf >= lo) & (rf < hi) if nm != "high" else (rf >= lo)
        mm = mf & tg[:, None]
        a0 = (base[mm] == Fo["truth"][mm]).mean(); a1 = (Hf[mm] == Fo["truth"][mm]).mean()
        print(f"  {name:13s} {nm:4s} r: {tg.sum():3d} targets, pool cells {mm.sum():6d}, analog {a0:.4f} -> stz {a1:.4f} "
              f"({a1-a0:+.4f})")
