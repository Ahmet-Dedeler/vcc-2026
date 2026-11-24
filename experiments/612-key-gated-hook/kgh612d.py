"""612d (= 612c with the CLEAN4 tau as TAU612 env; separates coverage from source): delivery-faithful fold check (base + s*add, as kgh612b) for FOUR hook forms: none, plain CLEAN4, key-gated CLEAN4,
key-gated IDENTITY (622's S = z(P_C4-0.5) + 0.5 z(I), Flex refs = the other Arc Flex profiles (A/B/C [+h1]) via hyb622's
profiles/identity; tau 0.05 on CLEAN4 where it exists, identity alone where CLEAN4 is undefined). Also the EV-vs-tail table."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
for p in ("experiments/480-nmae-instrument", "experiments/537-arm-pipeline", "experiments/622-identity-track"):
    sys.path.insert(0, str(ROOT / p))
import fold480 as F
from common480 import boot_delta, load_donor_rows, nmae_rows, load_h1
import foldsign537 as FS
import hyb622 as H
C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4B = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
h = load_h1(); G = np.isfinite(h["T"])
L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"], float))
c3 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, h["targets"], h["genes"], pc_cpm=1.0)[0] for d in F.DON]), 0))
BETA = float(np.median(np.abs(L0[G])) / np.median(np.abs(c3[G])))
prof = H.profiles()
import os
TAU = float(os.environ.get("TAU612", "0.05"))
print(f"CLEAN4 tau {TAU}")
res = {}
for name in C4B:
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], [str(g) for g in S["genes"]], S["targets"]
    cons0 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, np.array(genes), pc_cpm=1.0)[0] for d in S["donors"]]), 0))
    P4 = FS.clean4_P(np.array(genes), C4B[name])
    s4 = np.where(np.isfinite(P4) & (np.abs(P4 - 0.5) >= TAU), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    zc = np.where(np.isfinite(P4), (P4 - 0.5) / np.nanstd(P4 - 0.5), np.nan)
    zi = H.identity(prof, name, [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != name], genes)
    sc_low = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan_to_num(zi))
    sI = np.where(np.isfinite(sc_low) & (sc_low != 0), np.where(sc_low > 0, -1.0, 1.0), 0.0)
    S4, SI = np.broadcast_to(s4, T.shape), np.broadcast_to(sI, T.shape)
    sc = np.sign(cons0); add = F.mde_add(cnt)[None, :]
    weak = np.zeros(T.shape, bool)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i]); nz = a > 0
        weak[i] = (~nz) | (a < (np.quantile(a[nz], 0.5) if nz.any() else 0.0))
    base = BETA * cons0
    P = {"nohook": base + sc * add, "plain": base + np.where(S4 != 0, S4, sc) * add,
         "kgC4": base + np.where(weak & (S4 != 0), S4, sc) * add, "kgID": base + np.where(weak & (SI != 0), SI, sc) * add}
    rows = {k: nmae_rows(v, T) for k, v in P.items()}
    res[name] = {k: boot_delta(rows[k], rows["nohook"]) for k in ("plain", "kgC4", "kgID")}
    res[name]["kgID_vs_kgC4"] = boot_delta(rows["kgID"], rows["kgC4"])
    print(f"{name:14s} vs no hook (Δnmae, + = better): " + "  ".join(f"{k} {v[0]:+.5f} [{v[1]:+.5f},{v[2]:+.5f}]" for k, v in res[name].items()))
# EV vs tail: friendly = mean over the two Orion folds; hostile = k562; neutral = h1
fr = {k: np.mean([res[s][k][0] for s in ("orion_hct116", "orion_hek293t")]) for k in ("plain", "kgC4", "kgID")}
ho = {k: res["k562_gwps"][k][0] for k in ("plain", "kgC4", "kgID")}
print("\nEV per context (raw nmae gain vs no hook) as a function of P(hostile line); friendly = Orion mean, hostile = K562 fold:")
for ph in (0.0, 0.1, 0.2, 1 / 3, 0.5):
    print(f"  P(hostile) {ph:.2f}: " + "  ".join(f"{k} {(1 - ph) * fr[k] + ph * ho[k]:+.5f}" for k in fr))
for a, b in (("plain", "kgC4"), ("plain", "kgID")):
    df, dh = fr[a] - fr[b], ho[b] - ho[a]
    print(f"  break-even P(hostile) for {b} over {a}: {df / (df + dh):.3f} (friendly edge of {a} {df:+.5f}, hostile edge of {b} {dh:+.5f})")
