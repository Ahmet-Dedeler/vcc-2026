"""517: does a same-lab programme transfer better on genes whose control expression is SIMILAR between the source line and
the recipient? (A filter computable for A/B/C: H1 controls vs the context's.) labhook517's instrument, Xaira pair, the
partner's full-programme P_down, complement (hook territory) and block, stratified by |ln frac_L - ln frac_partner| tercile."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from xline517 import BANKS, ESS, THR, NpyRows
gl = {k: [l.strip() for l in open(d / "genes.txt")] for k, (_, d) in BANKS.items()}
common = sorted(set.intersection(*map(set, gl.values()))); G = len(common)
R, fr = {}, {}
for k, (_, d) in BANKS.items():
    pos = {g: i for i, g in enumerate(gl[k])}; col = np.array([pos[g] for g in common])
    e = np.expm1(np.load(d / "b_ctrl.npy")[col]); fr[k] = e / e.sum()
    R[k] = (NpyRows(d / "delta.npy"), [l.strip() for l in open(d / "perts.txt")], col)
def read(k, idx):
    dl, perts, col = R[k]; out = np.empty((len(idx), G), np.float32); order = np.argsort(idx); srt = [idx[i] for i in order]
    for s in range(0, len(srt), 1000):
        out[order[s:s + 1000]] = np.asarray(dl.rows(srt[s:s + 1000]), dtype=np.float32)[:, col]
    return out
rng = np.random.default_rng(5181)
weissman = [k for k in BANKS if BANKS[k][0] == "Weissman"]
for L, partner in (("HCT116", "HEK293T"), ("HEK293T", "HCT116")):
    Xp = read(partner, list(range(len(R[partner][1]))))
    nz = (Xp != 0).sum(0); dn = (Xp < 0).sum(0); Pp = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan); del Xp
    z = 0.5 - Pp; a = np.abs(z); z = np.where(a >= np.nanquantile(a, 0.37), z, np.nan)       # the most confident 63 %
    dexp = np.abs(np.log(np.maximum(fr[L], 1e-9)) - np.log(np.maximum(fr[partner], 1e-9)))
    terc = np.digitize(dexp, np.quantile(dexp, [1 / 3, 2 / 3]))
    lp = R[L][1]
    tg = [i for i, p in enumerate(lp) if p not in ESS and any(p in R[k][1] for k in weissman)]
    tg = sorted(rng.choice(tg, min(2000, len(tg)), replace=False).tolist())
    XL = read(L, tg)
    D = np.zeros((len(tg), G), np.float32); cnt = np.zeros(len(tg))
    for k in weissman:
        idx = {p: i for i, p in enumerate(R[k][1])}
        have = [(j, idx[lp[t]]) for j, t in enumerate(tg) if lp[t] in idx]
        if have:
            Y = read(k, [h for _, h in have])
            for (j, _), y in zip(have, Y): D[j] += y; cnt[j] += 1
    D /= np.maximum(cnt, 1)[:, None]
    blk = np.zeros(G, bool); blk[np.argsort(-fr[L])[:int(round(0.22 * G))]] = True
    acc = {}
    for j in range(len(tg)):
        head = np.zeros(G, bool); head[np.argsort(-np.abs(D[j]))[:int(round(0.04 * G))]] = True
        gate = np.abs(XL[j]) >= THR
        for reg, m0 in (("complement", ~head & ~blk & gate), ("block", blk & ~head & gate)):
            for tt in range(3):
                m = m0 & (terc == tt)
                for src, v in (("programme", z), ("design", np.where(D[j] != 0, D[j], np.nan))):
                    mm = m & np.isfinite(v)
                    k_ = (reg, tt, src); h, n = acc.get(k_, (0, 0))
                    acc[k_] = (h + int((np.sign(v[mm]) == np.sign(XL[j][mm])).sum()), n + int(mm.sum()))
    print(f"\n{L} (programme from {partner}); |Δ ln expression| terciles (similar → dissimilar)")
    for reg in ("complement", "block"):
        print(f"  {reg:10s} " + " | ".join(f"T{tt}: programme {acc[(reg, tt, 'programme')][0] / max(acc[(reg, tt, 'programme')][1], 1):.4f} "
                                          f"design {acc[(reg, tt, 'design')][0] / max(acc[(reg, tt, 'design')][1], 1):.4f} (n {acc[(reg, tt, 'programme')][1]})"
                                          for tt in range(3)))
