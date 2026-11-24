"""280 -- the reliability-vs-depth curve for GWCD4i, and the centring robustness check.

(a) DEPTH: the guide split (donors identical on both sides) and the donor split, both thinned to a ladder of
    counting budgets, so "would more cells reach F1's Flex-400 regime?" is answered by measurement plus a
    Spearman-Brown extrapolation rather than by assumption.
(b) CENTRING: recompute the primary across-donor cosine with the residual centred on 500 random NON-PANEL
    targets (a genome-wide sample) instead of the 297 panel targets.
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from lib280 import DORDER, SEED, Bank, boot_median, cosines, thin  # noqa: E402

P = lambda *a: print(*a, flush=True)  # noqa: E731
COND = "Stim48hr"
bank = Bank(".")
SC, CC = bank.scored, bank.cols_cd4
o = bank.obs
tg = o[(o.gtype == "targeting") & o.keep_for_DE]
IDX = {(c, d, t): g.sort_values("guide") for (c, d, t), g in tg.groupby(["cond", "D", "pert"])}
ALLT = sorted(set(tg.pert))

# ---------------------------------------------------------------- (b) the centring sample
zx = np.load("obs280.npz", allow_pickle=True)
xrows = np.load("extra280_rows.npy")
XX = np.load("extra280_X.npy", mmap_mode="r")
xp = zx["pert"][xrows].astype(str)
xd = pd.Series(zx["donor"][xrows]).astype(str).map({"CE0008162": "D1", "CE0010866": "D2",
                                                    "CE0008678": "D3", "CE0006864": "D4"}).values
xkeep = zx["keep_for_DE"][xrows]


def extra_delta(donors, targets):
    b = bank.ctrl_profile(COND, donors)[CC]
    out = []
    for t in targets:
        m = (xp == t) & np.isin(xd, list(donors)) & xkeep
        if m.sum() == 0:
            continue
        S = np.asarray(XX[np.where(m)[0]], dtype=np.float64).sum(0)
        out.append(bank.profile(S)[CC] - b)
    return np.asarray(out)


# ---------------------------------------------------------------- (a) the depth ladder
ok = [t for t in ALLT if all((COND, d, t) in IDX and len(IDX[(COND, d, t)]) >= 2 for d in DORDER)]
allfour = [t for t in ALLT if all((COND, d, t) in IDX for d in DORDER)]
P(f"{COND}: {len(ok)} targets with 2 guides in all 4 donors; {len(allfour)} present in all 4 donors")
rows = []
for frac in (0.0625, 0.125, 0.25, 0.5, 1.0):
    rng = np.random.default_rng(SEED)

    def build(spec, targets):
        donors = sorted({d for d, _ in spec}) if spec and spec[0][1] is not None else sorted({d for d, _ in spec})
        b = bank.ctrl_profile(COND, donors)[CC]
        M, nc = np.zeros((len(targets), len(CC))), np.zeros(len(targets))
        for k, t in enumerate(targets):
            S = np.zeros(bank.X.shape[1])
            for d, r in spec:
                g = IDX[(COND, d, t)]
                idx = g.i.values if r is None else [int(g.i.values[r])]
                S = S + bank.counts(idx)
                nc[k] += float(g.n_cells.sum() if r is None else g.n_cells.values[r])
            if frac < 1.0:
                S = thin(S, frac * S.sum(), rng)
            M[k] = bank.profile(S)[CC] - b
            nc[k] *= frac
        R = M - M.mean(0, keepdims=True)
        R[:, ~SC] = 0.0
        return R[:, SC], nc

    A, ca = build([(d, 0) for d in DORDER], ok)
    B, cb = build([(d, 1) for d in DORDER], ok)
    g = cosines(A, B)
    mg, glo, ghi, ng = boot_median(g)
    dacc = []
    for pr in (("D1", "D2"), ("D1", "D3"), ("D1", "D4")):
        oth = tuple(d for d in DORDER if d not in pr)
        A2, c2 = build([(d, None) for d in pr], allfour)
        B2, c3 = build([(d, None) for d in oth], allfour)
        dacc.append(cosines(A2, B2))
    dd = np.nanmean(np.vstack(dacc), 0)
    md, dlo, dhi, nd = boot_median(dd)
    eff = float(np.median((ca + cb) / 2))
    P(f"  frac {frac:<6.4f} (~{eff:6.1f} cell-equivalents/side): GUIDE split median {mg:+.4f} [{glo:+.4f},{ghi:+.4f}]"
      f"  |  DONOR split median {md:+.4f} [{dlo:+.4f},{dhi:+.4f}]  (<=0 {np.nanmean(dd<=0):.3f})")
    rows.append(dict(frac=frac, cells_per_side=eff, guide_median=mg, guide_lo=glo, guide_hi=ghi, n_guide=ng,
                     donor_median=md, donor_lo=dlo, donor_hi=dhi, n_donor=nd,
                     donor_frac_le0=float(np.nanmean(dd <= 0))))
df = pd.DataFrame(rows)
df.to_csv("depth280_curve.csv", index=False)

# Spearman-Brown from the measured full-depth donor split to a 400-cell-per-side reference
r = df[df.frac == 1.0].iloc[0]
c1 = float(r.donor_median)
n1 = float(r.cells_per_side) * 2       # a donor-split side is 2 donors of the 4 -> ~half the target's cells
for tgt in (400, 800):
    k = tgt / (float(r.cells_per_side))
    sb = k * c1 / (1 + (k - 1) * c1)
    P(f"  Spearman-Brown: donor-split reliability {c1:+.4f} at {r.cells_per_side:.0f} cells/side "
      f"-> {sb:+.4f} at {tgt} cells/side")

# ---------------------------------------------------------------- (b) centring robustness
P("\n=== centring robustness: panel-centred vs genome-wide-centred (Stim48hr, donor halves) ===")
tt = allfour
for pr in (("D1", "D2"), ("D1", "D3"), ("D1", "D4")):
    oth = tuple(d for d in DORDER if d not in pr)
    out = {}
    for lab, donors in (("A", pr), ("B", oth)):
        b = bank.ctrl_profile(COND, donors)[CC]
        M = np.zeros((len(tt), len(CC)))
        for k, t in enumerate(tt):
            M[k] = bank.profile(bank.counts(np.concatenate([IDX[(COND, d, t)].i.values for d in donors])))[CC] - b
        ex = extra_delta(donors, sorted(set(xp)))
        out[lab] = (M, ex)
    res = {}
    for cname, mu in (("panel", None), ("genome-wide", "ex")):
        cs = []
        for lab in ("A", "B"):
            M, ex = out[lab]
            m = M.mean(0, keepdims=True) if mu is None else ex.mean(0, keepdims=True)
            R = M - m
            R[:, ~SC] = 0.0
            cs.append(R[:, SC])
        res[cname] = cosines(cs[0], cs[1])
    P(f"  {'+'.join(pr)}|{'+'.join(oth)}: panel-centred median {np.nanmedian(res['panel']):+.4f}"
      f"   genome-wide-centred median {np.nanmedian(res['genome-wide']):+.4f}"
      f"   (n {len(tt)}, genome sample {len(out['A'][1])} targets)")
