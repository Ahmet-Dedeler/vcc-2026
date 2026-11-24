"""516 step 2 (H1): unsupervised combinations of same-target sources as sign sources, priced against production.

Sources (panel-available: they cover A/B/C targets): K562-GWPS (K), Orion HCT116 (H), Orion HEK293T (E), CD4 Rest /
Stim8hr / Stim48hr (R, S8, S48; GWCD4i DESeq2), VIPerturb K562 (V). Noise model from sources516 (per-gene mad_g,
per-row tau_p, bank programme med_g). Constructions:
  maj(set)        sign majority vote, override only where |votes margin| >= m
  stz(set)        Stouffer of uncentred standardised scores x / (mad_g tau_p), capped at |c|; override where |Z| >= thr
  stzc(set)       same on centred scores (x - med_g): the target-specific deviation from the bank's programme
  cd4z            Stouffer of the three CD4 Wald z
Head protection: production keeps its sign on its own top-N genes per target (reach lives there; 513).
Every row is priced on all 85 targets and on the 64 targets no essential screen covers (430: essential donors cover
0/300 A/B/C targets; here no essential bank is used as a source either way).
    PYTHONDONTWRITEBYTECODE=1 uv run python vote516.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

from common516 import WORK, known_answers, load_h1, price, summary_line

d = load_h1()
got, ok, _ = known_answers(d)
print(f"KNOWN ANSWER {'PASS' if ok else 'FAIL'}: fid {got['fid']!r} jac {got['jac']!r} nmae {got['nmae']!r} "
      f"reach {got['reach']!r}", flush=True)
assert ok
T, G = d["targets"], d["genes"]
z = np.load(WORK / "src-h1.npz", allow_pickle=True)
assert (z["targets"] == T).all() and (z["genes"] == G).all()
ess = np.zeros(len(T), bool)
for bn in ("k562_essential", "rpe1", "nadig_jurkat"):
    ess |= z[f"{bn}__has"]
SUBS = {"all85": np.ones(len(T), bool), "no_ess": ~ess}

# production's rank position of each gene in its own order (0 = first)
pos = np.empty((len(T), len(G)), np.int32)
for k in range(len(T)):
    pos[k, d["orders"][k]] = np.arange(len(G))


def std_score(name, centred=False, cap=None):
    """x / (mad_g tau_p) (or centred), NaN where absent."""
    if name.startswith("cd4_"):
        x, med, mad, tau = z[f"{name}__log_fc"], z[f"{name}__med"], z[f"{name}__mad"], z[f"{name}__tau"]
    else:
        x, med, mad, tau = z[f"{name}__x"], z[f"{name}__med"], z[f"{name}__mad"], z[f"{name}__tau"]
    x = x.astype(float)
    if centred:
        x = x - med[None, :]
    s = x / (mad[None, :] * tau[:, None])
    if cap is not None:
        s = np.clip(s, -cap, cap)
    return s


SRC = {"K": "k562_gwps", "H": "orion_hct116", "E": "orion_hek293t", "V": "viperturb_k562",
       "R": "cd4_Rest", "S8": "cd4_Stim8hr", "S48": "cd4_Stim48hr"}


def maj(keys, margin):
    S = np.stack([np.sign(np.nan_to_num(z[f"{SRC[k]}__x"] if not SRC[k].startswith("cd4_")
                                        else z[f"{SRC[k]}__log_fc"])) for k in keys])
    v = S.sum(0)
    return np.where(np.abs(v) >= margin, np.sign(v), 0.0)


def stouffer(keys, centred=False, cap=3.0, weights=None):
    Ss = [std_score(SRC[k], centred, cap) for k in keys]
    w = np.ones(len(keys)) if weights is None else np.asarray(weights, float)
    num = sum(wi * np.nan_to_num(s) for wi, s in zip(w, Ss))
    den = np.sqrt(sum((wi ** 2) * np.isfinite(s) for wi, s in zip(w, Ss)))
    return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)


def gate(Z, thr, head=0):
    C = np.where(np.nan_to_num(np.abs(Z)) >= thr, np.sign(np.nan_to_num(Z)), 0.0)
    if head:
        C = np.where(pos < head, 0.0, C)
    return C


rows = []


def run(C, label):
    r = price(d, C, label, subsets=SUBS, verbose=False)
    for sn in SUBS:
        line = summary_line(r[sn])
        print(line, flush=True)
        rows.append(dict(label=label, subset=sn, acc=r[sn]["acc"], acc_tmean=r[sn]["acc_tmean"],
                         acc_w=r[sn]["acc_w"], members=r[sn]["members"], coverage=r[sn]["coverage_pool"]))


if __name__ == "__main__":
    six = ["K", "H", "E", "R", "S8", "S48"]
    # sanity: production itself (all deltas must be exactly 0)
    run(np.zeros((len(T), len(G))), "KA: empty candidate = production")
    for m in (1, 2, 3):
        run(maj(six, m), f"maj6 margin>={m}")
    for keys, lab in ((["K", "H", "E"], "maj panel3"), (["R", "S8", "S48"], "maj cd4x3"), (["K", "R", "S8", "S48"],
                      "maj K+cd4x3"), (["K", "H", "E", "S48"], "maj panel3+S48"), (six + ["V"], "maj7 (+V)")):
        run(maj(keys, 1), f"{lab} margin>=1")
    for cap in (3.0,):
        for thr in (0.0, 0.5, 1.0, 1.5, 2.0):
            run(gate(stouffer(six, False, cap), thr), f"stz6 uncentred cap{cap} |Z|>={thr}")
        for thr in (0.0, 1.0, 2.0):
            run(gate(stouffer(six, True, cap), thr), f"stz6 centred cap{cap} |Z|>={thr}")
    for head in (10, 25, 50, 100):
        run(gate(stouffer(six, False, 3.0), 0.0, head), f"stz6 uncentred |Z|>=0 head{head} protected")
        run(np.where(pos < head, 0.0, maj(six, 1)), f"maj6 margin>=1 head{head} protected")
    json.dump(rows, open(WORK / "vote516.json", "w"), indent=1)
