"""516 step 1: inventory of same-target sign sources on H1's DE pool.

For every donor bank (never h1_2025 / h1_2025_full: they cover h1-test's targets), the sign of its same-target row
(log2 pseudobulk ratio, pc 1 CPM; CD4 = DESeq2 log2FC per condition) is compared with production's emitted sign on
exactly the pool cells the donor covers, per pool band, with a target-bootstrap CI of the difference. Also: the
donor consensus, the CLEAN4 programme, and the essential-screen coverage flag (430: essential screens cover 0/300 A/B/C
targets, so their H1 value cannot transfer).
    PYTHONDONTWRITEBYTECODE=1 uv run python inv516.py
"""
from __future__ import annotations

import json

import numpy as np

from common516 import (ESS_DONORS, PANEL_DONORS, ROOT, WORK, cd4_rows, donor_rows, fmt_acc, known_answers,
                       load_h1, sign_table)

d = load_h1()
got, ok, _ = known_answers(d)
print(f"KNOWN ANSWER {'PASS' if ok else 'FAIL'}: fid {got['fid']!r} jac {got['jac']!r} nmae {got['nmae']!r} "
      f"reach {got['reach']!r}", flush=True)
assert ok
T, G = d["targets"], d["genes"]
pool, band, truth, prod = d["pool"], d["band"], d["truth"], d["prod"]
print(f"H1: {len(T)} targets, pool cells {int(pool.sum())}, bands (targets) "
      f"{[int((band == b).sum()) for b in range(3)]}, pool cells per band "
      f"{[int(pool[band == b].sum()) for b in range(3)]}")

rows, has = {}, {}
for dn in PANEL_DONORS + ESS_DONORS + ["viperturb_k562", "194-debias-orion_hct116", "194-debias-orion_hek293t"]:
    L2, h, _ = donor_rows(dn, T, G)
    rows[dn], has[dn] = L2, h
cd4 = cd4_rows(T, G, layers=("log_fc", "zscore"))
for c, r in cd4.items():
    rows[f"cd4_{c}"] = r["log_fc"]
    has[f"cd4_{c}"] = np.isfinite(r["log_fc"]).any(1)
ess_cov = np.zeros(len(T), bool)
for dn in ESS_DONORS:
    ess_cov |= has[dn]
print(f"essential-screen donors cover {ess_cov.sum()} of {len(T)} H1 targets; they hold "
      f"{pool[ess_cov].sum() / pool.sum():.3f} of the pool cells")

z = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
pidx = {g: i for i, g in enumerate(z["genes"].astype(str))}
Pd = np.array([z["P_down"][pidx[g]] if g in pidx else np.nan for g in G])
rows["clean4_programme"] = np.tile(np.where(np.isfinite(Pd), 0.5 - Pd, np.nan), (len(T), 1))
has["clean4_programme"] = np.ones(len(T), bool)

pan = np.stack([rows[k] for k in PANEL_DONORS])
rows["cons_panel3"] = np.where(np.isfinite(pan).any(0), np.nanmean(pan, 0), np.nan)
has["cons_panel3"] = np.isfinite(rows["cons_panel3"]).any(1)
cd4s = np.stack([rows[k] for k in rows if k.startswith("cd4_")])
rows["cd4_mean3"] = np.where(np.isfinite(cd4s).any(0), np.nanmean(cd4s, 0), np.nan)
has["cd4_mean3"] = np.isfinite(rows["cd4_mean3"]).any(1)

out = {"ka": got, "ess_cov_targets": [str(t) for t in T[ess_cov]]}
subsets = {"all85": np.ones(len(T), bool), "no_ess": ~ess_cov}
print(f"\n{'source':28s} {'tg':>4s} {'poolcov':>7s}  pooled sign acc on covered pool cells: production -> source, "
      f"delta [95% target bootstrap] n")
for name, R in rows.items():
    S = np.sign(np.nan_to_num(R))
    cov = np.isfinite(R) & (S != 0)
    out[name] = {}
    for sn, sub in subsets.items():
        Ssub = np.where(sub[:, None], S, 0.0)
        H = np.where(cov & sub[:, None], S, prod)
        tab = sign_table(truth, pool & sub[:, None], band, prod, H, cover=cov & sub[:, None])
        out[name][sn] = {b: list(map(float, v)) for b, v in tab.items()}
        pc = (cov & pool & sub[:, None]).sum() / max(1, (pool & sub[:, None]).sum())
        print(f"{name + ' [' + sn + ']':28s} {int((has[name] & sub).sum()):4d} {pc:7.3f}  " +
              " | ".join(f"{b}: {fmt_acc(tab[b])}" for b in ("all", "<20", "20-100", ">=100")), flush=True)
json.dump(out, open(WORK / "inv516.json", "w"), indent=1)
