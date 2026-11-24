"""307 probe (SEQUENCE SIDE ONLY): why does a random hCRISPRi-v2 spacer look 5x more promiscuous
than a 2025-library spacer?  If the 2025 library took only the top-ranked (low predicted off-target)
sgRNAs, the null pool must be restricted the same way -- that IS part of "the rule that selected
them".  No expression value is read here."""
from __future__ import annotations

import collections
import csv
import gzip
import json
import pickle
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E301 = ROOT / "experiments/301-metadata-acquisitions/work"
E298 = ROOT / "experiments/298-reviewer-checks/work"
HERE = Path(__file__).resolve().parent
LS = (12, 11, 10, 9)

D = pickle.load(open(HERE / "poolhits307.pkl", "rb"))
pool, hits = D["pool"], D["hits"]
sp2ix = collections.defaultdict(list)
for i, (sp, g) in enumerate(pool):
    sp2ix[sp].append(i)

hv = pickle.load(open(E301 / "horlbeck_supp3.pkl", "rb"))["hCRISPRi-v2"]
ix = {h: i for i, h in enumerate(hv[0])}


def clean(v) -> str:
    v = str(v)
    return v[2:-1] if v.startswith('="') and v.endswith('"') else v


meta_rows = [r for r in hv[1]
             if len(clean(r[ix["protospacer sequence"]])) >= 12
             and not (set(clean(r[ix["protospacer sequence"]]).upper()) - set("ACGT"))]
assert len(meta_rows) == len(pool), (len(meta_rows), len(pool))


def num(v):
    try:
        return float(v)
    except Exception:
        return float("nan")


rank = np.asarray([num(clean(r[ix["selection rank"]])) for r in meta_rows])
stringency = np.asarray([clean(r[ix["off-target stringency"]]) for r in meta_rows])
sublib = np.asarray([clean(r[ix["Sublibrary"]]) for r in meta_rows])

spac = {}
for r in csv.DictReader(open(E301 / "vcc2025_spacers.csv")):
    if r["route"] == "UNRESOLVED":
        continue
    spac[r["guide_id"]] = r
arms = list(csv.DictReader(open(E298 / "seed298.csv")))
real = []
for a in arms:
    r = spac[a["construct"] + "_A|" + a["construct"] + "_B"] if (
        a["construct"] + "_A|" + a["construct"] + "_B") in spac else None
    if r is None:
        for k, v in spac.items():
            if k.startswith(a["construct"] + "_A|"):
                r = v
                break
    for w in ("A", "B"):
        real.append((a["construct"], w, r["spacer_" + w], r["sgid_" + w], a["target"]))
ntc_real = [(g, w, spac[g]["spacer_" + w], spac[g]["sgid_" + w], "")
            for g in sorted(spac) if g.startswith("non-targeting") for w in ("A", "B")]
print(f"real 2025 spacers: {len(real)} targeting + {len(ntc_real)} non-targeting")

# where do the real spacers sit in the pool?
found = 0
ranks, strs_ = [], []
for c, w, sp, sgid, tgt in real:
    js = sp2ix.get(sp.upper(), [])
    if js:
        found += 1
        ranks.append(float(np.nanmin(rank[js])))
        strs_.append(stringency[js[0]])
print(f"real targeting spacers found verbatim in hCRISPRi-v2: {found}/{len(real)}")
print(f"  their 'selection rank': {collections.Counter(int(r) for r in ranks if r == r)}")
print(f"  their 'off-target stringency': {collections.Counter(strs_)}")
print(f"pool 'selection rank' distribution: "
      f"{dict(sorted(collections.Counter(int(r) for r in rank if r == r).items())[:8])}")
print(f"pool 'off-target stringency': {collections.Counter(stringency.tolist())}")

# per-spacer off-target hit rate: real vs pool, overall and stratified by selection rank
def rate(idxs, L):
    return float(np.mean([1.0 if hits[L][i] else 0.0 for i in idxs]))


real_ix = [sp2ix[sp.upper()][0] for c, w, sp, sgid, tgt in real if sp.upper() in sp2ix]
print("\nper-spacer P(>=1 off-target hit), by seed length:")
print(f"  {'set':34s} " + "  ".join(f"L={L}" for L in LS))
print(f"  {'real 2025 targeting (in pool)':34s} "
      + "  ".join(f"{rate(real_ix, L):.3f}" for L in LS))
print(f"  {'whole hCRISPRi-v2 pool':34s} "
      + "  ".join(f"{rate(range(len(pool)), L):.3f}" for L in LS))
for rk in (1, 2, 3, 4, 5):
    idx = np.flatnonzero(rank == rk)
    if idx.size:
        print(f"  {'pool, selection rank ' + str(rk):34s} "
              + "  ".join(f"{rate(idx, L):.3f}" for L in LS) + f"   (n={idx.size:,})")
for s in sorted(set(stringency.tolist())):
    idx = np.flatnonzero(stringency == s)
    print(f"  {'pool, stringency ' + s:34s} "
          + "  ".join(f"{rate(idx, L):.3f}" for L in LS) + f"   (n={idx.size:,})")

# does the 2025 library's own sublibrary/rank explain it?
out = dict(n_real=len(real), found_in_pool=found,
           real_ranks=collections.Counter(int(r) for r in ranks if r == r),
           real_stringency=collections.Counter(strs_),
           real_rate={str(L): rate(real_ix, L) for L in LS},
           pool_rate={str(L): rate(range(len(pool)), L) for L in LS},
           rank_rate={str(rk): {str(L): rate(np.flatnonzero(rank == rk).tolist(), L) for L in LS}
                      for rk in (1, 2, 3, 4, 5) if (rank == rk).any()},
           stringency_rate={s: {str(L): rate(np.flatnonzero(stringency == s).tolist(), L) for L in LS}
                            for s in sorted(set(stringency.tolist()))})
(HERE / "probe307.json").write_text(json.dumps(out, indent=1, default=str))
print("\nwrote probe307.json")
