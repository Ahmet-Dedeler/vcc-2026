"""307 step 2 (SEQUENCE SIDE ONLY -- no expression value is read here).

The OFF-TARGET selection propensity, measured the way 301 measured it.

`seq307.py`'s `pi_L` counted every hCRISPRi-v2 spacer's hit on its OWN target promoter, which 301's
scan excludes (`if sym == tgt: continue`).  This script rescans the whole hCRISPRi-v2 pool through
301's exact rule -- own-target exclusion and the 10 kb distal filter included -- and reports

  * `pi_off_L(g)`  = P(a spacer drawn uniformly from hCRISPRi-v2 hits gene g at seed length L,
                      OFF TARGET, in the CRISPRi active window, after the distal filter);
  * the empirical distribution of (L*, k) for a pseudo-construct = a random PAIR of pool spacers,
    which is the size/level distribution the real 34 constructs must be compared against;
  * `poolhits307.pkl` -- every pool spacer's off-target hit set at each L, for the secondary
    real-spacer null.

    PYTHONDONTWRITEBYTECODE=1 uv run python pool307.py
"""
from __future__ import annotations

import collections
import csv
import gzip
import json
import pickle
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E301 = ROOT / "experiments/301-metadata-acquisitions/work"
HERE = Path(__file__).resolve().parent
LS = (12, 11, 10, 9)
WIN = (-50, 300)
DISTAL_KB = 10
COMP = str.maketrans("ACGT", "TGCA")
SEED307 = 20260907


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def main() -> None:
    t0 = time.time()
    log = []

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        log.append(s)

    # ---------------- promoter corpus + active-window PAM sites (identical to seq307.py)
    prom = []
    with gzip.open(E301 / "promoters_cage_tss_hg19.fa.gz", "rt") as fh:
        hdr = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                hdr = line[1:].split("|")
            else:
                prom.append((hdr[0], hdr[1], hdr[2], int(hdr[3]), int(hdr[4]),
                             int(hdr[5]), int(hdr[6]), line.upper()))
    wins = collections.defaultdict(list)
    for g, t, c, s, tss, lo, hi, seq in prom:
        wins[g].append((c, lo, hi))

    _dcache: dict = {}

    def distal(tgt: str, cg: str) -> bool:
        if not tgt:
            return True
        key = (tgt, cg)
        v = _dcache.get(key)
        if v is None:
            v = True
            for c1, l1, h1 in wins.get(tgt, ()):
                for c2, l2, h2 in wins.get(cg, ()):
                    if c1 == c2 and max(l2 - h1, l1 - h2, 0) < DISTAL_KB * 1000:
                        v = False
                        break
                if not v:
                    break
            _dcache[key] = v
        return v

    site_kmer, site_gene = [], []
    for sym, eid, chrom, strand, tss, lo, hi, seq in prom:
        for orient, s in (("+", seq), ("-", rc(seq))):
            start = 0
            while True:
                j = s.find("GG", start)
                if j < 0:
                    break
                start = j + 1
                pam0 = j - 1
                if pam0 < 12:
                    continue
                gpos = (lo + pam0 - 1) if orient == "+" else (hi - (pam0 - 1))
                off = (gpos - tss) if strand == 1 else (tss - gpos)
                if not (WIN[0] <= off <= WIN[1]):
                    continue
                k = s[pam0 - 12:pam0]
                if len(k) != 12 or set(k) - set("ACGT"):
                    continue
                site_kmer.append(k)
                site_gene.append(sym)
    P(f"active-window PAM 12-mer sites {len(site_kmer):,}")

    lmer_genes = {}
    for L in LS:
        d = collections.defaultdict(set)
        for k, g in zip(site_kmer, site_gene):
            d[k[-L:]].add(g)
        lmer_genes[L] = {m: tuple(sorted(v)) for m, v in d.items()}
        P(f"  L={L}: {len(d):,} distinct seeds")

    # ---------------- hCRISPRi-v2 pool
    hv = pickle.load(open(E301 / "horlbeck_supp3.pkl", "rb"))["hCRISPRi-v2"]
    ix = {h: i for i, h in enumerate(hv[0])}

    def clean(v) -> str:
        """Horlbeck's supplementary sheet Excel-escapes text cells as `="SYMBOL"` --
        301's `fetch_promoters_hg19.clean`, verbatim.  Without this the own-target
        exclusion never fires and every spacer's ON-target hit is miscounted as an
        off-target (measured: P(>=1 L=12 hit) 0.597 vs 0.107 corrected)."""
        v = str(v)
        return v[2:-1] if v.startswith('="') and v.endswith('"') else v

    pool = []
    for r in hv[1]:
        sp = clean(r[ix["protospacer sequence"]]).upper()
        gene = clean(r[ix["gene"]])
        if len(sp) < 12 or set(sp) - set("ACGT"):
            continue
        pool.append((sp, "" if gene == "negative_control" else gene))
    n_neg = sum(1 for _, g in pool if not g)
    P(f"hCRISPRi-v2 pool {len(pool):,} spacers ({n_neg:,} negative controls, "
      f"{len({g for _, g in pool if g}):,} target labels)")

    # ---------------- off-target hit set of every pool spacer, at every L
    hits = {L: [] for L in LS}          # index-aligned with pool
    for sp, tgt in pool:
        for L in LS:
            hs = lmer_genes[L].get(sp[-L:], ())
            if hs:
                hs = tuple(g for g in hs if g != tgt and distal(tgt, g))
            hits[L].append(hs)
    for L in LS:
        tot = sum(len(h) for h in hits[L])
        P(f"  L={L}: {tot:,} off-target hits, mean {tot / len(pool):.4f} genes per spacer, "
          f"{sum(1 for h in hits[L] if h) / len(pool):.4f} of spacers with >=1")

    pi_off = {}
    for L in LS:
        cnt = collections.Counter()
        for h in hits[L]:
            cnt.update(h)
        pi_off[L] = {g: c / len(pool) for g, c in cnt.items()}
        v = np.asarray(list(pi_off[L].values()))
        P(f"  pi_off_{L}: {len(v):,} reachable genes, sum {v.sum():.4f}, "
          f"median {np.median(v):.3e}, max {v.max():.3e}")

    # ---------------- pseudo-construct = random PAIR of pool spacers
    rng = np.random.default_rng(SEED307)
    NP = 200_000
    ia = rng.integers(0, len(pool), NP)
    ib = rng.integers(0, len(pool), NP)
    lstar, ksize = np.zeros(NP, dtype=np.int8), np.zeros(NP, dtype=np.int32)
    for i in range(NP):
        for L in LS:
            u = set(hits[L][ia[i]]) | set(hits[L][ib[i]])
            if u:
                lstar[i], ksize[i] = L, len(u)
                break
    P("pseudo-construct (random hCRISPRi-v2 spacer PAIR, 301's rule) over "
      f"{NP:,} draws:")
    for L in (12, 11, 10, 9, 0):
        m = lstar == L
        lab = f"L*={L}" if L else "no hit at any L"
        P(f"  {lab:16s} {m.sum():>7,} ({m.mean():.4f})"
          + (f"   k: mean {ksize[m].mean():.2f} median {np.median(ksize[m]):.0f} "
             f"max {ksize[m].max()}" if m.any() and L else ""))

    with open(HERE / "poolhits307.pkl", "wb") as fh:
        pickle.dump(dict(pool=pool, hits=hits, pi_off=pi_off), fh)
    meta = dict(
        n_pool=len(pool), n_neg_controls=n_neg,
        mean_offtarget_hits_per_spacer={str(L): sum(len(h) for h in hits[L]) / len(pool) for L in LS},
        frac_spacers_with_hit={str(L): sum(1 for h in hits[L] if h) / len(pool) for L in LS},
        pi_off_sum={str(L): float(sum(pi_off[L].values())) for L in LS},
        pseudo_construct={f"L*={L}": dict(n=int((lstar == L).sum()),
                                          frac=float((lstar == L).mean()),
                                          k_mean=float(ksize[lstar == L].mean()) if (lstar == L).any() else None,
                                          k_max=int(ksize[lstar == L].max()) if (lstar == L).any() else None)
                          for L in (12, 11, 10, 9)},
        pseudo_no_hit=int((lstar == 0).sum()), n_pseudo=NP, seed=SEED307,
        seconds=round(time.time() - t0, 1))
    (HERE / "pool307_meta.json").write_text(json.dumps(meta, indent=1))
    (HERE / "pool307.log").write_text("\n".join(log) + "\n")
    P(f"done in {meta['seconds']}s")


if __name__ == "__main__":
    main()
