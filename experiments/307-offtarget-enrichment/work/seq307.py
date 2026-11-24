"""307 step 1 (SEQUENCE SIDE ONLY -- no expression value is read here).

Reproduces 301's seed-scan selection rule exactly from its saved promoter corpus and spacer table,
then builds, for EVERY gene measured in the 2025 h5ads, the sequence-side covariates the joint null
needs:

  * `sites_L`      -- the number of PAM-anchored L-mer sites of the gene inside the CRISPRi active
                      window (TSS -50 .. +300), i.e. the number of chances the rule has to pick it;
  * `pi_L`         -- the EXACT marginal probability that a spacer drawn uniformly from the
                      hCRISPRi-v2 library (the library the 2025 guides come from) hits the gene at
                      seed length L under the identical rule.  = sum over the gene's distinct
                      active-window L-mers of the library's empirical frequency of that L-mer;
  * `gc`           -- GC fraction of the active window;
  * `n_windows`    -- number of CAGE promoter windows (P1/P2/P1P2 transcripts) = TSS multiplicity;
  * `tss_density`  -- number of OTHER genes with a promoter window within 10 kb on the same chromosome;
  * `map_mult`     -- mappability proxy: mean number of DISTINCT OTHER genes sharing each of the
                      gene's active-window 12-mers inside the promoter corpus (1.0 = fully unique);
  * `lowcomp`      -- fraction of active-window sites whose 12-mer has <= 2 distinct bases or a
                      >= 6-nt homopolymer run (low-complexity proxy).

Also emits the frozen candidate table (301's rule) and L*(c) per construct.

    PYTHONDONTWRITEBYTECODE=1 uv run python seq307.py
"""
from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import json
import pickle
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E301 = ROOT / "experiments/301-metadata-acquisitions/work"
E298 = ROOT / "experiments/298-reviewer-checks/work"
E305 = ROOT / "experiments/305-t7-expression-test/work"
HERE = Path(__file__).resolve().parent

LS = (12, 11, 10, 9)          # 301's seed lengths, longest first
WIN = (-50, 300)              # 301's CRISPRi active window
DISTAL_KB = 10                # 301's self-locus filter
COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    t0 = time.time()
    log = []

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        log.append(s)

    # ------------------------------------------------------------------ promoter corpus
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
    by_sym = collections.defaultdict(list)
    for p in prom:
        by_sym[p[0]].append(p)
    P(f"promoter corpus: {len(prom)} windows over {len(by_sym)} genes, "
      f"{sum(len(p[7]) for p in prom):,} bp")

    wins = collections.defaultdict(list)
    for g, t, c, s, tss, lo, hi, seq in prom:
        wins[g].append((c, lo, hi))

    def distal(tgt: str, cg: str, kb: int = DISTAL_KB) -> bool:
        """301's self-locus filter, verbatim."""
        if not tgt:
            return True
        for c1, l1, h1 in wins.get(tgt, []):
            for c2, l2, h2 in wins.get(cg, []):
                if c1 == c2 and max(l2 - h1, l1 - h2, 0) < kb * 1000:
                    return False
        return True

    # ------------------------------------------- PAM-anchored 12-mer sites, active window only
    # 301's scan: enumerate 'GG', pam0 = j-1, require pam0 >= 12, seed = s[pam0-L:pam0].
    # Storing the 12-mer means every shorter seed is its suffix -- exactly faithful.
    site_gene: list[str] = []
    site_kmer: list[str] = []
    site_off: list[int] = []
    for sym, eid, chrom, strand, tss, lo, hi, seq in prom:
        for orient, s in (("+", seq), ("-", rc(seq))):
            start = 0
            while True:
                j = s.find("GG", start)
                if j < 0:
                    break
                start = j + 1
                pam0 = j - 1
                if pam0 < LS[0]:
                    continue
                gpos = (lo + pam0 - 1) if orient == "+" else (hi - (pam0 - 1))
                off = (gpos - tss) if strand == 1 else (tss - gpos)
                if not (WIN[0] <= off <= WIN[1]):
                    continue
                k = s[pam0 - 12:pam0]
                if len(k) != 12 or set(k) - set("ACGT"):
                    continue
                site_gene.append(sym)
                site_kmer.append(k)
                site_off.append(off)
    P(f"active-window PAM-anchored 12-mer sites: {len(site_gene):,}")

    site_gene_a = np.asarray(site_gene)
    site_off_a = np.asarray(site_off, dtype=np.int32)

    # 12-mer -> gene list (for the mappability proxy and for L=12 lookups)
    k12_genes: dict[str, set] = collections.defaultdict(set)
    for k, g in zip(site_kmer, site_gene):
        k12_genes[k].add(g)
    P(f"distinct active-window 12-mers: {len(k12_genes):,}")

    # per-L: L-mer -> {gene: [offsets]}
    lmer_genes = {}
    for L in LS:
        d: dict[str, dict[str, list]] = collections.defaultdict(lambda: collections.defaultdict(list))
        for k, g, o in zip(site_kmer, site_gene, site_off):
            d[k[-L:]][g].append(o)
        lmer_genes[L] = d
        P(f"  L={L}: distinct L-mers {len(d):,}")

    # ------------------------------------------------------------------ the 2025 library
    spac = {}
    for r in csv.DictReader(open(E301 / "vcc2025_spacers.csv")):
        if r["route"] == "UNRESOLVED":
            continue
        spac[r["guide_id"]] = r
        if not r["guide_id"].startswith("non-targeting"):
            spac[r["guide_id"].split("|")[0][:-2]] = r
    arms = list(csv.DictReader(open(E298 / "seed298.csv")))
    fail = [a["construct"] for a in arms if a["failing"] == "True"]
    work = [a["construct"] for a in arms if a["failing"] != "True"]
    ntc = sorted(g for g in spac if g.startswith("non-targeting"))
    target_of = {a["construct"]: a["target"] for a in arms}
    split_of = {a["construct"]: a["split"] for a in arms}
    kd_of = {a["construct"]: float(a["kd_lfc"]) for a in arms}
    cells_of = {a["construct"]: int(a["cells"]) for a in arms}
    assert (len(fail), len(work), len(ntc)) == (13, 21, 31), (len(fail), len(work), len(ntc))
    P(f"2025 library: {len(fail)} failing + {len(work)} working targeting constructs, {len(ntc)} NTC")

    # genes any 2025 construct intends to target (exclusion set for the null)
    lib_targets = set()
    for s in ("train", "validation", "test"):
        m = json.load(open(E305 / f"sums305_{s}.json"))
        lib_targets |= {v for v in m["target_of_construct"].values() if v and v != "non-targeting"}
    P(f"genes targeted by SOME construct in the 2025 library (all splits): {len(lib_targets)}")

    # ------------------------------------------------------------------ the scan (301's rule)
    def scan(spacer_a: str, spacer_b: str, own_target: str):
        """Returns {L: {gene: [offsets]}} after own-target exclusion and the distal filter."""
        out = {}
        for L in LS:
            g2o: dict[str, list] = collections.defaultdict(list)
            for sp in (spacer_a, spacer_b):
                for g, offs in lmer_genes[L].get(sp[-L:], {}).items():
                    if own_target and g == own_target:
                        continue
                    g2o[g].extend(offs)
            out[L] = {g: sorted(o) for g, o in g2o.items()
                      if distal(own_target, g)}
        return out

    cand_L = {}
    for c in fail + work:
        r = spac[c]
        cand_L[c] = scan(r["spacer_A"], r["spacer_B"], target_of[c])
    for c in ntc:
        r = spac[c]
        cand_L[c] = scan(r["spacer_A"], r["spacer_B"], "")

    counts = {
        "failing": [sum(1 for c in fail if cand_L[c][12]), len(fail)],
        "working": [sum(1 for c in work if cand_L[c][12]), len(work)],
        "ntc": [sum(1 for c in ntc if cand_L[c][12]), len(ntc)],
    }
    P("KA-2 301 primary reading (L=12, CRISPRi window, distal): "
      + "  ".join(f"{k} {a}/{b}" for k, (a, b) in counts.items()))
    ka2 = (counts["failing"] == [4, 13] and counts["working"] == [4, 21] and counts["ntc"] == [4, 31])
    P(f"      expected failing 4/13, working 4/21, ntc 4/31 -> {'PASS' if ka2 else 'FAIL'}")

    # L*(c) = longest L with >= 1 candidate after the active-window + distal filters
    rows = []
    for grp, members in (("failing", fail), ("working", work), ("ntc", ntc)):
        for c in members:
            Lstar = next((L for L in LS if cand_L[c][L]), None)
            genes = sorted(cand_L[c][Lstar]) if Lstar else []
            rows.append(dict(
                construct=c, group=grp, unit="targeting" if grp != "ntc" else "non-targeting",
                intended_target=target_of.get(c, ""), split=split_of.get(c, "ALL"),
                kd_lfc=kd_of.get(c, ""), cells_full=cells_of.get(c, ""),
                L_star=Lstar, k=len(genes), candidates=",".join(genes),
                n12=len(cand_L[c][12]), n11=len(cand_L[c][11]),
                n10=len(cand_L[c][10]), n9=len(cand_L[c][9]),
                offsets=";".join(f"{g}:{','.join(str(o) for o in cand_L[c][Lstar][g])}" for g in genes),
            ))
    with open(HERE / "constructs307.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    n_targ12 = sum(1 for r in rows if r["unit"] == "targeting" and r["n12"] > 0)
    P(f"targeting constructs with >=1 L=12 candidate: {n_targ12}  (coordinator says 8)")
    P("L* distribution, targeting: "
      + str(dict(collections.Counter(r["L_star"] for r in rows if r["unit"] == "targeting"))))
    P("L* distribution, NTC      : "
      + str(dict(collections.Counter(r["L_star"] for r in rows if r["unit"] == "non-targeting"))))
    P("k distribution, targeting : "
      + str(dict(sorted(collections.Counter(r["k"] for r in rows if r["unit"] == "targeting").items()))))

    # ------------------------------------------------------------------ hCRISPRi-v2 spacer pool
    hv = pickle.load(open(E301 / "horlbeck_supp3.pkl", "rb"))["hCRISPRi-v2"]
    hdrs, hrows = hv[0], hv[1]
    ix = {h: i for i, h in enumerate(hdrs)}
    pool = []
    for r in hrows:
        sp = str(r[ix["protospacer sequence"]]).upper()
        gene = str(r[ix["gene"]])
        if len(sp) < 12 or set(sp) - set("ACGT"):
            continue
        pool.append((sp, gene))
    P(f"hCRISPRi-v2 spacer pool: {len(pool):,} usable protospacers over "
      f"{len({g for _, g in pool}):,} gene labels")

    # exact marginal hit probability pi_L(g) under a uniformly drawn library spacer
    pi = {}
    for L in LS:
        freq = collections.Counter(sp[-L:] for sp, _ in pool)
        n = len(pool)
        acc = collections.defaultdict(float)
        for m, gs in lmer_genes[L].items():
            f = freq.get(m, 0)
            if f:
                for g in gs:
                    acc[g] += f / n
        pi[L] = acc
        P(f"  pi_{L}: {len(acc):,} genes reachable; "
          f"mean {np.mean(list(acc.values())):.3e}, max {max(acc.values()):.3e}")

    # ------------------------------------------------------------------ per-gene covariates
    genes_all = sorted(by_sym)
    # GC + window geometry
    gc, nwin, awlen = {}, {}, {}
    for g in genes_all:
        tot = gcn = 0
        for sym, eid, chrom, strand, tss, lo, hi, seq in by_sym[g]:
            # active window in genomic coordinates
            if strand == 1:
                a, b = tss + WIN[0], tss + WIN[1]
            else:
                a, b = tss - WIN[1], tss - WIN[0]
            i0, i1 = max(0, a - lo), min(len(seq), b - lo + 1)
            sub = seq[i0:i1]
            tot += len(sub)
            gcn += sum(1 for ch in sub if ch in "GC")
        gc[g] = gcn / tot if tot else 0.0
        nwin[g] = len(by_sym[g])
        awlen[g] = tot

    # TSS density: other genes with a window within 10 kb on the same chromosome
    by_chrom = collections.defaultdict(list)
    for g, t, c, s, tss, lo, hi, seq in prom:
        by_chrom[c].append((tss, g))
    for c in by_chrom:
        by_chrom[c].sort()
    tssden = {}
    for g in genes_all:
        n = set()
        for sym, eid, chrom, strand, tss, lo, hi, seq in by_sym[g]:
            arr = by_chrom[chrom]
            pos = np.asarray([x[0] for x in arr])
            i0, i1 = np.searchsorted(pos, tss - 10000), np.searchsorted(pos, tss + 10000)
            n |= {arr[i][1] for i in range(i0, i1)}
        tssden[g] = len(n - {g})

    # mappability proxy + low-complexity, over active-window 12-mer sites
    mapmult = collections.defaultdict(list)
    lowcomp = collections.defaultdict(list)
    for k, g in zip(site_kmer, site_gene):
        mapmult[g].append(len(k12_genes[k]))
        homop = any(b * 6 in k for b in "ACGT")
        lowcomp[g].append(1.0 if (len(set(k)) <= 2 or homop) else 0.0)

    sites = {L: collections.defaultdict(int) for L in LS}
    for k, g in zip(site_kmer, site_gene):
        for L in LS:
            sites[L][g] += 1        # site count is L-independent (same anchors); kept per L for clarity
    mean_off = {g: float(site_off_a[site_gene_a == g].mean()) if (site_gene_a == g).any() else 0.0
                for g in genes_all}

    with open(HERE / "genes307.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "n_sites", "gc", "n_windows", "aw_bp", "tss_density_10kb",
                    "map_mult", "lowcomp", "mean_offset",
                    "pi12", "pi11", "pi10", "pi9", "is_library_target"])
        for g in genes_all:
            w.writerow([g, sites[12][g], round(gc[g], 5), nwin[g], awlen[g], tssden[g],
                        round(float(np.mean(mapmult[g])) if mapmult[g] else 0.0, 5),
                        round(float(np.mean(lowcomp[g])) if lowcomp[g] else 0.0, 5),
                        round(mean_off[g], 2),
                        f"{pi[12].get(g, 0.0):.6e}", f"{pi[11].get(g, 0.0):.6e}",
                        f"{pi[10].get(g, 0.0):.6e}", f"{pi[9].get(g, 0.0):.6e}",
                        int(g in lib_targets)])
    P(f"wrote genes307.csv ({len(genes_all)} genes)")

    # eligible-pool sizes at each L (pi>0, not a library target)
    for L in LS:
        n = sum(1 for g in genes_all if pi[L].get(g, 0) > 0 and g not in lib_targets)
        P(f"  eligible pool at L={L} (pi>0, not a 2025 library target): {n:,} genes")

    meta = dict(
        ka2_301_counts=counts, ka2_pass=bool(ka2),
        n_targeting_with_L12=n_targ12,
        n_windows=len(prom), n_genes_corpus=len(by_sym),
        n_active_sites=len(site_gene), n_distinct_12mers=len(k12_genes),
        n_pool_spacers=len(pool), n_library_targets=len(lib_targets),
        library_targets=sorted(lib_targets),
        Lstar_targeting=dict(collections.Counter(
            r["L_star"] for r in rows if r["unit"] == "targeting")),
        k_targeting=dict(sorted(collections.Counter(
            r["k"] for r in rows if r["unit"] == "targeting").items())),
        eligible_pool={str(L): sum(1 for g in genes_all
                                   if pi[L].get(g, 0) > 0 and g not in lib_targets) for L in LS},
        inputs={str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha256_file(p))
                for p in (E301 / "promoters_cage_tss_hg19.fa.gz", E301 / "vcc2025_spacers.csv",
                          E301 / "horlbeck_supp3.pkl", E298 / "seed298.csv")},
        seconds=round(time.time() - t0, 1),
    )
    (HERE / "seq307_meta.json").write_text(json.dumps(meta, indent=1))
    (HERE / "seq307.log").write_text("\n".join(log) + "\n")
    # pickled lookups for the null builder (sequence side only)
    with open(HERE / "seq307_lookup.pkl", "wb") as fh:
        pickle.dump(dict(lmer_genes={L: {m: dict(v) for m, v in lmer_genes[L].items()}
                                     for L in (12,)},
                         pi=dict(pi), lib_targets=sorted(lib_targets),
                         pool=pool, cand_L={c: {L: dict(v) for L, v in d.items()}
                                            for c, d in cand_L.items()}), fh)
    P(f"done in {meta['seconds']}s")


if __name__ == "__main__":
    main()
