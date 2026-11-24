"""301 item 4 (T7): seed-match scan of the 2025 library's two-construct arms against promoter windows.

Recipe = experiments/298-reviewer-checks/work/seed298_plan.md steps 1-7 (steps 2 and 4 unblocked by
301 items 1 and 2).

For every arm of the 17 two-construct targets (seed298.csv) take its two protospacers
(`vcc2025_spacers.csv`) and look for exact matches of the PAM-proximal seed (the 3' L nt of the
protospacer, L in {12,11,10,9}) inside a TSS +/- 500 bp window (`promoters_tss1000_grch38.fa.gz`),
on either strand, with a canonical NGG PAM immediately 3' of the matched site in protospacer
orientation. Hits in the arm's own target gene's promoter are excluded from the candidate list.
The scan is PAM-anchored: enumerate NGG sites, then test the L nt 5' of each against the seed set.

Known-answer row (KA-1): the same machinery, WITHOUT the own-target exclusion, must recover each
spacer's intended on-target site -- a full-length match of the 3' 19 nt of the protospacer (the 5'
base of a CRISPRi-v2 spacer is a cloning G that need not be genomic) with an NGG PAM, inside its own
target gene's TSS +/- 500 window. This validates the spacer table, the promoter windows, the strand
handling and the PAM rule at once.
Known-answer row (KA-2): the failing/working split read back from seed298.csv must be 13 / 21.
Background group: the 31 non-targeting constructs' 62 spacers scanned identically.
"""
import csv, gzip, collections, json, sys

FASTA = sys.argv[1] if len(sys.argv) > 1 else 'promoters_cage_tss_hg19.fa.gz'
TAG = sys.argv[2] if len(sys.argv) > 2 else 'cage'

COMP = str.maketrans('ACGT', 'TGCA')
def rc(s): return s.translate(COMP)[::-1]
LS = (12, 11, 10, 9)

prom = []
with gzip.open(FASTA, 'rt') as fh:
    hdr = None
    for line in fh:
        line = line.rstrip('\n')
        if line.startswith('>'): hdr = line[1:].split('|')
        else: prom.append((hdr[0], hdr[1], hdr[2], int(hdr[3]), int(hdr[4]), int(hdr[5]),
                           int(hdr[6]), line.upper()))
by_sym = collections.defaultdict(list)
for p in prom: by_sym[p[0]].append(p)
print(f'{FASTA}: promoter windows {len(prom)}  genes {len(by_sym)}  '
      f'bp {sum(len(p[7]) for p in prom)}', flush=True)

spac = {}
for r in csv.DictReader(open('vcc2025_spacers.csv')):
    if r['route'] == 'UNRESOLVED': continue
    spac[r['guide_id']] = r
    if not r['guide_id'].startswith('non-targeting'):
        spac[r['guide_id'].split('|')[0][:-2]] = r
arms = list(csv.DictReader(open('~/Code/virtual-cell-challenge-2026/experiments/'
                                '298-reviewer-checks/work/seed298.csv')))
n_fail = sum(1 for a in arms if a['failing'] == 'True')
print(f'KA-2 failing/working arms from seed298.csv: {n_fail} / {len(arms) - n_fail}'
      f'   (expected 13 / 21)  ->  {"PASS" if (n_fail, len(arms) - n_fail) == (13, 21) else "FAIL"}',
      flush=True)

queries = []       # (group, construct, own_target, failing, which, spacer)
for a in arms:
    r = spac[a['construct']]
    for w in ('A', 'B'):
        queries.append(('arm', a['construct'], a['target'], a['failing'] == 'True', w, r['spacer_' + w]))
ntcs = sorted(g for g in spac if g.startswith('non-targeting'))
for g in ntcs:
    r = spac[g]
    for w in ('A', 'B'):
        queries.append(('ntc', g, None, False, w, r['spacer_' + w]))
print(f'arms: {len(arms)}  ntc constructs: {len(ntcs)}  spacers scanned: {len(queries)}', flush=True)

# ---- KA-1: on-target recovery, own promoter, 19 nt + NGG ------------------
ka = {'checked': 0, 'recovered': 0, 'no_promoter': 0, 'missed': []}
for group, construct, tgt, failing, w, sp in queries:
    if group != 'arm': continue
    ka['checked'] += 1
    ps = by_sym.get(tgt)
    if not ps: ka['no_promoter'] += 1; continue
    core = sp[-19:]; ok = False
    for p in ps:
        for s in (p[7], rc(p[7])):
            i = s.find(core)
            while i >= 0:
                if s[i + 20:i + 22] == 'GG': ok = True; break
                i = s.find(core, i + 1)
            if ok: break
        if ok: break
    ka['recovered'] += ok
    if not ok: ka['missed'].append(f'{construct}_{w}')
print(f'KA-1 on-target recovery: {ka["recovered"]}/{ka["checked"]} arm spacers matched their own '
      f'target promoter at 19 nt + NGG (no promoter window for {ka["no_promoter"]})', flush=True)
print(f'   misses: {ka["missed"]}', flush=True)

# ---- the scan ------------------------------------------------------------
seeds = {L: collections.defaultdict(list) for L in LS}
for q in queries:
    for L in LS: seeds[L][q[5][-L:]].append(q)
hits = []
for sym, eid, chrom, strand, tss, lo, hi, seq in prom:
    for orient, s in (('+', seq), ('-', rc(seq))):
        start = 0
        while True:
            j = s.find('GG', start)
            if j < 0: break
            start = j + 1
            pam0 = j - 1
            if pam0 < LS[0]: continue
            for L in LS:
                for group, construct, tgt, failing, w, sp in seeds[L].get(s[pam0 - L:pam0], ()):
                    if sym == tgt: continue
                    gpos = (lo + pam0 - 1) if orient == '+' else (hi - (pam0 - 1))
                    off = (gpos - tss) if strand == 1 else (tss - gpos)
                    hits.append(dict(group=group, construct=construct, target=tgt or '', failing=failing,
                                     spacer=w, L=L, seed=s[pam0 - L:pam0], gene=sym, ensembl=eid,
                                     orient=orient, pam=s[pam0:pam0 + 3], offset_from_tss=off))
print(f'raw seed+PAM hits (own target excluded): {len(hits)}', flush=True)
cols = ['group', 'construct', 'target', 'failing', 'spacer', 'L', 'seed', 'gene', 'ensembl',
        'orient', 'pam', 'offset_from_tss']
with open(f'seed_hits301_{TAG}.csv', 'w', newline='') as fh:
    w_ = csv.DictWriter(fh, fieldnames=cols); w_.writeheader(); w_.writerows(hits)

# candidate = a gene hit at the LONGEST L that yields any hit for that construct (plan step 5)
best = {}
for h in hits: best[h['construct']] = max(best.get(h['construct'], 0), h['L'])
cand = collections.defaultdict(set)
for h in hits:
    if h['L'] == best[h['construct']]: cand[h['construct']].add(h['gene'])
fail = [a['construct'] for a in arms if a['failing'] == 'True']
work = [a['construct'] for a in arms if a['failing'] != 'True']
def rate(cs): return sum(1 for c in cs if cand.get(c)), len(cs)
# fixed-L candidate sets: a gene hit at exactly seed length L (cleaner than the plan's
# "longest L with any hit", whose denominator moves with the arm)
candL = {L: collections.defaultdict(set) for L in LS}
for h in hits: candL[h['L']][h['construct']].add(h['gene'])
def rateL(cs, L): return sum(1 for c in cs if candL[L].get(c)), len(cs)
fixed = {L: dict(failing=rateL(fail, L), working=rateL(work, L), ntc=rateL(ntcs, L),
                 mean_cand_failing=round(sum(len(candL[L].get(c, [])) for c in fail) / len(fail), 3),
                 mean_cand_working=round(sum(len(candL[L].get(c, [])) for c in work) / len(work), 3),
                 mean_cand_ntc=round(sum(len(candL[L].get(c, [])) for c in ntcs) / len(ntcs), 3))
         for L in LS}
summ = dict(ka1=ka, ka2=[n_fail, len(arms) - n_fail],
            failing=rate(fail), working=rate(work), ntc=rate(ntcs),
            hits_by_L=dict(collections.Counter(h['L'] for h in hits)),
            n_cand_failing={c: len(cand.get(c, [])) for c in fail},
            n_cand_working={c: len(cand.get(c, [])) for c in work},
            n_cand_ntc={c: len(cand.get(c, [])) for c in ntcs},
            fixed_L=fixed, longest_L=best,
            candidates={c: sorted(v) for c, v in cand.items()},
            candidates_L12={c: sorted(v) for c, v in candL[12].items()})
summ['fasta'] = FASTA
json.dump(summ, open(f'seed_scan301_{TAG}.json', 'w'), indent=1)
print(f'\nfailing arms with >=1 candidate : {rate(fail)[0]} / {rate(fail)[1]}')
print(f'working arms with >=1 candidate : {rate(work)[0]} / {rate(work)[1]}')
print(f'non-targeting  with >=1 candidate: {rate(ntcs)[0]} / {rate(ntcs)[1]}')
print('hits by seed length L:', dict(collections.Counter(h['L'] for h in hits)))
print('\nfixed-L candidate rates (constructs with >=1 candidate / n; mean candidates per construct):')
for L in LS:
    f = fixed[L]
    print(f'  L={L}: failing {f["failing"][0]}/{f["failing"][1]} ({f["mean_cand_failing"]}) | '
          f'working {f["working"][0]}/{f["working"][1]} ({f["mean_cand_working"]}) | '
          f'non-targeting {f["ntc"][0]}/{f["ntc"][1]} ({f["mean_cand_ntc"]})')
print()
for c in fail: print(f'  FAIL {c:32s} L*={best.get(c, "-")} n_cand={len(cand.get(c, []))} {sorted(cand.get(c, []))[:10]}')
for c in work: print(f'  work {c:32s} L*={best.get(c, "-")} n_cand={len(cand.get(c, []))} {sorted(cand.get(c, []))[:10]}')
