"""301 item 4, final reading: failing-vs-working contrast with (a) the CRISPRi active window
(TSS -50 .. +300) and (b) a DISTAL filter that drops hits whose candidate promoter window lies
within 10 kb of any window of the arm's own target gene -- those are the arm's own locus re-found
under a neighbouring gene's symbol (ZNF581 -> ZNF580 overlap at 0 bp), not off-targets."""
import csv, gzip, collections, json
from scipy.stats import fisher_exact, mannwhitneyu

wins = collections.defaultdict(list)
with gzip.open('promoters_cage_tss_hg19.fa.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            g, t, c, s, tss, lo, hi, n = line[1:].strip().split('|')
            wins[g].append((c, int(lo), int(hi)))
def distal(tgt, cg, kb=10):
    if not tgt: return True
    for c1, l1, h1 in wins.get(tgt, []):
        for c2, l2, h2 in wins.get(cg, []):
            if c1 == c2 and max(l2 - h1, l1 - h2, 0) < kb * 1000: return False
    return True

hits = list(csv.DictReader(open('seed_hits301_cage.csv')))
for h in hits: h['L'] = int(h['L']); h['off'] = int(h['offset_from_tss'])
arms = list(csv.DictReader(open('~/Code/virtual-cell-challenge-2026/experiments/'
                                '298-reviewer-checks/work/seed298.csv')))
fail = [a['construct'] for a in arms if a['failing'] == 'True']
work = [a['construct'] for a in arms if a['failing'] != 'True']
ntc = [r['guide_id'] for r in csv.DictReader(open('vcc2025_spacers.csv'))
       if r['guide_id'].startswith('non-targeting') and r['route'] != 'UNRESOLVED']
dropped = [h for h in hits if not distal(h['target'], h['gene'])]
print(f'hits dropped by the 10 kb self-locus filter: {len(dropped)} '
      f'({len({(h["construct"], h["gene"]) for h in dropped})} construct-gene pairs, '
      f'{sorted({h["gene"] for h in dropped})[:8]})')
hits = [h for h in hits if distal(h['target'], h['gene'])]

out = {}
for label, filt in [('whole_window_pm500', lambda h: True),
                    ('crispri_window_-50_+300', lambda h: -50 <= h['off'] <= 300)]:
    sub = [h for h in hits if filt(h)]
    print(f'\n=== {label} (distal only): {len(sub)} hits')
    res = {}
    for L in (12, 11, 10, 9):
        cand = collections.defaultdict(set)
        for h in sub:
            if h['L'] == L: cand[h['construct']].add(h['gene'])
        a = sum(1 for c in fail if cand.get(c)); b = len(fail) - a
        c_ = sum(1 for c in work if cand.get(c)); d = len(work) - c_
        orr, p = fisher_exact([[a, b], [c_, d]])
        nf = [len(cand.get(x, [])) for x in fail]; nw = [len(cand.get(x, [])) for x in work]
        nn = [len(cand.get(x, [])) for x in ntc]
        try: _, pu = mannwhitneyu(nf, nw, alternative='two-sided')
        except ValueError: pu = float('nan')
        print(f'  L={L:2d}: failing {a:2d}/{len(fail)}  working {c_:2d}/{len(work)}  '
              f'non-targeting {sum(1 for x in nn if x):2d}/{len(ntc)} | '
              f'Fisher OR={orr:.2f} p={p:.3f} | mean candidates F={sum(nf)/len(nf):.2f} '
              f'W={sum(nw)/len(nw):.2f} N={sum(nn)/len(nn):.2f} | MWU p={pu:.3f}')
        res[L] = dict(failing=[a, len(fail)], working=[c_, len(work)],
                      ntc=[sum(1 for x in nn if x), len(ntc)],
                      fisher_or=None if orr != orr else round(float(orr), 3),
                      fisher_p=round(float(p), 4), mwu_p=round(float(pu), 4),
                      mean_cand=dict(failing=round(sum(nf)/len(nf), 3),
                                     working=round(sum(nw)/len(nw), 3),
                                     ntc=round(sum(nn)/len(nn), 3)),
                      failing_candidates={x: sorted(cand.get(x, [])) for x in fail},
                      working_candidates={x: sorted(cand.get(x, [])) for x in work})
    out[label] = res
out['self_locus_hits_dropped'] = [[h['construct'], h['gene'], h['L'], h['off']] for h in dropped]
json.dump(out, open('seed_contrast301.json', 'w'), indent=1)
print('\nPRIMARY READING -- L=12 seed, CRISPRi active window, distal:')
for grp, cs in (('FAIL', fail), ('work', work)):
    for c in cs:
        g = out['crispri_window_-50_+300'][12][
            'failing_candidates' if grp == 'FAIL' else 'working_candidates'].get(c, [])
        if g: print(f'  {grp} {c:30s} -> {g}')
