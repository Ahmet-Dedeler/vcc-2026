"""301 item 1: resolve the VCC 2025 obs/guide_id names to protospacer sequences.
Source: Replogle et al. 2022 Cell Table S1 (mmc1.xlsx), tab TabA_K562_day8_library -- the
genome-scale dual-guide CRISPRi library whose 'protospacer sequences and cloning strategy'
Arc reused for the 2025 challenge library (Arc, 'Behind the Data of the Virtual Cell Challenge')."""
import csv, pickle, openpyxl, collections

wb = openpyxl.load_workbook('replogle2022_mmc1.xlsx', read_only=True)
rows = list(wb['TabA_K562_day8_library'].iter_rows(values_only=True))
ix = {h: i for i, h in enumerate(rows[0])}
by_key, by_sgid_pair = {}, {}
for r in rows[1:]:
    rec = dict(pair_id=r[ix['unique sgRNA pair ID']], gene=r[ix['gene']], transcript=r[ix['transcript']],
               ensembl=r[ix['ensembl gene id']], sgID_A=r[ix['sgID_A']], seq_A=r[ix['targeting sequence A']],
               sgID_B=r[ix['sgID_B']], seq_B=r[ix['targeting sequence B']])
    by_key.setdefault(f"{rec['gene']}_{rec['transcript']}", rec)
    by_sgid_pair[f"{rec['sgID_A']}|{rec['sgID_B']}"] = rec

guides = [l.strip() for l in open(
    '/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/'
    '42edd049-4f44-4b3c-b669-9213db727327/scratchpad/guide_ids.txt') if l.strip()]
out, stats = [], collections.Counter()
for g in guides:
    a, b = g.split('|')
    rec = by_sgid_pair.get(g)
    route = 'sgID_pair'
    if rec is None and a.endswith('_A') and b.endswith('_B') and a[:-2] == b[:-2]:
        rec = by_key.get(a[:-2]); route = 'gene_transcript'
    if rec is None:
        stats['UNRESOLVED'] += 1
        out.append(dict(guide_id=g, sgid_A=a, sgid_B=b, gene='', transcript='', ensembl='',
                        spacer_A='', spacer_B='', route='UNRESOLVED', pair_id=''))
        continue
    stats[route] += 1
    out.append(dict(guide_id=g, sgid_A=rec['sgID_A'], sgid_B=rec['sgID_B'], gene=rec['gene'],
                    transcript=rec['transcript'], ensembl=rec['ensembl'], spacer_A=rec['seq_A'],
                    spacer_B=rec['seq_B'], route=route, pair_id=rec['pair_id']))
with open('vcc2025_spacers.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
res = [o for o in out if o['route'] != 'UNRESOLVED']
targ = [o for o in out if not o['guide_id'].startswith('non-targeting')]
ntc = [o for o in out if o['guide_id'].startswith('non-targeting')]
print('guide_ids:', len(out), '| resolved:', len(res), '| routes:', dict(stats))
print('  targeting  :', sum(1 for o in targ if o['route'] != 'UNRESOLVED'), '/', len(targ))
print('  non-target :', sum(1 for o in ntc if o['route'] != 'UNRESOLVED'), '/', len(ntc))
lens = collections.Counter(len(o['spacer_A']) for o in res) + collections.Counter(len(o['spacer_B']) for o in res)
print('  spacer lengths:', dict(lens))
bad = [o for o in res if set(o['spacer_A']) - set('ACGT') or set(o['spacer_B']) - set('ACGT')]
print('  non-ACGT spacers:', len(bad))
print('  unresolved      :', [o['guide_id'] for o in out if o['route'] == 'UNRESOLVED'])
# the 17 two-construct targets of 293/298
arms = [l.split(',')[3] for l in open('~/Code/virtual-cell-challenge-2026/experiments/'
                                      '298-reviewer-checks/work/seed298.csv').readlines()[1:]]
have = {o['guide_id'].split('|')[0][:-2]: o for o in res}
print('  seed298 arms with spacers:', sum(1 for a in arms if a in have), '/', len(arms))
print('  missing arms:', [a for a in arms if a not in have])
