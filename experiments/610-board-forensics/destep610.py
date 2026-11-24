# 610: catalogue of DE-ONLY board steps (consecutive entries of one team with raw pds and mse unchanged < 2e-4), all teams,
# |d overall| >= 0.002, with member signature and the arriving entry's description (mechanism, where disclosed).
import json, sys, collections
es = json.load(open(sys.argv[1]))
K = dict(pds='pds_cosine', mse='expr_mse_unbiased_capped_norm', nmae='de_wilcoxon_lfc_nmae',
         fid='de_wilcoxon_direction_fidelity_yield_raw', reach='de_wilcoxon_direction_reach_raw', jac='de_wilcoxon_sig_jaccard')
SC = dict(nmae='score_nmae', fid='score_fid', reach='score_reach', jac='score_jac')
by = collections.defaultdict(list)
for e in es:
    if e.get('score_avg') is not None and e.get(K['fid']) is not None: by[e['team_name']].append(e)
rows = []
for t, L in by.items():
    L.sort(key=lambda e: e['submission_date'])
    for a, b in zip(L, L[1:]):
        if abs(a[K['pds']] - b[K['pds']]) < 2e-4 and abs(a[K['mse']] - b[K['mse']]) < 2e-4:
            dov = b['score_avg'] - a['score_avg']
            if abs(dov) < 0.002: continue
            d = {m: b[K[m]] - a[K[m]] for m in ('nmae', 'fid', 'reach', 'jac')}
            ds = {m: (b[SC[m]] - a[SC[m]]) / 6 for m in SC}   # overall-equivalent per member
            main = max(ds, key=lambda m: abs(ds[m]))
            rows.append((dov, t, a['submission_date'][:10], b['submission_date'][:10], d, ds, main, b.get('model_name') or '', (b.get('description') or '').replace('\n', ' '), a['score_avg']))
rows.sort(key=lambda r: -r[0])
print('DE-only steps with |d overall| >= 0.002: %d (positive %d), teams %d' % (len(rows), sum(r[0] > 0 for r in rows), len(set(r[1] for r in rows))))
print('\n## largest POSITIVE DE-only steps (overall-equivalent per member in brackets)')
for r in [x for x in rows if x[0] > 0][:45]:
    dov, t, d0, d1, d, ds, main, mn, desc, base = r
    print('%+.4f %-18s %s (base %.3f) main=%-5s | nmae %+.4f fid %+.4f reach %+.4f jac %+.4f [%s] | %s: %s' % (
        dov, t[:18], d1, base, main, d['nmae'], d['fid'], d['reach'], d['jac'],
        ' '.join('%s%+.4f' % (m[0], ds[m]) for m in ('nmae', 'fid', 'reach', 'jac')), mn[:28], desc[:330]))
sig = collections.Counter(r[6] for r in rows if r[0] > 0)
print('\nmain member of positive DE-only steps:', dict(sig))
for m in ('nmae', 'fid', 'reach', 'jac'):
    xs = sorted([r[5][m] for r in rows if r[0] > 0 and r[6] == m], reverse=True)
    if xs: print('  %-5s n=%3d median %+.4f  top5 %s' % (m, len(xs), xs[len(xs)//2], ' '.join('%+.4f' % x for x in xs[:5])))
