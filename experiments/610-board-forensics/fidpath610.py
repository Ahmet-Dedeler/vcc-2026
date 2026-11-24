# 610: fid trajectory of the fid-attractor teams (first entry, peak, largest single-step jumps, with descriptions).
import json, sys, collections
es=json.load(open(sys.argv[1]))
F='de_wilcoxon_direction_fidelity_yield_raw'; J='de_wilcoxon_sig_jaccard'; R='de_wilcoxon_direction_reach_raw'; N='de_wilcoxon_lfc_nmae'; PD='pds_cosine'
by=collections.defaultdict(list)
for e in es:
    if e.get(F) is not None: by[e['team_name']].append(e)
cellpred=sorted(t for t,L in by.items() if any((x.get('model_name') or '').startswith('cellpred') for x in L))
teams=cellpred+['Testh','sunny','Pythia','The_Way','CellWM','chill','cqawesome','nailong king']
print('cellpred-model teams:', cellpred)
for t in teams:
    L=sorted(by[t], key=lambda e:e['submission_date'])
    if not L: continue
    pk=max(L,key=lambda e:e[F])
    first_hi=next((e for e in L if e[F]>=0.525), None)
    print('\n== %s: %d entries; first %s fid %.4f (%s); peak %.4f %s; first >=0.525: %s'%(t,len(L),L[0]['submission_date'][:10],L[0][F],(L[0].get('model_name') or '')[:22],pk[F],pk['submission_date'][:10],(first_hi['submission_date'][:16]+' '+str(first_hi.get('model_name'))[:24]) if first_hi else 'never'))
    steps=[(b[F]-a[F],a,b) for a,b in zip(L,L[1:])]
    for dF,a,b in sorted(steps,key=lambda s:-s[0])[:3]:
        if dF<0.003: continue
        print('   step %s->%s dfid %+.4f (fid %.4f->%.4f) djac %+.4f dreach %+.4f dnmae %+.4f dpds %+.4f | %s'%(a['submission_date'][5:16],b['submission_date'][5:16],dF,a[F],b[F],b[J]-a[J],b[R]-a[R],b[N]-a[N],b[PD]-a[PD],((b.get('model_name') or '')+': '+(b.get('description') or ''))[:260].replace('\n',' ')))
