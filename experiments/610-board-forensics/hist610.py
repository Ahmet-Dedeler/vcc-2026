import json, glob, os, sys
# 610 (board forensics, 2026-10-02): every distinct leaderboard entry across data/leaderboard/*.json + a live snapshot,
# per team, with scaled-member deltas between consecutive entries.  Usage: python3 hist610.py <dir-with-lb_now.json> "Team A,Team B"
# Writes <dir>/all_entries.json (all distinct entries).  Read-only on the repo.
S=sys.argv[1]
files=sorted(glob.glob('data/leaderboard/*.json'))+[S+'/lb_now.json']
keys=['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard']
sk=['score_pds','score_mse','score_nmae','score_fid','score_reach','score_jac']
seen={}
for f in files:
    try: d=json.load(open(f))
    except Exception as ex: continue
    ents=d['entries'] if isinstance(d,dict) and 'entries' in d else d
    for e in ents:
        if not isinstance(e,dict) or e.get('score_avg') is None: continue
        k=e['id']
        if k not in seen: seen[k]=dict(e); seen[k]['_first_seen']=os.path.basename(f)
json.dump(list(seen.values()), open(S+'/all_entries.json','w'))
teams={}
for e in seen.values(): teams.setdefault(e['team_name'],[]).append(e)
print('distinct entries', len(seen), 'teams', len(teams))
want=sys.argv[2].split(',')
for t in want:
    es=sorted(teams.get(t,[]), key=lambda e:e['submission_date'])
    print('\n=== %s (%d entries seen)'%(t,len(es)))
    prev=None
    for e in es:
        raw=' '.join('%.4f'%e[k] for k in keys)
        d=''
        if prev:
            d=' Δsc '+' '.join('%+.3f'%(e[k]-prev[k]) for k in sk)
        print('%s %.4f %-24s | %s%s | %s'%(e['submission_date'][:16], e['score_avg'], (e.get('model_name') or '')[:24], raw, d, (e.get('description') or '')[:150].replace('\n',' ')))
        prev=e
