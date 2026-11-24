import json, glob, collections, math, os
BASE='~/Code/virtual-cell-challenge-2026/'
FILES = sorted(glob.glob(BASE+'data/leaderboard/lb_*.json'))
SNAPS=[]
for f in FILES:
    tag=os.path.basename(f).replace('lb_','').replace('.json','')
    d=json.load(open(f))['entries']
    SNAPS.append((tag,{x['team_name']:x for x in d}))
SC=['score_pds','score_mse','score_nmae','score_fid','score_reach','score_jac']
RAW=['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard']
CUR=SNAPS[-1][1]
def g(e,k): 
    v=e.get(k); return 0.0 if v is None else float(v)
