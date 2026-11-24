"""181 §board — recompute the 2026-09-02 board decomposition from the raw snapshot.

Known-answer rows (rule 10): leader cqawesome overall 0.21687, our a23 entry
0.1733597079 must both appear verbatim in the file; the six scaled members of
each entry must average to its own score_avg.
"""
import json, numpy as np

d = json.load(open('data/leaderboard/2026-09-02T0630Z.json'))
E = d['entries']
MEM = ['score_pds','score_mse','score_nmae','score_fid','score_reach','score_jac']
RAW = {'score_pds':'pds_cosine','score_mse':'expr_mse_unbiased_capped_norm',
       'score_nmae':'de_wilcoxon_lfc_nmae','score_fid':'de_wilcoxon_direction_fidelity_yield_raw',
       'score_reach':'de_wilcoxon_direction_reach_raw','score_jac':'de_wilcoxon_sig_jaccard'}

# --- KA gate: members average to overall
bad = 0
for e in E:
    m = [e[k] for k in MEM]
    if any(v is None for v in m):
        bad += 1; continue
    if abs(np.mean(m) - e['score_avg']) > 2e-9:
        bad += 1
print(f'KA rule10: entries whose 6 scaled members do NOT average to score_avg: {bad} / {len(E)}')

# --- find us
us = [e for e in E if 'Reviewer 2' in (e['model_name'] or '') or abs(e['score_avg']-0.1733597079)<1e-9]
print('our entries found:', [(e['team_name'], e['model_name'], e['score_avg']) for e in us])
me = us[0]
lead = max(E, key=lambda e: e['score_avg'])
print(f"leader: {lead['team_name']} {lead['model_name']} {lead['score_avg']:.10f}")
print(f"us:     {me['team_name']} {me['model_name']} {me['score_avg']:.10f}")
print(f"our overall rank: {1+sum(1 for e in E if e['score_avg']>me['score_avg'])} of {len(E)}")

print()
print(f"{'member':6} {'us':>9} {'leader':>9} {'edge/6':>10} {'share':>8} {'ourrank':>8} {'fieldbest':>10} {'headroom':>9}")
gap = lead['score_avg']-me['score_avg']
rows=[]
for k in MEM:
    vals = np.array([e[k] for e in E if e[k] is not None])
    rank = 1+sum(1 for v in vals if v > me[k])
    best = vals.max()
    edge = (lead[k]-me[k])/6
    rows.append((k, me[k], lead[k], edge, edge/gap, rank, best, (best-me[k])/6))
    print(f"{k[6:]:6} {me[k]:9.4f} {lead[k]:9.4f} {edge:10.5f} {edge/gap:8.1%} {rank:8d} {best:10.4f} {(best-me[k])/6:9.4f}")
print(f"gap to leader: {gap:.5f}   sum of edges: {sum(r[3] for r in rows):.5f}")

# --- field regression raw reach vs raw pds
ok = [e for e in E if e['pds_cosine'] is not None and e['de_wilcoxon_direction_reach_raw'] is not None]
x = np.array([e['pds_cosine'] for e in ok]); y = np.array([e['de_wilcoxon_direction_reach_raw'] for e in ok])
b1, b0 = np.polyfit(x, y, 1)
resid = y - (b0+b1*x)
print()
print(f'n={len(ok)}  corr(raw_reach, raw_pds) = {np.corrcoef(x,y)[0,1]:.4f}')
print(f'raw_reach = {b0:.4f} + {b1:.4f}*raw_pds ; residual sd {resid.std(ddof=2):.4f}')
mx, my = me['pds_cosine'], me['de_wilcoxon_direction_reach_raw']
pred = b0+b1*mx
print(f'us: raw pds {mx:.4f}, raw reach {my:.4f}; field predicts {pred:.4f}; residual {(my-pred):+.4f} = {(my-pred)/resid.std(ddof=2):+.2f} sd')
below = sum(1 for r in resid if r < (my-pred))
print(f'teams further below the line: {below} of {len(ok)} = {below/len(ok):.1%}')
print(f'closing to the field line on reach = {(pred-my)*0.895/6:+.5f} overall  (span 0.895)')

# --- teams with raw pds >= ours
print()
print('teams with raw pds >= ours:')
for e in sorted([e for e in ok if e['pds_cosine']>=mx], key=lambda e:-e['pds_cosine']):
    print(f"  {e['team_name'][:22]:22} pds {e['pds_cosine']:.4f} reach {e['de_wilcoxon_direction_reach_raw']:.4f} mse {e['expr_mse_unbiased_capped_norm']:.3f} overall {e['score_avg']:.4f}")

# --- call-set vs profile partition: is nmae with the profile group?
print()
print('rank distribution check — the claimed 5/23/68 vs 176/169/192 partition:')
for k in MEM:
    vals = np.array([e[k] for e in E if e[k] is not None])
    print(f'  {k[6:]:6} rank {1+sum(1 for v in vals if v>me[k]):3d}  pctile {(1-(sum(1 for v in vals if v>me[k]))/len(vals)):.3f}')
