"""386 audit: board rows recomputed from on-disk JSON with own code (rev 2)."""
import json, glob, os

latest = max(glob.glob('data/leaderboard/*.json'), key=os.path.getmtime)
assert latest.endswith('2026-09-11T1910Z.json'), latest
print('latest by mtime:', latest)
d = json.load(open(latest))
ents = d['entries']
print('n entries:', len(ents))
MEM = ['score_pds','score_mse','score_nmae','score_fid','score_reach','score_jac']

def mean6(e):
    return sum(e[k] for k in MEM)/6

# snapshot is one-row-per-team: check
print('distinct team_id:', len(set(e['team_id'] for e in ents)))

# Illumina + cqawesome from snapshot
ill = sorted([e for e in ents if e.get('team_name')=='Illumina AI'], key=lambda e:-e['score_avg'])[0]
cq = sorted([e for e in ents if e.get('team_name')=='cqawesome'], key=lambda e:-e['score_avg'])[0]
print('\nILL id=%s model=%s avg=%.7f (PI 0.2872)' % (ill['id'], ill.get('model_name'), ill['score_avg']))
print('CQ  id=%s model=%s avg=%.7f (PI 0.2506)' % (cq['id'], cq.get('model_name'), cq['score_avg']))
print('ill mean6-check diff=%.2e cq mean6-check diff=%.2e' % (mean6(ill)-ill['score_avg'], mean6(cq)-cq['score_avg']))
print('ill reach raw=%.4f (PI 0.367)' % ill['de_wilcoxon_direction_reach_raw'])

# ours from board-entries (not in snapshot: submitted after snapshot's last row)
p2 = json.load(open('data/board-entries/lqWNXX5YaYArdglu6Stn.json'))['entry']
p = json.load(open('data/board-entries/SvzrTOQJAAaNjawzSyrC.json'))['entry']
print('\n319p2 avg=%.7f (PI 0.2201159) mean6-check diff=%.2e' % (p2['score_avg'], mean6(p2)-p2['score_avg']))
print('319p  avg=%.7f (PI 0.2198190) mean6-check diff=%.2e' % (p['score_avg'], mean6(p)-p['score_avg']))
print('delta=%.6f (PI +0.000297)' % (p2['score_avg']-p['score_avg']))
for k in MEM:
    print('  d%-10s %+.6f' % (k, p2[k]-p[k]))

# gap decomp ill - p2
print('\nGAP (ill-p2)/6 vs PI [nmae .0224 reach .0170 mse .0156 fid .0079 pds .0067 jac -.0023]:')
for k in MEM:
    print('  %-10s %+.6f' % (k, (ill[k]-p2[k])/6))
print('  sum-check: %.7f vs avg-diff %.7f' % (sum((ill[k]-p2[k])/6 for k in MEM), ill['score_avg']-p2['score_avg']))

# rank plausibility
above = [e for e in ents if e['score_avg'] > p2['score_avg']]
print('\nsnapshot rows above 319p2: %d -> implied rank %d (PI rank 15 of 893)' % (len(above), len(above)+1))

# Illumina step 09-09 -> 09-11 (same team_id, different entry: rename builders->Illumina)
d0 = json.load(open('data/leaderboard/2026-09-09T2137Z.json'))
b0 = sorted([e for e in d0['entries'] if e['team_id']=='0b01OX4HQx2BB1tRiny3'], key=lambda e:-e['score_avg'])[0]
print('\n09-09 team-best: id=%s name=%s model=%s avg=%.7f' % (b0['id'], repr(b0.get('team_name')), b0.get('model_name'), b0['score_avg']))
print('09-11 team-best: id=%s name=%s model=%s avg=%.7f' % (ill['id'], repr(ill.get('team_name')), ill.get('model_name'), ill['score_avg']))
print('STEP (PI: mse +0.151 scaled, rest within +-0.011):')
for k in MEM + ['score_avg']:
    print('  %-10s 09-09=%+.6f 09-11=%+.6f step=%+.6f' % (k, b0[k], ill[k], ill[k]-b0[k]))
# team_id sharing
n0 = set(e.get('team_name') for e in d0['entries'] if e['team_id']=='0b01OX4HQx2BB1tRiny3')
n1 = set(e.get('team_name') for e in ents if e['team_id']=='0b01OX4HQx2BB1tRiny3')
print('team_id 0b01... names 09-09:', n0, ' 09-11:', n1)
