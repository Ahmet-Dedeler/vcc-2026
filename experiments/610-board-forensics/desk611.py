# 611 desk check (council 600, 2026-10-02): board vs 537's H1 regime prices for the board DE pairs that map onto 537 modes
# without leakage. Reads only data/board-entries/*.json and experiments/537-arm-pipeline/work/price537.json.
import json
RAW = dict(fid='de_wilcoxon_direction_fidelity_yield_raw', jac='de_wilcoxon_sig_jaccard', reach='de_wilcoxon_direction_reach_raw', nmae='de_wilcoxon_lfc_nmae')
def ent(i): return json.load(open(f'data/board-entries/{i}.json'))['entry']
P = json.load(open('experiments/537-arm-pipeline/work/price537.json'))['comparisons']
pairs = [('cp block (446bj2zcp - 446bj2z)', 'K1Q0y2OFYjTNRqh9jKfk', 'be90smZUDcIV7j5MgQ03', 'REF cp (528, kept) on 446bj2z')]
for name, a, b, key in pairs:
    A, B = ent(a), ent(b)
    bd = {m: A[RAW[m]] - B[RAW[m]] for m in RAW}
    print(name)
    print('  board raw delta: ' + '  '.join('%s %+.5f' % (m, bd[m]) for m in RAW))
    for reg in ('1.0', '0.5', '0.35', 'cmp'):
        h = P[key][reg]['d_mean']
        rat = {m: (bd[m] / h[m] if abs(h[m]) > 1e-9 else float('nan')) for m in RAW}
        sgn = sum(1 for m in RAW if abs(h[m]) > 1e-4 and (bd[m] > 0) == (h[m] > 0))
        print('  H1 f=%-4s delta: %s | board/H1: %s | sign-agree %d/%d' % (reg, '  '.join('%s %+.5f' % (m, h[m]) for m in RAW),
              ' '.join('%s %+.2f' % (m, rat[m]) for m in RAW), sgn, sum(1 for m in RAW if abs(h[m]) > 1e-4)))
print('zs post-pass (446bj2z - 446bj2): board nmae %+.6f, other members exactly 0; H1 (435) predicted nmae -0.0013 ... -0.0029 -> ratio ~0.75 ... 1.7, inside' % (ent('be90smZUDcIV7j5MgQ03')[RAW['nmae']] - ent('YvVru3DAF89Jiyr34olR')[RAW['nmae']]))
