"""Board-level natural experiments: consecutive submissions by the SAME team.

The board is the only surface validated for `reach`.  Every team that submits
twice runs a paired experiment on it for us.  This script reconstructs each
team's submission history from all 14 leaderboard snapshots in
`data/leaderboard/` and extracts the pairs that isolate one factor:

  * AMPLITUDE-ONLY-ish : raw mse moves a lot, raw pds barely moves.
    -> reads d(reach)/d(log mse) with the direction held fixed.
  * ORDERING-ONLY-ish  : pds and mse both barely move, reach moves a lot.
    -> proves reach is reachable without touching the profile, and names who did it.

Known-answer rows: our own a15->abc-knee pair must show d reach = +0.0031 and
d pds = +0.0001 (BOARD-RESULT.md), and jw_13288's `pm-r14-pds3o` must show the
mse-101 / reach-0.328 blowup recorded in board/trajectories.md.
"""
import glob
import json
import os

import numpy as np

ROOT = '../../data/leaderboard'
K = dict(pds='pds_cosine', mse='expr_mse_unbiased_capped_norm',
         nmae='de_wilcoxon_lfc_nmae', fid='de_wilcoxon_direction_fidelity_yield_raw',
         reach='de_wilcoxon_direction_reach_raw', jac='de_wilcoxon_sig_jaccard')
MEMBERS = list(K)

files = sorted(glob.glob(os.path.join(ROOT, '*.json')))
files.append('board/leaderboard_2026-09-03.json')

hist = {}
for f in files:
    d = json.load(open(f))
    E = d['entries'] if isinstance(d, dict) else d
    for e in E:
        if e.get('score_avg') is None:
            continue
        key = (e['team_name'], e['submission_date'][:19], e['submission_count'])
        hist.setdefault(e['team_name'], {})[key] = e

print(f'{len(files)} snapshots, {len(hist)} teams')

# one row per distinct (team, submission_date, submission_count)
pairs = []
for team, subs in hist.items():
    rows = sorted(subs.values(), key=lambda e: (e['submission_count'], e['submission_date']))
    for a, b in zip(rows, rows[1:]):
        if a['submission_date'] == b['submission_date']:
            continue
        d = {m: b[K[m]] - a[K[m]] for m in MEMBERS}
        pairs.append(dict(team=team, a=a, b=b, d=d))
print(f'{len(pairs)} consecutive submission pairs\n')


def fmt(p):
    a, b, d = p['a'], p['b'], p['d']
    return (f"{p['team'][:20]:<20} {str(a['model_name'])[:16]:<16}->"
            f"{str(b['model_name'])[:16]:<16} "
            f"pds {a[K['pds']]:.4f}->{b[K['pds']]:.4f} ({d['pds']:+.4f})  "
            f"mse {a[K['mse']]:6.3f}->{b[K['mse']]:6.3f}  "
            f"reach {a[K['reach']]:.4f}->{b[K['reach']]:.4f} ({d['reach']:+.4f})  "
            f"jac {d['jac']:+.4f} fid {d['fid']:+.4f} nmae {d['nmae']:+.4f}")


# ---------------------------------------------------------------- known answers
print('KNOWN-ANSWER ROWS')
ours = [p for p in pairs if p['team'] == 'Ahmet Dedeler'
        and p['b']['model_name'].startswith('We regret') and p['b']['submission_count'] == 25]
for p in ours:
    print('  ' + fmt(p))
    assert abs(p['d']['reach'] - 0.0031) < 6e-4, p['d']['reach']
    assert abs(p['d']['pds'] - 0.0001) < 6e-4, p['d']['pds']
print('  -> abc-knee reproduces BOARD-RESULT.md (+0.0031 reach, +0.0001 pds).  OK')
blow = [p for p in pairs if p['team'] == 'jw_13288' and 'pds3o' in str(p['b']['model_name'])]
for p in blow:
    print('  ' + fmt(p))
print()

# --------------------------------------------- AMPLITUDE-ONLY natural experiments
print('=' * 118)
print('AMPLITUDE ARMS: |d pds| <= 0.010 and mse ratio >= 1.8 either way  '
      '(direction ~fixed, amplitude moved)')
print('=' * 118)
amp = []
for p in pairs:
    a, b = p['a'], p['b']
    if a[K['mse']] <= 0 or b[K['mse']] <= 0:
        continue
    ratio = b[K['mse']] / a[K['mse']]
    if abs(p['d']['pds']) <= 0.010 and (ratio >= 1.8 or ratio <= 1 / 1.8) \
            and max(a[K['pds']], b[K['pds']]) > 0.70:
        p['ratio'] = ratio
        amp.append(p)
amp.sort(key=lambda p: p['ratio'])
for p in amp:
    print(f"  x{p['ratio']:5.2f} mse | " + fmt(p))
if amp:
    lr = np.array([np.log(p['ratio']) for p in amp])
    dr = np.array([p['d']['reach'] for p in amp])
    dj = np.array([p['d']['jac'] for p in amp])
    df = np.array([p['d']['fid'] for p in amp])
    A = np.vstack([lr, np.ones_like(lr)]).T
    for name, y in (('reach', dr), ('jac', dj), ('fid', df)):
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ beta
        print(f'\n  d {name:<5} = {beta[0]:+.4f} * d log(mse) {beta[1]:+.4f}   '
              f'(n={len(amp)}, resid sd {r.std():.4f}, '
              f'corr {np.corrcoef(lr, y)[0,1]:+.3f})')
    print(f'  -> halving raw mse (d log mse = -0.69) costs '
          f'{-0.69*np.linalg.lstsq(A, dr, rcond=None)[0][0]:+.4f} raw reach')

# ------------------------------------------ ORDERING-ONLY natural experiments
print()
print('=' * 118)
print('ORDERING ARMS: |d pds| <= 0.006 and |d log mse| <= 0.12 and |d reach| >= 0.020')
print('  (profile essentially unchanged, reach moved anyway -> the p-ORDER moved)')
print('=' * 118)
ordr = []
for p in pairs:
    a, b = p['a'], p['b']
    if a[K['mse']] <= 0 or b[K['mse']] <= 0:
        continue
    if abs(p['d']['pds']) <= 0.006 and abs(np.log(b[K['mse']] / a[K['mse']])) <= 0.12 \
            and abs(p['d']['reach']) >= 0.020 and max(a[K['pds']], b[K['pds']]) > 0.70:
        ordr.append(p)
ordr.sort(key=lambda p: -p['d']['reach'])
for p in ordr:
    print('  ' + fmt(p))

print()
print('=' * 118)
print('BIGGEST reach GAINS anywhere on the board, any pair, pds>0.75 after')
print('=' * 118)
big = sorted([p for p in pairs if p['b'][K['pds']] > 0.75],
             key=lambda p: -p['d']['reach'])[:20]
for p in big:
    print('  ' + fmt(p))
