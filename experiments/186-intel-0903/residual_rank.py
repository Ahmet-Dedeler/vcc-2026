"""The ordering-channel residual, on the coverage-saturated cohort.

Restrict to teams whose `fid` sits in [0.503, 0.512]: there `fid` is pure
directional precision because the yield cap min(1, n_pred/n_real) is saturated,
so differences in `jac` are call-set differences and not precision differences.
Then regress raw reach on everything else the board shows, and read the residual
— what is left is the p-value ORDERING channel.

Known-answer row: our own residual must come out at -0.0411 raw (-1.92 sd), the
number quoted in RESULTS.md section 2 / M1(iii).
"""
import json

import numpy as np

K = dict(pds='pds_cosine', mse='expr_mse_unbiased_capped_norm',
         nmae='de_wilcoxon_lfc_nmae',
         fid='de_wilcoxon_direction_fidelity_yield_raw',
         reach='de_wilcoxon_direction_reach_raw', jac='de_wilcoxon_sig_jaccard')

D = [e for e in json.load(open('board/leaderboard_2026-09-03.json'))['entries']
     if e.get('score_avg') is not None]
S = [e for e in D if 0.503 <= e[K['fid']] <= 0.512 and e[K['pds']] > 0.79]

r = np.array([e[K['reach']] for e in S])
p = np.array([e[K['pds']] for e in S])
n = np.array([e[K['nmae']] for e in S])
m = np.log(np.array([e[K['mse']] for e in S]))
j = np.array([e[K['jac']] for e in S])

print(f'coverage-saturated cohort: n = {len(S)}')
for name, x in (('jac', j), ('log mse', m), ('pds', p), ('nmae', n)):
    print(f'  corr(reach, {name:<7}) = {np.corrcoef(r, x)[0,1]:+.3f}')

A = np.vstack([np.ones(len(S)), p, n, m, j]).T
beta, *_ = np.linalg.lstsq(A, r, rcond=None)
res = r - A @ beta
print(f'\nreach ~ 1 + pds + nmae + log mse + jac :  pds {beta[1]:+.3f}  '
      f'nmae {beta[2]:+.3f}  log mse {beta[3]:+.4f}  jac {beta[4]:+.2f}')
print(f'  R2 = {1 - res.var()/r.var():.3f}   resid sd = {res.std():.4f}')

i_us = [k for k, e in enumerate(S) if e['team_name'] == 'Ahmet Dedeler'][0]
print(f'\nKNOWN ANSWER — our residual: {res[i_us]:+.4f} raw reach '
      f'({res[i_us]/res.std():+.2f} sd) = {res[i_us]*1.1292/6:+.5f} overall')
assert abs(res[i_us] + 0.0411) < 5e-4, res[i_us]
print('  matches RESULTS.md M1(iii).  OK')

o = np.argsort(res)
print(f'\nMOST NEGATIVE (ordering deficit)          rank of {len(S)}')
for k, i in enumerate(o[:5], 1):
    star = '  <== us' if S[i]['team_name'] == 'Ahmet Dedeler' else ''
    print(f'  {k}. {S[i]["team_name"][:22]:<22} `{str(S[i]["model_name"])[:24]:<24}` '
          f'resid {res[i]:+.4f} ({res[i]/res.std():+.2f} sd)  reach {r[i]:.4f}  '
          f'ovr {S[i]["score_avg"]:.4f}{star}')
print('MOST POSITIVE (ordering surplus)')
for k, i in enumerate(o[::-1][:5], 1):
    print(f'  {k}. {S[i]["team_name"][:22]:<22} `{str(S[i]["model_name"])[:24]:<24}` '
          f'resid {res[i]:+.4f} ({res[i]/res.std():+.2f} sd)  reach {r[i]:.4f}  '
          f'ovr {S[i]["score_avg"]:.4f}')

print('\nEXCHANGE RATE along the field\'s own jac<->reach line '
      '(anchors: reach span 0.8857, jac span 0.3668)')
d_jac = 0.005
print(f'  give up {d_jac} raw jac    = {-d_jac*2.7263/6:+.5f} overall')
print(f'  buys {-beta[4]*d_jac:+.4f} raw reach = {-beta[4]*d_jac*1.1292/6:+.5f} overall')
print(f'  NET                      = '
      f'{(-beta[4]*d_jac*1.1292 - d_jac*2.7263)/6:+.5f} overall  '
      '-> the trade line is worth ~nothing; the prize is the RESIDUAL')
