"""What moved on the board between 2026-08-28 13:25 UTC and 2026-09-03 00:30 UTC,
and whether the reach spread inside the shrunk cluster is explained by sign
accuracy `q` (the nmae channel) or by the ordering channel.

nmae identity (docs/CURRENT_STATE.md 2026-09-03 02:15):
    nmae(k, q) = q|k-1| + (1-q)(k+1)
    k <= 1  ->  q = (1 + k - nmae) / (2k)
    k >= 1  ->  q = (k + 1 - nmae) / 2
Known-answer row: our raw nmae 0.9432 at k = 0.75 must give q = 0.537.
"""
import json
import numpy as np

K = dict(pds='pds_cosine', mse='expr_mse_unbiased_capped_norm',
         nmae='de_wilcoxon_lfc_nmae', fid='de_wilcoxon_direction_fidelity_yield_raw',
         reach='de_wilcoxon_direction_reach_raw', jac='de_wilcoxon_sig_jaccard')

old = json.load(open('../../data/leaderboard/lb_20260828T1325.json'))
old = old['entries'] if isinstance(old, dict) else old
new = json.load(open('board/leaderboard_2026-09-03.json'))['entries']
old = [e for e in old if e.get('score_avg') is not None]
new = [e for e in new if e.get('score_avg') is not None]
om = {e['team_name']: e for e in old}

print(f'teams: {len(old)} (Aug 28 13:25) -> {len(new)} (Sep 3 00:30), '
      f'+{len(new)-len(old)}')
so = sorted([e['score_avg'] for e in old])[::-1]
sn = sorted([e['score_avg'] for e in new])[::-1]
for lbl, i in (('1st', 0), ('5th', 4), ('10th', 9), ('20th', 19), ('50th', 49)):
    print(f'  {lbl:>4}: {so[i]:.5f} -> {sn[i]:.5f}  ({sn[i]-so[i]:+.5f})')
print(f'  teams over 0.170: {sum(1 for s in so if s>0.17)} -> {sum(1 for s in sn if s>0.17)}')
print(f'  teams with raw reach >= 0.25: '
      f"{sum(1 for e in old if e[K['reach']]>=0.25)} -> "
      f"{sum(1 for e in new if e[K['reach']]>=0.25)}")
print(f'  teams with raw mse < 1.0: '
      f"{sum(1 for e in old if e[K['mse']]<1.0)} -> "
      f"{sum(1 for e in new if e[K['mse']]<1.0)}")
print(f'  teams with raw pds > 0.82: '
      f"{sum(1 for e in old if e[K['pds']]>0.82)} -> "
      f"{sum(1 for e in new if e[K['pds']]>0.82)}")

print('\nTOP 20 NOW, and where each was on Aug 28')
print(f"{'#':>3} {'team':<22} {'now':>7} {'aug28':>8} {'delta':>8} "
      f"{'d pds':>7} {'d mse':>8} {'d reach':>8} {'d nmae':>8}")
for i, e in enumerate(sorted(new, key=lambda e: -e['score_avg'])[:22], 1):
    o = om.get(e['team_name'])
    if o is None:
        print(f'{i:>3} {e["team_name"][:22]:<22} {e["score_avg"]:7.4f} '
              f'{"NEW":>8}')
        continue
    print(f'{i:>3} {e["team_name"][:22]:<22} {e["score_avg"]:7.4f} '
          f'{o["score_avg"]:8.4f} {e["score_avg"]-o["score_avg"]:+8.4f} '
          f'{e[K["pds"]]-o[K["pds"]]:+7.4f} {e[K["mse"]]-o[K["mse"]]:+8.3f} '
          f'{e[K["reach"]]-o[K["reach"]]:+8.4f} {e[K["nmae"]]-o[K["nmae"]]:+8.4f}')

# ------------------------------------------------------------------ q inversion
print('\n' + '=' * 100)
print('IS THE REACH SPREAD IN THE SHRUNK CLUSTER EXPLAINED BY SIGN ACCURACY q?')
print('=' * 100)


def q_of(nmae, k):
    return (1 + k - nmae) / (2 * k) if k <= 1 else (k + 1 - nmae) / 2


ours = [e for e in new if e['team_name'] == 'Ahmet Dedeler'][0]
assert abs(q_of(ours[K['nmae']], 0.75) - 0.537) < 1.5e-3, q_of(ours[K['nmae']], 0.75)
print(f"known-answer: our nmae {ours[K['nmae']]:.4f} at k=0.75 -> "
      f'q = {q_of(ours[K["nmae"]], 0.75):.4f}  (CURRENT_STATE says 0.537)  OK\n')

C = sorted([e for e in new if e[K['mse']] < 1.2 and e[K['pds']] > 0.80],
           key=lambda e: -e[K['reach']])
print(f"{'team':<22} {'reach':>7} {'nmae':>7} | q at k=  {'0.6':>6} {'0.75':>6} "
      f"{'0.9':>6} {'1.0':>6} {'1.1':>6}")
for e in C:
    qs = '  '.join(f'{q_of(e[K["nmae"]], k):6.4f}' for k in (0.6, 0.75, 0.9, 1.0, 1.1))
    star = '  <== us' if e['team_name'] == 'Ahmet Dedeler' else ''
    print(f'{e["team_name"][:22]:<22} {e[K["reach"]]:7.4f} {e[K["nmae"]]:7.4f} |'
          f'          {qs}{star}')

print('\nThe decisive pairs (same k assumption cancels: whoever has the LOWER nmae'
      '\nhas the HIGHER q at EVERY k, because q is monotone decreasing in nmae):')
hi = [e for e in C if e[K['reach']] >= 0.245]
lo = [e for e in C if e[K['reach']] <= 0.19]
inv = 0
tot = 0
for a in hi:
    for b in lo:
        tot += 1
        if b[K['nmae']] < a[K['nmae']]:   # b has BETTER q but WORSE reach
            inv += 1
            print(f'  {b["team_name"][:20]:<20} reach {b[K["reach"]]:.4f} '
                  f'nmae {b[K["nmae"]]:.4f}   BEATS ON q   '
                  f'{a["team_name"][:20]:<20} reach {a[K["reach"]]:.4f} '
                  f'nmae {a[K["nmae"]]:.4f}')
print(f'\n  {inv}/{tot} high-reach x low-reach pairs are SIGN-INVERTED on q.')
print('  A pure sign-accuracy explanation of reach requires 0 inversions at any '
      'common k.')
