"""Anatomy of raw `reach` on the 2026-09-03 board.

Question: what gives raw reach >= 0.25 at raw pds >= 0.84?
`reach` reads only (a) per-gene sign accuracy on the reference's DE pool and
(b) the order our emitted p-values impose on that pool.  Everything a team does
that touches either reaches the board through the six numbers we can see, so
this script asks which of the visible numbers actually carry the reach signal.

Known-answer row: our own entry (Ahmet Dedeler) must come out at the raw values
recorded in the mandate (pds 0.8368, reach 0.1773, mse 0.874) and its reach
residual must reproduce the mandate's -1.04 sigma.
"""
import json
import numpy as np

D = json.load(open('board/leaderboard_2026-09-03.json'))['entries']

K = dict(
    pds='pds_cosine',
    mse='expr_mse_unbiased_capped_norm',
    nmae='de_wilcoxon_lfc_nmae',
    fid='de_wilcoxon_direction_fidelity_yield_raw',
    reach='de_wilcoxon_direction_reach_raw',
    jac='de_wilcoxon_sig_jaccard',
)


def get(e, k):
    return e[K[k]]


# ---------------------------------------------------------------- known answer
ours = [e for e in D if e['team_name'] == 'Ahmet Dedeler'][0]
print('KNOWN-ANSWER ROW (ours):',
      f"pds {get(ours,'pds'):.4f} reach {get(ours,'reach'):.4f} "
      f"mse {get(ours,'mse'):.3f} fid {get(ours,'fid'):.4f} jac {get(ours,'jac'):.4f} "
      f"nmae {get(ours,'nmae'):.4f} overall {ours['score_avg']:.5f}")
assert abs(get(ours, 'pds') - 0.8368) < 5e-4
assert abs(get(ours, 'reach') - 0.1773) < 5e-4
assert abs(get(ours, 'mse') - 0.874) < 5e-3
print('  -> matches the mandate table.  OK\n')

# reference scaling: recover per-member (b, r) from raw + scaled, so we can
# convert raw deltas into overall points without assuming the doc's ranges.
print('IMPLIED REFERENCE ANCHORS (fit raw->scaled over the whole board)')
anchors = {}
for m in K:
    raw = np.array([get(e, m) for e in D], float)
    sc = np.array([e['score_' + m] for e in D], float)
    ok = np.isfinite(raw) & np.isfinite(sc)
    # scaled = (raw - b)/(r - b): fit a line, but skip clamped members
    if m == 'mse':
        ok &= (sc > 1e-9) & (sc < 1 - 1e-9)
    A = np.vstack([raw[ok], np.ones(ok.sum())]).T
    (slope, icept), *_ = np.linalg.lstsq(A, sc[ok], rcond=None)
    resid = sc[ok] - A @ np.array([slope, icept])
    b = -icept / slope
    span = 1.0 / slope
    anchors[m] = (b, span)
    print(f'  {m:<5} scaled = (raw - {b:+.4f}) / {span:+.4f}   '
          f'(fit rms {resid.std():.2e}, n={ok.sum()})')
print()

# ------------------------------------------------- the two high-reach families
print('EVERY TEAM WITH raw reach >= 0.24 (sorted by overall)')
hdr = (f"{'team':<22} {'model':<26} {'ovr':>7} {'pds':>7} {'mse':>7} {'nmae':>7} "
       f"{'fid':>7} {'reach':>7} {'jac':>7}")
print(hdr)
rows = sorted([e for e in D if get(e, 'reach') >= 0.24],
              key=lambda e: -e['score_avg'])
for e in rows:
    print(f"{e['team_name'][:22]:<22} {str(e['model_name'])[:26]:<26} "
          f"{e['score_avg']:7.4f} {get(e,'pds'):7.4f} {get(e,'mse'):7.3f} "
          f"{get(e,'nmae'):7.4f} {get(e,'fid'):7.4f} {get(e,'reach'):7.4f} "
          f"{get(e,'jac'):7.4f}")
print()

lo = [e for e in rows if get(e, 'mse') < 1.2]
hi = [e for e in rows if get(e, 'mse') >= 1.2]
print(f'  shrunk family  (raw mse < 1.2): n={len(lo)}  '
      f"median jac {np.median([get(e,'jac') for e in lo]):.4f}  "
      f"median fid {np.median([get(e,'fid') for e in lo]):.4f}  "
      f"median reach {np.median([get(e,'reach') for e in lo]):.4f}")
print(f'  loud family    (raw mse >=1.2): n={len(hi)}  '
      f"median jac {np.median([get(e,'jac') for e in hi]):.4f}  "
      f"median fid {np.median([get(e,'fid') for e in hi]):.4f}  "
      f"median reach {np.median([get(e,'reach') for e in hi]):.4f}")
print()

# ---------------------------------------------------- what predicts reach?
# Restrict to serious entries so we are not fitting the dead zone.
S = [e for e in D if get(e, 'pds') > 0.70 and np.isfinite(get(e, 'reach'))]
print(f'REGRESSIONS on the pds>0.70 cohort, n={len(S)}')
y = np.array([get(e, 'reach') for e in S])


def fit(names, label):
    X = [np.ones(len(S))]
    for n in names:
        if n == 'logmse':
            X.append(np.log(np.array([get(e, 'mse') for e in S])))
        else:
            X.append(np.array([get(e, n) for e in S]))
    X = np.vstack(X).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    r = y - pred
    r2 = 1 - r.var() / y.var()
    terms = '  '.join(f'{n}={b:+.4f}' for n, b in zip(names, beta[1:]))
    print(f'  {label:<34} R2={r2:.3f}  resid sd={r.std():.4f}   {terms}')
    return beta, r


_, r_pds = fit(['pds'], 'reach ~ pds')
fit(['pds', 'logmse'], 'reach ~ pds + log mse')
fit(['pds', 'logmse', 'jac'], 'reach ~ pds + log mse + jac')
fit(['pds', 'logmse', 'jac', 'fid'], 'reach ~ pds + log mse + jac + fid')
_, r_full = fit(['pds', 'logmse', 'jac', 'fid', 'nmae'], 'reach ~ all five others')

i_us = [i for i, e in enumerate(S) if e['team_name'] == 'Ahmet Dedeler'][0]
print(f'\n  our reach residual vs pds only : {r_pds[i_us]/r_pds.std():+.2f} sigma'
      '   (mandate says -1.04)')
print(f'  our reach residual vs all five : {r_full[i_us]/r_full.std():+.2f} sigma')
print()

# --------------------------------------- jac is the discriminator: show it flat
print('WITHIN THE SHRUNK CLUSTER (raw mse < 1.2, pds > 0.80): reach vs jac')
C = sorted([e for e in D if get(e, 'mse') < 1.2 and get(e, 'pds') > 0.80],
           key=lambda e: -get(e, 'reach'))
print(f"{'team':<22} {'model':<24} {'pds':>7} {'mse':>6} {'reach':>7} {'jac':>7} "
      f"{'fid':>7} {'nmae':>7} {'ovr':>7}")
for e in C:
    star = ' <== us' if e['team_name'] == 'Ahmet Dedeler' else ''
    print(f"{e['team_name'][:22]:<22} {str(e['model_name'])[:24]:<24} "
          f"{get(e,'pds'):7.4f} {get(e,'mse'):6.3f} {get(e,'reach'):7.4f} "
          f"{get(e,'jac'):7.4f} {get(e,'fid'):7.4f} {get(e,'nmae'):7.4f} "
          f"{e['score_avg']:7.4f}{star}")
jj = np.array([get(e, 'jac') for e in C])
rr = np.array([get(e, 'reach') for e in C])
print(f'\n  n={len(C)}  corr(reach, jac) = {np.corrcoef(rr, jj)[0,1]:+.3f}')
pp = np.array([get(e, 'pds') for e in C])
print(f'            corr(reach, pds) = {np.corrcoef(rr, pp)[0,1]:+.3f}')
nn = np.array([get(e, 'nmae') for e in C])
print(f'            corr(reach, nmae)= {np.corrcoef(rr, nn)[0,1]:+.3f}')
mm = np.log(np.array([get(e, 'mse') for e in C]))
print(f'            corr(reach, lmse)= {np.corrcoef(rr, mm)[0,1]:+.3f}')
ff = np.array([get(e, 'fid') for e in C])
print(f'            corr(reach, fid) = {np.corrcoef(rr, ff)[0,1]:+.3f}')

# partial: reach ~ jac controlling pds
A = np.vstack([np.ones(len(C)), pp, jj]).T
beta, *_ = np.linalg.lstsq(A, rr, rcond=None)
print(f'\n  reach ~ 1 + pds + jac  ->  d reach / d jac = {beta[2]:+.3f} '
      f'per unit jac  ({beta[2]*0.005:+.4f} raw reach per 0.005 jac)')
print(f'  our jac {get(ours,"jac"):.4f}; shrunk-cluster reach>=0.245 median jac '
      f'{np.median([get(e,"jac") for e in C if get(e,"reach")>=0.245]):.4f}')

# ----------------------------------------------------- overall-point arithmetic
b_r, span_r = anchors['reach']
b_j, span_j = anchors['jac']
print('\nPOINT ARITHMETIC (overall = mean of six scaled members)')
print(f'  +0.01 raw reach  = {0.01/span_r/6:+.5f} overall')
print(f'  -0.005 raw jac   = {-0.005/span_j/6:+.5f} overall')
print(f'  closing our reach gap to cqawesome (+{get(rows[0],"reach")-get(ours,"reach"):.4f} raw)'
      f' = {(get(rows[0],"reach")-get(ours,"reach"))/span_r/6:+.5f} overall')
