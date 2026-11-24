"""181 — THE known-answer row 177 and 179 never ran: the board has already scored the
resample-tail emission law.  a21 = a15 + `resample_x 0.14` (keep the top 14% of genes on the
parametric law, RESAMPLE THE REST) — the same mechanism as `hy_kt8000`, at a larger dose.
Its meta (submissions/a21-resample-x014.h5ad.meta.json) differs from a15's in `resample_x`
alone: alpha 0.982, keep_top 8000, keep_spread 0.25, tail 1.0, phi 4, min_cpm 10, cd4_target_f 0.1.

Compare 177's h1-test prediction for the same family against a21's realised board deltas.
Board deltas are SCALED (docs/LOG.md 2026-09-01); h1-test deltas are raw and are converted
with the per-member 1/span fitted on the 581-entry board snapshot (matches AGENTS.md S4's
published spans on all four members).
"""
import numpy as np, sys, json
sys.path.insert(0,'experiments/153-de-bench')
from debench import Surface
P='experiments/177-resample-emitter/de177.%s.npz'

# 1/span, fitted on data/leaderboard/2026-09-02T0630Z.json (see board181.py); AGENTS.md S4 spans
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)
SPAN_PUB = dict(reach='0.86-0.93 -> 1.075-1.163', fid='0.28-0.33 -> 3.03-3.57',
                jac='0.34-0.39 -> 2.56-2.94', nmae='0.57-0.63 -> 1.59-1.75')
for k,v in INV.items(): print(f'  1/span[{k}] fitted {v:.3f}   AGENTS.md S4 {SPAN_PUB[k]}')

# a21 board deltas vs a15, docs/LOG.md 2026-09-01 (scaled)
A21 = dict(pds=0.0, mse=0.0, nmae=-0.064214, fid=-0.029464, reach=-0.018482, jac=-0.016583)
A22 = dict(pds=0.0, mse=0.0, nmae=+0.001520, fid=-0.004447, reach=-0.010656, jac=-0.001523)

def rows(a, b='base_mf'):
    sA=Surface(P%a); sB=Surface(P%b)
    def prec(s):
        o=[]
        for p in range(85):
            i=np.flatnonzero(s.pred_sig[p]&s.adj[p])
            if i.size: o.append((np.sign(s.pred_lfc[p,i])==np.sign(s.real_lfc[p,i])).mean())
        return float(np.mean(o))
    return dict(reach=sA.reach()-sB.reach(), fid=sA.fid()-sB.fid(), jac=sA.jac()-sB.jac(),
                nmae=-(sA.nmae()-sB.nmae()), prec=prec(sA)-prec(sB))   # nmae: lower raw is better

print()
print('a21 = resample 86% of genes (the tail).  177 s closest measured arm is `rs_all` (100%).')
d = rows('rs_all')
print(f"{'member':7} {'h1-test raw d':>14} {'-> scaled':>11} {'a21 BOARD scaled':>18} {'sign':>6} {'ratio':>8}")
for m in ('reach','fid','jac','nmae'):
    hs = d[m]*INV[m]
    ok = 'OK' if np.sign(hs)==np.sign(A21[m]) else 'WRONG'
    print(f'{m:7} {d[m]:+14.5f} {hs:+11.5f} {A21[m]:+18.6f} {ok:>6} {hs/A21[m]:8.2f}')
print(f"  h1-test precision delta {d['prec']:+.5f} (rises, as 164 predicts); "
      f"board a21 raw fid fell by {(A21['fid'])/INV['fid']:+.5f} despite n_pred 4,420 -> 1,141")

print()
print('177 s section 6 table predicts, for `rs_all` at the settled n_real ~240:  +0.0211 overall')
real = sum(A21[m] for m in ('nmae','fid','reach','jac'))/6
print(f'a21 realised, the same four DE members in overall units:                {real:+.4f} overall')
print(f'   -> the table is off by {0.0211-real:+.4f} overall and its SIGN is inverted on the one')
print('      board point its own family has.')

print()
print('the same for a22 (resample the top 10% — the opposite assignment; 177 does cite this one):')
d2 = rows('hy_kt2000')   # nearest available: keeps 2000 parametric, resamples the rest (not the same shape)
print('   (no h1-test arm in 177 resamples only the head, so a22 has no matched row here)')
print(f'   board a22: reach {A22["reach"]:+.6f}  nmae {A22["nmae"]:+.6f}  fid {A22["fid"]:+.6f}  jac {A22["jac"]:+.6f}')
print()
print('Rule-13 toll check: the -0.0146 emitted-profile/emission-law reach toll is exactly')
print(f'   mean(a21, a22) reach = {(A21["reach"]+A22["reach"])/2:+.6f}.  Both board points, no seed noise')
print('   (2026-09-01: the board scorer has zero run-to-run variance).  177 s counter-argument')
print('   ("rule 13 s reach was fitted on arms inside the seed spread") addresses PLOLO s 2/6')
print('   signs, not this toll, whose two calibration points are exact board measurements of')
print('   THE SAME EMISSION LAW 177 proposes.')

# ---- addendum: a21's kept fraction, and the matched h1-test bracket -------------------
print()
print('ADDENDUM — a21 keeps k = round(0.14 * 18,533) = %d genes parametric by |effect| over the'
      % round(0.14*18533))
print('FULL gene axis; `profile_from_effect` zeroes the effect outside the pb-CPM-10 gate, so all')
print('of them are gated genes.  The challenge gate holds ~9,135 / 8,668 / 9,127 genes, so a21 keeps')
print('~28%% parametric and resamples ~72%% of the gated axis  ->  it is `hy_kt~2600`, which sits')
print('between 177 s hy_kt2000 (18.6%% kept on h1-test) and hy_kt4000 (37%%) — the region where 177')
print('measures the LARGEST reach gain.')
print()
print(f"{'analogue':12} {'member':7} {'h1-test raw d':>14} {'-> scaled':>11} {'a21 BOARD':>12} {'sign':>6}")
for a in ('hy_kt2000','hy_kt4000'):
    d = rows(a)
    for m in ('reach','fid','jac','nmae'):
        hs = d[m]*INV[m]
        print(f'{a:12} {m:7} {d[m]:+14.5f} {hs:+11.5f} {A21[m]:+12.6f} '
              f'{"OK" if np.sign(hs)==np.sign(A21[m]) else "WRONG":>6}')
print()
print('CONTRACT CHECK (the one test that could overturn this): at commit 9e521f5 — the state that')
print('built a21/a22/a23 — the builder sets  resample_ctrl = control_matrix  (the raw context h5ad')
print('CSR, never rescaled),  libraries = control_matrix.sum(axis=1)  (the cells own libraries) and')
print('resample_profile = totals/totals.sum()  (the POOLED control profile).  Both of 177 section 7 s')
print('call-site contracts were honoured by a21.  a21 is a valid board test of the resample law.')
