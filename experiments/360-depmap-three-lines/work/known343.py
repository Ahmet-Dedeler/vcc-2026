import numpy as np
S=np.load('~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work/stats-cap15.npz')
R=np.load('~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work/reference.npz')
pred_lfc, ref_lfc, ref_q = S['lfc'], R['lfc'], R['q']
elig = (ref_q < 0.05).sum(axis=1) > 0
print('eligible targets:', elig.sum())
nmae=[]
for i in np.where(elig)[0]:
    m = ref_q[i] < 0.05
    nmae.append(np.abs(pred_lfc[i][m]-ref_lfc[i][m]).sum()/np.abs(ref_lfc[i][m]).sum())
nmae=np.array(nmae)
print('mean NMAE: %.6f (record 0.950836)' % nmae.mean())
print('n targets:', len(nmae))
