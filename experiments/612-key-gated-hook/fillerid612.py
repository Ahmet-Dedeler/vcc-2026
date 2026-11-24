# 612 side test: does the IDENTITY-augmented sign (622's S, Flex refs, lam 0.5) predict the sign of NON-DE genes on the
# CPM block better than CLEAN4? (fid's filler precision; d4cd5f's exchange rate: +0.01 fid-like ~ +0.0055 overall).
# h1-test truth (213's real DE table; 400 cells x 85 targets vs 18,400 controls = the board shape). Analysis only.
import numpy as np
d = np.load('experiments/213-transport-law/work/de_real_h1.npz', allow_pickle=True)
genes = list(d['genes']); perts = list(d['perts']); padj = d['real_padj'].copy(); lfc = d['real_lfc'].copy()
P, G = lfc.shape; gi = {g: i for i, g in enumerate(genes)}
for p, t in enumerate(perts):
    if t in gi: padj[p, gi[t]] = 1.0; lfc[p, gi[t]] = np.nan
sig = padj < 0.05; adj = np.isfinite(lfc) & (lfc != 0)
z = np.load('experiments/610-board-forensics/work/nullskew_h1test.npz', allow_pickle=True)
assert list(z['genes']) == genes
cpm = z['cpm']; order = np.argsort(-cpm)
idz = np.load('experiments/612-key-gated-hook/work/idsign612.npz', allow_pickle=True)
ipos = {g: i for i, g in enumerate(idz['genes'].astype(str))}
s_id = np.array([idz['sign'][ipos[g]] if g in ipos else 0.0 for g in genes])
c4 = np.load('experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz', allow_pickle=True)
c4g = {g: i for i, g in enumerate(c4['genes'])}
pd4 = np.array([c4['P_down'][c4g[g]] if g in c4g else np.nan for g in genes])
s_c4 = np.where(np.isfinite(pd4), np.where(pd4 > 0.5, -1.0, 1.0), 0.0)
s_dn = -np.ones(G)
def prec(s, M):   # precision of sign s on NON-DE adjudicable genes of mask M (genes where s != 0), mean over targets
    v = []
    for p in range(P):
        m = M & adj[p] & ~sig[p] & (s != 0)
        if m.sum(): v.append((np.sign(lfc[p, m]) == s[m]).mean())
    return np.mean(v), int((M & (s != 0)).sum())
def fid(s, M):
    F = []
    for p in range(P):
        c = M & adj[p] & (s != 0); npred = int((M & (s != 0)).sum()); nreal = int(sig[p].sum())
        k = int((np.sign(lfc[p, c]) == s[c]).sum())
        if nreal == 0 and npred == 0: continue
        F.append(k / npred if nreal == 0 else k / max(npred, nreal))
    return np.mean(F)
print('h1-test, block = top-K genes by control CPM; precision on NON-DE genes [n genes with a sign] | fid')
for K in (1000, 2200, 4000):
    M = np.zeros(G, bool); M[order[:K]] = True
    both = M & (s_c4 != 0) & (s_id != 0)
    r = {n: prec(s, both) for n, s in (('alldown', s_dn), ('CLEAN4', s_c4), ('identity622', s_id))}
    f = {n: fid(s, both) for n, s in (('alldown', s_dn), ('CLEAN4', s_c4), ('identity622', s_id))}
    print('K %5d (both signed: %d) | ' % (K, both.sum()) + '  '.join('%s %.4f' % (n, r[n][0]) for n in r) + ' | fid ' + '  '.join('%s %.4f' % (n, f[n]) for n in f))
    dis = both & (s_c4 != s_id)
    rd = {n: prec(s, dis)[0] for n, s in (('CLEAN4', s_c4), ('identity622', s_id))}
    print('        where they disagree (%d genes): CLEAN4 %.4f identity %.4f' % (dis.sum(), rd['CLEAN4'], rd['identity622']))
