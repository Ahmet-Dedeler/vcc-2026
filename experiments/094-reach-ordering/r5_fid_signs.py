"""094 step 5 — is there a sign bias in the reference we fail to match?

`fid`'s numerator k counts sign matches over OUR CALL SET intersected with the
adjudicable genes, which is dominated by genes the reference did NOT call
significant. Their reference lfc is small, but small is not the same as
sign-random. If the reference drifts systematically one way and we do not, that
is free `fid` -- and for genes outside the reference's DE set it costs nothing on
`nmae` (computed only over the reference's DE set) or `reach` (same pool).
"""
import os
import numpy as np

E092 = "~/Code/virtual-cell-challenge-2026/experiments/092-call-set-truncation"
d = np.load(os.path.join(E092, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pq, pl = d["pred_padj"], d["pred_lfc"]
G, P = genes.size, perts.size
ALPHA = 0.05
gpos = {g: i for i, g in enumerate(genes)}
B_FID, R_FID = 0.51239, 0.80964

pos_all, pos_de, pos_nonde, acc_call, frac_up_call = [], [], [], [], []
fid_prod, fid_allup, fid_alldn, fid_match = [], [], [], []
for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA)
    ps = ok & (pq[i] < ALPHA)
    ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])
    scored = ps & ind
    if scored.sum() == 0: continue
    n_real = int(rs.sum())
    rsign = np.sign(rl[i]); psign = np.sign(np.nan_to_num(pl[i]))
    pos_all.append((rsign[ind] > 0).mean())
    if (rs & ind).sum(): pos_de.append((rsign[rs & ind] > 0).mean())
    nd = ind & ~rs
    if nd.sum(): pos_nonde.append((rsign[nd] > 0).mean())
    frac_up_call.append((psign[scored] > 0).mean())
    m = int((scored & (psign == rsign)).sum())
    den = max(int(scored.sum()), n_real)
    fid_prod.append(m / den)
    acc_call.append(m / int(scored.sum()))
    # counterfactual: same call set, but every call declared UP / DOWN
    fid_allup.append(int((scored & (rsign > 0)).sum()) / den)
    fid_alldn.append(int((scored & (rsign < 0)).sum()) / den)

def sc(v): return (v - B_FID) / (R_FID - B_FID)
print("=== reference lfc sign, by gene class (mean over perturbations) ===")
print(f"  all adjudicable genes           : {np.mean(pos_all):.4f} positive")
print(f"  reference-SIGNIFICANT genes     : {np.mean(pos_de):.4f} positive")
print(f"  reference-NONsignificant genes  : {np.mean(pos_nonde):.4f} positive")
print(f"\n  our call set: fraction we predict UP = {np.mean(frac_up_call):.4f}")
print(f"  our sign accuracy on our call set    = {np.mean(acc_call):.4f}")
print("\n=== fid counterfactuals, same call set, board anchors ===")
F = np.mean(fid_prod)
for name, v in (("production", F), ("declare every call UP", np.mean(fid_allup)),
                ("declare every call DOWN", np.mean(fid_alldn))):
    print(f"  {name:<26} raw fid {v:.4f}   scaled {sc(v):+.4f}   Δoverall {(sc(v)-sc(F))/6:+.5f}")
print("\n(the UP/DOWN rows are NOT buildable as stated -- they would move nmae and pds.")
print(" They bound how much sign bias is on the table if it exists at all.)")
