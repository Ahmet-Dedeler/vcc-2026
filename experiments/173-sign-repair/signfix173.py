"""173 -- does correcting whitening-flipped signs help, and what does it cost on pds?

On the h1-test reference DE pool, whitened H disagrees in sign with the un-whitened consensus
on 15.7% of pairs; H is right 45.6% there, cons 54.4%.  Arm: H' = H with sign set to cons's
where they disagree (magnitude kept).  Price pds/cos in the delivered space exactly as 168d did.
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments/168-sharpening"))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from h1_tiebreak import group_profiles, H1, PRED, CTRL, TS            # noqa: E402
from massfix import GATE_CPM, mass_fix                                 # noqa: E402
from screen import endpoints, render                                   # noqa: E402
from debench import Surface                                            # noqa: E402

# ---- truth, exactly as 168d built it -------------------------------------------------
tg, genes_t, acc_t = group_profiles(H1)
ci = int(np.flatnonzero(tg == CTRL)[0])
prof_t = np.log1p(TS * acc_t / acc_t.sum(1, keepdims=True))
b_ctrl = prof_t[ci]
keep_t = np.asarray([i for i in range(len(tg)) if i != ci])
T_all = prof_t[keep_t] - b_ctrl[None, :]; targets = tg[keep_t]

# ---- design effect H (169's recovery from the same artifact) --------------------------
E = np.load(os.path.join(REPO, "experiments/169-sharpening-six-member/work/h1test-effect.npz"), allow_pickle=True)
assert np.array_equal(np.asarray([str(g) for g in E["genes"]]), np.asarray([str(g) for g in genes_t])), "gene axes"
tpos = {str(t): i for i, t in enumerate(E["targets"])}
have = np.asarray([i for i, t in enumerate(targets) if str(t) in tpos])
T = T_all[have]; targets = targets[have]
H = E["H"][np.asarray([tpos[str(t)] for t in targets])].astype(np.float64)
gate = E["gate"].astype(bool)
expr = np.expm1(b_ctrl); cpm = 1e6 * expr / expr.sum()
gpos = {str(g): i for i, g in enumerate(genes_t)}
drop = {gpos[str(t)] for t in targets if str(t) in gpos}
mask = np.asarray([i for i in range(len(genes_t)) if i not in drop])

# ---- un-whitened consensus from 172, on the DE-universe axis -> full axis --------------
F = np.load(os.path.join(REPO, "experiments/172-sign-confidence-key/work/h1test-feats.npz"), allow_pickle=True)
S = Surface(os.path.join(REPO, "experiments/164-callset-artefact/de164.a15.npz"))
assert list(S.perts) == [str(t) for t in E["targets"]] or set(S.perts) == set(str(t) for t in E["targets"])
de_gi = np.asarray([gpos[g] for g in S.genes])                    # DE-universe gene -> full axis
de_ti = np.asarray([list(S.perts).index(str(t)) for t in targets])
cons = np.zeros_like(H); cons[:, de_gi] = F["cons"][de_ti].astype(np.float64)

sH, sC = np.sign(H), np.sign(cons)
flip = (sH != 0) & (sC != 0) & (sH == -sC)
print(f"h1-test: {len(targets)} targets, {int(gate.sum())} gated, {len(mask)} scored genes")
print(f"whitening flipped the sign on {flip[:, gate].mean():.1%} of gated (target,gene) entries; "
      f"|H| share of energy in flipped entries: {(H[flip]**2).sum()/(H**2).sum():.1%}")

def price(Heff, label):
    eff = mass_fix(np.where(gate[None, :], Heff, 0.0), b_ctrl, gate)
    P = render(eff, b_ctrl)[:, mask]; Q = T[:, mask]
    e = endpoints(P, Q)
    Pc = P - P.mean(0, keepdims=True); Qc = Q - Q.mean(0, keepdims=True)
    ec = endpoints(Pc, Qc)
    # per-gene sign accuracy of the DELIVERED effect against truth, on the h1 DE pool
    print(f"  {label:<44} pds {e['pds']:.5f}  cos {e['cos']:.5f}  | mean-removed pds {ec['pds']:.5f} cos {ec['cos']:.5f}")
    return e, ec

print("\n=== delivered-space pricing (168d's exact pipeline; a15-era base) ===")
b0, c0 = price(H, "BASE  H as shipped")
H_fix = np.where(flip, -H, H)
b1, c1 = price(H_fix, "ARM   flip H's sign to cons where they disagree")
H_zero = np.where(flip, 0.0, H)
b2, c2 = price(H_zero, "CTRL  zero the disagreeing entries instead")
rng = np.random.default_rng(0)
perm = np.zeros_like(flip)
for i in range(flip.shape[0]):
    k = int(flip[i].sum()); idx = np.flatnonzero(gate); perm[i, rng.choice(idx, k, replace=False)] = True
b3, c3 = price(np.where(perm, -H, H), "NULL  flip the SAME COUNT of random gated entries")
for lab, (b, c) in [("ARM", (b1, c1)), ("CTRL zero", (b2, c2)), ("NULL", (b3, c3))]:
    print(f"  d {lab:<10} pds {b['pds']-b0['pds']:+.5f}  cos {b['cos']-b0['cos']:+.5f}   mean-removed pds {c['pds']-c0['pds']:+.5f} cos {c['cos']-c0['cos']:+.5f}")

# ---- and the DE side: sign accuracy on the reference pool, key-counterfactual reach ----
pool = S.real_sig & S.adj
Hde = H[:, de_gi][np.argsort(de_ti)] if False else np.zeros(S.real_lfc.shape)
Hde = np.zeros(S.real_lfc.shape); Hde[de_ti] = H[:, de_gi]
Hfde = np.zeros(S.real_lfc.shape); Hfde[de_ti] = H_fix[:, de_gi]
sT = np.sign(S.real_lfc)
print(f"\nsign accuracy on the reference DE pool:  H {np.mean(np.sign(Hde)[pool]==sT[pool]):.4f}   H_fix {np.mean(np.sign(Hfde)[pool]==sT[pool]):.4f}")
np.savez_compressed(os.path.join(HERE, "signfix173.npz"), flip=flip, targets=targets, genes=genes_t)
