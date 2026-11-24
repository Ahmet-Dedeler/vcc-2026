"""224 Phase 1, rung 5 — the factorised bilinear / MLP over (target embedding x recipient
basal), fitted directly on the macro-cosine loss.

    P_hat[t,:] = P[t,:] * ( g0 + sum_r a_r(e_t) * c_r(x_g) )

  * `g0` is a free per-gene gain, so the model CONTAINS rung 4's diagonal exactly;
  * `e_t` is the PCA-50 target embedding (DepMap co-essentiality + STRING);
  * `x_g` is the per-gene recipient context [b_A, b_D, b_A - b_D], standardised — so the
    modulation is a function of the RECIPIENT, which is what makes the wrong-recipient
    null (Phase-2 falsifier 1, the one that killed 184) meaningful for it;
  * an optional additive low-rank term `sum_r u_r(e_t) v_r(x_g)` is fitted alongside.

Trained with Adam on the negative mean cosine of the inner leave-one-training-LINE-out
pairs, early-stopped on a holdout that is BOTH a different training line AND a disjoint
target block — the same honest selection `core224.select_operator` uses.

    uv run python mlp224.py --surface S7
"""
from __future__ import annotations

import argparse
import itertools
import json
import time

import numpy as np
import torch

from core224 import (WORK, decoy_cos, donor_mean, inner_pairs, load_surface,  # noqa: E402
                     macro_cos, per_target_cos, pds, unit)
from features224 import pca, target_features  # noqa: E402

SEED = 20260904
DEV = torch.device("mps" if torch.backends.mps.is_available() else "cpu")


class Modulator(torch.nn.Module):
    def __init__(self, G: int, d_e: int, d_g: int, rank: int, hid: int = 64,
                 additive: bool = True):
        super().__init__()
        self.g0 = torch.nn.Parameter(torch.ones(G))
        self.a = torch.nn.Sequential(torch.nn.Linear(d_e, hid), torch.nn.GELU(),
                                     torch.nn.Linear(hid, rank))
        self.c = torch.nn.Sequential(torch.nn.Linear(d_g, hid), torch.nn.GELU(),
                                     torch.nn.Linear(hid, rank))
        self.additive = additive
        if additive:
            self.u = torch.nn.Sequential(torch.nn.Linear(d_e, hid), torch.nn.GELU(),
                                         torch.nn.Linear(hid, rank))
            self.v = torch.nn.Sequential(torch.nn.Linear(d_g, hid), torch.nn.GELU(),
                                         torch.nn.Linear(hid, rank))
        for m in self.modules():                       # start at the identity
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.zeros_(m.bias)
        for net in ([self.a, self.c] + ([self.u, self.v] if additive else [])):
            torch.nn.init.zeros_(net[-1].weight)

    def forward(self, P, E, X):
        """P (n,G) donor mean; E (n,d_e) target embeddings; X (G,d_g) recipient context."""
        mod = self.g0[None, :] + self.a(E) @ self.c(X).T
        out = P * mod
        if self.additive:
            out = out + (self.u(E) @ self.v(X).T) * P.abs().mean()
        return out


def cos_loss(P, Y):
    return -(torch.nn.functional.normalize(P, dim=1)
             * torch.nn.functional.normalize(Y, dim=1)).sum(1).mean()


def gene_ctx(basal, A, D):
    bA = basal[A].astype(np.float64)
    bD = basal[list(D)].mean(0).astype(np.float64)
    X = np.stack([bA, bD, bA - bD], 1)
    return (X - X.mean(0, keepdims=True)) / np.maximum(X.std(0, keepdims=True), 1e-8)


def fit_one(C, basal, D, A, tsel, tval, tfit, teval, E, rank, epochs, lr, wd, rng,
            wrong_recipient=None):
    G = C.shape[2]
    hold = int(rng.integers(len(D)))
    Ps, Ys, Fs = inner_pairs(C, D, tsel)
    Pv, Yv, Fv = inner_pairs(C, D, tval)
    Es = np.tile(E[tsel], (len(D), 1))
    Ev = np.tile(E[tval], (len(D), 1))
    ms, mv = Fs != hold, Fv == hold

    # each inner recipient gets its OWN gene context, exactly as at deployment
    def ctx_stack(rows, Fmask, tset):
        out = []
        for j, B in enumerate(D):
            if not Fmask[j]:
                continue
            out.append(gene_ctx(basal, B, [x for x in D if x != B]))
        return out

    ctxs = {j: gene_ctx(basal, B, [x for x in D if x != B]) for j, B in enumerate(D)}
    t = lambda x: torch.as_tensor(np.ascontiguousarray(x), dtype=torch.float32, device=DEV)

    model = Modulator(G, E.shape[1], 3, rank).to(DEV)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    tr = [(t(Ps[(Fs == j) & ms]), t(Ys[(Fs == j) & ms]), t(Es[(Fs == j) & ms]), t(ctxs[j]))
          for j in range(len(D)) if (Fs == j).any() and j != hold]
    va = [(t(Pv[Fv == hold]), t(Yv[Fv == hold]), t(Ev[Fv == hold]), t(ctxs[hold]))]

    best, best_state, patience = np.inf, None, 0
    for ep in range(epochs):
        model.train()
        opt.zero_grad()
        loss = sum(cos_loss(model(P, Ee, X), Y) for P, Y, Ee, X in tr) / len(tr)
        loss.backward()
        opt.step()
        if ep % 5 == 0 or ep == epochs - 1:
            model.eval()
            with torch.no_grad():
                v = float(sum(cos_loss(model(P, Ee, X), Y) for P, Y, Ee, X in va) / len(va))
            if v < best - 1e-6:
                best, patience = v, 0
                best_state = {k: p.detach().clone() for k, p in model.state_dict().items()}
            else:
                patience += 1
                if patience >= 8:
                    break
    if best_state is not None:
        model.load_state_dict(best_state)

    # deployment: refit is NOT done (early stopping already used the val fold); predict
    model.eval()
    Pheld = donor_mean(C, D, teval)
    Theld = C[A][teval].astype(np.float64)
    rec = A if wrong_recipient is None else wrong_recipient
    Xg = gene_ctx(basal, rec, D)
    with torch.no_grad():
        out = model(t(Pheld), t(E[teval]), t(Xg)).cpu().numpy().astype(np.float64)
    v = per_target_cos(out, Theld)
    dec = decoy_cos(out, Theld)
    return dict(line=float(v.mean()), pds=pds(out, Theld), decoy=dec,
                excess=float(v.mean()) - dec, val=float(-best), epochs=int(ep + 1),
                rank=rank), v.astype(np.float32)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--surface", default="S7")
    ap.add_argument("--ranks", default="4,16,64")
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--lr", type=float, default=0.02)
    ap.add_argument("--wd", type=float, default=1e-4)
    ap.add_argument("--max-subsets", type=int, default=2)
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    C, basal, members, groups, targets, genes, dropped = load_surface(args.surface)
    L, T = len(members), C.shape[1]
    perm = rng.permutation(T)
    tfit = np.sort(perm[: T // 2])
    teval = np.sort(perm[T // 2:])
    q = rng.permutation(len(tfit))
    tsel = np.sort(tfit[q[: int(0.6 * len(tfit))]])
    tval = np.sort(tfit[q[int(0.6 * len(tfit)):]])
    Xf, have, fmeta = target_features(targets)
    E = pca(Xf, np.asarray(tfit), 50)
    print(f"surface {args.surface} {members}  device={DEV}  {T}x{C.shape[2]}")

    rows, t0 = [], time.time()
    for A in range(L):
        others = [x for x in range(L) if x != A]
        for k in (len(others), max(2, len(others) - 2), 2):
            combos = list(itertools.combinations(others, k))
            if len(combos) > args.max_subsets:
                combos = [combos[i] for i in
                          rng.choice(len(combos), args.max_subsets, replace=False)]
            for D in combos:
                base = macro_cos(donor_mean(C, list(D), teval), C[A][teval].astype(np.float64))
                for rank in [int(x) for x in args.ranks.split(",")]:
                    r, _ = fit_one(C, basal, list(D), A, tsel, tval, tfit, teval, E,
                                   rank, args.epochs, args.lr, args.wd, rng)
                    wrong = [x for x in range(L) if x != A and x not in D]
                    wr = None
                    if wrong:
                        wr, _ = fit_one(C, basal, list(D), A, tsel, tval, tfit, teval, E,
                                        rank, args.epochs, args.lr, args.wd,
                                        np.random.default_rng(SEED),
                                        wrong_recipient=int(wrong[0]))
                    rows.append(dict(held=members[A], k=k, donors=[members[x] for x in D],
                                     R0=base, R5=r,
                                     R5_wrong_recipient=wr,
                                     wrong=None if not wrong else members[wrong[0]]))
                    print(f"  [{time.time() - t0:5.0f}s] {members[A]:14s} k={k} r={rank:3d} "
                          f"R0 {base:.5f} R5 {r['line']:.5f} "
                          f"({100 * (r['line'] - base) / base:+.2f}%) dpds "
                          f"{r['pds'] - pds(donor_mean(C, list(D), teval), C[A][teval].astype(np.float64)):+.4f}"
                          f"  wrong-recipient "
                          f"{'-' if wr is None else f'{wr[chr(108)+chr(105)+chr(110)+chr(101)]:.5f}'}",
                          flush=True)

    (WORK / f"mlp224-{args.surface}.json").write_text(json.dumps(
        dict(surface=args.surface, members=members, rows=rows, device=str(DEV),
             epochs=args.epochs, lr=args.lr, wd=args.wd), indent=1))
    print(f"\nwrote mlp224-{args.surface}.json  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
