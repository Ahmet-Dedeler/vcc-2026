"""188 Phase 3 — resolve 169 §6.4's `pds` sign contradiction, by decomposition not by assertion.

THE CONTRADICTION
-----------------
For the SAME lfc-space sharpening arm at q = 1.15, 169 recorded

    PLOLO, challenge frame, delivered log-pseudobulk effect, 269-294 decoys   **+0.017631** scaled
    h1-test DE table, log2 fold change, 85 decoys                            **-0.0228**

169 trusted PLOLO (correctly: it is the board-validated instrument, 6/6 signs, and it computes
`pds_cosine` in the scorer's own space) but flagged the disagreement as something to resolve
before the family is revisited.  Three things differ at once -- the SPACE the cosine is taken
in, the DECOY POOL size, and the GENE SET.  This script varies them one at a time.

  1. On h1-test, convert BOTH sides of the same DE table into the scorer's `h` space using the
     harness's own control CPM (`h = log1p(TS*1e-6*cpm*2^lfc) - log1p(TS*1e-6*cpm)`, exact and
     invertible) and recompute `pds` on the SAME 85 decoys and the SAME genes.  Only the space
     changes.
  2. On PLOLO, recompute the same arm's `pds` in `lfc` space, and separately with the decoy
     pool subsampled to 85.  Only the space, then only the pool, changes.
  3. Measure the mechanism: per-gene sign accuracy against CONTROL EXPRESSION decile on both
     surfaces.  `lfc` up-weights low-expression genes relative to `h` by roughly
     `expm1(b)/log1p(...)`; if sign accuracy falls with expression, the two statistics must
     disagree on any arm that moves magnitude between expression bands -- which is exactly what
     `sharpen_lfc` does (169 §6.4).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
sys.path.insert(0, HERE)
sys.path.insert(0, E169)
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface                                   # noqa: E402
from budget188 import FOLDS, PDS_SPAN, deliver, load, pds_of   # noqa: E402
from mapfit import to_lfc                                      # noqa: E402
from sharp import sharpen, sharpen_lfc                          # noqa: E402

TS = 50_000.0
RNG = np.random.default_rng(188_000)
NSUB = 40                      # decoy-subsample repeats


def lfc_to_h(L, cpm):
    """Exact: the log-pseudobulk effect a log2 fold change implies, given the control CPM."""
    base = TS * 1e-6 * np.asarray(cpm, dtype=np.float64)
    return np.log1p(base[None, :] * np.exp2(np.nan_to_num(L))) - np.log1p(base)[None, :]


def pds_pool(P, T, npool, reps=NSUB):
    """`pds` with the decoy pool subsampled to `npool` targets, averaged over `reps` draws."""
    n = P.shape[0]
    if npool >= n:
        return pds_of(P, T), 0.0
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    own = np.diag(S).copy()
    vals = []
    for _ in range(reps):
        sel = RNG.choice(n, size=npool, replace=False)
        sub = S[np.ix_(sel, sel)]
        o = np.diag(sub).copy()
        beaten = (sub < o[:, None]).sum(1)
        ties = (sub == o[:, None]).sum(1) - 1
        vals.append(float(np.mean((beaten + 0.5 * ties) / (npool - 1))))
    return float(np.mean(vals)), float(np.std(vals))


def h1_block():
    """h1-test: the same DE tables 169 scored, plus the control CPM needed for the `h` space."""
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    egenes = np.asarray([str(g) for g in E["genes"]])
    ecpm = E["cpm"].astype(np.float64)
    pos = {g: i for i, g in enumerate(egenes)}
    arms = {}
    for tag, fn in (("q1.00", "de169.q1.00_mf.npz"), ("q1.15", "de169.q1.15_mf.npz"),
                    ("q1.15L", "de169.q1.15L_mf.npz"), ("q1.25", "de169.q1.25_mf.npz"),
                    ("q1.25L", "de169.q1.25L_mf.npz")):
        p = os.path.join(E169, fn)
        if os.path.exists(p):
            arms[tag] = p
    s0 = Surface(arms["q1.00"])
    cpm = np.array([ecpm[pos[g]] for g in s0.genes])
    return arms, cpm, s0


PART = os.environ.get("PART", "all")       # "h1" / "plolo" / "all" — RSS is shared, 4 agents


def main() -> None:
    print("=" * 100)
    print("1.  h1-test: ONE table, TWO spaces, SAME 85 decoys, SAME genes")
    arms, cpm, s0 = h1_block()
    ctrl_h = np.log1p(TS * 1e-6 * cpm)
    S = {t: Surface(p) for t, p in arms.items()}
    base = S["q1.00"]
    # the panel-target columns pds_cosine drops
    sidx = set(int(i) for i in base.self_idx if i >= 0)
    good = np.array([i for i in range(len(base.genes)) if i not in sidx])
    RT_l = np.nan_to_num(base.real_lfc)[:, good]
    RT_h = lfc_to_h(np.nan_to_num(base.real_lfc), cpm)[:, good]
    print(f"    {len(base.perts)} perturbations, {good.size} scored genes, "
          f"control CPM median {np.median(cpm):.2f}")
    print(f"    {'arm':<8}{'pds lfc':>10}{'d lfc':>10}{'pds h':>10}{'d h':>10}"
          f"{'d lfc scaled':>14}{'d h scaled':>12}")
    rows = {}
    b_l = pds_of(np.nan_to_num(base.pred_lfc)[:, good], RT_l)
    b_h = pds_of(lfc_to_h(np.nan_to_num(base.pred_lfc), cpm)[:, good], RT_h)
    for t in arms:
        s = S[t]
        pl = pds_of(np.nan_to_num(s.pred_lfc)[:, good], RT_l)
        ph = pds_of(lfc_to_h(np.nan_to_num(s.pred_lfc), cpm)[:, good], RT_h)
        rows[t] = dict(pds_lfc=pl, pds_h=ph, d_lfc=pl - b_l, d_h=ph - b_h)
        print(f"    {t:<8}{pl:>10.5f}{pl-b_l:>+10.5f}{ph:>10.5f}{ph-b_h:>+10.5f}"
              f"{(pl-b_l)/PDS_SPAN:>+14.5f}{(ph-b_h)/PDS_SPAN:>+12.5f}")
    print(f"    (169 §6.4 quotes the lfc-space `q1.15L` delta as -0.0228 scaled; "
          f"debench's own `Surface.pds` is reproduced below)")
    print(f"    debench Surface.pds direct: " + "  ".join(
        f"{t}={S[t].pds():.5f}" for t in arms))

    if PART == "h1":
        with open(os.path.join(WORK, "phase3-h1.json"), "w") as fh:
            json.dump(dict(h1=rows), fh, indent=1, default=float)
        print(f"\nwrote {WORK}/phase3-h1.json (PART=h1)")
        return

    print("\n" + "=" * 100)
    print("2.  PLOLO: the SAME arms, the space and the decoy pool varied one at a time")
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    bA, gA = z["bA"], z["gA"]
    F = {f: load(f) for f in FOLDS}
    out = {}
    for held in FOLDS:
        P0 = F[held]["pred"]
        sc = F[held]["scored"]
        Tr = F[held]["truth"][:, sc]
        Tl = to_lfc(F[held]["truth"], F[held]["b"], gA)[:, sc]
        for name, Q in (("q1.00", np.where(gA[None, :], P0, 0.0)),
                        ("q1.15", sharpen(P0, 1.15, gA, mode="row")),
                        ("q1.15L", sharpen_lfc(P0, 1.15, bA, gA, mode="row")),
                        ("q1.25", sharpen(P0, 1.25, gA, mode="row")),
                        ("q1.25L", sharpen_lfc(P0, 1.25, bA, gA, mode="row"))):
            D = deliver(Q, bA, gA)
            Ph = D[:, sc]
            Pl = to_lfc(D, bA, gA)[:, sc]
            p85, sd85 = pds_pool(Ph, Tr, 85)
            out[(held, name)] = dict(h=pds_of(Ph, Tr), lfc=pds_of(Pl, Tl),
                                     h85=p85, h85sd=sd85)
        print(f"    [{held}] done", flush=True)
    print(f"    {'arm':<8}{'d pds h':>10}{'scaled':>10}{'d pds lfc':>12}{'scaled':>10}"
          f"{'d pds h, 85 decoys':>21}{'scaled':>10}")
    for name in ("q1.15", "q1.15L", "q1.25", "q1.25L"):
        dh = float(np.mean([out[(f, name)]["h"] - out[(f, "q1.00")]["h"] for f in FOLDS]))
        dl = float(np.mean([out[(f, name)]["lfc"] - out[(f, "q1.00")]["lfc"] for f in FOLDS]))
        d8 = float(np.mean([out[(f, name)]["h85"] - out[(f, "q1.00")]["h85"] for f in FOLDS]))
        print(f"    {name:<8}{dh:>+10.5f}{dh/PDS_SPAN:>+10.5f}{dl:>+12.5f}"
              f"{dl/PDS_SPAN:>+10.5f}{d8:>+21.5f}{d8/PDS_SPAN:>+10.5f}")

    print("\n" + "=" * 100)
    print("3.  THE MECHANISM: sign accuracy against CONTROL EXPRESSION, both surfaces")
    print("    h1-test, reference-DE genes only (the pool `nmae`/`reach` are paid on)")
    ex = np.argsort(cpm)
    q = np.array_split(ex, 8)
    print(f"      {'octile':>7}{'ctrl CPM med':>14}{'n pairs':>10}{'P(sign ok)':>12}"
          f"{'mean |lfc_ref|':>15}{'mean |h_ref|':>14}")
    R = np.nan_to_num(base.real_lfc)
    Pd = np.nan_to_num(base.pred_lfc)
    sig = base.real_sig
    for i, sel in enumerate(q):
        m = sig[:, sel]
        if m.sum() == 0:
            continue
        ok = (np.sign(Pd[:, sel]) == np.sign(R[:, sel])) & m
        hr = lfc_to_h(R, cpm)[:, sel]
        print(f"      {i+1:>7}{np.median(cpm[sel]):>14.2f}{int(m.sum()):>10}"
              f"{ok.sum()/m.sum():>12.4f}{np.abs(R[:, sel])[m].mean():>15.4f}"
              f"{np.abs(hr)[m].mean():>14.4f}")

    print("    PLOLO (4 folds pooled), all scored+gated genes")
    print(f"      {'octile':>7}{'ctrl CPM med':>14}{'n pairs':>10}{'P(sign ok)':>12}"
          f"{'mean |r| h':>13}{'mean |r| lfc':>14}")
    accs, cpms = [], []
    for held in FOLDS:
        b = F[held]["b"]
        c = F[held]["scored"]
        c = c[gA[c]]
        cp = 1e6 * np.expm1(b[c]) / float(np.expm1(b).sum())
        P = F[held]["pred"][:, c]
        Rt = F[held]["truth"][:, c]
        accs.append((np.sign(P) == np.sign(Rt)) & (P != 0) & (Rt != 0))
        cpms.append(np.broadcast_to(cp, P.shape))
    A = np.concatenate([a.ravel() for a in accs])
    C = np.concatenate([c.ravel() for c in cpms])
    L = np.concatenate([np.abs(to_lfc(F[f]["truth"], F[f]["b"], gA)[:, F[f]["scored"][
        gA[F[f]["scored"]]]]).ravel() for f in FOLDS])
    H = np.concatenate([np.abs(F[f]["truth"][:, F[f]["scored"][gA[F[f]["scored"]]]]).ravel()
                        for f in FOLDS])
    ordc = np.argsort(C)
    for i, sel in enumerate(np.array_split(ordc, 8)):
        print(f"      {i+1:>7}{np.median(C[sel]):>14.2f}{sel.size:>10}"
              f"{A[sel].mean():>12.4f}{H[sel].mean():>13.4f}{L[sel].mean():>14.4f}")

    with open(os.path.join(WORK, "phase3.json"), "w") as fh:
        json.dump(dict(h1=rows, plolo={f"{k[0]}|{k[1]}": v for k, v in out.items()}),
                  fh, indent=1, default=float)
    print(f"\nwrote {WORK}/phase3.json")


if __name__ == "__main__":
    sys.exit(main())
