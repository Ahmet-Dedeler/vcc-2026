"""193 Phase 3 — the four DE members and a second `pds` surface, on h1-test.

WHY THE ARM CAN BE EXPRESSED HERE AT ALL
----------------------------------------
Two operators built on the same frozen eigenbasis differ by a single linear map:

    W_alt W_prod^-1 = V diag( f_alt(ev) / f_prod(ev) ) V^T

so an already-whitened effect matrix can be converted into the one the alternative operator
would have produced WITHOUT inverting anything ill-conditioned -- the ratio is bounded and the
eigenbasis is shared.  Production then restores the pre-whitening median row norm, and since
both arms restore it to the SAME target, re-normalising the converted block to the incumbent's
median row norm reproduces the alternative arm exactly.  (AGENTS.md 11.1 rule 11: this is only
legitimate against the operator that did the whitening, so the h1-test artifact -- an a15-era
build -- is converted with a15's own `union_pruned`, sha `e2f65ecf`, asserted below.)

TWO APPROXIMATIONS, STATED
--------------------------
1. `H` is RECOVERED from the shipped artifact's emitted counts (169 `probe_effect.py`), so it
   is the DELIVERED effect, after the render.  Applying the operator ratio to it is exact only
   up to the render's nonlinearity.  169/183 make the same substitution; the `pre`/`post`
   placements bracket it and both are reported.
2. `H` also contains the post-whitening CD4 residual (a15's `f = 0.10`), which production adds
   AFTER the operator and which this transform therefore also touches.  ~10 % of each row's
   norm is transformed that should not be.  Both approximations inflate the arm's apparent
   effect if anything, so a NULL result here is the conservative reading.

MEMBERS
-------
`pds`   -- delivered log-pseudobulk effect, own-target removed.  A SECOND SURFACE only: 188
           Phase 3 forbids `debench.Surface.pds()` from sizing a `pds` arm, and 183 records
           that h1-test's 85-row panel makes absolute values non-transferable.  Sign only.
`nmae`  -- exactly, from 188 §2.4's identity `nmae = 1 + (U - 2M)/D` on the reference's own
           BH-significant pool (median ~89 genes/row), which is the significance-selected pool
           188 requires; a top-N-by-magnitude pool inverts the sign of the deltas.
`reach` -- `debench.Surface.reach` with the delivered `|lfc|` as the ordering key.
`fid`   -- read as PRECISION with the call set HELD at the incumbent's (rule 9: `n_real` is
           unreadable on this surface), so it moves only through sign flips.
`jac`   -- same held call set, so it too moves only through sign flips.
"""
from __future__ import annotations

import gc
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in ("experiments/168-sharpening", "experiments/156-delivered-geometry",
          "experiments/153-de-bench", "src"):
    sys.path.insert(0, str(ROOT / p))
sys.path.insert(0, str(HERE))

from h1_tiebreak import group_profiles, H1, CTRL, TS      # noqa: E402
from massfix import mass_fix                              # noqa: E402
from screen import render                                 # noqa: E402
from debench import Surface                               # noqa: E402
import maps193                                            # noqa: E402
from screen193 import sha256_file                         # noqa: E402

OP = ROOT / "data/donors/whitener_union_pruned.npz"
OP_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
N_BASIS = 29394
GAMMA = 0.65
BASE_PDS_168D = 0.86134          # 173's reproduced known answer, via 183
CANDIDATES = ("prod(lam=0)", "ridge1e-06", "ridge1e-05", "ridge0.0001", "ridge0.007",
              "floorrel1e-06", "floorq0.1", "cliprank12000", "clipMP", "rie-LP-eta0.5",
              "nls-LW2020")
CANDIDATES_B = ("tail^1.2", "tail^0.8", "head^1.2")


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def pds_of(P, T):
    S = unit(P) @ unit(T).T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


def h_to_lfc(h, cpm):
    """Exact inverse of 188's `lfc_to_h`: the log2 fold change a log-space effect implies."""
    base = TS * 1e-6 * np.asarray(cpm, dtype=np.float64)
    lin = np.expm1(np.asarray(h) + np.log1p(base)[None, :])
    return np.log2(np.maximum(lin, 1e-300) / np.maximum(base, 1e-300)[None, :])


def main() -> None:
    assert sha256_file(OP) == OP_SHA, "operator sha mismatch — rule 11"

    tg, genes_t, acc_t = group_profiles(H1)
    ci = int(np.flatnonzero(tg == CTRL)[0])
    prof_t = np.log1p(TS * acc_t / acc_t.sum(1, keepdims=True))
    b_ctrl = prof_t[ci]
    keep_t = np.asarray([i for i in range(len(tg)) if i != ci])
    T_all = prof_t[keep_t] - b_ctrl[None, :]
    targets_all = tg[keep_t]

    E = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                allow_pickle=True)
    assert np.array_equal(np.asarray([str(g) for g in E["genes"]]),
                          np.asarray([str(g) for g in genes_t])), "gene axes disagree"
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    have = np.asarray([i for i, t in enumerate(targets_all) if str(t) in tpos])
    T = T_all[have]
    targets = targets_all[have]
    H = E["H"][np.asarray([tpos[str(t)] for t in targets])].astype(np.float64)
    gate = E["gate"].astype(bool)
    cpm_h1 = E["cpm"].astype(np.float64)
    gpos = {str(g): i for i, g in enumerate(genes_t)}
    drop = {gpos[str(t)] for t in targets if str(t) in gpos}
    mask = np.asarray([i for i in range(len(genes_t)) if i not in drop])

    # ---- operator columns, official axis -> h1-test axis -------------------------------
    official = [ln.strip() for ln in
                (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
                if ln.strip()]
    panel = {ln.split(",")[0] for ln in
             (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
             if ln.strip()}
    wsyms = sorted({g for g in official} - panel)          # the builder's union-axis rule
    with np.load(OP) as z:
        evals = np.asarray(z["evals"], float)
        V = np.asarray(z["evecs"], float)
    p = len(evals)
    # The operator axis is the DONOR union minus the panel, not all official genes, so it
    # cannot be re-derived from the bundle alone; recover it from the fold cache 162 built,
    # which asserted it against the operator itself.
    wcols = np.load(ROOT / "experiments/162-retrodicting-validator/work/fold-kolf2.npz",
                    allow_pickle=True)["wcols"]
    assert len(wcols) == p, f"{len(wcols)} operator columns vs {p} eigenvalues"
    del wsyms
    op_syms = [official[i] for i in wcols]
    h1_of = np.asarray([gpos.get(s, -1) for s in op_syms])
    ok = h1_of >= 0
    print(f"operator {OP.name} sha {sha256_file(OP)[:8]} p={p}; "
          f"{int(ok.sum())}/{p} operator genes present on h1-test's axis")
    rows_ok = np.flatnonzero(ok)
    cols_h1 = h1_of[ok]

    fmaps = maps193.build_maps(evals, N_BASIS, gamma=GAMMA)
    fmaps.update(maps193.build_maps_b(evals, N_BASIS, gamma=GAMMA))
    f_prod = fmaps["prod(lam=0)"]

    S = Surface(str(ROOT / "experiments/164-callset-artefact/de164.a15.npz"))
    de_gi = np.asarray([gpos[g] for g in S.genes])
    de_ti = np.asarray([list(S.perts).index(str(t)) for t in targets])
    pool = S.real_sig & S.adj
    call_held = S.pred_sig.copy()                     # rule 9: hold the call set
    sT = np.sign(S.real_lfc)
    absR = np.abs(S.real_lfc)
    print(f"h1-test reference pool: {int(pool.sum())} pairs, "
          f"median {np.median(pool.sum(1)):.0f} genes/row over {len(targets)} targets")

    Hg = np.where(gate[None, :], H, 0.0)
    block0 = Hg[:, cols_h1]
    norm0 = float(np.median(np.linalg.norm(block0, axis=1))) or 1.0

    out = {}
    for name in CANDIDATES + CANDIDATES_B:
        ratio = fmaps[name] / f_prod
        Y = np.zeros((Hg.shape[0], p), dtype=np.float64)
        Y[:, rows_ok] = block0
        Y = ((Y @ V) * ratio) @ V.T
        blk = Y[:, rows_ok]
        blk = blk * (norm0 / (float(np.median(np.linalg.norm(blk, axis=1))) or 1.0))
        Hn = Hg.copy()
        Hn[:, cols_h1] = blk
        Hn = np.where(gate[None, :], Hn, 0.0)
        del Y, blk
        gc.collect()

        D = render(mass_fix(Hn, b_ctrl, gate), b_ctrl)
        rec = dict(map=name, pds=pds_of(D[:, mask], T[:, mask]))
        L = h_to_lfc(D, cpm_h1)
        P = np.zeros(S.real_lfc.shape)
        P[de_ti] = L[:, de_gi]
        S.pred_lfc = P
        S.pred_sig = call_held
        # --- 188 §2.4's exact nmae identity on the reference's own BH pool ---------------
        U = np.abs(P) * pool
        M = np.minimum(np.abs(P), absR) * (np.sign(P) == sT) * pool
        Dn = absR * pool
        per = 1.0 + (U.sum(1) - 2 * M.sum(1)) / np.maximum(Dn.sum(1), 1e-300)
        rec["nmae_identity"] = float(np.mean(per))
        rec["U_over_D"] = float(U.sum() / Dn.sum())
        rec["M_over_D"] = float(M.sum() / Dn.sum())
        rec["nmae_debench"] = S.nmae()
        rec["q_sign"] = float(np.mean(np.sign(P)[pool] == sT[pool]))
        rec["k_mag"] = float(np.median(np.abs(P)[pool] / np.maximum(absR[pool], 1e-12)))
        rec["reach"] = S.reach(key=np.abs(P), call=call_held)
        rec["fid_precision"] = S.fid(call=call_held)
        rec["jac"] = S.jac(call=call_held)
        out[name] = rec
        print(f"  {name:16s} pds {rec['pds']:.5f}  nmae {rec['nmae_identity']:.5f}  "
              f"q {rec['q_sign']:.5f}  k {rec['k_mag']:.4f}  reach {rec['reach']:.5f}  "
              f"fid {rec['fid_precision']:.5f}  jac {rec['jac']:.5f}", flush=True)
        del D, L, P, Hn
        gc.collect()

    base = out["prod(lam=0)"]
    d = abs(base["pds"] - BASE_PDS_168D)
    print(f"\nGATE 5 (identity map reproduces 168d/173's h1-test base pds {BASE_PDS_168D}): "
          f"{base['pds']:.5f}, |d| = {d:.2e} ({'PASS' if d < 5e-4 else 'FAIL'})")
    print(f"\n{'map':16s}{'d pds':>10}{'d nmae':>10}{'d reach':>10}{'d fid':>10}"
          f"{'d jac':>10}{'d q':>10}")
    for name, r in out.items():
        print(f"{name:16s}{r['pds']-base['pds']:>+10.5f}"
              f"{r['nmae_identity']-base['nmae_identity']:>+10.5f}"
              f"{r['reach']-base['reach']:>+10.5f}"
              f"{r['fid_precision']-base['fid_precision']:>+10.5f}"
              f"{r['jac']-base['jac']:>+10.5f}{r['q_sign']-base['q_sign']:>+10.5f}")
    (HERE / "work/h1six.json").write_text(json.dumps(
        dict(operator=str(OP.relative_to(ROOT)), op_sha256=sha256_file(OP),
             gate5=dict(base_pds=base["pds"], expected=BASE_PDS_168D, diff=d),
             n_targets=len(targets), pool_pairs=int(pool.sum()), arms=out), indent=1))
    print("\nwrote work/h1six.json")


if __name__ == "__main__":
    main()
