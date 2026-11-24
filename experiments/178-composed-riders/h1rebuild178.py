"""178 — the h1-test effect REBUILT FROM THE DONORS, so gamma / donor weights / sign repair
are all expressible on one surface, and priced in the DELIVERED space (168d's pipeline).

WHY A REBUILD AND NOT 173's RECOVERY
------------------------------------
173 priced sign repair by transforming the effect `H` recovered from the frozen 095 h1-test
artifact. That works for a sign flip, which is a function of `H` and the un-whitened consensus.
It does NOT work for `WHITEN_GAMMA` (the exponent has to be re-applied by the operator that
whitened it, rule 11) or for donor weights (the artifact used two donors; the weighting is a
three-donor object). So the consensus is rebuilt here from
`data/donors/{k562_gwps, orion_hct116, orion_hek293t}` on h1-test's own 18,080-gene axis with
`MultiDonorModel.effect`'s exact arithmetic, and whitened by the operator production would use
for a line that is neither Orion nor K562:

    data/donors/whitener_union_pruned.npz   sha256 e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123

which is the SHIPPED operator (a15/a23's `operator.sha256`), re-exponentiated at each arm's
gamma. h1-test is H1 hESC: it contributed no row to that operator's basis, and
`data/donors/h1_2025` (which covers 85/85 of the targets and is pure oracle leakage) never
enters. The operator's columns are the challenge official-axis genes minus the 300 panel
targets, mapped onto h1-test's axis by SYMBOL; genes the operator does not carry pass through
un-whitened, exactly as production leaves off-axis coordinates alone.

FOUR KNOWN-ANSWER ROWS (rule 10's corollary)
--------------------------------------------
  KA1  `signfix173.py` reproduces its own base `pds` 0.86134 — run it separately; this script
       re-reads the same recovered `H` and re-prices it, which must give the same number.
  KA2  175's per-donor sign accuracies on the h1-test reference DE pool, unscaled and
       un-whitened: K562 0.5670 / HCT116 0.5375 / HEK293T 0.5308, (1,1,1) 0.5697, (2,1,1)
       0.5740, per-target oracle 0.6032. This is the row that validates the donor machinery
       AND the weighting switch, and it is measured on the exact code path 175 used.
  KA3  the REBUILD-vs-RECOVERY reconstruction check (097 §9b's discipline): per-row cosine of
       the rebuilt gamma-0.65 equal-weight effect against the artifact's own `H`, plus the
       delivered `pds`/`cos` of both. Absolute levels are NOT comparable across the two
       (different donor set, different operator, different alpha); every arm below is a RATIO
       inside the rebuild, where one pipeline produced every row.
  KA4  173's flip statistics on the rebuild: whitening must flip a similar fraction of gated
       entries (13.9% measured on the recovery) carrying a similar sliver of energy (1.6%).

CONTROL: the random-flip null (`permflip`), which flips the same number of gated entries per
row at random. If sign repair were a generic perturbation of a whitened matrix rather than a
targeted repair, the null would read the same way. 173 measured it at `pds` -0.034 / cos -0.010
against the arm's -0.004 / +0.003.

No cell-eval2, no build, no submission. Peak RSS ~4 GB.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
for p in ("src", "experiments/168-sharpening", "experiments/156-delivered-geometry",
          "experiments/153-de-bench", "experiments/172-sign-confidence-key"):
    sys.path.insert(0, str(ROOT / p))

from h1_tiebreak import CTRL, H1, TS, group_profiles          # noqa: E402
from massfix import mass_fix                                   # noqa: E402
from screen import endpoints, render                           # noqa: E402
from debench import Surface                                    # noqa: E402
from donors172 import DonorHandle                              # noqa: E402
from vcc2026.whiten import Whitener                            # noqa: E402

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")        # PRODUCTION_DONORS order
OP_PATH = ROOT / "data/donors/whitener_union_pruned.npz"
OP_META = ROOT / "data/donors/whitener_union_pruned.json"
OP_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
EFF_NPZ = ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz"
DE_NPZ = ROOT / "experiments/164-callset-artefact/de164.a15.npz"
ALPHA, EFFECT_SCALE = 0.982, 6.817                             # a15/a23 as shipped
CHUNK = 2000

# arm -> (weights in DONORS order, gamma, sign repair)
ARMS: dict[str, tuple[tuple[float, float, float], float, bool]] = {
    "a15":        ((1.0, 1.0, 1.0), 0.65, False),
    "g045":       ((1.0, 1.0, 1.0), 0.45, False),
    "g030":       ((1.0, 1.0, 1.0), 0.30, False),
    "flip":       ((1.0, 1.0, 1.0), 0.65, True),
    "w211":       ((2.0, 1.0, 1.0), 0.65, False),
    "flip_w211":  ((2.0, 1.0, 1.0), 0.65, True),
    "g045_flip":  ((1.0, 1.0, 1.0), 0.45, True),
    "g045_w211":  ((2.0, 1.0, 1.0), 0.45, False),
    "all3_g045":  ((2.0, 1.0, 1.0), 0.45, True),
    "all3_g030":  ((2.0, 1.0, 1.0), 0.30, True),
}


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_builder():
    """The production module, for `sign_repair_to_consensus` and `mass_consistent_effect`."""
    spec = importlib.util.spec_from_file_location(
        "prod_builder_178_h1", ROOT / "scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def shared_scale(handle: DonorHandle, dcols: np.ndarray) -> float:
    """`MultiDonorModel._scale[n]`: median over ALL the donor's perturbations of the norm of
    its mean-centred row restricted to the columns SHARED with the recipient axis.

    `DonorHandle.scale` is the same median over the donor's FULL axis, which is a different
    number; production renormalises on the shared columns, so that is what is used here.
    """
    cache = WORK / f"scale_shared.{handle.name}.{len(dcols)}.json"
    if cache.exists():
        return float(json.loads(cache.read_text())["scale"])
    n = handle.delta.shape[0]
    norms = np.empty(n, dtype=np.float64)
    m = handle.m.astype(np.float32)[dcols]
    for i in range(0, n, CHUNK):
        block = np.asarray(handle.delta[i:i + CHUNK][:, dcols], dtype=np.float32) - m
        norms[i:i + block.shape[0]] = np.linalg.norm(block, axis=1)
    value = float(np.median(norms))
    cache.write_text(json.dumps({"scale": value, "n_cols": int(len(dcols))}))
    return value


def main() -> None:
    t0 = time.time()
    M = load_builder()

    # ---------------------------------------------------------------- truth and axes
    tg, genes_t, acc_t = group_profiles(H1)
    ci = int(np.flatnonzero(tg == CTRL)[0])
    prof_t = np.log1p(TS * acc_t / acc_t.sum(1, keepdims=True))
    b_ctrl = prof_t[ci]
    keep_t = np.asarray([i for i in range(len(tg)) if i != ci])
    T_all = prof_t[keep_t] - b_ctrl[None, :]
    targets_all = tg[keep_t]

    E = np.load(EFF_NPZ, allow_pickle=True)
    assert np.array_equal(np.asarray([str(g) for g in E["genes"]]),
                          np.asarray([str(g) for g in genes_t])), "gene axes"
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    have = np.asarray([i for i, t in enumerate(targets_all) if str(t) in tpos])
    T = T_all[have]
    targets = [str(t) for t in targets_all[have]]
    H_rec = E["H"][np.asarray([tpos[t] for t in targets])].astype(np.float64)
    gate = E["gate"].astype(bool)
    gpos = {str(g): i for i, g in enumerate(genes_t)}
    drop = {gpos[t] for t in targets if t in gpos}
    mask = np.asarray([i for i in range(len(genes_t)) if i not in drop])
    print(f"h1-test: {len(targets)} targets, {len(genes_t)} genes, {int(gate.sum())} gated, "
          f"{len(mask)} scored  [{time.time()-t0:.0f}s]", flush=True)

    # ---------------------------------------------------------------- donor consensus
    S: dict[str, np.ndarray] = {}
    carried: dict[str, np.ndarray] = {}
    scale: dict[str, float] = {}
    for n in DONORS:
        h = DonorHandle(n)
        tcols, dcols = h.project(genes_t)
        rows, car = h.rows(targets, genes_t, tcols, dcols)
        S[n] = rows.astype(np.float64)
        carried[n] = car
        scale[n] = shared_scale(h, dcols)
        print(f"  {n:<15} covers {len(tcols):>6}/{len(genes_t)} genes, carries "
              f"{int(car.sum()):>3}/{len(targets)} targets, _scale {scale[n]:.5f}", flush=True)
        del h
    ref = float(np.median(list(scale.values())))
    # `MultiDonorModel`: renorm[n] = ref/scale[n], _global = effect_scale/ref, so the effect is
    #   alpha * effect_scale * sum_n (w_n/tot) * s_n / scale[n]
    print(f"  ref (median _scale) {ref:.5f}; production renorm "
          f"{ {n: round(ref/scale[n], 4) for n in DONORS} }", flush=True)

    def consensus(weights: tuple[float, float, float]) -> np.ndarray:
        acc = np.zeros((len(targets), len(genes_t)), dtype=np.float64)
        tot = np.zeros(len(targets), dtype=np.float64)
        for w, n in zip(weights, DONORS):
            tot += w * carried[n]
        assert (tot > 0).all(), "a target with no carrying donor would need the mean fallback"
        for w, n in zip(weights, DONORS):
            if w == 0.0:
                continue
            acc += (w / tot)[:, None] * (S[n] / scale[n])
        return ALPHA * EFFECT_SCALE * acc

    # ---------------------------------------------------------------- the operator
    assert sha_file(OP_PATH) == OP_SHA, "operator sha does not match the shipped one"
    op_meta = json.loads(OP_META.read_text())
    official = [line.strip() for line in
                (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
                if line.strip()]
    assert len(official) == 18533
    op_cols = np.load(WORK / "opcolumns.npy")          # built by effects178.py stage `pre`
    assert len(op_cols) == int(op_meta["n_cols"]) == 17942
    op_genes = [official[int(j)] for j in op_cols]
    h_idx, o_idx = [], []
    for k, g in enumerate(op_genes):
        j = gpos.get(g)
        if j is not None:
            h_idx.append(j)
            o_idx.append(k)
    h_idx = np.asarray(h_idx, dtype=np.int64)
    o_idx = np.asarray(o_idx, dtype=np.int64)
    print(f"  operator {OP_PATH.name} sha {OP_SHA[:8]} n_cols {len(op_cols)}; "
          f"{len(h_idx)} of them are on h1-test's axis "
          f"({100*len(h_idx)/len(op_cols):.1f}%)", flush=True)
    W = Whitener.load(OP_PATH)
    assert W.evecs.shape == (len(op_cols), len(op_cols))
    print(f"  whitener loaded {W.evecs.shape} [{time.time()-t0:.0f}s]", flush=True)

    def whiten(pre: np.ndarray, gamma: float) -> np.ndarray:
        """`apply_v9_whitening`'s arithmetic, on the operator columns present on this axis."""
        out = pre.copy()
        if gamma == 0.0:
            return out
        u = np.zeros((pre.shape[0], len(op_cols)), dtype=np.float64)
        u[:, o_idx] = pre[:, h_idx]
        before = float(np.median(np.linalg.norm(u, axis=1))) or 1.0
        t = W.apply(u, gamma, 0.0, renorm=False)
        after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
        out[:, h_idx] = t[:, o_idx] * (before / after)
        return out

    # ---------------------------------------------------------------- pricing helpers
    Q = T[:, mask]

    def price(Heff: np.ndarray) -> dict:
        eff = mass_fix(np.where(gate[None, :], Heff, 0.0), b_ctrl, gate)
        P = render(eff, b_ctrl)[:, mask]
        e = endpoints(P, Q)
        Pc, Qc = P - P.mean(0, keepdims=True), Q - Q.mean(0, keepdims=True)
        ec = endpoints(Pc, Qc)
        return dict(pds=e["pds"], cos=e["cos"], sse=e["sse"], gstar=e["gstar"], rr=e["rr"],
                    pds_c=ec["pds"], cos_c=ec["cos"], sse_c=ec["sse"],
                    rownorm=float(np.median(np.linalg.norm(eff, axis=1))))

    # DE-pool sign accuracy, on the reference's own DE calls (164's a15 table)
    Sf = Surface(str(DE_NPZ))
    pool = Sf.real_sig & Sf.adj
    sT = np.sign(Sf.real_lfc)
    de_gi = np.asarray([gpos[str(g)] for g in Sf.genes])
    de_ti = np.asarray([list(Sf.perts).index(t) for t in targets])
    pool_r = pool[de_ti]
    sT_r = sT[de_ti]

    def sign_acc(Heff: np.ndarray) -> float:
        key = np.sign(Heff[:, de_gi])
        m = pool_r & (key != 0)
        return float((key[m] == sT_r[m]).mean())

    # ================================================================ KA2 — 175's rows
    print(f"\n=== KA2: 175's per-donor sign accuracy (unscaled, un-whitened) ===", flush=True)
    K, Hc, He = (S[n][:, de_gi] for n in DONORS)
    allc = np.all([carried[n] for n in DONORS], axis=0)
    m3 = pool_r & allc[:, None] & (np.sign(K) != 0) & (np.sign(Hc) != 0) & (np.sign(He) != 0)
    # 175 reports TWO blocks and RESULTS.md's table is the second one. Both are checked:
    #   own-support   = pool & this donor carries & its entry is non-zero (n differs per donor)
    #   all-three     = pool & all three carry & all three non-zero (n = 39,945, the table)
    want_own = {"k562_gwps": 0.5668, "orion_hct116": 0.5366, "orion_hek293t": 0.5292}
    want_m3 = {"k562_gwps": 0.5670, "orion_hct116": 0.5375, "orion_hek293t": 0.5308}
    for n, arr in zip(DONORS, (K, Hc, He)):
        mm = pool_r & carried[n][:, None] & (np.sign(arr) != 0)
        got = float((np.sign(arr)[mm] == sT_r[mm]).mean())
        got3 = float((np.sign(arr)[m3] == sT_r[m3]).mean())
        print(f"  {n:<15} own-support {got:.4f} (want {want_own[n]:.4f}, n={int(mm.sum())}) "
              f"{'PASS' if abs(got - want_own[n]) < 2e-4 else 'FAIL'} | "
              f"all-three {got3:.4f} (want {want_m3[n]:.4f}) "
              f"{'PASS' if abs(got3 - want_m3[n]) < 2e-4 else 'FAIL'}")
    ka2 = {}
    for w, wnt in (((1, 1, 1), 0.5697), ((2, 1, 1), 0.5740), ((3, 1, 1), 0.5742)):
        c = (w[0] * K + w[1] * Hc + w[2] * He) / sum(w)
        got = float((np.sign(c)[m3] == sT_r[m3]).mean())
        ka2[str(w)] = got
        print(f"  weights {str(w):<9} {got:.4f}  (want {wnt:.4f}, n={int(m3.sum())})  "
              f"{'PASS' if abs(got - wnt) < 2e-4 else 'FAIL'}")

    # ================================================================ KA3 — rebuild check
    print(f"\n=== KA3: rebuild vs the 095 artifact's recovered H (ratios only) ===", flush=True)
    pre111 = consensus((1.0, 1.0, 1.0))
    reb = whiten(pre111, 0.65)
    a = reb[:, mask]
    b = H_rec[:, mask]
    rowcos = np.einsum("ij,ij->i", a, b) / np.maximum(
        np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-300)
    rec = price(H_rec)
    base = price(reb)
    print(f"  per-row cosine rebuild vs recovered H: median {np.median(rowcos):.5f}  "
          f"min {rowcos.min():.5f}  p05 {np.quantile(rowcos, 0.05):.5f}")
    print(f"  recovered H (095, 2 donors, 095 operator, alpha 1.08): "
          f"pds {rec['pds']:.5f} cos {rec['cos']:.5f} sign {sign_acc(H_rec):.4f} "
          f"rownorm {rec['rownorm']:.4f}   [signfix173 reads pds 0.86134 cos 0.03531]")
    print(f"  REBUILD   (3 donors, union_pruned, alpha 0.982):        "
          f"pds {base['pds']:.5f} cos {base['cos']:.5f} sign {sign_acc(reb):.4f} "
          f"rownorm {base['rownorm']:.4f}")

    # ================================================================ KA4 — flip stats
    flip_mask = (np.sign(reb[:, h_idx]) * np.sign(pre111[:, h_idx])) < 0.0
    on_gate = gate[h_idx]
    print(f"\n=== KA4: what whitening flips on the rebuild ===")
    print(f"  flipped {flip_mask[:, on_gate].mean():.1%} of gated on-axis entries "
          f"(173 on the recovery: 13.9%); energy share "
          f"{np.square(reb[:, h_idx][flip_mask]).sum() / np.square(reb).sum():.2%} "
          f"(173: 1.6%)")
    fm_full = np.zeros_like(reb, dtype=bool)
    fm_full[:, h_idx] = flip_mask
    p_flip = pool_r & fm_full[:, de_gi] & (np.sign(reb[:, de_gi]) != 0)
    if p_flip.sum():
        kr = np.sign(reb[:, de_gi])
        kc = np.sign(pre111[:, de_gi])
        print(f"  on the reference DE pool, {p_flip.sum() / max(pool_r.sum(), 1):.1%} of pairs "
              f"are flipped; whitened sign right {(kr[p_flip] == sT_r[p_flip]).mean():.4f}, "
              f"consensus right {(kc[p_flip] == sT_r[p_flip]).mean():.4f} "
              f"(173: 0.4556 / 0.5444)")

    # ================================================================ arms
    print(f"\n=== arms (delivered space, 168d's pipeline; deltas vs the a15 arm) ===",
          flush=True)
    pre = {(1.0, 1.0, 1.0): pre111, (2.0, 1.0, 1.0): consensus((2.0, 1.0, 1.0))}
    out: dict[str, dict] = {}
    keep: dict[str, np.ndarray] = {}
    for arm, (w, gamma, repair) in ARMS.items():
        Heff = whiten(pre[w], gamma)
        if repair:
            Heff, sr = M.sign_repair_to_consensus(Heff, pre[w], h_idx)
        else:
            sr = {"flipped_entries": 0}
        keep[arm] = Heff.astype(np.float64)
        row = price(Heff)
        row.update(arm=arm, weights=w, gamma=gamma, repair=repair,
                   sign=sign_acc(Heff), flips=int(sr["flipped_entries"]))
        out[arm] = row
        print(f"  {arm:<11} w{''.join(str(int(x)) for x in w)} g{gamma:<5} r{int(repair)}  "
              f"pds {row['pds']:.5f}  cos {row['cos']:.5f}  sign {row['sign']:.4f}  "
              f"pds_c {row['pds_c']:.5f} cos_c {row['cos_c']:.5f}", flush=True)

    # control: the random-flip null at the same per-row count as the `flip` arm
    rng = np.random.default_rng(178)
    perm = np.zeros_like(flip_mask)
    gated_axis = np.flatnonzero(on_gate)
    for i in range(flip_mask.shape[0]):
        k = int(flip_mask[i, on_gate].sum())
        perm[i, rng.choice(gated_axis, k, replace=False)] = True
    Hnull = reb.copy()
    Hnull[:, h_idx] = np.where(perm, -reb[:, h_idx], reb[:, h_idx])
    row = price(Hnull)
    row.update(arm="permflip_NULL", weights=(1.0, 1.0, 1.0), gamma=0.65, repair=False,
               sign=sign_acc(Hnull), flips=int(perm.sum()))
    out["permflip_NULL"] = row
    print(f"  {'permflip':<11} NULL, same count of RANDOM gated flips       "
          f"pds {row['pds']:.5f}  cos {row['cos']:.5f}  sign {row['sign']:.4f}")

    # ---------------------------------------------------------------- deltas + additivity
    b0 = out["a15"]
    print(f"\n{'arm':<14}{'d pds':>10}{'d cos':>10}{'d cos %':>9}{'d sign':>9}"
          f"{'d pds_c':>10}{'d cos_c':>10}")
    for arm, row in out.items():
        print(f"{arm:<14}{row['pds']-b0['pds']:>+10.5f}{row['cos']-b0['cos']:>+10.5f}"
              f"{100*(row['cos']/b0['cos']-1):>+9.2f}{row['sign']-b0['sign']:>+9.4f}"
              f"{row['pds_c']-b0['pds_c']:>+10.5f}{row['cos_c']-b0['cos_c']:>+10.5f}")

    print(f"\n=== additivity: composite vs the sum of its parts (h1-test delivered) ===")
    parts = {"g045": "g045", "flip": "flip", "w211": "w211"}
    for metric in ("pds", "cos", "sign", "cos_c"):
        s_all = sum(out[p][metric] - b0[metric] for p in parts.values())
        for comp, members in (("g045_flip", ("g045", "flip")),
                              ("g045_w211", ("g045", "w211")),
                              ("all3_g045", ("g045", "flip", "w211"))):
            s = sum(out[p][metric] - b0[metric] for p in members)
            d = out[comp][metric] - b0[metric]
            ratio = d / s if abs(s) > 1e-12 else float("nan")
            print(f"  {metric:<6} {comp:<11} composite {d:+.5f}  sum of parts {s:+.5f}  "
                  f"ratio {ratio:.3f}  "
                  f"{'super-additive' if ratio > 1.05 else 'sub-additive' if ratio < 0.95 else 'additive'}")
        del s_all

    # the rebuilt effect matrices, for `realize178.py`'s DE pass (same axis, same targets)
    np.savez_compressed(WORK / "h1arms178.npz", targets=np.asarray(targets),
                        genes=np.asarray([str(g) for g in genes_t]), gate=gate,
                        **{f"H_{k}": v for k, v in keep.items()})
    (WORK / "h1rebuild178.json").write_text(json.dumps(
        {"arms": {k: {kk: (list(vv) if isinstance(vv, tuple) else vv)
                      for kk, vv in v.items()} for k, v in out.items()},
         "ka2": ka2,
         "recovered": {k: v for k, v in rec.items()},
         "rowcos_median": float(np.median(rowcos)),
         "operator": {"path": str(OP_PATH.relative_to(ROOT)), "sha256": OP_SHA,
                      "n_cols": int(len(op_cols)), "on_h1_axis": int(len(h_idx))}},
        indent=1))
    print(f"\nwrote {WORK/'h1rebuild178.json'} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
