"""182 -- score pretrained ST-HVG against our production consensus on per-gene SIGN.

Rows reported:
  KNOWN-ANSWER  K562 (a line INSIDE the checkpoint's training set), truth = the measured
                Replogle K562 pseudobulk delta. If ST does not beat chance here by a lot,
                the inference path is broken and nothing below is readable.
  HELD          H1 (a line the checkpoint never saw), truth = the h1-test reference DE
                pool, exactly the object production's 0.537 is measured on.

Every accuracy is reported against its own context-mean null (the same prediction for
every target), because a common stress axis can carry a lot of apparent sign accuracy
without any target-specific information -- that is the `fid` baseline of AGENTS.md §4.
"""

import json
import pickle
import sys

import anndata as ad
import numpy as np

sys.path.insert(0, "experiments/153-de-bench")
sys.path.insert(0, "experiments/172-sign-confidence-key")
from debench import Surface  # noqa: E402
from donors172 import DonorHandle  # noqa: E402

MODEL_DIR = "data/external/state-pretrained/ST-HVG-Replogle/zeroshot/hepg2"
WORK = "experiments/182-pretrained-state/work"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")  # production's equal-weight trio


def state_delta(tag, targets):
    """Predicted effect per target = mean(predicted target cells) - mean(predicted controls)."""
    A = ad.read_h5ad(f"{WORK}/{tag}_pred.h5ad")
    H = np.asarray(A.obsm["X_hvg"], dtype=np.float64)
    lab = A.obs["gene"].astype(str).values
    ctrl = H[lab == "non-targeting"].mean(axis=0)
    return np.stack([H[lab == t].mean(axis=0) - ctrl for t in targets])


def acc(pred, truth, mask):
    """Sign accuracy over `mask`, ignoring entries where either side is exactly zero."""
    sp, st = np.sign(pred), np.sign(truth)
    m = mask & (sp != 0) & (st != 0)
    return float((sp[m] == st[m]).mean()), int(m.sum())


def boot_targets(P, T, mask, n_boot=2000, seed=20260903):
    """95% CI on sign accuracy, resampling TARGETS -- the axis the metric averages over."""
    rng = np.random.default_rng(seed)
    per = []
    for i in range(P.shape[0]):
        m = mask[i] & (np.sign(P[i]) != 0) & (np.sign(T[i]) != 0)
        per.append((int((np.sign(P[i][m]) == np.sign(T[i][m])).sum()), int(m.sum())))
    per = np.array(per, dtype=float)
    idx = rng.integers(0, len(per), size=(n_boot, len(per)))
    s = per[idx].sum(axis=1)
    d = s[:, 0] / np.maximum(s[:, 1], 1)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def with_null(P, T, mask, label, out):
    """Report a predictor and its own context-mean null on the same mask."""
    a, n = acc(P, T, mask)
    null = np.repeat(P.mean(axis=0, keepdims=True), P.shape[0], axis=0)
    a0, _ = acc(null, T, mask)
    lo, hi = boot_targets(P, T, mask)
    print(
        f"  {label:<34} {a:.4f} [{lo:.4f},{hi:.4f}]   "
        f"(context-mean null {a0:.4f}, excess {a - a0:+.4f}, n={n})"
    )
    out[label] = {
        "sign_acc": a, "ci95": [lo, hi], "context_mean_null": a0,
        "excess": a - a0, "n_pairs": n,
        "median_abs_effect": float(np.median(np.abs(P[mask]))),
    }
    return a


def discrimination(P, T):
    """Cosine rank of each target's own real effect among all targets (a pds analogue)."""
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    C = Pn @ Tn.T
    n = C.shape[0]
    # fraction of wrong targets that the correct one beats, averaged (0.5 = chance)
    return float(np.mean([(C[i, i] > np.delete(C[i], i)).mean() for i in range(n)]))


def main():
    res = {}
    vd = pickle.load(open(f"{MODEL_DIR}/var_dims.pkl", "rb"))
    mgenes = np.array([str(g) for g in vd["gene_names"]])
    targets = json.load(open(f"{WORK}/prep_meta.json"))["targets"]
    print(f"targets: {len(targets)}  model genes: {len(mgenes)}\n")

    # ---------------------------------------------------------------- KNOWN-ANSWER: K562
    print("KNOWN-ANSWER ROW -- K562, a line INSIDE the checkpoint's training set")
    h = DonorHandle("k562_essential")
    K, carried = h.rows(targets, mgenes)  # truth on the model's own gene axis
    # 3 targets were dropped from the donor build for having <20 cells; score the rest
    if not carried.all():
        print(f"  (K562 donor table lacks {list(np.array(targets)[~carried])}; scoring {int(carried.sum())})")
    K = K[carried]
    S = state_delta("k562", [t for t, c in zip(targets, carried) if c])
    # proxy DE pool: the 200 genes each target actually moves most in the real data
    pool_k = np.zeros(K.shape, dtype=bool)
    for i in range(K.shape[0]):
        pool_k[i, np.argsort(-np.abs(K[i]))[:200]] = True
    res["k562"] = {}
    with_null(S, K, pool_k, "ST-HVG zeroshot/hepg2 (K562 in-train)", res["k562"])
    res["k562"]["discrimination_cosine"] = discrimination(S, K)
    print(f"  {'discrimination (cosine, 20 targets)':<34} {res['k562']['discrimination_cosine']:.4f}   (chance 0.5)")

    # ------------------------------------------------------------------------- HELD: H1
    Sf = Surface("experiments/164-callset-artefact/de164.a15.npz")
    pool_full = Sf.real_sig & Sf.adj
    all_targets = [str(p) for p in Sf.perts]
    sT_full = np.sign(Sf.real_lfc)

    # known-answer row 2: reproduce 175's production consensus number on the FULL pool
    D = {n: DonorHandle(n).rows(all_targets, Sf.genes) for n in DONORS}
    cons_full = sum(v[0].astype(np.float64) for v in D.values()) / 3
    a_full, n_full = acc(cons_full, Sf.real_lfc, pool_full)
    # 175's headline table is on the subset where ALL THREE donors carry a nonzero entry
    # (its stated n = 39,945), not on the raw pool -- reproduce that exact mask.
    allc = np.all([v[1] for v in D.values()], axis=0)
    m3 = (
        pool_full
        & allc[:, None]
        & np.all([np.sign(v[0]) != 0 for v in D.values()], axis=0)
    )
    a_m3, n_m3 = acc(cons_full, Sf.real_lfc, m3)
    a_k, _ = acc(D["k562_gwps"][0].astype(np.float64), Sf.real_lfc, m3)
    print("\nKNOWN-ANSWER ROW -- production consensus on the h1-test pool (reproducing 175)")
    print(f"  equal-mean, raw pool               {a_full:.4f}   (n={n_full})")
    print(f"  equal-mean, all-three-carrying     {a_m3:.4f}   (n={n_m3})   175 measured 0.5697 at n=39,945")
    print(f"  K562 alone, all-three-carrying     {a_k:.4f}                 175 measured 0.5670")
    res["known_answer_175"] = {
        "equal_mean_raw_pool": a_full,
        "n_raw_pool": n_full,
        "equal_mean_all_three": a_m3,
        "n_all_three": n_m3,
        "expected_equal_mean_175": 0.5697,
        "k562_alone_all_three": a_k,
        "expected_k562_alone_175": 0.5670,
    }

    # restrict to the 20 expressible targets and the model's gene axis
    ti = [all_targets.index(t) for t in targets]
    gi = np.array([i for i, g in enumerate(Sf.genes) if g in set(mgenes)])
    gsel = {g: i for i, g in enumerate(mgenes)}
    mi = np.array([gsel[Sf.genes[i]] for i in gi])
    T = Sf.real_lfc[np.ix_(ti, gi)]
    pool = pool_full[np.ix_(ti, gi)]
    print(f"\nHELD ROW -- H1, a line the checkpoint NEVER saw")
    print(f"  restricted pool: {len(ti)} targets x {len(gi)} shared genes, {int(pool.sum())} pool pairs")
    res["h1"] = {"n_targets": len(ti), "n_genes": len(gi), "n_pool_pairs": int(pool.sum())}

    consR = cons_full[np.ix_(ti, gi)]
    with_null(consR, T, pool, "production consensus (K562+Orion)", res["h1"])
    with_null(D["k562_gwps"][0].astype(np.float64)[np.ix_(ti, gi)], T, pool, "K562 donor alone", res["h1"])
    for tag, label in (("h1", "ST-HVG zeroshot (own-median norm)"), ("h1_scalematched", "ST-HVG zeroshot (scale-matched)")):
        P = state_delta(tag, targets)[:, mi]
        with_null(P, T, pool, label, res["h1"])
        res["h1"][f"{label}_discrimination"] = discrimination(P, T)

    # same-pool discrimination head-to-head (a pds analogue over these 20 targets)
    print("\n  discrimination on the SAME 20 targets (cosine rank, chance 0.5):")
    res["h1"]["discrimination_consensus"] = discrimination(consR, T)
    print(f"    {'production consensus':<32} {res['h1']['discrimination_consensus']:.4f}")
    for label in ("ST-HVG zeroshot (own-median norm)", "ST-HVG zeroshot (scale-matched)"):
        print(f"    {label:<32} {res['h1'][label + '_discrimination']:.4f}")

    json.dump(res, open(f"{WORK}/scores.json", "w"), indent=2, default=float)


if __name__ == "__main__":
    main()
