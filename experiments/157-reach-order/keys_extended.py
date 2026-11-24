"""157 Q3 — confidence keys 154 did not screen, evaluated at power-matched n_real.

154 screened donor-derived keys only (concordance, agreement, geometric magnitude) times
`|pred_lfc|`.  Three families of BUILD-TIME signal were never tested and are all free:

  * the reference's OWN CONTROL CELLS, which the challenge hands us -- per-gene control CPM,
    control dispersion, and the control's detection rate.  These are legitimate model inputs,
    not reference-truth, and they are exactly what sets a gene's measurability;
  * the DONOR PANEL's marginal per-gene response strength (how much a gene moves across ALL
    donor perturbations), which is a gene-level prior independent of the target;
  * combinations of the above with `|pred_lfc|`.

Every key is scored with `debench.Surface.reach`, which is the scorer's own definition, at
several `n_real_target` values so the challenge's power regime is represented rather than
the harness's.  Deterministic.  No emitter here, no cell-eval2, no build, no submission.
"""
from __future__ import annotations
import os
import sys

import anndata as ad
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface, overall_delta  # noqa: E402

TABLE = os.environ.get("PRED_TABLE", f"{ROOT}/experiments/095-clean-harness/de_thinned.npz")
# The surface's OWN control cells, i.e. the deployment-time analogue of the challenge's
# `context_{A,B,C}.h5ad`. Never the reference's perturbed cells.
CTRLS = os.environ.get("CTRLS", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad")
# The donor that IS the evaluation surface must be excluded. `h1_2025` covers 85/85
# h1-test targets, so leaving it in makes every donor key an oracle (it read 0.404
# reach before this was caught). `orion_hct116` is likewise excluded on de_thinned.
DONORS = os.environ.get("DONORS", "k562_gwps,orion_hek293t,h1_2025").split(",")
NS = [None, 68, 240, 400]


def control_features(genes):
    """Per-gene statistics of the REAL challenge control cells (a legitimate model input)."""
    A = ad.read_h5ad(CTRLS)
    g = np.asarray([str(x) for x in A.var_names])
    pos = {x: i for i, x in enumerate(g)}
    idx = np.array([pos.get(x, -1) for x in genes])
    rows = np.arange(A.shape[0])
    if "target_gene" in A.obs:
        rows = np.flatnonzero(np.asarray(A.obs["target_gene"].astype(str)) == "non-targeting")
    n = rows.size
    s1 = np.zeros(g.size); s2 = np.zeros(g.size); nz = np.zeros(g.size)
    for i in range(0, n, 2000):
        m = A.X[rows[i:i + 2000]]
        m = np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)
        t = m.sum(1); t[t == 0] = 1.0
        c = m * (1e6 / t[:, None])
        s1 += c.sum(0); s2 += (c ** 2).sum(0); nz += (m > 0).sum(0)
    mu = s1 / n
    var = np.maximum(s2 / n - mu ** 2, 0.0)
    cv = np.sqrt(var) / np.maximum(mu, 1e-12)
    det = nz / n
    ok = idx >= 0
    out = {}
    for name, v in (("ctrl_cpm", mu), ("ctrl_cv", cv), ("ctrl_detect", det)):
        w = np.zeros(genes.size)
        w[ok] = v[idx[ok]]
        out[name] = w
    print(f"control features from {os.path.basename(CTRLS)}: {n} cells, "
          f"{int(ok.sum())}/{genes.size} genes matched")
    return out


def donor_features(genes, perts, pred_lfc):
    P, G = pred_lfc.shape
    mats, haves = [], []
    marg = np.zeros(G); nmarg = 0
    for name in DONORS:
        root = f"{ROOT}/data/donors/{name}"
        dg = np.array(open(f"{root}/genes.txt").read().split())
        dp = np.array(open(f"{root}/perts.txt").read().split())
        delta = np.load(f"{root}/delta.npy", mmap_mode="r")
        gi = {g: i for i, g in enumerate(dg)}
        cols = np.array([gi.get(g, -1) for g in genes])
        ok = cols >= 0
        pi = {p: i for i, p in enumerate(dp)}
        M = np.zeros((P, G)); have = np.zeros(P, bool)
        for j, p in enumerate(perts):
            r = pi.get(p)
            if r is None:
                continue
            M[j, ok] = np.asarray(delta[r], float)[cols[ok]]
            have[j] = True
        mats.append(M); haves.append(have)
        # gene-level marginal responsiveness over the WHOLE donor matrix, target-independent
        step = max(delta.shape[0] // 400, 1)
        sub = np.abs(np.asarray(delta[::step], dtype=np.float64))
        v = np.zeros(G); v[ok] = sub.mean(0)[cols[ok]]
        marg += v / max(v.mean(), 1e-12); nmarg += 1
    stack = np.stack(mats)
    mask = np.stack([np.broadcast_to(h[:, None], (P, G)) for h in haves])
    n_avail = mask.sum(0).astype(float)
    z = np.zeros((P, G))
    concord = np.divide(np.abs((np.sign(stack) * mask).sum(0)), n_avail, out=z.copy(),
                        where=n_avail > 0)
    absmean = np.divide((np.abs(stack) * mask).sum(0), n_avail, out=z.copy(), where=n_avail > 0)
    return concord, np.sqrt(np.maximum(absmean, 0.0)), marg / max(nmarg, 1)


def main() -> None:
    print(f"table: {TABLE}")
    s0 = Surface(TABLE)
    genes, perts = s0.genes, s0.perts
    a = np.abs(s0.pred_lfc)
    cf = control_features(genes)
    concord, geo, marg = donor_features(genes, perts, s0.pred_lfc)
    B = lambda v: np.broadcast_to(v, a.shape)          # noqa: E731  gene-level -> (P,G)

    KEYS = {
        "|pred_lfc| (154 baseline)": a,
        "geo x concord (154 best)": geo * concord,
        # --- never screened: control-cell features ---
        "ctrl_cpm": B(cf["ctrl_cpm"]),
        "ctrl_detect": B(cf["ctrl_detect"]),
        "-ctrl_cv": B(-cf["ctrl_cv"]),
        "|lfc| x log ctrl_cpm": a * B(np.log1p(cf["ctrl_cpm"])),
        "|lfc| x ctrl_detect": a * B(cf["ctrl_detect"]),
        "|lfc| / ctrl_cv": a / B(np.maximum(cf["ctrl_cv"], 1e-3)),
        # --- never screened: gene-level donor responsiveness prior ---
        "donor marginal |delta|": B(marg),
        "|lfc| x donor marginal": a * B(marg),
        "|lfc| / donor marginal": a / B(np.maximum(marg, 1e-3)),
        # --- best-of combinations ---
        "geo x concord x ctrl_detect": geo * concord * B(cf["ctrl_detect"]),
        "geo x concord x log cpm": geo * concord * B(np.log1p(cf["ctrl_cpm"])),
    }

    lab = {None: "native BH"}
    hdr = ["incumbent (emitted p-order)"] + list(KEYS) + ["ORACLE"]
    res = {}
    for N in NS:
        s = Surface(TABLE, n_real_target=N)
        orc = (np.sign(s.pred_lfc) == np.sign(s.real_lfc)).astype(float)
        res[N] = [s.reach()] + [s.reach(key=k) for k in KEYS.values()] + [s.reach(key=orc)]
    print(f"\n{'key':<30}" + "".join(f"{lab.get(N, f'N={N}'):>12}" for N in NS)
          + f"{'d_overall@240':>15}")
    print("-" * (30 + 12 * len(NS) + 15))
    i240 = NS.index(240)
    inc240 = res[240][0]
    for i, name in enumerate(hdr):
        d = overall_delta(d_reach=res[240][i] - inc240)
        print(f"{name:<30}" + "".join(f"{res[N][i]:>12.5f}" for N in NS) + f"{d:>+15.5f}")


if __name__ == "__main__":
    main()
