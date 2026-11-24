"""157 Q2 — does experiment 154's ordering closure survive power matching?

154 screened confidence keys on `de_thinned.npz` at that surface's NATIVE reference
power (BH-significant set, median n_real 68) and concluded the best buildable key
captures 3.4% of the oracle headroom, therefore `reach` is model-limited.

`reach = k*/n_real` and `k*` is bounded by the pool size, so a tiny pool cannot express
deep reach for ANY key.  The challenge's `n_real` is now CONFIRMED small (~240, a21 board
result, docs/LOG.md 2026-09-01) but that is still 3.5x this harness's 68.  This re-runs
154's headline comparison under `Surface(..., n_real_target=N)` power matching and reports
how incumbent / best-buildable / oracle each scale with N.

Deterministic.  No emitter, no seed, no cell-eval2, no build, no submission.
"""
from __future__ import annotations
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface, overall_delta  # noqa: E402

DE = f"{ROOT}/experiments/095-clean-harness/de_thinned.npz"
PRED_TABLE = os.environ.get("PRED_TABLE", DE)
HARNESS_DONOR = "orion_hct116"  # excluded: it IS the evaluation surface
# The donor that IS the evaluation surface must be excluded. `h1_2025` covers 85/85
# h1-test targets, so leaving it in makes every donor key an oracle (it read 0.404
# reach before this was caught). `orion_hct116` is likewise excluded on de_thinned.
DONORS = os.environ.get("DONORS", "k562_gwps,orion_hek293t,h1_2025").split(",")
NS = [None, 68, 132, 222, 240, 400, 800]


def load_keys(genes, perts, pred_lfc):
    """154's candidate keys, rebuilt on this gene axis. Higher = more confident."""
    P, G = pred_lfc.shape
    mats, haves = [], []
    for name in DONORS:
        root = f"{ROOT}/data/donors/{name}"
        try:
            dg = np.array(open(f"{root}/genes.txt").read().split())
            dp = np.array(open(f"{root}/perts.txt").read().split())
            delta = np.load(f"{root}/delta.npy", mmap_mode="r")
        except FileNotFoundError:
            print(f"  donor {name}: missing, skipped")
            continue
        gi = {g: i for i, g in enumerate(dg)}
        cols = np.array([gi.get(g, -1) for g in genes])
        pi = {p: i for i, p in enumerate(dp)}
        M = np.zeros((P, G))
        have = np.zeros(P, bool)
        ok = cols >= 0
        for j, p in enumerate(perts):
            r = pi.get(p)
            if r is None:
                continue
            M[j, ok] = np.asarray(delta[r], float)[cols[ok]]
            have[j] = True
        print(f"  donor {name}: {have.sum()}/{P} perturbations, {ok.sum()}/{G} genes")
        mats.append(M)
        haves.append(have)
    stack = np.stack(mats)
    mask = np.stack([np.broadcast_to(h[:, None], (P, G)) for h in haves])
    n_avail = mask.sum(0).astype(float)
    z = np.zeros((P, G))
    sgn = np.sign(stack) * mask
    concord = np.divide(np.abs(sgn.sum(0)), n_avail, out=z.copy(), where=n_avail > 0)
    absmean = np.divide((np.abs(stack) * mask).sum(0), n_avail, out=z.copy(), where=n_avail > 0)
    geo = np.sqrt(np.maximum(absmean, 0.0))
    agree = np.divide(((np.sign(stack) == np.sign(pred_lfc)[None]) * mask).sum(0), n_avail,
                      out=z.copy(), where=n_avail > 0)
    a = np.abs(pred_lfc)
    return {
        "|pred_lfc|": a,
        "donor concordance": concord,
        "concordance x |pred_lfc|": concord * a,
        "concordance^2 x |pred_lfc|": concord ** 2 * a,
        "agreement x |pred_lfc|": agree * a,
        "geo-magnitude x concordance": geo * concord,   # 154's winner
    }


def main() -> None:
    print(f"pred table : {PRED_TABLE}")
    base = Surface(PRED_TABLE)
    print(f"surface    : {base.perts.size} perturbations x {base.genes.size} gated genes")
    keys = load_keys(base.genes, base.perts, base.pred_lfc)
    print()

    rows = {}
    hdr = ["incumbent (emitted p-order)"] + list(keys) + ["ORACLE (correct first)"]
    for N in NS:
        s = Surface(PRED_TABLE, n_real_target=N)
        nr, npd = s.counts()
        oracle = (np.sign(s.pred_lfc) == np.sign(s.real_lfc)).astype(float)
        vals = [s.reach()] + [s.reach(key=k) for k in keys.values()] + [s.reach(key=oracle)]
        rows[N] = (vals, float(np.median(nr)), float(np.median(npd)),
                   s.sign_accuracy())
    lab = {None: "native BH"}
    print(f"{'key':<30}" + "".join(f"{lab.get(N, f'N={N}'):>12}" for N in NS))
    print("-" * (30 + 12 * len(NS)))
    for i, name in enumerate(hdr):
        print(f"{name:<30}" + "".join(f"{rows[N][0][i]:>12.5f}" for N in NS))
    print("-" * (30 + 12 * len(NS)))
    print(f"{'median n_real':<30}" + "".join(f"{rows[N][1]:>12.0f}" for N in NS))
    print(f"{'median n_pred':<30}" + "".join(f"{rows[N][2]:>12.0f}" for N in NS))
    print(f"{'pool sign accuracy':<30}" + "".join(f"{rows[N][3]:>12.4f}" for N in NS))

    print("\n--- headroom capture: (best buildable - incumbent) / (oracle - incumbent) ---")
    print(f"{'N':>10}{'incumbent':>12}{'best key':>12}{'oracle':>12}"
          f"{'capture%':>10}{'d_overall':>11}")
    for N in NS:
        v = rows[N][0]
        inc, orc = v[0], v[-1]
        bi = int(np.argmax(v[1:-1])) + 1
        best = v[bi]
        cap = (best - inc) / max(orc - inc, 1e-12) * 100
        print(f"{lab.get(N, N)!s:>10}{inc:>12.5f}{best:>12.5f}{orc:>12.5f}{cap:>10.1f}"
              f"{overall_delta(d_reach=best - inc):>+11.5f}  <- {hdr[bi]}")

    print("\n--- same, with 154's own baseline (|pred_lfc|) as the reference point ---")
    print(f"{'N':>10}{'|pred_lfc|':>12}{'best key':>12}{'oracle':>12}{'capture%':>10}")
    for N in NS:
        v = rows[N][0]
        b, orc = v[1], v[-1]
        best = max(v[1:-1])
        print(f"{lab.get(N, N)!s:>10}{b:>12.5f}{best:>12.5f}{orc:>12.5f}"
              f"{(best - b) / max(orc - b, 1e-12) * 100:>10.1f}")


if __name__ == "__main__":
    main()
