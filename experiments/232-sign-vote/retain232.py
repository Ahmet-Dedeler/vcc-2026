"""232 falsifier 3 — leave-top-3-out retention on COSINE and on `nmae` SEPARATELY, plus the
paired-target bootstrap.

222 §9 measured 14.6% retention on `pds` against 72.7% on cosine for the same arm: the
channels differ, so an arm whose gain is carried by three targets on one channel can be
distributed on another.  228's own table repeated the split (cosine 68.9%, `pds` -0.80).
The mandate makes the cosine/`nmae` pair mandatory here because `nmae` is the member §8.1d
says a sign arm should move the right way, so its concentration is the thing to check.

`pds` retention is carried too, since it is the largest positive term in the composite and
228 found it the least trustworthy one.

    uv run python retain232.py --arm 'production|vote_z@0'
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common232 import (ALL_FOLDS, POOL_N, WORK, FoldBundle, ctx_gate,  # noqa: E402
                       de_pools, deliver, flip_mask, impose_sign, load_sources,
                       stack_sources, to_lfc, umd, vote_bundle)

BOOT = 4000
SEED = 20260904


def per_target(final, bundle, bctx, gate, pool, LR):
    """Per-target x, a, d (the energy-weighted aggregate's terms), `pds`, and `nmae`."""
    sc = bundle.sc
    truth = bundle.truth
    n = final.shape[0]
    X = np.zeros(n); A = np.zeros(n); D = np.zeros(n); P = np.zeros(n)
    for c in ("A", "B", "C"):
        dl = deliver(final, bctx[c], gate[c], "minnorm")[:, sc]
        X += (dl * truth).sum(1); A += (dl * dl).sum(1); D += (truth * truth).sum(1)
        u = dl / np.maximum(np.linalg.norm(dl, axis=1, keepdims=True), 1e-300)
        v = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-300)
        S = u @ v.T
        own = np.diag(S).copy()
        P += ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    H_A = deliver(final, bctx["A"], gate["A"], None)
    LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
    U, M, Dn = umd(LP, LR, pool)
    nm = 1.0 + (U - 2.0 * M) / Dn
    return X / 3.0, A / 3.0, D / 3.0, P / 3.0, nm


def agg(x, a, d, keep=None):
    s = slice(None) if keep is None else keep
    return float(x[s].sum() / np.sqrt(max(a[s].sum(), 1e-300) * max(d[s].sum(), 1e-300)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["production|vote_z@0"])
    ap.add_argument("--folds", nargs="*", default=list(ALL_FOLDS))
    args = ap.parse_args()
    bctx, gate = {}, {}
    for c in ("A", "B", "C"):
        bctx[c], gate[c] = ctx_gate(c)

    out_path = WORK / "retain232.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    for held in args.folds:
        b = FoldBundle(held)
        src = load_sources(held)
        pools = de_pools(b)
        pool = pools["wald89"]
        LR = to_lfc(b.truth_h, b.fold["b_ctrl"], gate["A"])[:, b.sc]
        base = b.final()
        xb, ab, db, pb, nb = per_target(base, b, bctx, gate, pool, LR)

        for arm in args.arms:
            tag, _, spec = arm.partition("|")
            name, _, tt = spec.partition("@")
            names, E, avail, supp = stack_sources(src, use_extra=(tag != "production"))
            vb = vote_bundle(E, avail, supp)
            del E
            # `ctrl_mean` is the CONTROL for `vote_z`: the z-vote's own coordinate set,
            # with `sign(mean)` imposed there instead of the vote's sign.
            rule = "maj" if name.endswith("_maj") or name == "vote_maj" else "z"
            m, s = flip_mask(vb, rule, float(tt))
            if name.startswith("ctrl_mean"):
                s = vb["s_mean"]
            fin = impose_sign(base, m, s)
            xa, aa, da, pa, na = per_target(fin, b, bctx, gate, pool, LR)
            del fin, vb

            n = len(xb)
            full_b, full_a = agg(xb, ab, db), agg(xa, aa, da)
            d_full = full_a - full_b
            cb = xb / np.sqrt(np.maximum(ab * db, 1e-300))
            ca = xa / np.sqrt(np.maximum(aa * da, 1e-300))
            dcos, dpds, dnm = ca - cb, pa - pb, na - nb

            def retention(delta, better_is_negative=False):
                d = -delta if better_is_negative else delta
                top = np.argsort(-d)[:3]
                keep = np.ones(n, bool)
                keep[top] = False
                mean_all = float(delta.mean())
                if mean_all == 0:
                    return float("nan"), float("nan")
                return (float(delta[keep].mean() / mean_all),
                        float(delta[top].sum() / (n * mean_all)))

            top_c = np.argsort(-dcos)[:3]
            keep_c = np.ones(n, bool)
            keep_c[top_c] = False
            ret_cos = ((agg(xa, aa, da, keep_c) - agg(xb, ab, db, keep_c)) / d_full
                       if d_full != 0 else float("nan"))
            ret_pds, share_pds = retention(dpds)
            ret_nmae, share_nmae = retention(dnm, better_is_negative=True)

            # seeded per (fold, arm) so a row is reproducible on its own, independently of
            # which other arms happen to be in the same invocation
            rng = np.random.default_rng(SEED + zlib.crc32(f"{held}|{arm}".encode()) % 1_000_003)
            idx = rng.integers(0, n, size=(BOOT, n))
            bs = np.array([agg(xa[i], aa[i], da[i]) - agg(xb[i], ab[i], db[i]) for i in idx])
            bsp = dpds[idx].mean(1)
            bsn = dnm[idx].mean(1)

            r = dict(n_targets=int(n),
                     cos_agg_base=full_b, cos_agg_arm=full_a, d_cos_agg=float(d_full),
                     rel_cos_agg_pct=100.0 * d_full / abs(full_b),
                     d_pds=float(dpds.mean()), d_nmae=float(dnm.mean()),
                     nmae_base=float(nb.mean()), nmae_arm=float(na.mean()),
                     retention_cos_leave_top3=float(ret_cos),
                     retention_pds_leave_top3=float(ret_pds),
                     retention_nmae_leave_top3=float(ret_nmae),
                     top3_share_of_d_pds=float(share_pds),
                     top3_share_of_d_nmae=float(share_nmae),
                     boot_cos_ci=[float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))],
                     boot_cos_p_gt0=float((bs > 0).mean()),
                     boot_pds_ci=[float(np.quantile(bsp, 0.025)), float(np.quantile(bsp, 0.975))],
                     boot_pds_p_gt0=float((bsp > 0).mean()),
                     boot_nmae_ci=[float(np.quantile(bsn, 0.025)),
                                   float(np.quantile(bsn, 0.975))],
                     boot_nmae_p_lt0=float((bsn < 0).mean()))
            out.setdefault(arm, {})[held] = r
            out_path.write_text(json.dumps(out, indent=1))
            print(f"[{held}/{arm}] d_cos_agg {d_full:+.6f} ({r['rel_cos_agg_pct']:+.2f}%) "
                  f"ret_cos {ret_cos:.3f} | d_pds {r['d_pds']:+.6f} ret_pds {ret_pds:.3f} "
                  f"| d_nmae {r['d_nmae']:+.6f} ret_nmae {ret_nmae:.3f} | "
                  f"P(cos>0) {r['boot_cos_p_gt0']:.3f} P(pds>0) {r['boot_pds_p_gt0']:.3f} "
                  f"P(nmae<0) {r['boot_nmae_p_lt0']:.3f}", flush=True)
        del b, src, base
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
