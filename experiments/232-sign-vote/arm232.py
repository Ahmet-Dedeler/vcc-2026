"""232 Phases 2 and 3 — the sign arms, priced through all six members on PLOLO.

THE ARM, AND WHY IT IS A CLEAN TEST OF §8.1d
--------------------------------------------
Every arm here is built as

    final_arm = final_prod       everywhere except a mask M,
    final_arm[M] = |final_prod[M]| * s_target[M]

so `|effect|` is **elementwise bit-identical** to the incumbent's -- asserted at runtime by
`common232.impose_sign`, which raises if a single magnitude moves.  Only signs change, and
only on coordinates with **at least two covering donors** (a gene one donor carries has no
vote and is left alone).  That is what makes this the first arm in the programme that
§8.1d's rule classifies as a SIGN arm rather than a magnitude arm.

ARMS
  BASE                the shipped `h1union` stack, `t = 0`.  Known-answer row.
  vote_z@tau          the confidence-weighted sign vote (weights `|e_s|/sd_s(gene)`), imposed
                      where it DISAGREES with `sign(mean)` and its normalised margin >= tau.
  vote_maj@tau        the unweighted majority, same wiring.  At arity 2 its mask is EMPTY by
                      construction (every contested pair is a 1-1 tie), which is itself the
                      result; it is a real rule only on `kolf2` and the `+kolf2` voter set.
  ctrl_mean@tau       CONTROL: `sign(mean)` imposed on the IDENTICAL coordinate set.  It
                      absorbs everything the coordinate selection does on its own, so
                      `vote - ctrl_mean` is the vote and nothing else.
  null_rand@tau       FALSIFIER 1: the same NUMBER of signs per target flipped at random,
                      drawn from the coordinates with >= 2 covering donors.  Seed 20260904.
  null_wrong@tau      FALSIFIER 2: target pi(i)'s mask and vote signs applied to target i
                      under a derangement.  184 died exactly here.
  oracle@f            INSTRUMENT, not a shippable arm: a fraction `f` of the wrong signs
                      inside the reference's DE pool set to the TRUTH's sign.  It is the only
                      way to measure what a sign improvement is WORTH -- §8.1d predicts
                      `nmae` falls and `e = dM/dU` clears the 0.5 break-even.  With
                      `antioracle@f` it gives a two-sided dose response and the
                      **sign-accuracy -> overall-points exchange rate**, which nothing in the
                      programme has ever measured.

    uv run python arm232.py k562_gwps
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common232  # noqa: E402
from common232 import (ALL_FOLDS, N_REALS, POOL_N, WORK, FoldBundle,  # noqa: E402
                       block_digest, call_set, control_block, ctx_gate, de_pools,
                       deliver, flip_mask, impose_sign, key_precision, load_sources,
                       stack_sources, to_lfc, topn, umd, vote_bundle,
                       weighted_sign_accuracy)

SEED = 20260904
TAUS = (0.0, 0.1, 0.2, 0.3, 0.5, 0.7)
ORACLE_F = (0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.00)


def _mask_target_counts(mask, s_target, final):
    """How many coordinates actually CHANGE sign (a mask entry whose sign already agrees
    with the target is a no-op).  Reported per target so the nulls can match it exactly."""
    changed = mask & (np.sign(final) != s_target) & (s_target != 0)
    return changed, changed.sum(1)


def build_arm(name, tau, vb, final, bundle, pools, rng):
    """Return (mask, s_target, audit).  `mask` may be empty -- that is a result, not a bug."""
    sc = bundle.sc
    if name in ("vote_z", "vote_maj"):
        rule = "z" if name == "vote_z" else "maj"
        m, s = flip_mask(vb, rule, tau)
        return m, s, dict(rule=rule)
    if name in ("ctrl_mean", "ctrl_mean_maj"):
        m, _ = flip_mask(vb, "z" if name == "ctrl_mean" else "maj", tau)
        return m, vb["s_mean"], dict(rule="mean-on-the-vote's-own-coordinates")
    if name in ("null_rand", "null_rand_maj"):
        m0, s0 = flip_mask(vb, "z" if name == "null_rand" else "maj", tau)
        changed, per_target = _mask_target_counts(m0, s0, final)
        cand = vb["ncov"] >= 2
        m = np.zeros_like(m0)
        for i in range(final.shape[0]):
            idx = np.flatnonzero(cand[i])
            k = int(min(per_target[i], idx.size))
            if k:
                m[i, rng.choice(idx, size=k, replace=False)] = True
        return m, -np.sign(final).astype(np.int8), dict(
            rule="random flip, count-matched", n_matched=int(per_target.sum()))
    if name in ("null_wrong", "null_wrong_maj"):
        m0, s0 = flip_mask(vb, "z" if name == "null_wrong" else "maj", tau)
        n = final.shape[0]
        perm = rng.permutation(n)
        while (perm == np.arange(n)).any():          # a true derangement
            perm = rng.permutation(n)
        m = m0[perm] & (vb["ncov"] >= 2)
        return m, s0[perm], dict(rule="wrong recipient", derangement=perm.tolist())
    if name in ("oracle", "antioracle"):
        pool = pools["wald89"]
        st = np.sign(bundle.truth)
        sf = np.sign(final[:, sc])
        want = (st != 0) & ((sf != st) if name == "oracle" else (sf == st))
        cand = pool & want
        m_sc = np.zeros_like(cand)
        for i in range(cand.shape[0]):
            idx = np.flatnonzero(cand[i])
            k = int(round(tau * idx.size))
            if k:
                m_sc[i, rng.choice(idx, size=k, replace=False)] = True
        m = np.zeros(final.shape, dtype=bool)
        m[:, sc] = m_sc
        s = np.zeros(final.shape, dtype=np.int8)
        s[:, sc] = (st if name == "oracle" else -st).astype(np.int8)
        return m, s, dict(rule=name, f=tau, n_candidates=int(cand.sum()))
    raise ValueError(name)


def score_arm(bundle, final, pools, ctrl, gate, bctx) -> dict:
    sc, truth = bundle.sc, bundle.truth
    out = bundle.score(final)

    H_A = deliver(final, bctx["A"], gate["A"], None)         # a15 render, 188's frame
    LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
    b = bundle.fold["b_ctrl"]
    LR = to_lfc(bundle.truth_h, b, gate["A"])[:, sc]
    nm = {k: umd(LP, LR, p) for k, p in pools.items()}
    out["nmae_pools"] = {k: float(np.mean(1 + (u - 2 * m) / d)) for k, (u, m, d) in nm.items()}
    out["nmae_U"] = {k: float(np.mean(v[0] / v[2])) for k, v in nm.items()}
    out["nmae_M"] = {k: float(np.mean(v[1] / v[2])) for k, v in nm.items()}
    out["reach_prec89"] = float(key_precision(np.abs(H_A[:, sc]), H_A[:, sc], truth, n=POOL_N))
    out["q_w_pool"] = {k: weighted_sign_accuracy(np.sign(H_A[:, sc]).astype(np.int8),
                                                 truth, p)[0] for k, p in pools.items()}
    out["delivered_median_norm"] = float(np.median(np.linalg.norm(H_A[:, bundle.wcols], axis=1)))
    del H_A, LP

    real_de = {n: topn(np.abs(truth), n) for n in N_REALS}
    fid, jac, digests = [], {n: [] for n in N_REALS}, []
    for c in ("A", "B", "C"):
        keep, ref_lin = ctrl[c]
        ok = tot = 0
        calls = np.zeros((final.shape[0], len(sc)), dtype=bool)
        for i in range(final.shape[0]):
            p = bundle_profile(final[i], bctx[c], gate[c])
            digests.append(block_digest(p))
            h_row = np.log1p(50_000.0 * p) - bctx[c]
            call, _ = call_set(final[i], gate[c], keep, ref_lin)
            cs = call[sc]
            calls[i] = cs
            s_pred, s_true = np.sign(h_row[sc]), np.sign(truth[i])
            live = cs & (s_true != 0) & (s_pred != 0)
            ok += int((s_pred[live] == s_true[live]).sum())
            tot += int(live.sum())
        fid.append(ok / max(tot, 1))
        for n in N_REALS:
            inter = (calls & real_de[n]).sum(1)
            union = (calls | real_de[n]).sum(1)
            jac[n].append(float(np.mean(inter / np.maximum(union, 1))))
    out["fid_precision"] = float(np.mean(fid))
    out["jac"] = {str(n): float(np.mean(v)) for n, v in jac.items()}
    out["digest"] = hashlib.sha256("".join(digests).encode()).hexdigest()[:16]
    out["n_blocks"] = len(digests)
    out["_block_digests"] = digests
    return out


def bundle_profile(effect_row, b, g):
    return common232.plolo.profile_from_effect(effect_row, b, g, "minnorm")


def main() -> None:
    t0 = time.time()
    folds = [a for a in sys.argv[1:] if not a.startswith("-")] or list(ALL_FOLDS)
    out_path = WORK / "arm232.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    ctrl = {c: control_block(c) for c in ("A", "B", "C")}
    bctx, gate = {}, {}
    for c in ("A", "B", "C"):
        bctx[c], gate[c] = ctx_gate(c)

    for held in folds:
        b = FoldBundle(held)
        src = load_sources(held)
        pools = de_pools(b)
        base_final = b.final()
        dig_path = WORK / f"digest-BASE-{held}.json"
        base_digests = json.loads(dig_path.read_text()) if dig_path.exists() else None

        for voters in (("production", False), ("production+kolf2", True)):
            tag, use_extra = voters
            if use_extra and "extra_E" not in src:
                continue
            names, E, avail, supp = stack_sources(src, use_extra=use_extra)
            vb = vote_bundle(E, avail, supp)
            del E

            plan = []
            for tau in TAUS:
                for nm in ("vote_z", "vote_maj", "ctrl_mean", "null_rand", "null_wrong"):
                    plan.append((nm, tau))
                # the majority is a real rule only at arity >= 3; give it its OWN matched
                # control and nulls there (the ones above are count-matched to the z-vote)
                if len(names) >= 3:
                    for nm in ("ctrl_mean_maj", "null_rand_maj", "null_wrong_maj"):
                        plan.append((nm, tau))
            if not use_extra:
                plan += [("oracle", f) for f in ORACLE_F]
                plan += [("antioracle", f) for f in ORACLE_F]

            if not use_extra and (held, "BASE") not in done:
                e = score_arm(b, base_final, pools, ctrl, gate, bctx)
                base_digests = e.pop("_block_digests")
                dig_path.write_text(json.dumps(base_digests))
                rows.append(dict(fold=held, voters="-", arm="BASE", tau=0.0,
                                 n_flipped=0, n_mask=0, mag_invariant=True,
                                 n_targets=int(base_final.shape[0]), **e))
                out_path.write_text(json.dumps(rows, indent=1))
                print(f"[{held}/BASE] nmae {e['nmae_pools']['wald89']:.5f} "
                      f"pds {e['pds']:.6f} cos_agg {e['cos_agg']:.6f} "
                      f"q_w {e['q_w_pool']['wald89']:.5f} digest {e['digest']}", flush=True)

            for nm, tau in plan:
                arm = f"{tag}|{nm}@{tau:g}"
                if (held, arm) in done:
                    continue
                rng = np.random.default_rng(
                    SEED + zlib.crc32(f"{held}|{arm}".encode()) % 1_000_003)
                mask, s_target, aud = build_arm(nm, tau, vb, base_final, b, pools, rng)
                changed, per_t = _mask_target_counts(mask, s_target, base_final)
                n_ch = int(changed.sum())
                if n_ch == 0:
                    rows.append(dict(fold=held, voters=tag, arm=arm, tau=float(tau),
                                     n_flipped=0, n_mask=int(mask.sum()),
                                     mag_invariant=True, empty=True, audit=aud))
                    out_path.write_text(json.dumps(rows, indent=1))
                    print(f"[{held}/{arm}] EMPTY MASK (mask {int(mask.sum())}, "
                          f"0 sign changes) -- identical to BASE by construction", flush=True)
                    continue
                fin = impose_sign(base_final, mask, s_target)
                e = score_arm(b, fin, pools, ctrl, gate, bctx)
                dig = e.pop("_block_digests")
                nd = (sum(1 for x, y in zip(dig, base_digests) if x != y)
                      if base_digests else None)
                rows.append(dict(fold=held, voters=tag, arm=arm, tau=float(tau),
                                 n_flipped=n_ch, n_mask=int(mask.sum()),
                                 flips_per_target_mean=float(per_t.mean()),
                                 mag_invariant=True, n_blocks_differing_vs_BASE=nd,
                                 audit={k: v for k, v in aud.items()
                                        if k != "derangement"},
                                 n_targets=int(fin.shape[0]), **e))
                out_path.write_text(json.dumps(rows, indent=1))
                print(f"[{held}/{arm}] flips {n_ch:6d}  nmae {e['nmae_pools']['wald89']:.5f} "
                      f"pds {e['pds']:.6f} cos_agg {e['cos_agg']:.6f} "
                      f"q_w {e['q_w_pool']['wald89']:.5f} "
                      f"norm {e['delivered_median_norm']:.5f} ({time.time()-t0:.0f}s)",
                      flush=True)
                del fin
            del vb
        del b, src, base_final
    print(f"\n-> {out_path}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
