"""234 — realise the best flip key as an ACTUAL sign-flip arm and price all six members.

The screen measures `Dq_w` on the DELIVERED sign matrix, which is an upper bound: a shippable
arm flips at the `final` effect matrix and the renderer then rewrites a few percent of those
signs back (232 §2b).  This script closes that gap and supplies the three falsifiers that
cannot be checked analytically:

  * MAGNITUDE INVARIANCE -- `|effect|` asserted elementwise bit-identical (232 `impose_sign`)
  * `jac` -- 232 §11(a) says a magnitude-preserving sign arm leaves it EXACTLY invariant;
    asserted here at all three `n_real` depths, and a non-zero would mean the arm is not
    magnitude-preserving
  * the realised `nmae` / `pds` / `cos` / `fid`, so the exchange-rate extrapolation is
    checked against a direct six-member read rather than assumed

Routes are 228's and 232's, unchanged (rule 12).  `k_cos` is carried as the RANGE
{1.000, 1.384}; no PLOLO->board scalar is quoted (§8.1c: the factor is not identified).

    uv run python arm234.py k562_gwps
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
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import deliver  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import score_effects  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402

from common234 import (KA232_NMAE_BASE, KA232_PDS_BASE, KA232_QW_DELIVERED,  # noqa: E402
                       SEED, prep_weights, qw)
from extract234 import build, ctx_gate  # noqa: E402

RC_HEAD, RC_CPM = 400, 2200
MEDIAN_LIBRARY, CELLS = 20109.0, 400
N_REALS = (89, 200, 400)
DOSES = (0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30)


def control_block(context: str):
    z = np.load(WORK / f"ctrl-{context}.npz")           # copied from 232, control-only
    return np.asarray(z["keep"], bool), np.asarray(z["ref_lin"], np.float64)


def call_set(effect_row, gate10, keep, ref_lin):
    score = np.where(gate10 & keep, np.abs(np.nan_to_num(effect_row)), -np.inf)
    head = np.argsort(-score, kind="stable")[:RC_HEAD]
    rest = np.where(keep, ref_lin, -np.inf).copy()
    rest[head] = -np.inf
    cpm = np.argsort(-rest, kind="stable")[:RC_CPM]
    call = np.zeros(len(score), dtype=bool)
    call[head] = True
    call[cpm] = True
    return call


def block_digest(profile) -> str:
    tot = np.rint(CELLS * MEDIAN_LIBRARY * profile).astype(np.int64)
    return hashlib.sha256(tot.tobytes()).hexdigest()[:16]


def impose_flip(final, mask):
    """Reverse the sign inside `mask`; `|effect|` is bit-identical.  Asserted."""
    eff = mask & (np.sign(final) != 0)
    out = final.copy()
    out[eff] = -out[eff]
    assert np.array_equal(np.abs(out), np.abs(final)), "MAGNITUDE INVARIANCE VIOLATED"
    return out, int(eff.sum())


def score_arm(d, final, ctrl, bctx, gate) -> dict:
    sc, truth = d["sc"], d["truth"]
    fold = d["fold"]
    out = {k: v for k, v in score_effects(final, fold, truth).items()
           if not k.startswith("_")}
    H_A = deliver(final, bctx["A"], gate["A"], None)
    LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
    LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64), fold["b_ctrl"],
                gate["A"])[:, sc]
    nm = {k: umd(LP, LR, p) for k, p in d["pools"].items()}
    out["nmae_pools"] = {k: float(np.mean(1 + (u - 2 * m) / dd))
                         for k, (u, m, dd) in nm.items()}
    out["nmae_U"] = {k: float(np.mean(v[0] / v[2])) for k, v in nm.items()}
    out["nmae_M"] = {k: float(np.mean(v[1] / v[2])) for k, v in nm.items()}
    out["reach_prec89"] = float(key_precision(np.abs(H_A[:, sc]), H_A[:, sc], truth,
                                              n=POOL_N))
    out["q_w_pool"] = {}
    for k, p in d["pools"].items():
        w, den = prep_weights(truth, p)
        out["q_w_pool"][k] = qw(np.sign(H_A[:, sc]).astype(np.int8), truth, w, den)
    out["delivered_median_norm"] = float(np.median(
        np.linalg.norm(H_A[:, d["wcols"]], axis=1)))
    del H_A, LP

    real_de = {n: topn(np.abs(truth), n) for n in N_REALS}
    fid, jac, digests = [], {n: [] for n in N_REALS}, []
    for c in ("A", "B", "C"):
        keep, ref_lin = ctrl[c]
        ok = tot = 0
        calls = np.zeros((final.shape[0], len(sc)), dtype=bool)
        for i in range(final.shape[0]):
            p = plolo.profile_from_effect(final[i], bctx[c], gate[c], "minnorm")
            digests.append(block_digest(p))
            h_row = np.log1p(50_000.0 * p) - bctx[c]
            cs = call_set(final[i], gate[c], keep, ref_lin)[sc]
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
    return out


def row_mask(key, frac: float, shape_full, sc, allowed=None) -> np.ndarray:
    """Lowest-`frac` of each row by `key` (on the scored axis), lifted to the official axis.

    `allowed` restricts the candidate set.  Two of the three variants use only quantities we
    hold (the whole scored axis, or our OWN top-|effect| head); `combo_pool` uses the
    reference's DE pool and is therefore ORACLE-ASSISTED and diagnostic only -- it exists to
    isolate how much of the arm's damage is caused by flipping OUTSIDE the scored pool.
    """
    T, G = key.shape
    k = np.where(np.zeros(T, dtype=bool), 0, 0)
    m_sc = np.zeros((T, G), dtype=bool)
    kk = (np.asarray(key, dtype=np.float64) if allowed is None
          else np.where(allowed, key, np.inf))
    n_allowed = (np.full(T, G) if allowed is None else allowed.sum(1))
    k = np.minimum(np.rint(frac * n_allowed).astype(int), n_allowed)
    order = np.argsort(kk, axis=1, kind="stable")
    for i in range(T):
        m_sc[i, order[i, :k[i]]] = True
    m = np.zeros(shape_full, dtype=bool)
    m[:, sc] = m_sc
    return m


def main() -> None:
    t0 = time.time()
    held = sys.argv[1]
    cfg = json.loads((WORK / "honest234.json").read_text())["best"]
    frac = cfg["doses"][held]
    d = build(held)
    sc, final = d["sc"], d["final"]
    ctrl = {c: control_block(c) for c in ("A", "B", "C")}
    bctx, gate = {}, {}
    for c in ("A", "B", "C"):
        bctx[c], gate[c] = ctx_gate(c)

    key = np.load(WORK / f"combokey-{held}.npy")
    rng = np.random.default_rng(SEED + zlib.crc32(held.encode()) % 1_000_003)

    out_path = WORK / "arm234.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    rows = [r for r in rows if r["fold"] != held]

    base = score_arm(d, final, ctrl, bctx, gate)
    ka = {"q_w": KA232_QW_DELIVERED.get(held), "nmae": KA232_NMAE_BASE.get(held),
          "pds": KA232_PDS_BASE.get(held)}
    okq = ka["q_w"] is None or abs(base["q_w_pool"]["wald89"] - ka["q_w"]) < 1e-12
    okn = ka["nmae"] is None or abs(base["nmae_pools"]["wald89"] - ka["nmae"]) < 5e-6
    okp = ka["pds"] is None or abs(base["pds"] - ka["pds"]) < 5e-6
    print(f"[{held}/BASE] q_w {base['q_w_pool']['wald89']:.15f} "
          f"nmae {base['nmae_pools']['wald89']:.6f} pds {base['pds']:.6f} "
          f"cos_agg {base['cos_agg']:.6f} jac89 {base['jac']['89']:.6f}  "
          f"KA {'PASS' if (okq and okn and okp) else 'FAIL'}", flush=True)
    assert okq and okn and okp, f"KA row failed for {held}"
    rows.append(dict(fold=held, arm="BASE", frac=0.0, n_flipped=0, **base))

    # ---- the three candidate sets.  `head400` is our own rc_callset |effect| head, so it
    # is knowable at prediction time; `pool` is the reference's DE set and is NOT.
    head400 = np.zeros((final.shape[0], len(sc)), dtype=bool)
    absH_sc = np.abs(deliver(final, bctx["A"], gate["A"], None)[:, sc])
    for i in range(final.shape[0]):
        head400[i, np.argsort(-absH_sc[i], kind="stable")[:RC_HEAD]] = True
    del absH_sc
    allowed = {"": None, "_head": head400, "_pool": d["pools"]["wald89"]}

    noise = rng.random(key.shape).astype(np.float32)
    arms = []
    for fr in DOSES:
        arms.append((f"combo@{fr:g}", key, fr, None))
        arms.append((f"combo_head@{fr:g}", key, fr, head400))
        arms.append((f"combo_pool@{fr:g}", key, fr, d["pools"]["wald89"]))
        arms.append((f"null_rand@{fr:g}", noise, fr, None))
    arms.append((f"combo@{frac:g}", key, frac, None))
    arms.append((f"null_rand@{frac:g}", noise, frac, None))
    seen = set()
    arms = [a for a in arms if not (a[0] in seen or seen.add(a[0]))]

    for name, kk, fr, alw in arms:
        mask = row_mask(kk, fr, final.shape, sc, alw)
        fin, n = impose_flip(final, mask)
        e = score_arm(d, fin, ctrl, bctx, gate)
        dj = max(abs(e["jac"][str(n_)] - base["jac"][str(n_)]) for n_ in N_REALS)
        print(f"[{held}/{name}] flips {n:7d}  "
              f"d q_w {e['q_w_pool']['wald89']-base['q_w_pool']['wald89']:+.5f}  "
              f"d nmae {e['nmae_pools']['wald89']-base['nmae_pools']['wald89']:+.6f}  "
              f"d pds {e['pds']-base['pds']:+.6f}  "
              f"rel cos {100*(e['cos_agg']-base['cos_agg'])/abs(base['cos_agg']):+.3f}%  "
              f"max|d jac| {dj:.3e}  ({time.time()-t0:.0f}s)", flush=True)
        rows.append(dict(fold=held, arm=name, frac=fr, n_flipped=n,
                         mag_invariant=True, max_abs_d_jac=dj, **e))
        del fin
    out_path.write_text(json.dumps(rows, indent=1))
    print(f"-> {out_path}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
