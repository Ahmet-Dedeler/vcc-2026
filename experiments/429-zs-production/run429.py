#!/usr/bin/env python3
"""429 — 426's harness (copied verbatim: run426.py, k_emit426.py -> k_emit429.py) plus 423's exact-U down-gene pass
(`--zs dn_k250_zs`) applied after the ballast block and after the optional zero repair (428), i.e. in the builder's
run order. One process per (design, seed, MDE, repair): the final pre-pass blocks are scored as `<tag>` and the same
blocks after the pass as `<tag>_<variant>` (423's same-block design), with 423's verify record per target.
--zs-impl vendored (default) runs this directory's src_rank_mag.py (the production module); rp423 runs 423's file.
Known answers: `base S --bulk ps` == 416's base_s{S}_bps; `... --zs dn_k250_zs` == 423's dn_k250_zs_s{S}_bps;
`... --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3` == 426's base_s{S}_bps_mde0.1cap0.1min0.3 (stats bit-for-bit).

426's doc: 426 — MDE-structured magnitude on the delivered LFC (A-channel only), on 416's H1 harness. One variable vs 416's
`base --bulk ps` (the shipped 410b construction): the prior spec carries `mde_add` = c / sqrt(expected control count per
cell), applied in the sign of each steered gene's delivered target and mass-balanced (k_emit426.py). --mde-c 0 must
reproduce 416's base_s{seed}_bps bit-for-bit.

416's doc: 416 — DE-carrier menu on H1 under ballast, seed-replicated. Copy of 410's run410b.py (ballast, ratio mode, jackknife
boost) with the CARRIER design as the one variable and the emission seed as a replicate axis; 410's files are read-only.

Per target: 400 - K carrier cells emitted by the K-class block (k3c5bq80 unless the design says otherwise) from the
design's carrier effect; K ballast cells carry alloc(S, bulk composition) - carrier totals (ratio mode: the bulk profile
is applied as a ratio on the carriers' DELIVERED composition); pair-wise anti-correlated composition noise restores the
carriers' own jackknife. The bulk effect is the shipped cap15 matrix for every design (--bulk base) or 410b's key-scaled
bulk (--bulk ps, lo 0.3 / hi 3, energy-preserving: the shipped 410b construction).

usage: uv run python run416.py DESIGN SEED [--bulk base|ps] [--kfrom RUNTAG]
  DESIGN: base | g000 | g030 | g090 | gain | signrep | prior_off | q70 | q90 | conf4
          follow-ups: pre (410's de_unwhite: gamma 0 without CD4) | w065nocd4 | gm030 (prep416b.py)
  --kfrom: force each target's ballast count K to the one recorded in work/h1/diag-RUNTAG.json (K-matched control)
known answers: `x125 0` == 410's de_x1p25_ratio_jk (carriers x1.25, cap15 bulk, K 1-3, jk boost live);
               `base 0 --bulk ps` == 410's bal_lo03_hi3_ratio_jk (the shipped 410b construction).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "h1"
WORK.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import run400 as H  # noqa: E402  (read-only)
if os.environ.get("RUN426_EMITTER") == "prod":
    from vcc2026 import rc_emit_k as K  # noqa: E402  (the production emitter: 426's known answers)
else:
    import k_emit429 as K  # noqa: E402  (k_emit426.py verbatim: 416's emitter + optional mde_add)

H.WORK = WORK
R, E = H.R, H.E
LB = 999_000
TS = 50_000.0
OP_H1 = ROOT / "experiments/323-current-family-validation/work/whitener_h1_pruned.npz"
STAGE = ROOT / "experiments/323-current-family-validation/work/stage.npz"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CARRIERS = HERE / "work/carriers416.npz"
CLEAN4 = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")

PRIOR_MODE = "q"      # "q" (k3c5bq80 family) | "off" (k3plain: K-class, unsteered)
PROTECT_Q = 0.8
GATE_FULL = None      # full-axis bool: prior may flip only here (conf4)
MDE_ADD = None        # gated-axis array: 426's c / sqrt(count); None = off (416 verbatim)
MDE_SCOPE = "all"
MDE_CONC = False
MDE_SIGN = None       # gated-axis +-1/0: 426b's sign source for non-block genes
MDE_OVERFILL = None   # cap on each class's over-fill / sub-gate mass under the MDE dose
MDE_AGREE = None      # gated-axis +-1/0: 426e's agreement gate (covered CLEAN4 prior sign, 0 = silent)


# ---------------------------------------------------------------- 410's ballast machinery, verbatim (run410b.py)
def alloc(total: int, profile: np.ndarray) -> np.ndarray:
    raw = profile * total
    base = np.floor(raw).astype(np.int64)
    rem = int(total - base.sum())
    if rem > 0:
        idx = np.argsort(-(raw - base))[:rem]
        base[idx] += 1
    return base


def jk_sum(block):
    X = np.asarray(block, dtype=np.float64)
    P = X.sum(0); S = P.sum(); L = X.sum(1)
    V = np.log1p(TS * (P[None, :] - X) / (S - L)[:, None])
    n = X.shape[0]
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())


def _balance_pair_noise(eps):
    d = int(eps.sum())
    if d == 0:
        return eps
    idx = np.flatnonzero(eps > 0) if d > 0 else np.flatnonzero(eps < 0)
    order = idx[np.argsort(-np.abs(eps[idx]), kind="stable")]
    mags = np.abs(eps[order])
    cum = np.cumsum(mags)
    stop = int(np.searchsorted(cum, abs(d)))
    eps[order[:stop]] = 0
    left = abs(d) - (int(cum[stop - 1]) if stop > 0 else 0)
    if left:
        eps[order[stop]] -= left if d > 0 else -left
    if int(eps.sum()) != 0:
        raise AssertionError("pair noise balance failed")
    return eps


def jk_boost(block, n_car, K_, target, rng):
    if K_ < 2:
        return block, 0.0, jk_sum(block)
    base = block[n_car:].copy()
    G = base.shape[1]
    z = np.where(rng.random(G) < 0.5, -1, 1)

    def build(rho):
        out = base.copy()
        for j in range(0, K_ - 1, 2):
            a, b = out[j], out[j + 1]
            m = np.minimum(a, b)
            eps = np.floor(rho * m).astype(np.int64) * z
            eps = _balance_pair_noise(eps)
            out[j] = a + eps; out[j + 1] = b - eps
        return np.vstack([block[:n_car], out])
    lo, hi = 0.0, 0.95
    cur = jk_sum(block)
    if cur >= target:
        return block, 0.0, cur
    best = build(hi); v = jk_sum(best)
    if v < target:
        return best, hi, v
    for _ in range(12):
        mid = 0.5 * (lo + hi); cand = build(mid); v = jk_sum(cand)
        if v >= target:
            hi, best = mid, cand
        else:
            lo = mid
    return best, hi, jk_sum(best)


def cap_mass(Q, K_, cap=995_000):
    excess = int(Q.sum()) - K_ * cap
    if excess <= 0:
        return Q, 0
    take = np.floor(Q * (excess / Q.sum())).astype(np.int64)
    short = excess - int(take.sum())
    if short > 0:
        take[np.argsort(-(Q - take), kind="stable")[:short]] += 1
    return Q - take, excess


# ---------------------------------------------------------------- the one variable: the carrier block
def carrier_block(profile, libs, emit_rng, spread, key, i, emit_c):
    cs = H.rc_emit.CallSet(*H.CALLSET)
    kw = dict(spread=spread, overdispersion=4., control=emit_c, key=key, jitter_seed=1000 + i, call_set_first=True,
              call_set=cs, emit_counts=H.EMIT_COUNTS)
    if PRIOR_MODE == "off":
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng, classes=list(H.CLASSES3), lfc_target=None, **kw)
    gated = np.asarray(emit_c.keep, dtype=bool)
    spec = {"prior": np.ascontiguousarray(H.PRIOR_FULL[gated]), "shrink_floor": 0.002, "protect_q": PROTECT_Q}
    if GATE_FULL is not None:
        spec["gate"] = np.ascontiguousarray(GATE_FULL[gated])
    if MDE_ADD is not None:
        spec["mde_add"] = np.ascontiguousarray(MDE_ADD)
        spec["mde_scope"] = MDE_SCOPE
        if MDE_CONC:
            spec["mde_concordant"] = True
        if MDE_SIGN is not None:
            spec["mde_sign"] = np.ascontiguousarray(MDE_SIGN)
        if MDE_OVERFILL is not None:
            spec["mde_overfill"] = MDE_OVERFILL
        if MDE_AGREE is not None:
            spec["mde_agree"] = np.ascontiguousarray(MDE_AGREE)
    return K.rc_hybrid_block_k(profile, libs, .1, emit_rng, classes=list(H.CLASSES3), lfc_target=spec, **kw)


def ballast_block(profile_de, profile_bulk, libs, emit_rng, spread, key, i, emit_c, kmax, k_force=None):
    live = profile_de > 1e-9
    ratio = float(np.min(profile_bulk[live] / profile_de[live]))
    s_c_est = float(np.sum(libs[: 400 - kmax]))
    need = s_c_est * (1.0 / max(ratio, 1e-6) - 1.0)
    K_ = int(min(kmax, max(1, np.ceil(need / LB))))
    if k_force is not None:
        K_ = int(k_force)
    carriers, dg = carrier_block(profile_de, libs[: 400 - K_], emit_rng, spread, key, i, emit_c)
    carriers = np.asarray(carriers)
    T_c = carriers.sum(0).astype(np.int64)
    S = int(T_c.sum()) + K_ * LB
    # ratio mode (410b): the bulk/DE design ratio applied on the carriers' DELIVERED composition
    r = np.where(profile_de > 0, profile_bulk / np.maximum(profile_de, 1e-300), 1.0)
    tgt = (T_c / T_c.sum()) * r
    T_b = alloc(S, tgt / tgt.sum())
    Q = T_b - T_c
    clipped = Q < 0
    clip_mass = int(-Q[clipped].sum())
    Q = np.where(clipped, 0, Q)
    Q, trim = cap_mass(Q, K_)
    per = Q // K_
    ball = np.tile(per, (K_, 1))
    rem = Q - per * K_
    off = np.arange(Q.size) % K_
    for k in range(K_):
        ball[k] += (rem > ((k - off) % K_)).astype(np.int64)
    assert np.array_equal(ball.sum(0), Q) and ball.max() >= 0 and ball.sum(1).max() <= 1_000_000, ball.sum(1)
    block = np.vstack([carriers.astype(np.int64), ball]).astype(np.int64)
    jk_target = jk_sum(carriers)
    block, jk_rho, jk_after = jk_boost(block, carriers.shape[0], K_, jk_target, np.random.default_rng(1000 + i))
    assert np.array_equal(block.sum(0), np.vstack([carriers, ball]).sum(0)), "jk boost moved the totals"
    assert block.min() >= 0 and block.sum(1).max() <= 1_000_000
    dg = dict(dg)
    dg.update(jk_rho=jk_rho, jk_after=jk_after, jk_target=jk_target, sink_used=0)
    dg.update(ballast_K=K_, ballast_ratio=ratio, ballast_clip_genes=int(clipped.sum()), ballast_clip_mass=clip_mass,
              ballast_trim=trim, ballast_libs=[int(x) for x in ball.sum(1)], carrier_total=int(T_c.sum()))
    return block, dg


class H1Op:
    """run410h's minimal operator on the 323 H1 axis (read-only)."""
    def __init__(self):
        with np.load(OP_H1) as z:
            self.evals = np.asarray(z["evals"], float)
            self.evecs = np.asarray(z["evecs"], float)
        with np.load(STAGE, allow_pickle=False) as z:
            self.cols = np.asarray(z["columns"], dtype=np.int64)
            self.genes = [str(g) for g in z["genes"]]
        assert self.evecs.shape[0] == len(self.cols) == 18157

    def _check_axis(self, wcols):
        assert np.array_equal(np.asarray(wcols), self.cols)


def ps_bulk(eff, genes, targets):
    """410b's shipped bulk on H1: the key410 weights, clip(0.3, 3), energy-preserving over non-panel genes."""
    import key410  # noqa: PLC0415  (410's file, read-only)
    op = H1Op()
    key, carried, nd = key410.shared_energy_key(list(map(str, targets)), DONORS, op, op.cols, cols=None)
    w, med = key410.weights(key, carried, 0.3, 3.0)
    del op
    panel = set(map(str, targets))
    nonp = np.array([str(g) not in panel for g in genes])
    e2 = np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp])
    c = float(np.sqrt(e2.sum() / (w * w * e2).sum()))
    return eff * (w * c)[:, None]


def clean4_gate(genes_full):
    """conf4: True where the four CLEAN4 banks' P_down all sit on the same side of 0.5 (covered genes only)."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float)
    P = PB[[names.index(b) for b in CLEAN4]]
    unanimous = np.abs(np.sign(P - 0.5).sum(0)) == len(CLEAN4)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.zeros(len(genes_full), dtype=bool)
    out[col >= 0] = unanimous[col[col >= 0]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("design", choices=["base", "x125", "pre", "w065nocd4", "gm030", "g000", "g030", "g090", "gain", "signrep", "prior_off", "q70", "q90", "conf4"])
    ap.add_argument("seed", type=int)
    ap.add_argument("--bulk", choices=["base", "ps"], default="base")
    ap.add_argument("--kfrom", default=None)
    ap.add_argument("--mde-c", type=float, default=0.0)
    ap.add_argument("--mde-scope", choices=["all", "block", "nonblock"], default="all")
    ap.add_argument("--mde-cap", type=float, default=0.0, help="max add in log2 (0 = uncapped)")
    ap.add_argument("--mde-min-count", type=float, default=0.0, help="genes below this control count/cell get no add")
    ap.add_argument("--mde-concordant", action="store_true", help="add only where the LFC sign agrees with sign(key)")
    ap.add_argument("--mde-overfill", type=float, default=-1.0, help="cap each class's over-fill at this x its sub-gate mass (-1 off)")
    ap.add_argument("--mde-agree", action="store_true", help="426e: add only where our sign agrees with the covered CLEAN4 prior")
    ap.add_argument("--zero-repair", type=float, default=0.0, help="428: repair zero-total gated genes to this log2 floor (0 = off)")
    ap.add_argument("--conf-raise", type=float, default=0.0,
                    help="427: raise the CARRIER effect where >= 3 non-H1 donors agree with |d| >= 0.1 to kappa * mean|d| (0 = off)")
    ap.add_argument("--mde-prior-sign", type=float, default=-1.0,
                    help="426b: non-block genes take the CLEAN4 per-gene prior sign where covered and |P_down - 0.5| > this (-1 = off)")
    ap.add_argument("--zs", default="", help="429: 423's post-pass variant after ballast (+ zero repair), e.g. dn_k250_zs")
    ap.add_argument("--zs-impl", choices=["vendored", "rp423"], default="vendored",
                    help="vendored = this directory's src_rank_mag.py (the production module); rp423 = 423's file")
    ap.add_argument("--dump", type=int, default=0, help="save the first N final pre-pass blocks to work/blocks/ (then continue)")
    a = ap.parse_args()
    global PRIOR_MODE, PROTECT_Q, GATE_FULL, MDE_ADD, MDE_SCOPE, MDE_CONC, MDE_SIGN, MDE_OVERFILL, MDE_AGREE
    tag = f"{a.design}_s{a.seed}" + ("_bps" if a.bulk == "ps" else "") + (f"_K{a.kfrom}" if a.kfrom else "")
    if a.conf_raise > 0:
        tag += f"_conf{a.conf_raise:g}"
    if a.zero_repair < 0:
        tag += f"_zr{a.zero_repair:g}"
    if os.environ.get("RUN426_EMITTER") == "prod":
        tag += "_prodem"
    if a.mde_c > 0:
        tag += f"_mde{a.mde_c:g}{'' if a.mde_scope == 'all' else a.mde_scope}"
        tag += (f"cap{a.mde_cap:g}" if a.mde_cap > 0 else "") + (f"min{a.mde_min_count:g}" if a.mde_min_count > 0 else "")
        tag += "conc" if a.mde_concordant else ""
        tag += f"ps{a.mde_prior_sign:g}" if a.mde_prior_sign >= 0 else ""
        tag += f"of{a.mde_overfill:g}" if a.mde_overfill >= 0 else ""
        tag += "agr" if a.mde_agree else ""
    if (WORK / f"metrics-{tag}.json").exists():
        raise FileExistsError(tag)
    zs_spec = zs_tag = None
    if a.zs:
        sys.path.insert(0, str(ROOT / "experiments/423-rank-preserving-magnitude"))
        import rp423  # noqa: PLC0415  (423's file, read-only: verify / r1_tie, and the reference implementation)
        if a.zs_impl == "vendored":
            import src_rank_mag as ZS  # noqa: PLC0415
            zs_spec = ZS.parse_spec(a.zs)
        else:
            ZS = rp423
            zs_spec = rp423.parse_variant(a.zs)
        assert zs_spec["mode"] == "zs", zs_spec
        zs_tag = f"{tag}_{a.zs.replace('_', '')}" + ("" if a.zs_impl == "vendored" else "rp")
        if (WORK / f"metrics-{zs_tag}.json").exists():
            raise FileExistsError(zs_tag)
        acc_zs = dict(sums=[], jks=[], ps=[], qs=[], lfcs=[], diags=[])
    b = H.builder389()
    eff, genes, targets = H.native_cap15(b)
    eff_de_native = eff
    if a.design in ("g000", "g030", "g090", "gain", "signrep"):
        with np.load(CARRIERS, allow_pickle=False) as zc:
            assert (zc["genes"] == genes).all() and (zc["targets"] == targets).all()
            eff_de_native = np.asarray(zc[a.design], dtype=np.float64)
    elif a.design in ("w065nocd4", "gm030"):
        with np.load(HERE / "work/carriers416b.npz", allow_pickle=False) as zc:
            assert (zc["genes"] == genes).all() and (zc["targets"] == targets).all()
            eff_de_native = np.asarray(zc[a.design], dtype=np.float64)
    elif a.design == "pre":
        # 410's de_unwhite construction verbatim: 323's pre-whitening cap15 WITHOUT the CD4 residual, row-norm matched
        with np.load(STAGE, allow_pickle=False) as zs:
            assert (zs["targets"].astype(str) == np.asarray(targets).astype(str)).all() and (zs["genes"] == genes).all()
            pre = np.asarray(zs["cap15"], dtype=np.float64)
        eff_de_native = pre * (np.linalg.norm(eff, axis=1) / np.maximum(np.linalg.norm(pre, axis=1), 1e-12))[:, None]
    elif a.design == "x125":
        eff_de_native = eff * 1.25          # known-answer only: 410's de_x1p25_ratio_jk (carrier change + K > 1 + jk path)
    elif a.design == "prior_off":
        PRIOR_MODE = "off"
    elif a.design in ("q70", "q90"):
        PROTECT_Q = {"q70": 0.7, "q90": 0.9}[a.design]
    eff_bulk = eff.copy() if a.bulk == "base" else ps_bulk(eff, genes, targets)
    if a.conf_raise > 0:
        # 427: carrier-profile raise on donor-confident cells; the bulk (eff_bulk) is untouched (ballast pins it)
        DON = ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"]
        tg = [str(t) for t in targets]; gn = [str(g) for g in genes]
        Dv = np.full((len(DON), len(tg), len(gn)), np.nan, dtype=np.float32)
        for di, dn in enumerate(DON):
            bd = ROOT / "data/donors" / dn
            Pn = [l.strip() for l in open(bd / "perts.txt")]; Gn = [l.strip() for l in open(bd / "genes.txt")]
            pi = {q: i for i, q in enumerate(Pn)}; gc = {q: i for i, q in enumerate(Gn)}
            Mx = np.load(bd / "delta.npy", mmap_mode="r")
            cols = np.array([gc.get(q, -1) for q in gn]); okc = cols >= 0
            for k, t in enumerate(tg):
                if t in pi:
                    Dv[di, k, okc] = np.asarray(Mx[pi[t]], dtype=np.float32)[cols[okc]]
        ncv = np.isfinite(Dv).sum(0); pos = (np.nan_to_num(Dv) > 0).sum(0); neg = (np.nan_to_num(Dv) < 0).sum(0)
        agr = np.where(pos >= neg, 1.0, -1.0)
        allag = (np.maximum(pos, neg) == ncv) & (ncv >= 3)
        mina = np.where(np.isfinite(Dv), np.abs(Dv), np.inf).min(0)
        with np.errstate(invalid="ignore"):
            mabs = np.nanmean(np.where(np.isfinite(Dv), np.abs(Dv), np.nan), 0)
        conf = allag & (mina >= 0.1)
        panel_cols = np.array([g in set(tg) for g in gn])
        conf &= ~panel_cols[None, :]
        keep_own = np.where(np.sign(eff_de_native) == agr, np.abs(eff_de_native), 0.0)
        raised = np.where(conf, agr * np.maximum(keep_own, a.conf_raise * np.nan_to_num(mabs)), eff_de_native)
        print(f"427 conf-raise kappa {a.conf_raise}: {int(conf.sum())} cells ({conf.sum() / len(tg):.1f}/target); "
              f"changed {int((raised != eff_de_native).sum())}; row-norm ratio median "
              f"{np.median(np.linalg.norm(raised, axis=1) / np.maximum(np.linalg.norm(eff_de_native, axis=1), 1e-12)):.3f}", flush=True)
        eff_de_native = raised
    k_force = None
    if a.kfrom:
        dk = json.loads((WORK / f"diag-{a.kfrom}.json").read_text())
        k_force = [int(d["ballast_K"]) for d in dk]
    r, cc, emit_c = E.load_reference()
    if a.mde_c > 0:
        MDE_SCOPE = a.mde_scope
        MDE_CONC = bool(a.mde_concordant)
        MDE_OVERFILL = a.mde_overfill if a.mde_overfill >= 0 else None
        if a.mde_prior_sign >= 0 or a.mde_agree:
            cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
            names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float)
            Pd = PB[[names.index(b) for b in CLEAN4]].mean(0)
            pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
            gk = np.asarray(r["genes"]).astype(str)[np.asarray(emit_c.keep, dtype=bool)]
            col = np.array([pos.get(g, -1) for g in gk])
            pdn = np.where(col >= 0, Pd[np.maximum(col, 0)], 0.5)
            psgn = np.where((col >= 0) & (np.abs(pdn - 0.5) > max(a.mde_prior_sign, 0.0)), np.where(pdn > 0.5, -1.0, 1.0), 0.0)
            if a.mde_prior_sign >= 0:
                MDE_SIGN = psgn
            if a.mde_agree:
                MDE_AGREE = np.where(col >= 0, np.where(pdn > 0.5, -1.0, 1.0), 0.0)
            print(f"mde prior: {int((psgn != 0).sum())}/{psgn.size} gated genes carry a prior sign", flush=True)
        cnt = np.asarray(cc.ref_lin, dtype=np.float64) * float(np.median(r["libraries"])) / 1e6
        MDE_ADD = a.mde_c / np.sqrt(np.maximum(cnt, 1e-3))
        if a.mde_cap > 0:
            MDE_ADD = np.minimum(MDE_ADD, a.mde_cap)
        if a.mde_min_count > 0:
            MDE_ADD = np.where(cnt >= a.mde_min_count, MDE_ADD, 0.0)
        assert MDE_ADD.shape == (int(np.asarray(emit_c.keep).sum()),), MDE_ADD.shape
        print(f"mde: c {a.mde_c} scope {a.mde_scope}; add median {np.median(MDE_ADD):.4f} (q10 {np.quantile(MDE_ADD, .1):.4f}, "
              f"q90 {np.quantile(MDE_ADD, .9):.4f})", flush=True)
    H.PRIOR_FULL, cov = H.prior_signs_full(r["genes"])
    if a.design == "conf4":
        GATE_FULL = clean4_gate(r["genes"])
        print(f"conf4 gate: {GATE_FULL.mean():.4f} of the full axis, {GATE_FULL[np.asarray(emit_c.keep, bool)].mean():.4f} "
              f"of the DE universe", flush=True)
    eff_de = R.map_axis(eff_de_native, genes, r["genes"])
    eff_bk = R.map_axis(eff_bulk, genes, r["genes"])
    pick_rng, emit_rng = b.make_rngs(a.seed)
    sums, jks, ps, qs, lfcs, diags = [], [], [], [], [], []
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        prof_de = b.profile_from_effect(eff_de[i], r["b_ctrl"], min_cpm=10.)
        prof_bk = b.profile_from_effect(eff_bk[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(
            eff_de[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(eff_de[i], r["b_ctrl"], min_cpm=10.)[0]
        block, dg = ballast_block(prof_de, prof_bk, r["libraries"][picks], emit_rng, spread,
                                  np.where(gate, eff_de[i], 0.), i, emit_c, 16,
                                  k_force=None if k_force is None else k_force[i])
        if K.LAST_PRIOR_REC:
            dg = {**dg, **K.LAST_PRIOR_REC.pop()}
        if a.zero_repair < 0:
            sys.path.insert(0, str(ROOT / "experiments/428-zero-gene-audit"))
            import repair428  # noqa: PLC0415
            gfull = [str(g) for g in r["genes"]]
            tcol = gfull.index(str(target)) if str(target) in gfull else -1
            block, zd = repair428.zero_repair(np.asarray(block), 400 - int(dg["ballast_K"]), np.asarray(emit_c.keep, dtype=bool),
                                              np.asarray(cc.ref_lin, dtype=np.float64), tcol, np.random.default_rng(5000 + i),
                                              floor_log2=float(a.zero_repair))
            dg.update(zd)
        dg = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in dg.items()}
        dg["target"] = str(target)
        diags.append(dg)
        n_car = 400 - int(dg["ballast_K"])
        if i < a.dump:
            import scipy.sparse as sp  # noqa: PLC0415
            Xd = sp.csr_matrix(np.asarray(block).astype(np.int32))
            (HERE / "work/blocks").mkdir(parents=True, exist_ok=True)
            np.savez_compressed(HERE / f"work/blocks/{tag}_{i:02d}.npz", data=Xd.data, indices=Xd.indices, indptr=Xd.indptr,
                                shape=np.array(Xd.shape), n_car=n_car, target=str(target))
        s, j, p, q, l = E.block_stats(block, cc)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        if zs_spec is not None:
            # 429: the builder's run order -- emitter (+ MDE) -> ballast -> zero repair -> zs pass, on the final block
            blk2, pdg = ZS.postpass(block, n_car, cc, l, zs_spec)
            if a.zs_impl == "vendored":
                # every block: the vendored production module must equal 423's reference implementation bit-for-bit
                blk_ref, _ = rp423.postpass(block, n_car, cc, l, rp423.parse_variant(a.zs))
                if not np.array_equal(blk_ref, blk2):
                    raise AssertionError(f"{target}: vendored zs pass differs from rp423")
                del blk_ref
                pdg["rp423_equal"] = True
            s2, j2, p2, q2, l2 = E.block_stats(blk2, cc)
            pdg.pop("_fr_by_gene")
            pdg.update(rp423.verify(block, blk2, n_car, cc, s, p, q, l, s2, p2, q2, l2, touched=pdg.pop("_touched")))
            pdg["ballast_equal"] = bool(np.array_equal(np.asarray(blk2)[n_car:], np.asarray(block)[n_car:]))
            pdg["target"] = str(target)
            for nm, val in (("sums", s2), ("jks", j2), ("ps", p2), ("qs", q2), ("lfcs", l2)):
                acc_zs[nm].append(val)
            acc_zs["diags"].append(pdg)
            del blk2
        if (i + 1) % 20 == 0:
            print(f"{tag} {i + 1}/85 ({time.monotonic() - t0:.1f}s) K={dg['ballast_K']} clip={dg['ballast_clip_genes']}", flush=True)
    raw, rows, legs = E.score(r, cc, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{tag}.npz", targets=r["targets"], sums=np.asarray(sums), jk=np.asarray(jks),
             p=np.asarray(ps), q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{tag}.csv")
    R.write_json(WORK / f"diag-{tag}.json", diags)
    out = {"tag": tag, "design": a.design, "seed": a.seed, "bulk": a.bulk, "kfrom": a.kfrom, "raw": raw, "mse_legs": legs,
           "seconds": time.monotonic() - t0,
           "K_hist": {int(k): int(v) for k, v in zip(*np.unique([d["ballast_K"] for d in diags], return_counts=True))},
           "clip_genes_median": float(np.median([d["ballast_clip_genes"] for d in diags])),
           "clip_mass_total": int(sum(d["ballast_clip_mass"] for d in diags)),
           "prior_lambda_median": (float(np.median([d["prior_lambda"] for d in diags])) if "prior_lambda" in diags[0] else None),
           "prior_flips_median": (float(np.median([d["prior_n_flip_target"] for d in diags])) if "prior_n_flip_target" in diags[0] else None),
           "mde_c": a.mde_c, "mde_scope": a.mde_scope, "mde_cap": a.mde_cap, "mde_min_count": a.mde_min_count,
           "mde_fu_median": (float(np.median([d["mde_fu"] for d in diags])) if "mde_fu" in diags[0] else None),
           "mde_fd_median": (float(np.median([d["mde_fd"] for d in diags])) if "mde_fd" in diags[0] else None),
           "mde_lambda_median": (float(np.median([d.get("mde_lambda", np.nan) for d in diags])) if "mde_fd" in diags[0] else None),
           "mde_overfill_max": (float(np.max([d.get("mde_overfill_frac", np.nan) for d in diags])) if "mde_fd" in diags[0] else None),
           "unsuppressible_mean": float(np.mean([d.get("call_set_unsuppressible", 0) for d in diags])),
           "clip_frac_mean": float(np.mean([d.get("clip_frac", 0) for d in diags])),
           "prior_lambda_mean": (float(np.mean([d["prior_lambda"] for d in diags])) if "prior_lambda" in diags[0] else None)}
    R.write_json(WORK / f"metrics-{tag}.json", out)
    print(tag, json.dumps(out), flush=True)
    if zs_spec is not None:
        A = acc_zs
        raw2, rows2, legs2 = E.score(r, cc, A["sums"], A["jks"], A["ps"], A["qs"], A["lfcs"])
        np.savez(WORK / f"stats-{zs_tag}.npz", targets=r["targets"], sums=np.asarray(A["sums"]), jk=np.asarray(A["jks"]),
                 p=np.asarray(A["ps"]), q=np.asarray(A["qs"]), lfc=np.asarray(A["lfcs"]))
        pl.DataFrame(rows2).write_csv(WORK / f"per-target-{zs_tag}.csv")
        R.write_json(WORK / f"diag-{zs_tag}.json", A["diags"])
        out2 = {"tag": zs_tag, "pre_tag": tag, "zs": a.zs, "zs_impl": a.zs_impl, "spec": zs_spec, "seed": a.seed,
                "bulk": a.bulk, "raw": raw2, "mse_legs": legs2, "seconds": time.monotonic() - t0,
                "summary": rp423.summarise(A["diags"])}
        R.write_json(WORK / f"metrics-{zs_tag}.json", out2)
        print(zs_tag, json.dumps({"raw": raw2, "summary": out2["summary"]}), flush=True)


if __name__ == "__main__":
    main()
