#!/usr/bin/env python3
"""490 — exact reconstruction of the production K-class emitter's per-gene INTENT on 435's H1 runs, compared with the
DELIVERED DE table of the full 400-cell block (carriers + ballast) that the scorer reads.

Why exact is possible without re-emitting: `emit_counts(mode="quasi")` fixes each gene's pass-1 total at
floor(n * L * p_g) (no rng), `pick_rng` is its own stream (default_rng(seed)), and the ballast K per target is in the
recorded diag. So A0, E0, the head/CPM-block membership, the q80 threshold, kf, the MDE factors and the coarse split are
deterministic functions of recorded inputs. Known answer per target: kf, mde_fu/fd, mde_n_up/dn and prior_n_flip_target
must equal the run's own diag, else the script aborts.

Inputs (read-only): 323's reference.npz (truth DE; the 1.6 GB control arrays are never loaded) and effects.npz (cap15
design), 371's CLEAN4 consensus prior, 435's diag/stats. Output: work/recon_<tag>.npz (per-target per-gene arrays) and a
printed summary. Peak RSS well under 2 GB.

usage: uv run python recon_h1.py [TAG] (default base_s0_bps_prodem_mde0.1cap0.1min0.3coarse)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(WORK / "numba-cache"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))

from vcc2026 import rc_emit_k as K  # noqa: E402  (production emitter, read-only)

import run323 as R  # noqa: E402  (read-only: model_env, load_module, map_axis)

H1 = ROOT / "experiments/435-zs-on-coarse/work/h1"
REF = ROOT / "experiments/323-current-family-validation/work/reference.npz"
EFF = ROOT / "experiments/323-current-family-validation/work/effects.npz"
PRIOR = ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz"
PRIOR_BANKS = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
CLASSES3 = ((1 / 3, 16078.75), (1 / 3, 20190.5), (1 / 3, 25330.5))     # run400.CLASSES3 (the H1 harness's classes)
N_HEAD, N_CPM = 400, 2200
MDE_C, MDE_CAP, MDE_MIN = 0.1, 0.1, 0.3
PROTECT_Q, FLOOR = 0.8, 0.002


def prior_signs_full(genes_full):
    """run400.prior_signs_full verbatim: CLEAN4 sign, uncovered -> -1."""
    cp = np.load(PRIOR, allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    PB = cp["P_banks"].astype(float)
    P = PB[[names.index(b) for b in PRIOR_BANKS]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.full(len(genes_full), -1.0)
    covered = col >= 0
    out[covered] = np.where(P[col[covered]] > 0.5, -1.0, 1.0)
    return out


def reconstruct(i, dg, eff_row, r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=False):
    """One block's intent. Mirrors rc_emit_k.rc_hybrid_block_k's prior path line for line (lines 936-1099)."""
    n_car = 400 - int(dg["ballast_K"])
    libs = np.asarray(r["libraries"], dtype=np.float64)[picks][:n_car]
    profile = b.profile_from_effect(eff_row, r["b_ctrl"], min_cpm=10.)
    L = int(min(np.median(libs), 1_000_000))
    totals = np.floor(n_car * L * profile).astype(np.int64)            # emit_counts quasi: exact per-gene totals
    grand = int(totals.sum())
    gate = b.emission_rank_axis(eff_row, r["b_ctrl"], min_cpm=10.)[0]
    key = np.where(gate, eff_row, 0.0)
    key_gated = key[keep]
    n_keep = int(keep.sum())
    Ms = totals[keep].astype(np.int64)
    cell_class, Lk, T_cells, nk = K.class_cells(n_car, CLASSES3, grand)
    A0 = 1e6 * Ms.astype(np.float64) / float(grand)
    with np.errstate(divide="ignore", invalid="ignore"):
        E0v = np.log2(np.maximum(A0, 1e-300) / np.maximum(ref_emit, 1e-300))
    kg = np.abs(key_gated)
    head = np.zeros(n_keep, dtype=bool)
    head[np.argsort(-kg, kind="stable")[:N_HEAD]] = True
    blk_all = np.zeros(n_keep, dtype=bool)
    blk_all[np.argsort(-np.where(head, -np.inf, ref_emit), kind="stable")[:N_CPM]] = True
    blk = blk_all & (Ms > 0) & np.isfinite(E0v)
    tgt = prior_g * np.abs(E0v)
    thr = float(np.quantile(np.abs(E0v[blk]), PROTECT_Q)) if blk.any() else 0.0
    protected = blk & (np.abs(E0v) > thr)
    tgt = np.where(blk & (np.abs(E0v) <= thr), tgt, E0v)
    side = blk & (tgt < 0)
    if fix:                     # 484's patch (475): protected genes out of the kf shrink
        side = side & ~protected

    def _t2(kf):
        return np.where(side, np.minimum(tgt * kf, -FLOOR), tgt)

    def _R(kf):
        t2 = _t2(kf)
        return float((A0[blk] * (np.power(2.0, t2[blk] - E0v[blk]) - 1.0)).sum())
    r_lo, r_hi = _R(0.0), _R(1.0)
    if r_lo * r_hi < 0:
        a_, b_ = 0.0, 1.0
        for _ in range(50):
            m_ = 0.5 * (a_ + b_)
            if _R(a_) * _R(m_) <= 0:
                b_ = m_
            else:
                a_ = m_
        kf = 0.5 * (a_ + b_)
    else:
        kf = 0.0 if abs(r_lo) < abs(r_hi) else 1.0
    move = np.where(blk, _t2(kf) - E0v, 0.0)
    base_t = np.where(blk, _t2(kf), E0v)
    sgn = np.sign(base_t)
    elig = (Ms > 0) & np.isfinite(E0v) & np.isfinite(madd) & (madd > 0) & (sgn != 0)
    up = elig & (sgn > 0)
    dn = elig & (sgn < 0)
    Bm = A0 * np.power(2.0, np.where(np.isfinite(base_t - E0v), base_t - E0v, 0.0))

    def _U(f):
        return float((Bm[up] * (np.power(2.0, f * madd[up]) - 1.0)).sum())

    def _D(f):
        return float((Bm[dn] * (1.0 - np.power(2.0, -f * madd[dn]))).sum())
    U1, D1 = _U(1.0), _D(1.0)
    fu = fd = 1.0
    if U1 > D1 and U1 > 0:
        a_, b_ = 0.0, 1.0
        for _ in range(60):
            m_ = 0.5 * (a_ + b_)
            if _U(m_) > D1:
                b_ = m_
            else:
                a_ = m_
        fu = 0.5 * (a_ + b_)
    elif D1 > U1 and D1 > 0:
        a_, b_ = 0.0, 1.0
        for _ in range(60):
            m_ = 0.5 * (a_ + b_)
            if _D(m_) > U1:
                b_ = m_
            else:
                a_ = m_
        fd = 0.5 * (a_ + b_)
    move_m = np.where(up, fu * madd, 0.0) - np.where(dn, fd * madd, 0.0)
    steer_set = blk | elig
    w_ = (nk * Lk).astype(np.float64)
    w_ /= w_.sum()
    base_ = K.split_massprop(Ms, w_)
    A_ms = K.arith_of_split(base_, Lk, n_car)
    lam = float(dg["prior_lambda"])
    lam_m = float(dg.get("mde_lambda", 0.0))
    # intended arithmetic mean per gene over the CARRIERS (the emitter's own frame), at the recorded doses
    A_int = A_ms.copy()
    if lam > 0:
        A_int[blk] = A0[blk] * np.power(2.0, lam * move[blk] + lam_m * move_m[blk])
    stn = steer_set & ~blk
    if lam_m > 0 and stn.any():
        _, A_n = K.split_coarse_extremes(base_[stn], A_ms[stn] * np.power(2.0, lam_m * move_m[stn]), Lk, n_car)
        A_int[stn] = A_n
    with np.errstate(divide="ignore", invalid="ignore"):
        lfc_int = np.log2((A_int + 1e-9) / (ref_emit + 1e-9))
        lfc_score_int = np.log2((A_int + 1e-9) / (ref_score + 1e-9))
    flip_n = int((blk & (np.sign(_t2(kf)) != np.sign(E0v))).sum())
    ka = dict(kf=(kf, dg["prior_shrink_k"]), mde_fu=(fu, dg.get("mde_fu")), mde_fd=(fd, dg.get("mde_fd")),
              mde_n_up=(int(up.sum()), dg.get("mde_n_up")), mde_n_dn=(int(dn.sum()), dg.get("mde_n_dn")),
              flip=(flip_n, dg["prior_n_flip_target"]), n_block=(int(blk.sum()), dg["prior_n_block"]))
    return dict(E0v=E0v, key=key_gated, head=head, blk=blk, protected=protected, prior=prior_g,
                flipped=blk & (np.sign(_t2(kf)) != np.sign(E0v)), lfc_int=lfc_int, lfc_score_int=lfc_score_int,
                A_int=A_int, Ms=Ms, madd=madd, elig=elig, move_m=move_m, kf=kf, lam=lam, lam_m=lam_m,
                n_car=n_car, ka=ka)


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse"
    seed = int(tag.split("_s")[1].split("_")[0])
    t0 = time.monotonic()
    diag = json.loads((H1 / f"diag-{tag}.json").read_text())
    st = np.load(H1 / f"stats-{tag}.npz", allow_pickle=False)
    stats = {k: st[k] for k in ("targets", "p", "q", "lfc", "sums")}
    zr = np.load(REF, allow_pickle=False)
    r = {k: zr[k] for k in ("genes", "targets", "p", "q", "lfc", "keep", "ref_lin", "emit_ref_lin", "libraries",
                            "b_ctrl")}
    zr.close()
    assert np.array_equal(stats["targets"], r["targets"])
    keep = np.asarray(r["keep"], dtype=bool)
    ref_emit = np.asarray(r["emit_ref_lin"], dtype=np.float64)
    ref_score = np.asarray(r["ref_lin"], dtype=np.float64)
    with R.model_env("uniform"):
        b = R.load_module("builder490", ROOT / "scripts/make_submission_pm_cd4.py")
    assert b.MASS_FIX is None
    with np.load(EFF, allow_pickle=False) as z:
        eff = R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], r["genes"])
        assert np.array_equal(z["targets"], r["targets"])
    prior_g = prior_signs_full(r["genes"])[keep]
    cnt = ref_score * float(np.median(r["libraries"])) / 1e6
    madd = np.minimum(MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), MDE_CAP)
    madd = np.where(cnt >= MDE_MIN, madd, 0.0)
    pick_rng, _ = b.make_rngs(seed)
    genes_g = np.asarray(r["genes"])[keep]
    out = {k: [] for k in ("E0v", "key", "head", "blk", "protected", "flipped", "lfc_int", "lfc_score_int", "Ms",
                           "elig")}
    kf_all, worst = [], {}
    for i, t in enumerate(r["targets"]):
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        rec = reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks)
        assert diag[i]["target"] == str(t)
        for name, (got, want) in rec["ka"].items():
            if want is None:
                continue
            d = abs(float(got) - float(want))
            worst[name] = max(worst.get(name, 0.0), d)
            if d > (1e-9 if name in ("kf", "mde_fu", "mde_fd") else 0):
                raise AssertionError(f"{t}: known answer {name} got {got} want {want}")
        for k in out:
            out[k].append(rec[k])
        kf_all.append(rec["kf"])
    np.savez_compressed(WORK / f"recon_{tag}.npz", targets=r["targets"], genes_gated=genes_g,
                        **{k: np.asarray(v) for k, v in out.items()},
                        lfc_del=stats["lfc"], q_del=stats["q"], p_del=stats["p"],
                        lfc_real=r["lfc"], q_real=r["q"], p_real=r["p"], prior=prior_g, madd=madd,
                        ref_emit=ref_emit, ref_score=ref_score)
    print(json.dumps({"tag": tag, "known_answer_max_abs_diff": worst, "kf_median": float(np.median(kf_all)),
                      "seconds": time.monotonic() - t0}, indent=1))


if __name__ == "__main__":
    main()
