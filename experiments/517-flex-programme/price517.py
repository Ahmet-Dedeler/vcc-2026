#!/usr/bin/env python3
"""517 — analytic H1 prices of target-independent sign sources in the DE channel, through the production emitter's exact
reconstruction (490's `recon_h1.py`), with calls / p / q / reach order held at the base run's delivered values.

Two kinds of substitution, each on one of three BASE runs (the delivered H1 stats + diag of a real emitter run):
  base  c4   435's shipped construction (CLEAN4 q80 block prior, coarse MDE)               = 446b / 446bj2 on the board
        g33  467's gated33 (Arc's programme on the design-strong 2/3, CLEAN4 on the weak 1/3) = 458g
        arc  462l's arcfull (Arc's programme on every target)                                = 458d
  (a) HOOK: the MDE add on complement genes (not head, not block: the uncalled genes) takes a source's sign where
      |P - 0.5| >= tau (496's `mde_add_sign`, 478's analytic route). fu/fd re-solved exactly as the emitter does.
  (b) BLOCK: the block prior (the sign the unprotected CPM-block genes are steered to) is replaced per target, optionally
      with another protect_q. The intent is reconstructed at the base's recorded doses (prior_lambda, mde_lambda) and the
      change of intent is added to the delivered means (478's mapping). Real runs re-roll p/q; this holds them.

Sources (P-like per gene, NaN = uncovered): CLEAN4 (the mean of the four banks' P_down, 478's prior_P), Arc = P_down over
data/donors/h1_2025_full with h1-test's 85 targets' rows removed (467/462l's rule), their z-blends, and anything passed in
by a caller (the ranking script adds more).

Known answers, every base and seed (the script aborts on any failure):
  1. recon's own: kf, mde_fu/fd, n_up/dn, flips, n_block == the base run's diag (per target);
  2. the delivered tables reproduce the base run's metrics json (<= 1e-9);
  3. the MDE with the shipped signs reproduces recon's A_int (<= 1e-9), so a no-change variant is Δ 0 exactly;
  4. (printed) hook c4 tau 0.05 on base c4 = 478 / 496's nmae (-0.00147 on seed 0 in 496's real run).

usage: uv run python experiments/517-flex-programme/price517.py --base c4 --seeds 0,1,2 --set hook   (< 2 GB, ~1-2 min/seed)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/490-de-emitter-review"))
import recon_h1 as RC  # noqa: E402  (read-only)
import score_de  # noqa: E402  (read-only)
from recon_h1 import K  # noqa: E402  (production rc_emit_k, read-only)

TAG = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
BASES = {"c4": ROOT / "experiments/435-zs-on-coarse/work/h1",
         "g33": ROOT / "experiments/467-strength-gated-programme/work/h1_gated33",
         "arc": ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull"}
SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}     # 467's H1-like spans (pds, mse pinned: Δ 0)
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
FOUR = ("CNN3", "COL6A1", "ZNF664", "AASS")


# ------------------------------------------------------------------------------------------------ sources (P-like)
def clean4_P(genes_full):
    """478's prior_P: the mean of the four CLEAN4 banks' P_down (371's consensus file); NaN where uncovered."""
    cp = np.load(RC.PRIOR, allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in RC.PRIOR_BANKS]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.full(len(genes_full), np.nan)
    out[col >= 0] = P[col[col >= 0]]
    return out


def bank_P(bank_dir, genes_full, exclude=(), rows_keep=None):
    """P_down over a bank's rows (non-zero entries), `exclude`d perturbations removed; NaN where the bank has no entry."""
    d = Path(bank_dir)
    delta = np.load(d / "delta.npy", mmap_mode="r")
    perts = [l.strip() for l in open(d / "perts.txt")]
    genes = [l.strip() for l in open(d / "genes.txt")]
    ex = set(map(str, exclude))
    rows = [i for i, p in enumerate(perts) if p not in ex and (rows_keep is None or p in rows_keep)]
    dn = np.zeros(len(genes)); nz = np.zeros(len(genes))
    for s in range(0, len(rows), 512):
        X = np.asarray(delta[rows[s:s + 512]], dtype=np.float64)
        dn += (X < 0).sum(0); nz += (X != 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: j for j, g in enumerate(genes)}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    return np.where(col >= 0, pd[np.maximum(col, 0)], np.nan), len(rows)


def sign_from_P(P, fill=-1.0):
    """The builder's block-sign rule: P_down > 0.5 -> -1, <= 0.5 -> +1; uncovered -> fill (run400: -1)."""
    return np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), fill)


def zblend(Pa, Pb, w):
    """sign source z = w z(Pa - .5) + (1-w) z(Pb - .5), sd over each one's covered genes; returned as a pseudo-P with the
    same |P - .5| scale as Pa (so tau thresholds stay comparable): Pa where only Pa covers, Pb's where only Pb covers."""
    za = np.where(np.isfinite(Pa), (0.5 - Pa) / np.nanstd(Pa), np.nan)
    zb = np.where(np.isfinite(Pb), (0.5 - Pb) / np.nanstd(Pb), np.nan)
    z = np.where(np.isfinite(za) & np.isfinite(zb), w * za + (1 - w) * zb, np.where(np.isfinite(za), za, zb))
    return np.where(np.isfinite(z), 0.5 - z * np.nanstd(Pa), np.nan)


# ------------------------------------------------------------------------------------------------ the MDE (478 verbatim)
def mde_intent(sgn, madd, Ms, A0, E0v, base_t, blk, Lk, nk, n_car, lam, lam_m, move):
    """478's mde_intent without payers: recon_h1's MDE path with an arbitrary sign vector (fu/fd re-solved)."""
    elig = (Ms > 0) & np.isfinite(E0v) & np.isfinite(madd) & (madd > 0) & (sgn != 0)
    up, dn = elig & (sgn > 0), elig & (sgn < 0)
    Bm = A0 * np.power(2.0, np.where(np.isfinite(base_t - E0v), base_t - E0v, 0.0))
    U = lambda f: float((Bm[up] * (np.power(2.0, f * madd[up]) - 1.0)).sum())       # noqa: E731
    D = lambda f: float((Bm[dn] * (1.0 - np.power(2.0, -f * madd[dn]))).sum())      # noqa: E731
    U1, D1 = U(1.0), D(1.0)
    fu = fd = 1.0
    if U1 > D1 and U1 > 0:
        a_, b_ = 0.0, 1.0
        for _ in range(60):
            m_ = 0.5 * (a_ + b_)
            if U(m_) > D1:
                b_ = m_
            else:
                a_ = m_
        fu = 0.5 * (a_ + b_)
    elif D1 > U1 and D1 > 0:
        a_, b_ = 0.0, 1.0
        for _ in range(60):
            m_ = 0.5 * (a_ + b_)
            if D(m_) > U1:
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
    A_int = A_ms.copy()
    if lam > 0:
        A_int[blk] = A0[blk] * np.power(2.0, lam * move[blk] + lam_m * move_m[blk])
    stn = steer_set & ~blk
    if lam_m > 0 and stn.any():
        _, A_n = K.split_coarse_extremes(base_[stn], A_ms[stn] * np.power(2.0, lam_m * move_m[stn]), Lk, n_car)
        A_int[stn] = A_n
    return A_int, dict(fu=fu, fd=fd, n_up=int(up.sum()), n_dn=int(dn.sum()))


def block_target(E0v, blk, prior_g, kf_fixed=None, protect_q=0.8):
    """recon's block path: (tgt after q-protection, kf-solved t2, move, base_t). kf is re-solved unless fixed."""
    thr = float(np.quantile(np.abs(E0v[blk]), protect_q)) if blk.any() else 0.0
    tgt = np.where(blk & (np.abs(E0v) <= thr), prior_g * np.abs(E0v), E0v)
    side = blk & (tgt < 0)
    return tgt, side, thr


class Ctx:
    """Everything the per-target reconstruction needs for one (base, seed), loaded once."""

    def __init__(self, base, seed, prior_for_target):
        self.base, self.seed = base, seed
        d = BASES[base]
        tag = TAG.format(s=seed)
        self.diag = json.loads((d / f"diag-{tag}.json").read_text())
        st = np.load(d / f"stats-{tag}.npz", allow_pickle=False)
        self.st = {k: st[k] for k in ("targets", "p", "q", "lfc")}
        zr = np.load(RC.REF, allow_pickle=False)
        self.r = {k: zr[k] for k in ("genes", "targets", "keep", "ref_lin", "emit_ref_lin", "libraries", "b_ctrl",
                                     "lfc", "q")}
        zr.close()
        assert np.array_equal(self.st["targets"], self.r["targets"])
        self.keep = np.asarray(self.r["keep"], dtype=bool)
        self.ref_emit = np.asarray(self.r["emit_ref_lin"], dtype=np.float64)
        self.ref_score = np.asarray(self.r["ref_lin"], dtype=np.float64)
        with RC.R.model_env("uniform"):
            self.b = RC.R.load_module("builder517", ROOT / "scripts/make_submission_pm_cd4.py")
        with np.load(RC.EFF, allow_pickle=False) as z:
            self.eff = RC.R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], self.r["genes"])
        cnt = self.ref_score * float(np.median(self.r["libraries"])) / 1e6
        self.madd = np.where(cnt >= RC.MDE_MIN, np.minimum(RC.MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), RC.MDE_CAP), 0.0)
        self.T = self.r["targets"].astype(str)
        self.G = np.asarray(self.r["genes"])[self.keep].astype(str)
        self.lfc_del = self.st["lfc"]
        self.m_del = (self.ref_score + 1e-9) * np.power(2.0, self.lfc_del) - 1e-9
        self.raw_o, self.per_o = score_de.score(self.T, self.G, self.st["p"], self.st["q"], self.lfc_del)
        want = json.loads((d / f"metrics-{tag}.json").read_text())["raw"]
        dev = max(abs(self.raw_o[m] - want[m]) for m in self.raw_o)
        assert dev <= 1e-9, f"KA2 harness tables {base} s{seed}: {dev}"
        # truth, on the keep axis: gate = q < .05 & finite nonzero lfc, the target's own gene out
        self.lfc_real = np.asarray(self.r["lfc"], dtype=np.float64)
        self.q_real = np.asarray(self.r["q"], dtype=np.float64)
        pick_rng, _ = self.b.make_rngs(seed)
        self.blocks = []
        worst = {}
        for i, t in enumerate(self.T):
            picks = pick_rng.choice(len(self.r["libraries"]), 400, replace=False)
            prior_g = prior_for_target(i, t)
            o = RC.reconstruct(i, self.diag[i], self.eff[i], self.r, self.keep, self.ref_emit, self.ref_score, prior_g,
                               self.madd, self.b, picks, fix=False)
            for name, (got, want_) in o["ka"].items():
                if want_ is None:
                    continue
                dd = abs(float(got) - float(want_))
                worst[name] = max(worst.get(name, 0.0), dd)
                if dd > (1e-9 if name.startswith(("kf", "mde_f")) else 0):
                    raise AssertionError(f"KA1 {base} s{seed} {t}: {name} {got} != {want_}")
            n_car = o["n_car"]
            libs = np.asarray(self.r["libraries"], dtype=np.float64)[picks][:n_car]
            profile = self.b.profile_from_effect(self.eff[i], self.r["b_ctrl"], min_cpm=10.)
            L = int(min(np.median(libs), 1_000_000))
            totals = np.floor(n_car * L * profile).astype(np.int64)
            grand = int(totals.sum())
            Ms = totals[self.keep].astype(np.int64)
            _, Lk, _, nk = K.class_cells(n_car, RC.CLASSES3, grand)
            A0 = 1e6 * Ms.astype(np.float64) / float(grand)
            E0v, kf, blk = o["E0v"], o["kf"], o["blk"]
            tgt, side, _ = block_target(E0v, blk, prior_g)
            t2 = np.where(side, np.minimum(tgt * kf, -RC.FLOOR), tgt)
            move = np.where(blk, t2 - E0v, 0.0)
            base_t = np.where(blk, t2, E0v)
            sgn0 = np.sign(base_t)
            A_chk, _ = mde_intent(sgn0, self.madd, Ms, A0, E0v, base_t, blk, Lk, nk, n_car, o["lam"], o["lam_m"], move)
            assert np.allclose(A_chk, o["A_int"], rtol=0, atol=1e-9), f"KA3 {base} s{seed} {t}"
            self.blocks.append(dict(i=i, t=t, o=o, Ms=Ms, A0=A0, E0v=E0v, base_t=base_t, blk=blk, Lk=Lk, nk=nk,
                                    n_car=n_car, move=move, sgn0=sgn0, comp=~o["head"] & ~blk, prior=prior_g,
                                    picks=picks))
        self.ka_worst = worst

    # ---------------------------------------------------------------------------------------- substitutions
    def hook(self, sign_for_target):
        """(a): per target a complement sign vector (+-1, 0 = keep the delivered sign) -> new lfc table."""
        lfc_new = self.lfc_del.copy()
        rec = dict(flips=[], cover=[], fd=[], fu=[])
        for B in self.blocks:
            asg = sign_for_target(B["i"], B["t"])
            use = B["comp"] & (asg != 0)
            sgn = np.where(use, asg, B["sgn0"])
            A_new, inf = mde_intent(sgn, self.madd, B["Ms"], B["A0"], B["E0v"], B["base_t"], B["blk"], B["Lk"], B["nk"],
                                    B["n_car"], B["o"]["lam"], B["o"]["lam_m"], B["move"])
            m_new = self.m_del[B["i"]] + (B["n_car"] / 400.0) * (A_new - B["o"]["A_int"])
            lfc_new[B["i"]] = np.log2((np.maximum(m_new, 0.0) + 1e-9) / (self.ref_score + 1e-9))
            rec["flips"].append(int((use & (asg != B["sgn0"])).sum())); rec["cover"].append(int(use.sum()))
            rec["fd"].append(inf["fd"]); rec["fu"].append(inf["fu"])
        return np.where(np.isfinite(self.lfc_del), lfc_new, self.lfc_del), rec

    def block(self, prior_for_target, protect_q=0.8, hook_for_target=None):
        """(b): per target a new block prior (and protect_q); kf and fu/fd re-solved; optionally with a hook on top."""
        lfc_new = self.lfc_del.copy()
        rec = dict(flips_vs_base=[], blk_flipped=[], kf=[])
        for B in self.blocks:
            o, E0v, blk = B["o"], B["E0v"], B["blk"]
            pg = prior_for_target(B["i"], B["t"])
            tgt, side, _ = block_target(E0v, blk, pg, protect_q=protect_q)
            A0 = B["A0"]

            def _t2(kf):
                return np.where(side, np.minimum(tgt * kf, -RC.FLOOR), tgt)

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
            t2 = _t2(kf)
            move = np.where(blk, t2 - E0v, 0.0)
            base_t = np.where(blk, t2, E0v)
            sgn = np.sign(base_t)
            if hook_for_target is not None:
                asg = hook_for_target(B["i"], B["t"])
                use = B["comp"] & (asg != 0)
                sgn = np.where(use, asg, sgn)
            A_new, _ = mde_intent(sgn, self.madd, B["Ms"], A0, E0v, base_t, blk, B["Lk"], B["nk"], B["n_car"],
                                  o["lam"], o["lam_m"], move)
            m_new = self.m_del[B["i"]] + (B["n_car"] / 400.0) * (A_new - o["A_int"])
            lfc_new[B["i"]] = np.log2((np.maximum(m_new, 0.0) + 1e-9) / (self.ref_score + 1e-9))
            rec["flips_vs_base"].append(int((blk & (np.sign(base_t) != np.sign(B["base_t"]))).sum()))
            rec["blk_flipped"].append(int((blk & (np.sign(t2) != np.sign(E0v))).sum()))
            rec["kf"].append(kf)
        return np.where(np.isfinite(self.lfc_del), lfc_new, self.lfc_del), rec

    # ---------------------------------------------------------------------------------------- scoring
    def score(self, lfc_new, rng, n_boot=2000):
        raw_n, per_n = score_de.score(self.T, self.G, self.st["p"], self.st["q"], lfc_new)
        D = {m: np.array([per_n[m].get(t, np.nan) - self.per_o[m].get(t, np.nan) for t in self.T]) for m in SP}
        ov = lambda ii: sum(np.nanmean(D[m][ii]) / SP[m] for m in SP) / 6         # noqa: E731
        allix = np.arange(len(self.T))
        ne = np.array([t not in ESS for t in self.T])
        boots = np.array([ov(rng.integers(0, len(self.T), len(self.T))) for _ in range(n_boot)])
        nix = np.flatnonzero(ne)
        boots_ne = np.array([ov(nix[rng.integers(0, nix.size, nix.size)]) for _ in range(n_boot)])
        out = {m: float(np.nanmean(D[m])) for m in SP}
        out.update(overall_equiv=float(ov(allix)), ci95=[float(np.quantile(boots, .025)), float(np.quantile(boots, .975))],
                   p_le_0=float((boots <= 0).mean()),
                   ne_overall_equiv=float(ov(nix)), ne_ci95=[float(np.quantile(boots_ne, .025)), float(np.quantile(boots_ne, .975))],
                   ne_nmae=float(np.nanmean(D["nmae"][nix])), ne_fid=float(np.nanmean(D["fid"][nix])),
                   ne_reach=float(np.nanmean(D["reach"][nix])),
                   without_four=float(ov(np.flatnonzero(~np.isin(self.T, FOUR)))))
        return out, D

    def gate_accuracy(self, sign_for_target, where="comp"):
        """Sign accuracy vs the truth on gate genes of the chosen subset (complement / block) where the source gives a
        sign; also the delivered base sign's accuracy on the same genes. Pooled over (target, gene)."""
        hit = tot = hit0 = 0
        for B in self.blocks:
            i = B["i"]
            asg = sign_for_target(i, B["t"])
            lr, qr = self.lfc_real[i], self.q_real[i]
            gate = (qr < 0.05) & np.isfinite(lr) & (lr != 0) & (self.G != B["t"])
            sub = B["comp"] if where == "comp" else B["blk"]
            m = gate & sub & (asg != 0)
            tot += int(m.sum())
            hit += int((np.sign(asg[m]) == np.sign(lr[m])).sum())
            d0 = np.sign(self.lfc_del[i][m])
            hit0 += int((d0 == np.sign(lr[m])).sum())
        return dict(n=tot, acc=hit / max(tot, 1), acc_delivered=hit0 / max(tot, 1))


def load_sources(genes_full, keep):
    """P-like sources on the keep axis."""
    h1t = set(np.load(RC.REF, allow_pickle=False)["targets"].astype(str))
    Pc4 = clean4_P(genes_full)[keep]
    Parc, n_arc = bank_P(ROOT / "data/donors/h1_2025_full", genes_full, exclude=h1t)
    return {"c4": Pc4, "arc": Parc[keep]}, {"arc_rows": n_arc}


def strong_set(T):
    """467's q33 gate on H1: strength n_gt_0.2 >= its 33rd percentile (= 1.72 -> strength >= 2); 57 strong, 28 weak."""
    st = json.loads((ROOT / "experiments/467-strength-gated-programme/work/strength467.json").read_text())
    meas = st["chosen_measure"]
    order = st["targets_in_order"]
    vals = np.array([st["per_target"][str(i)][meas] for i in range(len(st["per_target"]))], dtype=float)
    thr = float(np.quantile(vals, 0.33))
    strong_names = {order[i] for i in np.flatnonzero(vals >= thr)}
    return np.array([t in strong_names for t in T]), thr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", choices=list(BASES), default="c4")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--set", default="hook", choices=["ka", "hook", "block", "q"])
    a = ap.parse_args()
    t0 = time.monotonic()
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    zr.close()
    SRC, info = load_sources(genes_full, keep)
    strong, thr = strong_set(T)
    c4_block = RC.prior_signs_full(genes_full)[keep]
    arc_block = sign_from_P(SRC["arc"])
    base_prior = {"c4": lambda i, t: c4_block,
                  "g33": lambda i, t: arc_block if strong[i] else c4_block,
                  "arc": lambda i, t: arc_block}[a.base]
    print(f"[517] base {a.base}; Arc rows {info['arc_rows']}; gate thr {thr:.3f}: strong {int(strong.sum())}/{len(T)}; "
          f"Arc covers {np.isfinite(SRC['arc']).mean():.3f} / CLEAN4 {np.isfinite(SRC['c4']).mean():.3f} of the gated axis",
          flush=True)

    def hk(P, tau):
        s = np.where(np.isfinite(P) & (np.abs(P - 0.5) >= tau), np.where(P > 0.5, -1.0, 1.0), 0.0)
        return lambda i, t: s

    def hk_gate(Pstrong, tau_s, Pweak, tau_w):
        ss = hk(Pstrong, tau_s)(0, "")
        sw = hk(Pweak, tau_w)(0, "") if Pweak is not None else np.zeros_like(ss)
        return lambda i, t: ss if strong[i] else sw

    variants = {}
    if a.set in ("ka", "hook"):
        variants["none"] = ("hook", lambda i, t: np.zeros(int(keep.sum())))
        variants["c4_t0.05"] = ("hook", hk(SRC["c4"], 0.05))
    if a.set == "hook":
        for tau in (0.0, 0.1):
            variants[f"c4_t{tau:g}"] = ("hook", hk(SRC["c4"], tau))
        for tau in (0.0, 0.02, 0.05, 0.1, 0.15):
            variants[f"arc_t{tau:g}"] = ("hook", hk(SRC["arc"], tau))
        for w in (0.5,):
            Pb = zblend(SRC["arc"], SRC["c4"], w)
            for tau in (0.0, 0.05, 0.1):
                variants[f"blend{w:g}_t{tau:g}"] = ("hook", hk(Pb, tau))
        for tau in (0.0, 0.05, 0.1):
            variants[f"gate_arc{tau:g}_c4.05"] = ("hook", hk_gate(SRC["arc"], tau, SRC["c4"], 0.05))
            variants[f"gate_arc{tau:g}_none"] = ("hook", hk_gate(SRC["arc"], tau, None, 0.0))
        # agreement: the sign only where both sources are confident and agree
        for tau in (0.0, 0.05):
            sa, sc = hk(SRC["arc"], tau)(0, ""), hk(SRC["c4"], tau)(0, "")
            ag = np.where((sa != 0) & (sa == sc), sa, 0.0)
            variants[f"agree_t{tau:g}"] = ("hook", lambda i, t, ag=ag: ag)
    seeds = [int(s) for s in a.seeds.split(",")]
    out = {"base": a.base, "info": info, "gate_thr": thr, "n_strong": int(strong.sum()), "seeds": {}}
    for seed in seeds:
        ctx = Ctx(a.base, seed, base_prior)
        print(f"[517] {a.base} seed {seed}: KA1-3 PASS (worst {json.dumps({k: float(v) for k, v in ctx.ka_worst.items()})}); "
              f"base raw {json.dumps({k: round(v, 6) for k, v in ctx.raw_o.items()})} ({time.monotonic() - t0:.0f}s)",
              flush=True)
        rng = np.random.default_rng(517 + seed)
        res = {}
        for name, (kind, fn) in variants.items():
            if kind == "hook":
                lfc_new, rec = ctx.hook(fn)
                acc = ctx.gate_accuracy(fn, "comp")
            sc, _ = ctx.score(lfc_new, rng)
            sc.update(flips_median=float(np.median(rec["flips"])), cover_median=float(np.median(rec["cover"])),
                      fd_median=float(np.median(rec["fd"])), fu_median=float(np.median(rec["fu"])), comp_gate_acc=acc)
            res[name] = sc
            if name == "none":
                assert all(abs(sc[m]) <= 1e-12 for m in SP), f"KA3b none variant moved: {sc}"
            print(f"  s{seed} {name:22s} nmae {sc['nmae']:+.5f} fid {sc['fid']:+.5f} reach {sc['reach']:+.5f} jac {sc['jac']:+.5f}"
                  f" | ov {sc['overall_equiv']:+.5f} [{sc['ci95'][0]:+.4f},{sc['ci95'][1]:+.4f}] ne {sc['ne_overall_equiv']:+.5f}"
                  f" | flips {sc['flips_median']:.0f} fd {sc['fd_median']:.3f} | uncalled-gate acc {acc['acc']:.4f} "
                  f"(delivered {acc['acc_delivered']:.4f}, n {acc['n']})", flush=True)
        out["seeds"][str(seed)] = res
        del ctx
    out["seconds"] = time.monotonic() - t0
    (WORK / f"price517_{a.set}_{a.base}.json").write_text(json.dumps(out, indent=1))
    print(f"[517] done {out['seconds']:.0f}s -> work/price517_{a.set}_{a.base}.json", flush=True)


if __name__ == "__main__":
    main()
