#!/usr/bin/env python3
"""478 (third session, at the coordinator's request) — can the coarse MDE carry BETTER SIGNS on uncalled genes?

The shipped coarse MDE adds magnitude madd (0.1/sqrt(count), cap 0.1) in the sign of each eligible gene's base target;
on the complement (not head, not block: uncalled genes) that sign is the design's (E0). 480 puts 65 % of the nmae sign
pot on uncalled genes; 441: CLEAN4's sign there is 0.568 vs the design's 0.536 (H1). 426 found prior signs over-fill the
classes. Here, analytically and exactly on the lfc channel (490's reconstruction, calls / p / order held):
  sign sources on complement genes covered by CLEAN4: design (shipped), CLEAN4 where |P-0.5| >= tau (tau 0, 0.05, 0.1, 0.2)
  delivery: the MDE's own up/down balance re-solved (fu, fd) with the new signs: an A-balanced coarse move set is
  count-balanced per class (426h), so balance == class capacity; optionally PAYERS: extra up-mass on sub-gate genes
  (below the 5-CPM filter, invisible to every DE member; pds/mse pinned by per-gene totals), capacity = a fraction
  rho of the block's sub-gate counts in the high-library class moved to the low class.
Known answers: the reconstruction's own (kf, mde_fu/fd, n_up/dn, flips, n_block == diag), the delivered tables ==
harness metrics, and the design-sign / no-payer arm == delivered exactly (Δ 0).

    uv run python experiments/478-uncalled-sign-delivery/price478.py        (< 2 GB, ~1 min)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/490-de-emitter-review"))
import recon_h1 as RC  # noqa: E402
import score_de  # noqa: E402
from recon_h1 import K  # noqa: E402  (production rc_emit_k)

SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}
TAUS = (0.0, 0.05, 0.1, 0.2)
RHOS = (0.0, 0.25, 0.5)
SOURCE = "recip" if "--recip" in sys.argv else "clean4"


def prior_P(genes_full):
    cp = np.load(RC.PRIOR, allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in RC.PRIOR_BANKS]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.full(len(genes_full), np.nan)
    out[col >= 0] = P[col[col >= 0]]
    return out


def recip_P(r, keep, Pk):
    """491c's best recipient-side programme on the keep axis, expressed as a pseudo-P (0.5 - 0.5*tanh(z/2)) so the same
    |P-0.5| confidence thresholds apply: z = z(CLEAN4) + z(proj_k2 of CLEAN4 on H1's control PCs) + z(-identity)."""
    import h5py, csv
    genes = np.asarray(r["genes"]).astype(str)
    kg = genes[keep]
    f = h5py.File(ROOT / "data/harness/h1-test.h5ad", "r")
    cats = f["obs"]["target_gene"]["categories"][:].astype(str)
    codes = f["obs"]["target_gene"]["codes"][:]
    rows = np.where(codes == int(np.where(cats == "non-targeting")[0][0]))[0]
    indptr = f["X"]["indptr"][:]
    kept_col = np.full(genes.size, -1, np.int64)
    kept_col[np.where(keep)[0]] = np.arange(int(keep.sum()))
    X = np.zeros((rows.size, int(keep.sum())), np.float32)
    i = 0
    while i < rows.size:
        a = rows[i]; j = i
        while j + 1 < rows.size and rows[j + 1] < a + 4000:
            j += 1
        lo, hi = indptr[a], indptr[rows[j] + 1]
        dat = f["X"]["data"][lo:hi].astype(np.float64); ind = f["X"]["indices"][lo:hi]
        for k in range(i, j + 1):
            rr = rows[k]; s_, t_ = indptr[rr] - lo, indptr[rr + 1] - lo
            d, ix = dat[s_:t_], ind[s_:t_]; c = kept_col[ix]; m = c >= 0
            X[k, c[m]] = (1e6 * d[m] / d.sum()).astype(np.float32)
        i = j + 1
    f.close()
    m_all = X.mean(0).astype(np.float64)
    np.log1p(X, out=X)
    X -= X.mean(0)
    mC = np.where(np.isfinite(Pk), 0.5 - Pk, 0.0)
    rg = np.random.default_rng(0)
    Qm = X @ rg.standard_normal((X.shape[1], 22)).astype(np.float32)
    for _ in range(2):
        Qm, _ = np.linalg.qr(Qm)
        Qm = X @ (X.T @ Qm)
    Qm, _ = np.linalg.qr(Qm)
    _, _, Vt = np.linalg.svd(Qm.T @ X, full_matrices=False)
    V = Vt[:2].T.astype(np.float64)
    del X
    proj = V @ (V.T @ mC)
    offi = [row["gene_name"] for row in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    opos = {g: i for i, g in enumerate(offi)}
    abc = []
    for ctx in "ABC":
        mc = np.load(ROOT / f"experiments/491-own-control-programme/work/ctrl_meancpm_{ctx}.npy")
        abc.append(np.array([mc[opos[g]] if g in opos else np.nan for g in kg]))
    with np.errstate(invalid="ignore"):
        abc = np.nanmean(np.vstack([np.log(np.maximum(v, 1e-3)) for v in abc]), 0)
    ident = np.nan_to_num(np.log(np.maximum(m_all, 1e-3)) - abc)
    # block = the 2,200 top control-CPM kept genes is where 491c standardised; use every kept gene's sd here
    zs = lambda v: v / (np.std(v[v != 0]) + 1e-12)
    z = zs(np.where(np.isfinite(Pk), mC, 0.0)) + zs(proj) + zs(-ident)
    return 0.5 - 0.5 * np.tanh(z / 2.0)       # pseudo-P: z > 0 (up) -> P < 0.5


def mde_intent(o, sgn, madd, Ms, A0, E0v, base_t, blk, Lk, nk, n_car, lam, lam_m, move, payer_A=0.0):
    """recon_h1's MDE path with an arbitrary sign vector (fu/fd re-solved; payer up-capacity in the U side)."""
    elig = (Ms > 0) & np.isfinite(E0v) & np.isfinite(madd) & (madd > 0) & (sgn != 0)
    up, dn = elig & (sgn > 0), elig & (sgn < 0)
    Bm = A0 * np.power(2.0, np.where(np.isfinite(base_t - E0v), base_t - E0v, 0.0))
    U = lambda f: float((Bm[up] * (np.power(2.0, f * madd[up]) - 1.0)).sum())
    D = lambda f: float((Bm[dn] * (1.0 - np.power(2.0, -f * madd[dn]))).sum())
    U1, D1 = U(1.0) + payer_A, D(1.0)
    fu = fd = 1.0
    if U1 > D1 and U1 > 0:
        a_, b_ = 0.0, 1.0
        for _ in range(60):
            m_ = 0.5 * (a_ + b_)
            (b_ := m_) if U(m_) + payer_A * m_ > D1 else (a_ := m_)
        fu = 0.5 * (a_ + b_)
    elif D1 > U1 and D1 > 0:
        a_, b_ = 0.0, 1.0
        for _ in range(60):
            m_ = 0.5 * (a_ + b_)
            (b_ := m_) if D(m_) > U1 else (a_ := m_)
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
    return A_int, dict(fu=fu, fd=fd, n_up=int(up.sum()), n_dn=int(dn.sum()), U1=U1, D1=D1)


def run_seed(seed):
    tag = f"base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
    diag = json.loads((RC.H1 / f"diag-{tag}.json").read_text())
    st = np.load(RC.H1 / f"stats-{tag}.npz", allow_pickle=False)
    zr = np.load(RC.REF, allow_pickle=False)
    r = {k: zr[k] for k in ("genes", "targets", "keep", "ref_lin", "emit_ref_lin", "libraries", "b_ctrl")}
    zr.close()
    keep = np.asarray(r["keep"], dtype=bool)
    ref_emit = np.asarray(r["emit_ref_lin"], dtype=np.float64)
    ref_score = np.asarray(r["ref_lin"], dtype=np.float64)
    with RC.R.model_env("uniform"):
        b = RC.R.load_module("builder478", ROOT / "scripts/make_submission_pm_cd4.py")
    with np.load(RC.EFF, allow_pickle=False) as z:
        eff = RC.R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], r["genes"])
    prior_g = RC.prior_signs_full(r["genes"])[keep]
    Pk = prior_P(r["genes"])[keep]
    if SOURCE == "recip":
        Pk = recip_P(r, keep, Pk)
    cnt = ref_score * float(np.median(r["libraries"])) / 1e6
    madd = np.where(cnt >= RC.MDE_MIN, np.minimum(RC.MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), RC.MDE_CAP), 0.0)
    T = r["targets"].astype(str)
    G = np.asarray(r["genes"])[keep].astype(str)
    lfc_del = st["lfc"]
    m_del = (ref_score + 1e-9) * np.power(2.0, lfc_del) - 1e-9
    raw_o, per_o = score_de.score(T, G, st["p"], st["q"], lfc_del)
    want = json.loads((RC.H1 / f"metrics-{tag}.json").read_text())["raw"]
    assert max(abs(raw_o[m] - want[m]) for m in raw_o) <= 1e-9, "harness tables"
    pick_rng, _ = b.make_rngs(seed)
    blocks = []
    for i, t in enumerate(r["targets"]):
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        o = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks, fix=False)
        for name, (got, want_) in o["ka"].items():
            if want_ is not None and abs(float(got) - float(want_)) > (1e-9 if name.startswith(("kf", "mde_f")) else 0):
                raise AssertionError(f"seed {seed} {t}: known answer {name} {got} != {want_}")
        # the recon's internals needed to re-run the MDE with other signs
        n_car = o["n_car"]
        libs = np.asarray(r["libraries"], dtype=np.float64)[picks][:n_car]
        profile = b.profile_from_effect(eff[i], r["b_ctrl"], min_cpm=10.)
        L = int(min(np.median(libs), 1_000_000))
        totals = np.floor(n_car * L * profile).astype(np.int64)
        grand = int(totals.sum())
        Ms = totals[keep].astype(np.int64)
        sub_gate = int(totals[~keep].sum())
        _, Lk, _, nk = K.class_cells(n_car, RC.CLASSES3, grand)
        A0 = 1e6 * Ms.astype(np.float64) / float(grand)
        E0v = o["E0v"]
        kf = o["kf"]
        blk = o["blk"]
        # the block's target after the prior/kf (recon's base_t)
        thr = float(np.quantile(np.abs(E0v[blk]), RC.PROTECT_Q)) if blk.any() else 0.0
        tgt = np.where(blk & (np.abs(E0v) <= thr), prior_g * np.abs(E0v), E0v)
        side = blk & (tgt < 0)
        t2 = np.where(side, np.minimum(tgt * kf, -RC.FLOOR), tgt)
        move = np.where(blk, t2 - E0v, 0.0)
        base_t = np.where(blk, t2, E0v)
        sgn0 = np.sign(base_t)
        A_chk, info0 = mde_intent(o, sgn0, madd, Ms, A0, E0v, base_t, blk, Lk, nk, n_car, o["lam"], o["lam_m"], move)
        assert np.allclose(A_chk, o["A_int"], rtol=0, atol=1e-9), f"KA mde_intent != recon A_int ({t})"
        comp = ~o["head"] & ~blk
        lo_c, hi_c = int(np.argmin(Lk)), int(np.argmax(Lk))
        d_unit = (1e6 / n_car) * (1.0 / Lk[lo_c] - 1.0 / Lk[hi_c])
        hi_share = float(nk[hi_c] * Lk[hi_c] / (nk * Lk).sum())
        blocks.append(dict(i=i, o=o, Ms=Ms, A0=A0, E0v=E0v, base_t=base_t, blk=blk, Lk=Lk, nk=nk, n_car=n_car,
                           move=move, sgn0=sgn0, comp=comp, sub_gate=sub_gate, d_unit=d_unit, hi_share=hi_share))
    res = {"ka": "recon diag, harness tables, mde_intent == recon A_int: PASS"}
    rng = np.random.default_rng(478)
    for tau in TAUS:
        for rho in RHOS:
            lfc_new = lfc_del.copy()
            fds, flips = [], []
            for B in blocks:
                i = B["i"]
                cover = np.isfinite(Pk) & (np.abs(Pk - 0.5) >= tau) & B["comp"]
                sgn = B["sgn0"].copy()
                sgn[cover] = np.where(Pk[cover] > 0.5, -1.0, 1.0)
                payer = rho * B["hi_share"] * B["sub_gate"] * B["d_unit"]       # A-units of up-capacity from sub-gate genes
                A_new, inf = mde_intent(B["o"], sgn, madd, B["Ms"], B["A0"], B["E0v"], B["base_t"], B["blk"], B["Lk"],
                                        B["nk"], B["n_car"], B["o"]["lam"], B["o"]["lam_m"], B["move"], payer_A=payer)
                m_new = m_del[i] + (B["n_car"] / 400.0) * (A_new - B["o"]["A_int"])
                lfc_new[i] = np.log2((np.maximum(m_new, 0.0) + 1e-9) / (ref_score + 1e-9))
                fds.append(inf["fd"]); flips.append(int((sgn != B["sgn0"]).sum()))
            lfc_new = np.where(np.isfinite(lfc_del), lfc_new, lfc_del)
            raw_n, per_n = score_de.score(T, G, st["p"], st["q"], lfc_new)
            D = {m: np.array([per_n[m].get(t, np.nan) - per_o[m].get(t, np.nan) for t in T]) for m in SP}
            ov = lambda ii: sum(np.nanmean(D[m][ii]) / SP[m] for m in SP) / 6
            bo = np.array([ov(rng.integers(0, len(T), len(T))) for _ in range(2000)])
            key = f"tau{tau:g}_rho{rho:g}"
            res[key] = {**{m: round(float(np.nanmean(D[m])), 6) for m in SP}, "overall_equiv": round(float(ov(np.arange(len(T)))), 6),
                        "ci95": [round(float(np.quantile(bo, .025)), 6), round(float(np.quantile(bo, .975)), 6)],
                        "fd_median": float(np.median(fds)), "flips_median": float(np.median(flips))}
            if tau == 0.0 and rho == 0.0:
                pass
            print(f"seed {seed} tau {tau:<4} rho {rho:<4}: " + " ".join(f"{m} {res[key][m]:+.5f}" for m in SP)
                  + f" | ov {res[key]['overall_equiv']:+.5f} {res[key]['ci95']} | fd {np.median(fds):.3f} flips {np.median(flips):.0f}", flush=True)
    return res


def main():
    t0 = time.monotonic()
    out = {f"seed{s}": run_seed(s) for s in (0, 1, 2)}
    out["seconds"] = time.monotonic() - t0
    (HERE / "work" / ("price478_recip.json" if SOURCE == "recip" else "price478.json")).write_text(json.dumps(out, indent=1))
    print(f"done {out['seconds']:.0f}s")


if __name__ == "__main__":
    main()
