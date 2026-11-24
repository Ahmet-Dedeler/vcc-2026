"""540b (A/B/C) — the cap extension of 540 (caps 4 ... uncapped; same code, CLIPS changed).

540 (A/B/C) — the LODO key's clip (floor, cap) around 446bj2's (0.5, 3): render energies, level maps and ballast.

446bj2's bulk = the builder's apply_pert_scale of the stack dump's unkeyed `carriers` with 420's LODO artifact w =
clip(key / median_carried(key), 0.5, 3) on carried targets, 1 elsewhere (`pert_scale_415pp_lo0.5_hi3.npz`, 8c82bef3).
Each arm re-clips the SAME key vector at (lo, hi), re-keys the same carriers (the energy-preserving constant restores the
non-panel Frobenius norm), then abc492g's route: sharpen q 1.75 on each context's gate (non-panel, row norm restored),
x 446bj2's factors, render (minnorm), realise by 414's ballast emulation (KMAX 16) against the DE carriers' design.

Per arm: realised energy ratios vs 446bj2 (A_own/B_own, ntclevel's convention), per-row ratio extremes, ballast K, the
DE-displacement proxy of ballast492g (top-100 design genes per block), and a realised energy map at factor multipliers
0.92 ... 1.08 for the re-level.
KNOWN ANSWERS (abort): the artifact's key re-clipped at (0.5, 3) reproduces its w exactly; keying the dump's carriers with
it reproduces the dump's bulk (|d| 0); the base's emulated K median/max = 446bj2's meta (7 / 16).

    uv run python experiments/540-lodo-clip-scan/abc540b.py        (~1 GB, a few minutes)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"            # APFS clones of the cached challenge b_ctrl; nothing is written through plolo
from plolo import MIN_CPM, TS, profile_from_effect  # noqa: E402

W417 = ROOT / "experiments/417-bulk-sharpen/work"
DUMP = W417 / "dump_h1x2_cd4x15_lodo05.npz"
DUMP_SHA = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667"
CARR = W417 / "dump_ka417a_g065.npz"
CARR_SHA = "aed59b0bffa4fc6d364ca92b2a4e6d0814a5d26b2bd33591f4e2d6180d0c187f"
KEY = ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz"
KEY_SHA = "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b"
Q = 1.75
FAC = {"A": 1.0962, "B": 1.1461, "C": 1.1681}      # 446bj2's VCC_BULK_CTX_SCALE
S_C, LB, KMAX, NCELL, TOPN = 7_919_901.0, 995_000.0, 16, 400, 100
EPS, LFLOOR = 1e-15, -5.0
CTX = ("A", "B", "C")
CLIPS = [(0.5, 3.0), (0.5, 4.0), (0.5, 6.0), (0.5, 8.0), (0.5, 10.0), (0.5, 15.0), (0.5, 1e9), (0.3, 6.0), (0.4, 6.0),
         (0.3, 1e9)]
FM = (0.92, 0.95, 0.98, 1.0, 1.02, 1.05, 1.08)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tag(lo, hi):
    return f"V1_{lo:g}_{hi:g}"


def realise(p_de, p_bulk):
    """ballast414.realise (KMAX 16), returning the normalised realised profile, K and the block totals."""
    live = p_de > 1e-9
    r = float(np.min(p_bulk[live] / p_de[live]))
    K = int(min(KMAX, max(1, int(np.ceil(S_C * (1.0 / max(r, 1e-6) - 1.0) / LB)))))
    S = S_C + K * LB
    want = S * p_bulk
    floor = S_C * p_de
    tot = np.where(want < floor, floor, want)
    return tot / tot.sum(), K, tot


def lfc(m, pc):
    return np.maximum(np.log2((m + EPS) / (pc + EPS)), LFLOOR)


def main():
    t0 = time.time()
    log = lambda *a: print("[abc540]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    got = {p.name: sha(p) for p in (DUMP, CARR, KEY)}
    assert got == {DUMP.name: DUMP_SHA, CARR.name: CARR_SHA, KEY.name: KEY_SHA}, f"input shas: {got}"
    z = np.load(DUMP, allow_pickle=False)
    perts = [str(p) for p in z["perts"]]
    genes = [str(g) for g in z["genes"]]
    stack_car = np.asarray(z["carriers"], np.float64)
    bulk0 = np.asarray(z["bulk"], np.float64)
    zc = np.load(CARR, allow_pickle=False)
    assert [str(p) for p in zc["perts"]] == perts and [str(g) for g in zc["genes"]] == genes
    carriers = np.asarray(zc["carriers"], np.float64)
    zk = np.load(KEY, allow_pickle=False)
    assert [str(t) for t in zk["targets"]] == perts
    key, car = np.asarray(zk["key"], np.float64), np.asarray(zk["carried"], bool)
    med = float(np.median(key[car]))
    gpos = {g: i for i, g in enumerate(genes)}
    nonpanel = np.ones(len(genes), bool)
    nonpanel[[gpos[t] for t in perts]] = False
    own = np.zeros((len(perts), len(genes)), bool)
    own[np.arange(len(perts)), [gpos[t] for t in perts]] = True
    n = len(perts)

    def w_of(lo, hi):
        return np.where(car, np.clip(key / med, lo, hi), 1.0)

    def keyed(wv):
        e2 = np.einsum("ij,ij->i", stack_car[:, nonpanel], stack_car[:, nonpanel])
        c = float(np.sqrt(e2.sum() / float((wv * wv * e2).sum())))
        return stack_car * (wv * c)[:, None], c
    assert np.array_equal(w_of(0.5, 3.0), np.asarray(zk["w"], np.float64)), "KA: re-clip at (0.5, 3) != artifact w"
    kb, _ = keyed(w_of(0.5, 3.0))
    ka_key = float(np.max(np.abs(kb - bulk0)))
    assert ka_key < 1e-9, f"KA key scale {ka_key}"
    frames, P_DE, TOP, LDE = {}, {}, {}, {}
    for c in CTX:
        b = plolo.challenge_bctrl(c)
        pc = np.expm1(b) / np.expm1(b).sum()
        gate = (1e6 * pc > MIN_CPM) & (b > 0)
        frames[c] = (b, gate, pc)
        P_DE[c] = [profile_from_effect(carriers[i], b, gate, "minnorm") for i in range(n)]
        G = np.flatnonzero(gate & nonpanel)
        TOP[c] = np.stack([G[np.argsort(-np.abs(lfc(P_DE[c][i][G], pc[G])))[:TOPN]] for i in range(n)])
        LDE[c] = np.stack([lfc(P_DE[c][i][TOP[c][i]], pc[TOP[c][i]]) for i in range(n)])
    log(f"KA key PASS (re-clip exact; key scale |d| {ka_key:.1e}); carried {int(car.sum())}")

    def sharpen(B, gate):
        G = gate & nonpanel
        out = B.copy()
        H = B[:, G]
        S = np.sign(H) * np.abs(H) ** Q
        out[:, G] = S * (np.linalg.norm(H, axis=1) / np.maximum(np.linalg.norm(S, axis=1), 1e-300))[:, None]
        return out

    def render(B, fmult=1.0, full=True):
        out = {}
        wgt = 300 / 18532
        for c in CTX:
            b, gate, pc = frames[c]
            F = sharpen(B, gate) * FAC[c] * fmult
            rea = np.empty_like(F)
            K = np.empty(n, int)
            L = np.empty((n, TOPN))
            share = np.empty(n)
            for i in range(n):
                p = profile_from_effect(F[i], b, gate, "minnorm")
                pr, k, tot = realise(P_DE[c][i], p)
                rea[i] = np.log1p(TS * pr) - b
                K[i] = k
                share[i] = float(pr.max())
                if full:
                    Qb = tot - S_C * P_DE[c][i]
                    m = ((NCELL - k) * P_DE[c][i] + k * Qb / max(Qb.sum(), 1e-300)) / NCELL
                    L[i] = lfc(m[TOP[c][i]], pc[TOP[c][i]])
            Ho = np.where(own, 0.0, rea)
            Hoc = Ho - Ho.mean(0, keepdims=True)
            out[c] = dict(E=dict(A_own=float((Ho * Ho).sum(1).mean() * wgt), B_own=float((Hoc * Hoc).sum(1).mean() * wgt)),
                          row=np.einsum("ij,ij->i", rea[:, nonpanel], rea[:, nonpanel]), K=K, L=L, share=share)
        return out

    res = dict(inputs=got, factors=FAC, q=Q, arms={}, maps={})
    base = render(bulk0)
    Kb = np.concatenate([base[c]["K"] for c in CTX])
    assert np.median(Kb) == 7 and Kb.max() == 16, "KA-c: base ballast K"
    log(f"KA-c PASS: base K median {np.median(Kb):.0f} max {Kb.max()} mean {Kb.mean():.3f}")
    for lo, hi in CLIPS:
        wv = w_of(lo, hi)
        B, cst = keyed(wv)
        R = base if (lo, hi) == (0.5, 3.0) else render(B)
        Ka = np.concatenate([R[c]["K"] for c in CTX])
        d_from = float(np.mean(np.concatenate([(np.abs(R[c]["L"] - LDE[c]) - np.abs(base[c]["L"] - LDE[c])).ravel() for c in CTX])))
        rec = dict(lo=lo, hi=hi, constant=cst, frac_floor=float(np.mean(car & (wv <= lo + 1e-12))),
                   frac_cap=float(np.mean(car & (wv >= hi - 1e-12))), K_mean=float(Ka.mean()), dK_mean=float(Ka.mean() - Kb.mean()),
                   D_from_design=d_from)
        for c in CTX:
            ratio = R[c]["row"] / np.maximum(base[c]["row"], 1e-300)
            rec[c] = dict(A_ratio=R[c]["E"]["A_own"] / base[c]["E"]["A_own"], B_ratio=R[c]["E"]["B_own"] / base[c]["E"]["B_own"],
                          max_row_ratio=float(ratio.max()), min_row_ratio=float(ratio.min()),
                          max_mass_share=float(R[c]["share"].max()), base_max_mass_share=float(base[c]["share"].max()),
                          frac_capped=float(np.mean(R[c]["K"] >= KMAX)))
        rec["E"] = {c: R[c]["E"] for c in CTX}
        mp = {}
        for fm in FM:
            RR = R if fm == 1.0 else render(B, fmult=fm, full=False)
            mp[f"{fm:g}"] = {c: RR[c]["E"] for c in CTX}
        res["maps"][tag(lo, hi)] = mp
        res["arms"][tag(lo, hi)] = rec
        log(f"{tag(lo, hi):12s} c {cst:.4f} floor {rec['frac_floor']:.3f} cap {rec['frac_cap']:.3f} | E " +
            " ".join(f"{c} x{rec[c]['A_ratio']:.4f}" for c in CTX) + f" | row max x{max(rec[c]['max_row_ratio'] for c in CTX):.2f} "
            f"min x{min(rec[c]['min_row_ratio'] for c in CTX):.2f} | K {rec['K_mean']:.2f} ({rec['dK_mean']:+.2f}) D {d_from:+.5f}")
    (WORK / "abc540b.json").write_text(json.dumps(res, indent=1, default=float))
    log("done")


if __name__ == "__main__":
    main()
