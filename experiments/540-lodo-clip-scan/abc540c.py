"""540c (A/B/C) — energies of 479 x 495: the LODO key at (0.5, 6) with P2g / S_P2, relative to 446bj2 (cap 3, q 1.75).
abc492g.py with the arms replaced (its known answers are re-run first).

(abc492g's docstring follows.)
492g (A/B/C) — per-target sharpening re-priced on 446bj2's REAL bulk: the stack keyed by the LODO key (clip 0.5-3).

Same route as abc492.py (itself 473's abc473.py): the stack dump's keyed bulk, sharpened (sign|h|^q on each context's 10-CPM
gate, non-panel columns only, each row's gated norm restored), times the per-context factor, rendered by production's
`profile_from_effect` (minnorm mass fix), then REALISED by 414's ballast emulation (KMAX 16) against the carriers' design
profiles. What changes vs abc492: the bulk is `dump_h1x2_cd4x15_lodo05.npz` (446bj2's VCC_BULK_LOAD), the factors are
446bj2's VCC_BULK_CTX_SCALE, the tercile key is the LODO artifact's clipped w (floor 0.5) with the floor guard.

KNOWN ANSWERS (abort if outside tolerance):
  KA-a  input shas (LODO dump, DE carriers dump, LODO key artifact, K3 key artifact, ntclevel-446bj.json) are the pinned
        ones; the LODO dump's `carriers` equal the K3 dump's (the same unkeyed stack) and the builder's apply_pert_scale on
        them with the LODO artifact's w reproduces the dump's `bulk` exactly.
  KA-b  the LODO bulk rendered at 446b's factors (= 446bj, which ntclevel measured) has realised energies within the
        documented route offset of ntclevel-446bj.json's measured A/B (abc492: -1.1 ... -2.6 %; only ratios are used).
  KA-c  the emulated ballast K median/max at 446bj2's factors = 446bj2's build meta (rc_callset/ballast_K: 7 / 16).

Also computed here (no PLOLO needed): the emulated ballast K of 433k's, 446b's and 446bj2's bulks (the ballast nmae cost
model of k492.py, and its test on the board's key step 446b -> 446bj2, where nmae moved +0.00006 raw).

    uv run python experiments/540-lodo-clip-scan/abc540c.py        (~1 GB, ~2 min)
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

plolo.WORK = WORK / "plolo"            # APFS clones of 414's cached challenge b_ctrl; nothing is written through plolo
from plolo import MIN_CPM, TS, profile_from_effect  # noqa: E402

W417 = ROOT / "experiments/417-bulk-sharpen/work"
DUMP = W417 / "dump_h1x2_cd4x15_lodo05.npz"
DUMP_SHA = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667"
DUMP_K3 = W417 / "dump_h1x2_cd4x15.npz"
DUMP_K3_SHA = "bb19f0e58bd00aa70e7aeb28b749fb3eee26d69dd62668e036da6370b8b792f2"
CARR = W417 / "dump_ka417a_g065.npz"
CARR_SHA = "aed59b0bffa4fc6d364ca92b2a4e6d0814a5d26b2bd33591f4e2d6180d0c187f"
KEY = ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz"
KEY_SHA = "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b"
KEY3 = ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz"
KEY3_SHA = "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"
NTC = ROOT / "experiments/446-q175-level/work/ntclevel-446bj.json"
Q = 1.75
FAC = {"A": 1.0962, "B": 1.1461, "C": 1.1681}      # 446bj2's VCC_BULK_CTX_SCALE (meta config bulk_ctx_scale)
F446 = {"A": 1.1004, "B": 1.1517, "C": 1.1720}     # 446b's = 446bj's (the configuration ntclevel-446bj measured)
F433 = {"A": 1.0464, "B": 1.092, "C": 1.114}       # 433k's (k492's calibration row)
S_C, LB, KMAX = 7_919_901.0, 995_000.0, 16          # ballast414.realise's constants (KMAX 16 = 446b's)
CTX = ("A", "B", "C")
LODO_FLOOR, K3_FLOOR = 0.5, 0.3


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def realise(p_de, p_bulk):
    """ballast414.realise, verbatim arithmetic (KMAX 16)."""
    live = p_de > 1e-9
    r = float(np.min(p_bulk[live] / p_de[live]))
    K = int(min(KMAX, max(1, int(np.ceil(S_C * (1.0 / max(r, 1e-6) - 1.0) / LB)))))
    S = S_C + K * LB
    want = S * p_bulk
    floor = S_C * p_de
    clipped = want < floor
    tot = np.where(clipped, floor, want)
    return tot / tot.sum(), r, K


def tercile(v, spec, floor=None):
    """492's knob: cut points np.quantile(v, [1/3, 2/3]) + np.digitize (low < q1 <= mid < q2 <= high); with `floor`, every
    target at the key's floor goes to LOW (the October guard). Returns the q vector and the tercile labels."""
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    if floor is not None:
        t = np.where(v <= floor + 1e-12, 0, t)
    return np.asarray([spec[k] for k in t], float), t


def main():
    t0 = time.time()
    log = lambda *a: print("[abc492g]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    got = {p.name: sha(p) for p in (DUMP, DUMP_K3, CARR, KEY, KEY3)}
    want = {DUMP.name: DUMP_SHA, DUMP_K3.name: DUMP_K3_SHA, CARR.name: CARR_SHA, KEY.name: KEY_SHA, KEY3.name: KEY3_SHA}
    assert got == want, f"KA-a input shas: {got}"
    got[NTC.name] = sha(NTC)
    z = np.load(DUMP, allow_pickle=False)
    perts = [str(p) for p in z["perts"]]
    genes = [str(g) for g in z["genes"]]
    bulk = np.asarray(z["bulk"], np.float64)
    stack_car = np.asarray(z["carriers"], np.float64)
    z3 = np.load(DUMP_K3, allow_pickle=False)
    assert [str(p) for p in z3["perts"]] == perts and [str(g) for g in z3["genes"]] == genes
    assert np.array_equal(np.asarray(z3["carriers"], np.float64), stack_car), "KA-a: LODO and K3 dumps' carriers differ"
    bulk_k3 = np.asarray(z3["bulk"], np.float64)
    zc = np.load(CARR, allow_pickle=False)
    assert [str(p) for p in zc["perts"]] == perts and [str(g) for g in zc["genes"]] == genes
    carriers = np.asarray(zc["carriers"], np.float64)
    assert tuple(genes) == plolo.OFFICIAL, "official gene axis"
    gpos = {g: i for i, g in enumerate(genes)}
    pp = np.asarray([gpos[t] for t in perts])
    nonpanel = np.ones(len(genes), bool)
    nonpanel[pp] = False
    own = np.zeros((len(perts), len(genes)), bool)      # each row's own target (the scorer drops it)
    own[np.arange(len(perts)), pp] = True
    zk = np.load(KEY, allow_pickle=False)
    assert [str(t) for t in zk["targets"]] == perts, "LODO key target order"
    w = np.asarray(zk["w"], np.float64)
    zk3 = np.load(KEY3, allow_pickle=False)
    k3pos = {str(t): i for i, t in enumerate(zk3["targets"])}
    w3 = np.asarray([float(zk3["w"][k3pos[t]]) for t in perts])

    def keyed(wv):
        """the builder's apply_pert_scale: rows x w c, c restoring the non-panel Frobenius norm."""
        e2 = np.einsum("ij,ij->i", stack_car[:, nonpanel], stack_car[:, nonpanel])
        c = float(np.sqrt(e2.sum() / float((wv * wv * e2).sum())))
        return stack_car * (wv * c)[:, None]
    ka_key = float(np.max(np.abs(keyed(w) - bulk)))
    ka_key3 = float(np.max(np.abs(keyed(w3) - bulk_k3)))
    assert ka_key < 1e-9 and ka_key3 < 1e-9, f"KA-a key scale: {ka_key} {ka_key3}"
    ntc = json.loads(NTC.read_text())["per_context"]
    frames = {}
    for c in CTX:
        b = plolo.challenge_bctrl(c)
        cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
        frames[c] = (b, (cpm > MIN_CPM) & (b > 0))
    log(f"KA-a PASS (shas; carriers identical; key scale |d| {ka_key:.1e} / K3 {ka_key3:.1e}); gates "
        f"{[int(frames[c][1].sum()) for c in CTX]}; LODO floor share {np.mean(w <= LODO_FLOOR + 1e-12):.3f}, "
        f"K3 {np.mean(w3 <= K3_FLOOR + 1e-12):.3f}")

    # carriers' design profiles (the DE side is untouched by every arm)
    P_DE = {c: [profile_from_effect(carriers[i], frames[c][0], frames[c][1], "minnorm") for i in range(len(perts))]
            for c in CTX}

    def bsharpen(B, gate, qv):
        """the builder's bulk sharpening, per-target q: sign|h|^q on the gated NON-PANEL genes, that part's row norm
        restored, panel columns (incl. the row's knockdown) untouched."""
        G = gate & nonpanel
        out = B.copy()
        H = B[:, G]
        qq = np.asarray(qv, float)
        S = np.sign(H) * np.abs(H) ** qq[:, None]
        n0, n1 = np.linalg.norm(H, axis=1), np.linalg.norm(S, axis=1)
        out[:, G] = S * (n0 / np.maximum(n1, 1e-300))[:, None]
        return out

    def render(B, qv, fac, fmult=1.0):
        """per context: sharpen, factor, render, realise. Returns energies (design/realised; A, B with and without the own
        target), per-row realised energies, K, max single-gene mass share of the realised blocks."""
        out = {}
        for c in CTX:
            b, gate = frames[c]
            F = bsharpen(B, gate, qv) * fac[c] * fmult
            n = F.shape[0]
            des = np.empty_like(F)
            rea = np.empty_like(F)
            K = np.empty(n, int)
            share = np.empty(n)
            for i in range(n):
                p = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p) - b
                pr, r, k = realise(P_DE[c][i], p)
                rea[i] = np.log1p(TS * pr) - b
                K[i] = k
                share[i] = float(pr.max())
            E = {}
            wgt = 300 / 18532               # ntclevel's weights: sum_p |row_p|^2 / G_p, G_p = 18532
            for tag, H in (("design", des), ("realised", rea)):
                Hn = H[:, nonpanel]
                Hc = Hn - Hn.mean(0, keepdims=True)
                Ho = np.where(own, 0.0, H)                        # ntclevel's convention: every gene but the own target
                Hoc = Ho - Ho.mean(0, keepdims=True)
                E[tag] = dict(A=float((Hn * Hn).sum(1).mean() * wgt), B=float((Hc * Hc).sum(1).mean() * wgt),
                              A_own=float((Ho * Ho).sum(1).mean() * wgt), B_own=float((Hoc * Hoc).sum(1).mean() * wgt))
            out[c] = dict(E=E, row=np.einsum("ij,ij->i", rea[:, nonpanel], rea[:, nonpanel]), K=K, share=share)
        return out

    n = len(perts)
    q175 = np.full(n, Q)
    res = dict(ka={}, arms={}, inputs=got, factors=FAC, q=Q)
    # ---------------- KA-b: the LODO bulk at 446b's factors (= 446bj) vs ntclevel-446bj's measured energies
    bj = render(bulk, q175, F446)
    ka_b = {c: dict(A=bj[c]["E"]["realised"]["A_own"] / ntc[c]["A"] - 1, B=bj[c]["E"]["realised"]["B_own"] / ntc[c]["B"] - 1,
                    A_design=bj[c]["E"]["design"]["A_own"] / ntc[c]["A"] - 1) for c in CTX}
    for c in CTX:
        log(f"KA-b {c} (446bj): realised A_own {bj[c]['E']['realised']['A_own']:.5f} vs measured {ntc[c]['A']:.5f} "
            f"({100 * ka_b[c]['A']:+.2f} %), B {100 * ka_b[c]['B']:+.2f} %; design A {100 * ka_b[c]['A_design']:+.2f} %")
    worst = max(abs(v) for c in CTX for k, v in ka_b[c].items() if k in ("A", "B"))
    res["ka"]["b_rel_err_446bj"] = ka_b
    assert worst < 0.05, f"KA-b grossly off: {worst:.4f}"
    Kbj = np.concatenate([bj[c]["K"] for c in CTX])
    # ---------------- base = 446bj2 (LODO bulk, q 1.75, 446bj2's factors); KA-c
    base = render(bulk, q175, FAC)
    Kall = np.concatenate([base[c]["K"] for c in CTX])
    res["ka"]["c_K"] = dict(median=float(np.median(Kall)), max=int(Kall.max()), min=int(Kall.min()),
                            median_446bj=float(np.median(Kbj)), frac_capped=float(np.mean(Kall >= KMAX)))
    log(f"KA-c emulated ballast K at 446bj2's factors: median {np.median(Kall):.0f} max {Kall.max()} min {Kall.min()} "
        f"(meta 7 / 16 / 1); at 446bj's: median {np.median(Kbj):.0f} (meta 7)")
    assert np.median(Kall) == 7 and Kall.max() == 16, "KA-c FAILED"
    # energy of 446bj2 relative to 446bj in this route (anchors the cliff model at 446bj2 from ntclevel-446bj)
    res["bj2_over_bj"] = {c: dict(A_own=base[c]["E"]["realised"]["A_own"] / bj[c]["E"]["realised"]["A_own"],
                                  B_own=base[c]["E"]["realised"]["B_own"] / bj[c]["E"]["realised"]["B_own"]) for c in CTX}
    log("446bj2/446bj realised energy: " + "  ".join(f"{c} A x{res['bj2_over_bj'][c]['A_own']:.4f} B x{res['bj2_over_bj'][c]['B_own']:.4f}"
                                                     for c in CTX))
    log("KA-a, KA-c PASS; KA-b within the documented route offset: arms read as ratios")

    # ---------------- arms: cap 6 (the artifact's key re-clipped at (0.5, 6)) x the per-target q
    key = np.asarray(zk["key"], np.float64)
    car = np.asarray(zk["carried"], bool)
    med = float(np.median(key[car]))
    w6 = np.where(car, np.clip(key / med, 0.5, 6.0), 1.0)
    bulk6 = keyed(w6)
    strength = np.linalg.norm(stack_car[:, nonpanel], axis=1)
    P2 = (2.5, 2.0, 1.75)
    qp2g6, ter6 = tercile(w6, P2, LODO_FLOOR)
    qp2g3, ter3 = tercile(w, P2, LODO_FLOOR)
    qs, _ = tercile(strength, P2)
    res["tercile_same_cap3_cap6"] = bool(np.array_equal(ter3, ter6))
    log(f"P2g terciles at cap 6 {np.bincount(ter6, minlength=3).tolist()}, identical to cap 3's: {res['tercile_same_cap3_cap6']}")
    arms = {"cap3 + P2g": (bulk, qp2g3), "cap3 + S_P2": (bulk, qs), "cap6": (bulk6, q175), "cap6 + P2g": (bulk6, qp2g6),
            "cap6 + S_P2": (bulk6, qs)}
    for name, (B, qv) in arms.items():
        R = render(B, qv, FAC)
        rec = dict(E={c: R[c]["E"] for c in CTX}, E_base={c: base[c]["E"] for c in CTX})
        Karm = np.concatenate([R[c]["K"] for c in CTX])
        for c in CTX:
            ratio = R[c]["row"] / np.maximum(base[c]["row"], 1e-300)
            rec[c] = dict(energy_ratio_realised=R[c]["E"]["realised"]["A_own"] / base[c]["E"]["realised"]["A_own"],
                          B_ratio_realised=R[c]["E"]["realised"]["B_own"] / base[c]["E"]["realised"]["B_own"],
                          max_row_ratio=float(ratio.max()), min_row_ratio=float(ratio.min()))
        rec["K_mean"] = float(Karm.mean())
        res["arms"][name] = rec
        log(f"{name:14s} " + "  ".join(f"{c}: E x{rec[c]['energy_ratio_realised']:.4f} row max x{rec[c]['max_row_ratio']:.2f}" for c in CTX)
            + f" | K mean {rec['K_mean']:.2f}")
    res["base_K_mean"] = float(Kall.mean())
    maps = {}
    FM = (0.92, 0.95, 0.98, 1.0, 1.02, 1.05, 1.08)
    for tag, (B, qv) in (("cap6", (bulk6, q175)), ("cap6 + P2g", (bulk6, qp2g6)), ("cap6 + S_P2", (bulk6, qs))):
        maps[tag] = {f"{fm:g}": {c: R[c]["E"] for c, R in [(c, RR) for RR in [render(B, qv, FAC, fmult=fm)] for c in CTX]} for fm in FM}
    res["level_maps"] = dict(factors=FM, maps=maps)
    (WORK / "abc540c.json").write_text(json.dumps(res, indent=1, default=float))
    log("done")


if __name__ == "__main__":
    main()
