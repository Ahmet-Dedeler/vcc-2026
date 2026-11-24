"""518 (A/B/C) — the strength key on the real panel: render energies, level maps and ballast, for the 540 pricer.

abc540.py's route verbatim (446bj2's stack dump carriers keyed by the LODO artifact 8c82bef3's key vector, sharpen q 1.75
per context gate, x 446bj2's factors, render (minnorm), realise by 414's ballast emulation at KMAX 16 against the DE
carriers), with each arm's per-target weight = clip(key/median, 0.5, cap) x exp(theta z_p):
    z_p = centred percentile rank of the robust response strength (strength518) over
          'src'  the three production sources (K562-GWPS, Orion HCT116, Orion HEK293T), or
          'meta' those plus CD4 (DESeq2 nDE), VIPerturb-K562 and the KOLF2 atlas panel rows,
    the mean over the screens that carry target p (0 if none).
KNOWN ANSWERS (abort): the artifact's key re-clipped at (0.5, 3) reproduces its w; keying the dump's carriers reproduces the
dump's bulk (|d| < 1e-9); the base's emulated ballast K median / max = 446bj2's meta (7 / 16); cap 6's energy ratios and
dK_mean reproduce abc540.json's `V1_0.5_6` (1e-12).

    uv run python experiments/518-bulk-mse-gap/abc518.py        (~1 GB)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/540-lodo-clip-scan"))
sys.path.insert(0, str(HERE))
import abc540 as A5  # noqa: E402  (its constants and realise(); its main() is never called)
from lofo518 import crank  # noqa: E402

A5.plolo.WORK = ROOT / "experiments/540-lodo-clip-scan/work/plolo"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
SRC = ("k562_gwps", "orion_hct116", "orion_hek293t")
META = SRC + ("cd4_deseq", "viperturb_k562", "kolf2")
ARMS = [("cap3", 3.0, None, 0.0), ("cap6", 6.0, None, 0.0),
        ("meta th1.25 on cap3", 3.0, "meta", 1.25), ("meta th1 on cap3", 3.0, "meta", 1.0), ("src th1 on cap3", 3.0, "src", 1.0),
        ("meta th0.5 on cap6", 6.0, "meta", 0.5), ("meta th0.75 on cap6", 6.0, "meta", 0.75)]


def zprod(perts, kind):
    S = np.load(WORK / "strength518.npz")
    panel = [str(p) for p in S["panel"]]
    idx = np.asarray([panel.index(p) for p in perts])
    cols = []
    for k in (SRC if kind == "src" else META):
        v = np.asarray(S[k], float)[idx]
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    m = np.nanmean(np.vstack(cols), 0)
    return np.where(np.isfinite(m), m, 0.0)


def main():
    t0 = time.time()
    log = lambda *a: print("[abc518]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    got = {p.name: A5.sha(p) for p in (A5.DUMP, A5.CARR, A5.KEY)}
    assert got == {A5.DUMP.name: A5.DUMP_SHA, A5.CARR.name: A5.CARR_SHA, A5.KEY.name: A5.KEY_SHA}, got
    z = np.load(A5.DUMP, allow_pickle=False)
    perts = [str(p) for p in z["perts"]]
    genes = [str(g) for g in z["genes"]]
    stack_car = np.asarray(z["carriers"], np.float64)
    bulk0 = np.asarray(z["bulk"], np.float64)
    carriers = np.asarray(np.load(A5.CARR, allow_pickle=False)["carriers"], np.float64)
    zk = np.load(A5.KEY, allow_pickle=False)
    key, car = np.asarray(zk["key"], np.float64), np.asarray(zk["carried"], bool)
    med = float(np.median(key[car]))
    gpos = {g: i for i, g in enumerate(genes)}
    nonpanel = np.ones(len(genes), bool)
    nonpanel[[gpos[t] for t in perts]] = False
    own = np.zeros((len(perts), len(genes)), bool)
    own[np.arange(len(perts)), [gpos[t] for t in perts]] = True
    n = len(perts)
    w_of = lambda lo, hi: np.where(car, np.clip(key / med, lo, hi), 1.0)

    def keyed(wv):
        e2 = np.einsum("ij,ij->i", stack_car[:, nonpanel], stack_car[:, nonpanel])
        c = float(np.sqrt(e2.sum() / float((wv * wv * e2).sum())))
        return stack_car * (wv * c)[:, None], c
    assert np.array_equal(w_of(0.5, 3.0), np.asarray(zk["w"], np.float64)), "KA re-clip"
    kb, _ = keyed(w_of(0.5, 3.0))
    assert float(np.max(np.abs(kb - bulk0))) < 1e-9, "KA key scale"
    frames, P_DE, TOP, LDE = {}, {}, {}, {}
    for c in A5.CTX:
        b = A5.plolo.challenge_bctrl(c)
        pc = np.expm1(b) / np.expm1(b).sum()
        gate = (1e6 * pc > A5.MIN_CPM) & (b > 0)
        frames[c] = (b, gate, pc)
        P_DE[c] = [A5.profile_from_effect(carriers[i], b, gate, "minnorm") for i in range(n)]
        G = np.flatnonzero(gate & nonpanel)
        TOP[c] = np.stack([G[np.argsort(-np.abs(A5.lfc(P_DE[c][i][G], pc[G])))[:A5.TOPN]] for i in range(n)])
        LDE[c] = np.stack([A5.lfc(P_DE[c][i][TOP[c][i]], pc[TOP[c][i]]) for i in range(n)])

    def sharpen(B, gate):
        G = gate & nonpanel
        out = B.copy()
        H = B[:, G]
        S = np.sign(H) * np.abs(H) ** A5.Q
        out[:, G] = S * (np.linalg.norm(H, axis=1) / np.maximum(np.linalg.norm(S, axis=1), 1e-300))[:, None]
        return out

    def render(B, fmult=1.0, full=True):
        out = {}
        wgt = 300 / 18532
        for c in A5.CTX:
            b, gate, pc = frames[c]
            F = sharpen(B, gate) * A5.FAC[c] * fmult
            rea = np.empty_like(F)
            K = np.empty(n, int)
            L = np.empty((n, A5.TOPN))
            for i in range(n):
                p = A5.profile_from_effect(F[i], b, gate, "minnorm")
                pr, k, tot = A5.realise(P_DE[c][i], p)
                rea[i] = np.log1p(A5.TS * pr) - b
                K[i] = k
                if full:
                    Qb = tot - A5.S_C * P_DE[c][i]
                    m = ((A5.NCELL - k) * P_DE[c][i] + k * Qb / max(Qb.sum(), 1e-300)) / A5.NCELL
                    L[i] = A5.lfc(m[TOP[c][i]], pc[TOP[c][i]])
            Ho = np.where(own, 0.0, rea)
            Hoc = Ho - Ho.mean(0, keepdims=True)
            out[c] = dict(E=dict(A_own=float((Ho * Ho).sum(1).mean() * wgt), B_own=float((Hoc * Hoc).sum(1).mean() * wgt)),
                          row=np.einsum("ij,ij->i", rea[:, nonpanel], rea[:, nonpanel]), K=K, L=L)
        return out

    base = render(bulk0)
    Kb = np.concatenate([base[c]["K"] for c in A5.CTX])
    assert np.median(Kb) == 7 and Kb.max() == 16, "KA base ballast K"
    ref540 = json.loads((W540 / "abc540.json").read_text())
    res = dict(arms={}, maps={}, z={})
    for tag, cap, kind, th in ARMS:
        m = np.exp(th * zprod(perts, kind)) if kind else np.ones(n)
        if kind:
            res["z"][tag] = dict(m_quantiles=np.quantile(m, [0, .1, .25, .5, .75, .9, 1]).tolist())
        wv = w_of(0.5, cap) * m
        B, cst = keyed(wv)
        R = base if tag == "cap3" else render(B)
        Ka = np.concatenate([R[c]["K"] for c in A5.CTX])
        d_from = float(np.mean(np.concatenate([(np.abs(R[c]["L"] - LDE[c]) - np.abs(base[c]["L"] - LDE[c])).ravel() for c in A5.CTX])))
        rec = dict(cap=cap, kind=kind, theta=th, constant=cst, K_mean=float(Ka.mean()), dK_mean=float(Ka.mean() - Kb.mean()),
                   D_from_design=d_from)
        for c in A5.CTX:
            ratio = R[c]["row"] / np.maximum(base[c]["row"], 1e-300)
            rec[c] = dict(A_ratio=R[c]["E"]["A_own"] / base[c]["E"]["A_own"], B_ratio=R[c]["E"]["B_own"] / base[c]["E"]["B_own"],
                          max_row_ratio=float(ratio.max()), min_row_ratio=float(ratio.min()))
        if tag == "cap6":
            r5 = ref540["arms"]["V1_0.5_6"]
            dd = max(abs(rec[c][k] - r5[c][k]) for c in A5.CTX for k in ("A_ratio", "B_ratio")) + abs(rec["dK_mean"] - r5["dK_mean"])
            assert dd < 1e-12, f"KA cap 6 vs abc540: {dd}"
            log(f"KA cap 6 == abc540 V1_0.5_6 (|d| {dd:.1e})")
        mp = {}
        for fm in A5.FM:
            RR = R if fm == 1.0 else render(B, fmult=fm, full=False)
            mp[f"{fm:g}"] = {c: RR[c]["E"] for c in A5.CTX}
        res["maps"][tag] = mp
        res["arms"][tag] = rec
        log(f"{tag:22s} E " + " ".join(f"{c} x{rec[c]['A_ratio']:.4f}" for c in A5.CTX) +
            f" | row max x{max(rec[c]['max_row_ratio'] for c in A5.CTX):.2f} min x{min(rec[c]['min_row_ratio'] for c in A5.CTX):.2f}"
            f" | K {rec['K_mean']:.2f} ({rec['dK_mean']:+.2f}) D {d_from:+.5f}")
    (WORK / "abc518.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
