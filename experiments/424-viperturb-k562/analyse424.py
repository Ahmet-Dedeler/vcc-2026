"""424 analysis: every VIPerturb arm -> board prediction vs 417a's bulk (PP + 410 key clip(0.3,3) + sharpening q 1.25),
through 420's pricer UNCHANGED (analyse420.Pricer, primary route r417a: ballast-realised K16 rows on PLOLO for pds / cos,
the arm's realised A/B/C render energy for the cliff-aware mse model, overall = (2.2138 dpds_bd - 1.051 dmse_raw) / 6).

KNOWN ANSWERS (abort on failure):
  KA-P420  420's pricer on 420's own `PP|K3+sh1.5` row reproduces analyse420.json (k = 1 and 1.407, every fold).
  KA-BASE  my base row (comp424 `PP|K3+sh1.25`, abc424 energies) priced as an arm reads 0 exactly (d_rel, dpds, dmse).
The K562-GWPS fold holds K562 out: VIPerturb (K562) is excluded there, so every arm reads 0 on it by construction; the
table shows the three-fold cancer mean (420's convention) AND the two-fold mean over the folds where VIPerturb can act.

    uv run python analyse424.py        # -> work/analyse424.{json,log}
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
W420 = ROOT / "experiments/420-bulk-stack"
sys.path.insert(0, str(W420))
import analyse420 as A0  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = (*CANCER, "kolf2")
REF = "PP|K3+sh1.25"
K_KEY = A0.K_KEY
CTX = "ABC"
ARM_ORDER = ("i:extra|K3+sh1.25", "ii:avg|K3+sh1.25", "ii-:avgonly|K3+sh1.25", "ii':ivw|K3+sh1.25", "iii:fill|K3+sh1.25",
             "iii0:fillpure|K3+sh1.25", "CTRL:shuf|K3+sh1.25", "i:extra|K3refit+sh1.25", "PP|K3refit+sh1.25")
LABEL = {"i:extra|K3+sh1.25": "(i) extra donor, equal weight", "ii:avg|K3+sh1.25": "(ii) one K562 row: 0.5/0.5 avg + fill",
         "ii-:avgonly|K3+sh1.25": "(ii-) avg only, no fill", "ii':ivw|K3+sh1.25": "(ii') avg by cell counts + fill",
         "iii:fill|K3+sh1.25": "(iii) fill the 30 only",
         "iii0:fillpure|K3+sh1.25": "(iii0) fill, K562 median fixed", "CTRL:shuf|K3+sh1.25": "CONTROL (i) rows permuted",
         "i:extra|K3refit+sh1.25": "(i) + 410 key refit incl. VIP", "PP|K3refit+sh1.25": "base + refit key only"}


def main(tag="424", abc_file="abc424.json", comp_tag="comp424", arrays_tag="arrays424"):
    L = []
    say = lambda *a: (print(*a, flush=True), L.append(" ".join(str(x) for x in a)))
    P = A0.Pricer()
    R420 = A0.load_folds()
    abc = json.loads((WORK / abc_file).read_text())["designs"]
    out = dict(ka={}, rows={}, boot={})

    # ---------------------------------------------------------------- KA-P420
    a420 = json.loads((W420 / "work/analyse420.json").read_text())["rows"]["PP|K3+sh1.5"]
    worst = 0.0
    for f in FOLDS:
        for k, key in ((1.0, "r417a_k1"), (K_KEY, f"r417a_k{K_KEY:.3g}")):
            got = P.price(R420[f], f, "PP|K3+sh1.5", k)["ov"]
            worst = max(worst, abs(got - a420[key][f]["ov"]))
    say(f"KA-P420: 420's pricer on 420's PP|K3+sh1.5 reproduces analyse420.json: max |d ov| {worst:.1e}")
    assert worst < 1e-12, "KA-P420 FAILED"
    out["ka"]["p420"] = worst

    al = lambda n: n.replace(":", "~")
    mine = {f: json.loads((WORK / f"{comp_tag}-{f}.json").read_text()) for f in FOLDS}
    R = {}
    for f in FOLDS:
        R[f] = dict(R420[f])
        for n, e in mine[f]["arms"].items():
            if n == REF:
                assert e["realised"]["16"]["cos_agg"] == R420[f][REF]["realised"]["16"]["cos_agg"], "base != 420 ref"
                continue
            R[f][al(n)] = e        # 420's pricer reads ':' as its level-variant separator: price under an alias

    def E_of(n):
        return abc[n]["realised"]["16"]

    # the pricer looks every arm's energy up in its abc table: register mine (420's names are never overwritten)
    for n, e in abc.items():
        if n not in P.abc:
            P.abc[al(n)] = e
    P.abc["__base__"] = abc[REF]

    # ---------------------------------------------------------------- KA-BASE
    for f in FOLDS:
        base_entry = dict(R420[f][REF])
        R[f]["__base__"] = base_entry
        pr = P.price(R[f], f, "__base__", 1.0, E_arm=E_of(REF))
        assert pr["d_rel"] == 0 and pr["dpds_bd"] == 0 and abs(pr["dmse"]) < 1e-15, f"KA-BASE {f}: {pr}"
    say("KA-BASE: the base priced as an arm reads 0 on every fold")

    def price(n, f, k):
        if al(n) not in R[f]:      # k562_gwps fold: VIPerturb excluded -> the arm IS the base
            return dict(d_rel=0.0, dpds_bd=0.0, dpds_lit=0.0, dmse=0.0, ov=0.0, excluded=True)
        return P.price(R[f], f, al(n), k, E_arm=E_of(n))

    def boot(f, arm, n=2000, seed=424):
        Z = np.load(WORK / f"{arrays_tag}-{f}.npz")
        if f"{arm}@R16__x" not in Z.files:
            return None
        d = Z["d"]
        rng = np.random.default_rng(seed)
        g = lambda nm, fld: Z[f"{nm}@R16__{fld}"]
        xa, aa, pa, xr, ar, pr = g(arm, "x"), g(arm, "a"), g(arm, "pdsl"), g(REF, "x"), g(REF, "a"), g(REF, "pdsl")
        agg = lambda x, a, i: np.mean([abs(x[c, i].sum()) / np.sqrt(a[c, i].sum() * d[i].sum()) for c in range(3)])
        rel, dp = [], []
        for _ in range(n):
            i = rng.integers(0, d.size, d.size)
            rel.append(agg(xa, aa, i) / agg(xr, ar, i) - 1)
            dp.append(pa[:, i].mean() - pr[:, i].mean())
        return dict(p_rel=float(np.mean(np.asarray(rel) > 0)), p_pds=float(np.mean(np.asarray(dp) > 0)),
                    pds_ci=np.quantile(dp, [0.025, 0.975]).tolist())

    say(f"\n==== VIPerturb arms vs 417a ({REF}); route r417a (realised K16); folds k562 / hct116 / hek293t | kolf2 (no vote)")
    say("k562_gwps fold: VIPerturb excluded (same line as the held truth) -> 0 by construction. mean3 = 420's 3-fold cancer "
        "mean; mean2 = hct116 + hek293t only.")
    hdr = (f"{'arm':34s} {'d rel cos (hct/hek | kolf2)':>30s} {'dpds_bd (hct/hek | kolf2)':>30s} {'dmse k1 (hct/hek)':>19s} | "
           f"{'ov k=1 per fold':>26s} {'mean3':>8s} {'mean2':>8s} | {'ov k=1.407 per fold':>26s} {'mean3':>8s} {'mean2':>8s} | "
           f"{'kolf2 k1/kK':>15s} | {'A/B/C energy vs 417a':>20s} | boot P(rel>0)/P(pds>0) hct hek")
    say(hdr)
    for n in ARM_ORDER:
        if n not in abc:
            continue
        row = {f"k{k:.3g}": {f: price(n, f, k) for f in FOLDS} for k in (1.0, K_KEY)}
        out["rows"][n] = row
        bt = {f: boot(f, n) for f in ("orion_hct116", "orion_hek293t")}
        out["boot"][n] = bt
        k1, kK = row["k1"], row[f"k{K_KEY:.3g}"]
        ov = lambda r, fl: [r[f]["ov"] for f in fl]
        m3 = lambda r: float(np.mean(ov(r, CANCER)))
        m2 = lambda r: float(np.mean(ov(r, CANCER[1:])))
        er = [E_of(n)[c]["A"] / E_of(REF)[c]["A"] for c in CTX]
        flag = "HI" if max(er) > 1.03 else ("LO" if min(er) < 0.90 else "")
        dm = "/".join(f"{k1[f]['dmse']:+.4f}" for f in CANCER[1:])
        f2 = lambda key: "/".join(f"{k1[f][key]:+.4f}" for f in CANCER[1:]) + f" | {k1['kolf2'][key]:+.4f}"
        say(f"{LABEL.get(n, n):34s} {f2('d_rel'):>30s} {f2('dpds_bd'):>30s} "
            f"{dm:>19s} | "
            f"{A0.fmt_fold(ov(k1, CANCER)):>26s} {m3(k1):+8.5f} {m2(k1):+8.5f} | {A0.fmt_fold(ov(kK, CANCER)):>26s} "
            f"{m3(kK):+8.5f} {m2(kK):+8.5f} | {k1['kolf2']['ov']:+.4f}/{kK['kolf2']['ov']:+.4f} | "
            f"{'/'.join(f'{x:.3f}' for x in er):>20s} {flag:2s} | "
            + " ".join(f"{bt[f]['p_rel']:.2f}/{bt[f]['p_pds']:.2f}" if bt[f] else "-" for f in bt))
        out["rows"][n]["summary"] = dict(mean3_k1=m3(k1), mean2_k1=m2(k1), mean3_kK=m3(kK), mean2_kK=m2(kK),
                                         kolf2_k1=k1["kolf2"]["ov"], kolf2_kK=kK["kolf2"]["ov"], energy_vs_417a=er,
                                         flag=flag, n_pos_2fold_k1=int(sum(x > 0 for x in ov(k1, CANCER[1:]))))
    (WORK / f"analyse{tag}.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / f"analyse{tag}.log").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    import os
    v = os.environ.get("VCC424_VIP", "")
    tg = f"_{v}" if v else ""
    main(tag=f"424{tg}", abc_file=f"abc424{tg}.json", comp_tag=f"comp424{tg}", arrays_tag=f"arrays424{tg}")
