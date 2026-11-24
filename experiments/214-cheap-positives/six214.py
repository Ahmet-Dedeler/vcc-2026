"""214 — the four DE members for arms A, B and A∘B, marginal to `cd4f034`, on h1-test.

Three pools, exactly 196/208's constructions:

  native   the reference's own BH<0.05 pool and our own BH<0.05 call set -- what a
           submission would actually score.
  wald89   188's SIGNIFICANCE-selected pool (`|lfc| * sqrt(control cpm)` at the reference's
           own BH depth, median `n_real` 89).  The mandate specifies this one, because a
           top-N-by-MAGNITUDE pool inverts `nmae`'s sign.
  matched  the same reference pool with our call set truncated per row to the incumbent
           arm's own `n_pred`, so no arm can buy anything through call-set SIZE.

KNOWN-ANSWER ROWS (rule 10, abort):
  * `art_floor` scored beside 202's own `de202.rc_callset.npz` -- five members equal.
  * `art_ham`  must reproduce 208's rider-A row: `n_pred` 2629, `reach` 0.22340,
    `fid` 0.49429, `jac` 0.07302, `nmae` 0.96040.

BOARD CONVERSION.  `reach` x 1.290 (`rc_hybrid`'s own measured h1-test -> board factor for
this arm class), with 164's stale 0.4277 as the sensitivity.  `fid` as PRECISION (the
challenge call set is 2,599 >> `n_real` ~ 222, so the yield cap is saturated).  `jac` through
208's surviving S* convention with the board `NP` pinned at `cd4f034`'s own MEASURED 2,599
(210 falsifier 5), so board `jac` moves only through recall.  `nmae` h1-test raw / span,
scaled by the measured challenge/h1-test tilt ratio.

    uv run python six214.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
from debench import ALPHA, Surface                                       # noqa: E402

E185 = ROOT / "experiments/185-rank-controlled-emission"
E202 = ROOT / "experiments/202-callset-rider/work"
TAG = os.environ.get("TAG", "de214")
SEED_TAGS = ("", "_s2", "_s3")
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)   # 181's board 1/span
K_REACH, K_REACH_ALT = 1.290, 0.4277
BOARD_NP = 2599.0          # cd4f034's MEASURED challenge call set (210 falsifier 5)
BOARD_JAC = 0.04043        # rc_callset's board raw `jac`; cd4f034 is unscored, see RESULTS
KA208 = dict(n_pred=2629.0, reach=0.22340, fid=0.49429, jac=0.07302, nmae=0.96040)


class WaldSurface(Surface):
    """188's pool: the reference's top-`n_real` genes by `|lfc| * sqrt(control cpm)`."""

    def __init__(self, path, cpm):
        super().__init__(path)
        z = np.load(path, allow_pickle=True)
        base = (z["real_padj"] < ALPHA) & np.isfinite(z["real_padj"])
        key = np.abs(np.nan_to_num(z["real_lfc"])) * np.sqrt(np.maximum(cpm, 0.0))[None, :]
        sig = np.zeros(key.shape, bool)
        for i in range(key.shape[0]):
            k = int(base[i].sum())
            if k:
                sig[i, np.argpartition(-key[i], min(k, key.shape[1]) - 1)[:k]] = True
        self.real_sig = sig & np.isfinite(z["real_p"])
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.real_sig[p, j] = False


def truncated_call(s, sizes):
    call = np.zeros(s.pred_p.shape, bool)
    for p in range(call.shape[0]):
        k = min(int(sizes[p]), call.shape[1])
        if k <= 0:
            continue
        pv = np.where(np.isfinite(s.pred_p[p]), s.pred_p[p], np.inf)
        call[p, np.argsort(pv, kind="stable")[:k]] = True
    return call & s.adj


def recall_of(s, call=None):
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        b = s.pred_sig[p] if call is None else call[p]
        out.append(int((s.real_sig[p] & b).sum()) / nr)
    return float(np.mean(out))


def row(s, call=None):
    nr, npd = s.counts()
    sel = call if call is not None else (s.pred_sig & s.adj)
    prec = []
    for p in range(len(s.perts)):
        k = np.flatnonzero(sel[p])
        if k.size:
            prec.append(float((np.sign(s.pred_lfc[p, k]) ==
                               np.sign(s.real_lfc[p, k])).mean()))
    return dict(n_pred=float(np.median(sel.sum(1) if call is not None else npd)),
                reach=s.reach(call=call), fid=s.fid(call=call), jac=s.jac(call=call),
                nmae=s.nmae(), prec=float(np.mean(prec)), recall=recall_of(s, call=call))


def jac_model(R, NP, NR):
    hit = R * NR
    return hit / (NP + NR - hit)


def main() -> None:
    meta = np.load(E185 / "work/ctrl_meta.npz", allow_pickle=True)
    cpm = meta["cpm"][meta["keep"]]
    out, res, fails = {}, {}, []

    # ---------------------------------------------------------------- known-answer rows
    print("KNOWN-ANSWER ROWS (rule 10)\n")
    ka_path = WORK / f"{TAG}.art_floor.npz"
    s_ka, s_me = Surface(E202 / "de202.rc_callset.npz"), Surface(ka_path)
    r_ka, r_me = row(s_ka), row(s_me)
    print("  KA1 — `art_floor` IS `rc_callset`, scored beside 202's own table")
    for k in ("n_pred", "reach", "fid", "jac", "nmae"):
        ok = abs(r_ka[k] - r_me[k]) <= (0.5 if k == "n_pred" else 1e-12)
        print(f"       {k:<7} 202 {r_ka[k]:>10.6f}   here {r_me[k]:>10.6f}  "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"KA1.{k}")
    print("\n  KA2 — `art_ham` must reproduce 208's rider-A row on the SAME emitter")
    r_ham = row(Surface(WORK / f"{TAG}.art_ham.npz"))
    for k, want in KA208.items():
        tol = 1.0 if k == "n_pred" else 5e-5
        ok = abs(r_ham[k] - want) <= tol
        print(f"       {k:<7} 208 {want:>10.5f}   here {r_ham[k]:>10.5f}  "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"KA2.{k}")
    if fails and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {fails}")
    print("\n  ALL DE KNOWN-ANSWER ROWS PASS\n")
    out["ka"] = dict(rc_callset=r_ka, art_floor=r_me, art_ham=r_ham, ka208=KA208)

    # ---------------------------------------------------------------- the arm grid
    GEOM_ARMS = ("reb_pm50_floor", "reb_pm00_floor", "reb_pm50_ham", "reb_pm00_ham")
    LABEL = {"reb_pm50_floor": "incumbent  (PM 0.50, floor)",
             "reb_pm00_floor": "A   drop PerturbMatch",
             "reb_pm50_ham": "B   Hamilton rounding",
             "reb_pm00_ham": "A∘B both"}
    BASE = "reb_pm50_floor"

    NR = None
    for pool in ("native", "wald89", "matched"):
        print(f"=== POOL: {pool} — h1-test, {len(SEED_TAGS)} matched seeds ===")
        print(f"{'arm':<30}{'seed':>6}{'n_pred':>8}{'recall':>9}{'reach':>9}{'fid':>9}"
              f"{'jac':>9}{'nmae':>9}")
        per = {a: [] for a in GEOM_ARMS}
        for st in SEED_TAGS:
            sizes = None
            for a in GEOM_ARMS:
                p = WORK / f"{TAG}.{a}{st}.npz"
                s = WaldSurface(p, cpm) if pool == "wald89" else Surface(p)
                call = None
                if pool == "matched":
                    if sizes is None:
                        b = Surface(WORK / f"{TAG}.{BASE}{st}.npz")
                        sizes = (b.pred_sig & b.adj).sum(1)
                    call = truncated_call(s, sizes)
                r = row(s, call=call)
                per[a].append(r)
                print(f"{LABEL[a]:<30}{(st or '_s1'):>6}{r['n_pred']:>8.0f}"
                      f"{r['recall']:>9.5f}{r['reach']:>9.5f}{r['fid']:>9.5f}"
                      f"{r['jac']:>9.5f}{r['nmae']:>9.5f}")
        print()
        # seed-averaged rows and the paired per-seed deltas (the estimand)
        mean = {a: {k: float(np.mean([r[k] for r in per[a]]))
                    for k in per[a][0]} for a in GEOM_ARMS}
        sem = {}
        for a in GEOM_ARMS:
            sem[a] = {}
            for k in per[a][0]:
                d = np.asarray([per[a][i][k] - per[BASE][i][k] for i in range(len(SEED_TAGS))])
                sem[a][k] = float(d.std(ddof=1) / np.sqrt(d.size)) if d.size > 1 else 0.0
        if NR is None:
            NR = float(brentq(lambda nr: jac_model(mean[BASE]["recall"], BOARD_NP, nr)
                              - BOARD_JAC, 1.0, 50000.0, xtol=1e-9))
            print(f"  jac model calibration: NR = {NR:.1f} solved so the mean-field Jaccard "
                  f"reproduces the incumbent's own board raw `jac` {BOARD_JAC} at "
                  f"NP = {BOARD_NP:,.0f}\n")
        print(f"{'arm':<30}{'d reach':>11}{'sem':>9}{'d fid':>11}{'d jac':>11}"
              f"{'d nmae':>11}{'d recall':>11}")
        for a in GEOM_ARMS:
            print(f"{LABEL[a]:<30}"
                  f"{mean[a]['reach']-mean[BASE]['reach']:>+11.5f}{sem[a]['reach']:>9.5f}"
                  f"{mean[a]['fid']-mean[BASE]['fid']:>+11.5f}"
                  f"{mean[a]['jac']-mean[BASE]['jac']:>+11.5f}"
                  f"{mean[a]['nmae']-mean[BASE]['nmae']:>+11.5f}"
                  f"{mean[a]['recall']-mean[BASE]['recall']:>+11.5f}")
        print()
        print(f"{'arm':<30}{'reach':>11}{'fid':>11}{'jac':>11}{'nmae':>11}{'SUM':>12}"
              f"{'SUM@0.4277':>12}")
        b = mean[BASE]
        for a in GEOM_ARMS:
            r = mean[a]
            d_reach = (r["reach"] - b["reach"]) * K_REACH * INV["reach"] / 6
            d_alt = (r["reach"] - b["reach"]) * K_REACH_ALT * INV["reach"] / 6
            d_fid = (r["fid"] - b["fid"]) * INV["fid"] / 6
            d_jac = (jac_model(r["recall"], BOARD_NP, NR)
                     - jac_model(b["recall"], BOARD_NP, NR)) * INV["jac"] / 6
            d_nmae = -(r["nmae"] - b["nmae"]) * INV["nmae"] / 6
            tot = d_reach + d_fid + d_jac + d_nmae
            res[f"{pool}|{a}"] = dict(reach=d_reach, reach_alt=d_alt, fid=d_fid, jac=d_jac,
                                      nmae=d_nmae, sum=tot,
                                      sum_alt=d_alt + d_fid + d_jac + d_nmae,
                                      raw=r, sem=sem[a])
            print(f"{LABEL[a]:<30}{d_reach:>+11.6f}{d_fid:>+11.6f}{d_jac:>+11.6f}"
                  f"{d_nmae:>+11.6f}{tot:>+12.6f}{d_alt+d_fid+d_jac+d_nmae:>+12.6f}")
        A = res[f"{pool}|reb_pm00_floor"]["sum"]
        B = res[f"{pool}|reb_pm50_ham"]["sum"]
        AB = res[f"{pool}|reb_pm00_ham"]["sum"]
        print(f"\n  composition on the DE channel: A + B = {A+B:+.6f}   measured A∘B = "
              f"{AB:+.6f}   difference {AB-(A+B):+.6f}\n")

    # ---- the artifact family: 208's own rider row, for reference ------------------------
    print("=== the ARTIFACT family (a15's own delivered profile), for comparison with 208 ===")
    print(f"{'arm':<18}{'n_pred':>8}{'reach':>9}{'fid':>9}{'jac':>9}{'nmae':>9}")
    for a in ("art_floor", "art_ham"):
        r = row(Surface(WORK / f"{TAG}.{a}.npz"))
        print(f"{a:<18}{r['n_pred']:>8.0f}{r['reach']:>9.5f}{r['fid']:>9.5f}"
              f"{r['jac']:>9.5f}{r['nmae']:>9.5f}")
        out.setdefault("artifact", {})[a] = r

    out["board"] = res
    out["NR"] = NR
    out["board_np"] = BOARD_NP
    json.dump(out, open(WORK / "six214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'six214.json'}")


if __name__ == "__main__":
    main()
