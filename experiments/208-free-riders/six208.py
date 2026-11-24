"""208 — the four DE members for the 2x2 rider grid, marginal to `rc_callset`, on h1-test.

Three pools, exactly 196's constructions, because the two riders move different things:

  native      the reference's own BH<0.05 pool and OUR own BH<0.05 call set -- what a
              submission would actually score.  `floor_tier` IS `rc_callset` and its row is
              the known-answer gate against 202's own emitted table.
  wald89      188's SIGNIFICANCE-selected pool: the reference's top-`n_real` genes by
              `|lfc| * sqrt(control cpm)` at the reference's own BH depth (median 89).  The
              mandate specifies this one, because a top-N-by-MAGNITUDE pool inverts `nmae`.
  matched     the same reference pool with our call set truncated per row to `floor_tier`'s
              own `n_pred`, so neither rider can buy anything through call-set SIZE.

BOARD CONVERSION.  `reach` at the 1.290 factor `rc_hybrid` measured for this arm class (and
at 0.4277 as the rule-9 sensitivity).  `fid` as PRECISION (n_pred 2,599 >> the challenge's
n_real ~ 222, so the yield cap is saturated and `fid` IS precision).  `jac` through the
convention Phase 1 adjudicated -- the arm's OWN measured challenge `n_pred`, which for both
riders is `rc_callset`'s 2,599 UNCHANGED, because the call set is PRESCRIBED and neither
rider can move it.  So board `jac` moves only through RECALL.

    uv run python six208.py
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
from debench import ALPHA, Surface  # noqa: E402

E185 = ROOT / "experiments/185-rank-controlled-emission"
E202 = ROOT / "experiments/202-callset-rider/work"
TAG = os.environ.get("TAG", "de208")
ARMS = ("floor_tier", "ham_tier", "floor_flat", "ham_flat")
BASE = "floor_tier"
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)   # 181's board 1/span
K_REACH = 1.290                # rc_hybrid's own measured h1-test -> board factor
K_REACH_ALT = 0.4277           # 164's stale spread-policy factor, as the sensitivity
BOARD_NP = 2599.0              # rc_callset's MEASURED challenge call set (202 falsifier 2)
BOARD_JAC = 0.04043            # rc_callset's board raw `jac`


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
                nmae=s.nmae(), prec=float(np.mean(prec)),
                recall=recall_of(s, call=call))


def jac_model(R, NP, NR):
    hit = R * NR
    return hit / (NP + NR - hit)


def main() -> None:
    meta = np.load(E185 / "work/ctrl_meta.npz", allow_pickle=True)
    cpm = meta["cpm"][meta["keep"]]
    path = {a: WORK / f"{TAG}.{a}.npz" for a in ARMS}
    for a, p in path.items():
        if not p.exists():
            raise SystemExit(f"missing {p}")
    S = {a: Surface(path[a]) for a in ARMS}
    T = len(S[BASE].perts)
    print(f"h1-test, {T} perturbations, {S[BASE].genes.size} gated genes; median n_real "
          f"{np.median(S[BASE].real_sig.sum(1)):.0f}\n")

    # ---- known-answer gate: floor_tier IS rc_callset -----------------------------------
    ka = Surface(E202 / "de202.rc_callset.npz")
    r_ka, r_me = row(ka), row(S[BASE])
    print("KNOWN-ANSWER ROW (rule 10): `floor_tier` scored beside 202's own `rc_callset`")
    bad = []
    for k in ("n_pred", "reach", "fid", "jac", "nmae"):
        ok = abs(r_ka[k] - r_me[k]) <= (0.5 if k == "n_pred" else 1e-12)
        print(f"  {k:<7} 202 {r_ka[k]:>10.6f}   here {r_me[k]:>10.6f}  "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(k)
    if bad:
        raise SystemExit(f"ABORT (rule 10): {bad}")
    # a second, independent gate: 202's incumbent `rc_hybrid` row, as 196/191 recorded it
    rh = row(Surface(E202 / "de202.rc_hybrid.npz"))
    print(f"  rc_hybrid (context) reach {rh['reach']:.5f}  n_pred {rh['n_pred']:.0f}  "
          f"jac {rh['jac']:.5f}  nmae {rh['nmae']:.5f}   [191: 0.21681 / 6043 / 0.06827 / "
          f"0.96670]\n")

    res, out = {}, {}
    for pool in ("native", "wald89", "matched"):
        print(f"=== POOL: {pool} ===")
        print(f"{'arm':<12}{'n_pred':>8}{'recall':>9}{'reach':>9}{'fid':>9}{'jac':>9}"
              f"{'nmae':>9}{'d overall':>11}")
        sizes = None
        base_r = None
        for a in ARMS:
            s = Surface(path[a]) if pool != "wald89" else WaldSurface(path[a], cpm)
            call = None
            if pool == "matched":
                if sizes is None:
                    sizes = (S[BASE].pred_sig & S[BASE].adj).sum(1)
                call = truncated_call(s, sizes)
            r = row(s, call=call)
            if a == BASE:
                base_r = r
            d = ((r["reach"] - base_r["reach"]) * K_REACH * INV["reach"]
                 + (r["fid"] - base_r["fid"]) * INV["fid"]
                 + (r["jac"] - base_r["jac"]) * INV["jac"]
                 - (r["nmae"] - base_r["nmae"]) * INV["nmae"]) / 6
            res[(pool, a)] = dict(r, d_overall_h1units=d)
            print(f"{a:<12}{r['n_pred']:>8.0f}{r['recall']:>9.5f}{r['reach']:>9.5f}"
                  f"{r['fid']:>9.5f}{r['jac']:>9.5f}{r['nmae']:>9.5f}{d:>+11.5f}")
        print()

    # ---- the honest board conversion ---------------------------------------------------
    NR = float(brentq(lambda nr: jac_model(res[("native", BASE)]["recall"], BOARD_NP, nr)
                      - BOARD_JAC, 1.0, 50000.0, xtol=1e-9))
    print(f"=== BOARD CONVERSION (marginal to `rc_callset`) ===")
    print(f"  `reach`  x {K_REACH} (rc_hybrid's own measured factor) / span {1/INV['reach']:.5f}")
    print(f"  `fid`    precision, yield cap saturated at n_pred {BOARD_NP:,.0f} >> n_real")
    print(f"  `jac`    Phase 1's surviving convention: board NP = the arm's OWN measured "
          f"challenge n_pred,\n           which is {BOARD_NP:,.0f} for EVERY arm here (the call "
          f"set is prescribed and neither\n           rider can move it), so board `jac` moves "
          f"only through RECALL.  NR calibrated on\n           rc_callset's own board point: "
          f"NR = {NR:.1f}.")
    print(f"  `nmae`   h1-test raw, factor 1.0 (rule 13: 3/6 signs -- indicative, not priced)\n")
    for pool in ("native", "wald89", "matched"):
        print(f"--- pool {pool} ---")
        print(f"{'arm':<12}{'reach':>11}{'fid':>11}{'jac':>11}{'nmae':>11}{'SUM':>12}"
              f"{'SUM@0.4277':>12}")
        b = res[(pool, BASE)]
        for a in ARMS:
            r = res[(pool, a)]
            d_reach = (r["reach"] - b["reach"]) * K_REACH * INV["reach"] / 6
            d_reach_alt = (r["reach"] - b["reach"]) * K_REACH_ALT * INV["reach"] / 6
            d_fid = (r["fid"] - b["fid"]) * INV["fid"] / 6
            jb = jac_model(b["recall"], BOARD_NP, NR)
            ja = jac_model(r["recall"], BOARD_NP, NR)
            d_jac = (ja - jb) * INV["jac"] / 6
            d_nmae = -(r["nmae"] - b["nmae"]) * INV["nmae"] / 6
            tot = d_reach + d_fid + d_jac + d_nmae
            out[f"{pool}|{a}"] = dict(reach=d_reach, reach_alt=d_reach_alt, fid=d_fid,
                                      jac=d_jac, nmae=d_nmae, sum=tot,
                                      sum_alt=d_reach_alt + d_fid + d_jac + d_nmae,
                                      board_jac_pred=ja, recall=r["recall"])
            print(f"{a:<12}{d_reach:>+11.6f}{d_fid:>+11.6f}{d_jac:>+11.6f}"
                  f"{d_nmae:>+11.6f}{tot:>+12.6f}"
                  f"{d_reach_alt+d_fid+d_jac+d_nmae:>+12.6f}")
        print()

    json.dump(dict(rows={f"{k[0]}|{k[1]}": v for k, v in res.items()}, board=out,
                   NR=NR, board_np=BOARD_NP, k_reach=K_REACH, inv=INV),
              open(WORK / "six208.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
