"""220 Phase 2 — the six members for every DELIVERED ordering key, and rule 13's gate.

The arm class is `emission-law, profile-PRESERVING`: `rc_hybrid_block` asserts
`block.sum(0) == floor(x)` for every block and `x` is the artifact's own intended profile,
which no key touches.  `emit220.py` prints the per-gene panel-total digest for all 85 blocks
against `absH` (= `rc_callset`); the class test is that measurement, not this docstring.
Rule 13's table therefore charges NO `reach` toll and pins `pds` and `mse` at exactly zero
(board-verified to sixteen significant figures on `rc_hybrid` and again on `rc_callset`).
`nmae` is NOT pinned by the digest -- it reads `pred_lfc`, which depends on where inside the
block a gene's fixed total is placed -- so it is measured rather than assumed.

KNOWN-ANSWER ROWS (rule 10, abort):
  * arm `absH` scored beside 202's own `de202.rc_callset.npz`: five members equal to 1e-12.
  * `reach(production)` = 0.223341, `fid` 0.494362, `jac` 0.072975, `nmae` 0.966702.

BOARD CONVERSION — 185's board point, not 164's stale 0.4277:
    reach x 1.290 ; fid as precision ; jac through 208's S* with board NP pinned at 2,599 ;
    nmae raw / span.  1/span from 181: reach 1.12927, fid 3.32482, jac 2.72081, nmae 1.64787.

    uv run python six220.py
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
from debench import ALPHA, Surface                                        # noqa: E402

E185 = ROOT / "experiments/185-rank-controlled-emission"
E202 = ROOT / "experiments/202-callset-rider/work"
TAG = os.environ.get("TAG", "de220")
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)
K_REACH, K_REACH_ALT = 1.290, 0.4277
BOARD_NP = 2599.0
BOARD_JAC = 0.04043
GATE_SHAPE = 0.0019
KA = dict(reach=0.223341, fid=0.494362, jac=0.072975, nmae=0.966702)

THREE = ("absH", "k189_l050", "k190_km25")
SEED_TAGS = ("", "_s2", "_s3")


class WaldSurface(Surface):
    """188's significance-selected pool (median n_real 89), 214's `wald89`."""

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


def recall_of(s, call=None):
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        b = s.pred_sig[p] if call is None else call[p]
        out.append(int((s.real_sig[p] & b).sum()) / nr)
    return float(np.mean(out))


def row(s):
    sel = s.pred_sig & s.adj
    prec = []
    for p in range(len(s.perts)):
        k = np.flatnonzero(sel[p])
        if k.size:
            prec.append(float((np.sign(s.pred_lfc[p, k]) ==
                               np.sign(s.real_lfc[p, k])).mean()))
    return dict(n_pred=float(np.median(sel.sum(1))), reach=s.reach(), fid=s.fid(),
                jac=s.jac(), nmae=s.nmae(), prec=float(np.mean(prec)),
                recall=recall_of(s))


def jac_model(R, NP, NR):
    hit = R * NR
    return hit / (NP + NR - hit)


def main() -> None:
    meta = np.load(E185 / "work/ctrl_meta.npz", allow_pickle=True)
    cpm = meta["cpm"][meta["keep"]]
    diag = json.load(open(WORK / f"{TAG}_diag.json"))
    arms = diag["arms"]
    out, fails = {}, []

    print("KNOWN-ANSWER ROWS (rule 10)\n")
    s_ka = Surface(str(E202 / "de202.rc_callset.npz"))
    s_me = Surface(str(WORK / f"{TAG}.absH.npz"))
    r_ka, r_me = row(s_ka), row(s_me)
    print("  KA1 — arm `absH` IS `rc_callset`, scored beside 202's own table")
    for k in ("n_pred", "reach", "fid", "jac", "nmae"):
        ok = abs(r_ka[k] - r_me[k]) <= (0.5 if k == "n_pred" else 1e-12)
        print(f"       {k:<7} 202 {r_ka[k]:>11.6f}   here {r_me[k]:>11.6f}  "
              f"{'PASS' if ok else 'FAIL'}")
        fails += [] if ok else [f"KA1.{k}"]
    print("\n  KA2 — the published absolute rows")
    for k, want in KA.items():
        ok = abs(r_me[k] - want) < 5e-6
        print(f"       {k:<7} want {want:>11.6f}   here {r_me[k]:>11.6f}  "
              f"{'PASS' if ok else 'FAIL'}")
        fails += [] if ok else [f"KA2.{k}"]
    if fails and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {fails}")
    print("\n  ALL KNOWN-ANSWER ROWS PASS\n")
    out["ka"] = dict(rc_callset=r_ka, absH=r_me)

    # ---- the digest class test, read back from the emitter's own record ---------------
    dg = diag["digest"]
    base_dg = dg["absH"]
    print("PER-GENE PANEL-TOTAL DIGEST (rule 13's class test — measured, not asserted)")
    cls = {}
    for a in arms:
        same = sum(x == y for x, y in zip(dg[a], base_dg))
        cls[a] = same
    n_blocks = len(base_dg)
    allsame = all(v == n_blocks for v in cls.values())
    print(f"  every one of the {len(arms)} arms is identical to `absH` on "
          f"{min(cls.values())}-{max(cls.values())} of {n_blocks} blocks   "
          f"{'ALL PROFILE-PRESERVING' if allsame else 'SOME PROFILE-CHANGING'}")
    print(f"  => rule 13 class = emission-law, profile-PRESERVING; `reach` toll = 0; "
          f"`pds` and `mse` pinned at 0.000000")
    out["digest_class"] = dict(all_preserving=bool(allsame), per_arm=cls,
                               n_blocks=n_blocks)
    if not allsame:
        fails.append("digest")

    # ---- the arm grid -----------------------------------------------------------------
    res = {}
    NR = None
    for pool in ("native", "wald89"):
        print(f"\n=== POOL: {pool} — h1-test, DELIVERED through `rc_callset` ===")
        rows = {}
        for a in arms:
            p = WORK / f"{TAG}.{a}.npz"
            s = WaldSurface(str(p), cpm) if pool == "wald89" else Surface(str(p))
            rows[a] = row(s)
        if NR is None:
            NR = float(brentq(lambda nr: jac_model(rows["absH"]["recall"], BOARD_NP, nr)
                              - BOARD_JAC, 1.0, 50000.0, xtol=1e-9))
            print(f"  jac model calibration: NR = {NR:.1f} so the mean-field Jaccard "
                  f"reproduces `rc_callset`'s own board raw jac {BOARD_JAC}\n")
        # seed sem for the three-seed arms
        sem = {}
        for k in THREE:
            d = np.asarray([rows[k if s == 0 else f'{k}_s{s+1}']["reach"]
                            - rows["absH" if s == 0 else f"absH_s{s+1}"]["reach"]
                            for s in range(3)])
            sem[k] = dict(mean=float(d.mean()),
                          sem=float(d.std(ddof=1) / np.sqrt(3)), per_seed=list(map(float, d)))
        print(f"{'arm':<16}{'n_pred':>8}{'recall':>9}{'reach':>10}{'d reach':>10}"
              f"{'fid':>9}{'jac':>9}{'nmae':>10}{'d nmae':>10}")
        b = rows["absH"]
        for a in arms:
            r = rows[a]
            print(f"{a:<16}{r['n_pred']:>8.0f}{r['recall']:>9.5f}{r['reach']:>10.5f}"
                  f"{r['reach']-b['reach']:>+10.5f}{r['fid']:>9.5f}{r['jac']:>9.5f}"
                  f"{r['nmae']:>10.6f}{r['nmae']-b['nmae']:>+10.6f}")
        print(f"\n  three-seed paired deltas (the SEED axis, 214's 0.00009 instrument):")
        for k in THREE:
            print(f"    {k:<14} mean {sem[k]['mean']:+.6f}  seed sem "
                  f"{sem[k]['sem']:.6f}   per seed " +
                  " ".join(f"{x:+.6f}" for x in sem[k]["per_seed"]))

        print(f"\n{'arm':<16}{'reach':>11}{'fid':>11}{'jac':>11}{'nmae':>11}"
              f"{'SUM':>12}{'gate':>9}{'verdict':>10}")
        for a in arms:
            r = rows[a]
            d_reach = (r["reach"] - b["reach"]) * K_REACH * INV["reach"] / 6
            d_alt = (r["reach"] - b["reach"]) * K_REACH_ALT * INV["reach"] / 6
            d_fid = (r["fid"] - b["fid"]) * INV["fid"] / 6
            d_jac = (jac_model(r["recall"], BOARD_NP, NR)
                     - jac_model(b["recall"], BOARD_NP, NR)) * INV["jac"] / 6
            d_nmae = -(r["nmae"] - b["nmae"]) * INV["nmae"] / 6
            tot = d_reach + d_fid + d_jac + d_nmae
            res[f"{pool}|{a}"] = dict(reach=d_reach, reach_alt=d_alt, fid=d_fid, jac=d_jac,
                                      nmae=d_nmae, sum=tot, raw=r)
            print(f"{a:<16}{d_reach:>+11.6f}{d_fid:>+11.6f}{d_jac:>+11.6f}"
                  f"{d_nmae:>+11.6f}{tot:>+12.6f}{GATE_SHAPE:>9.4f}"
                  f"{'SHIP' if tot > GATE_SHAPE else 'VETO':>10}")
        out[f"sem_{pool}"] = sem

    out["board"] = res
    out["NR"] = NR
    json.dump(out, open(WORK / "six220.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'six220.json'}")
    if fails:
        raise SystemExit(f"FAILED: {fails}")


if __name__ == "__main__":
    main()
