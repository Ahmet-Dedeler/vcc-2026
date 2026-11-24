"""196 — read the four DE members off the h1-test depth tables, three ways.

Three pools / call sets, because the depth knob moves DE POWER and profile FIDELITY at once
and the two have to be separated:

  native      the reference's own BH<0.05 pool and OUR own BH<0.05 call set.  This is what a
              submission at that depth would actually score.  Its `m = 1` rows are the
              known-answer gate (185/191: a15 reach 0.18283 / n_pred 6034 / jac 0.06824 /
              nmae 0.96659; rc_hybrid 0.21681 / 6043 / 0.06827 / 0.96670).
  matched     the same reference pool, but OUR call set truncated per row to the `m = 1`
              arm's own `n_pred`, by our own p-value.  Depth then buys only profile
              fidelity: the call-set size is pinned, so `fid`'s yield term and `jac`'s union
              cannot move for that reason alone.
  wald89      188's significance-selected pool construction -- the reference's top-`n_real`
              genes by `|lfc| * sqrt(cpm)` instead of its BH set -- at the reference's own
              BH depth (h1-test median `n_real` = 90).  188 measured that a top-N-by-
              magnitude pool inverts `nmae`'s sign, and this is the construction that does
              not.

    uv run python six196.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
from debench import ALPHA, Surface                                     # noqa: E402

TAGS = tuple(os.environ.get("TAGS", "de196,rc196").split(","))
INV = dict(reach=1.12927, fid=3.32482, jac=2.72081, nmae=1.64787)      # 181's board 1/span
KA = {
    "quasi_m1": dict(reach=0.18283, n_pred=6034, jac=0.06824, nmae=0.96659),
    "rc_hybrid_m1": dict(reach=0.21681, n_pred=6043, jac=0.06827, nmae=0.96670),
}


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


ROW_ATTRS = ("perts", "real_p", "real_padj", "real_lfc", "pred_p", "pred_padj",
             "pred_lfc", "pred_z", "absH", "adj", "real_sig", "pred_sig", "self_idx")


def subset_rows(s, ok):
    for attr in ROW_ATTRS:
        if hasattr(s, attr):
            setattr(s, attr, getattr(s, attr)[ok])
    return s


def truncated_call(s, sizes):
    """Our call set cut to `sizes[p]` genes per row by our own predicted p-value."""
    call = np.zeros(s.pred_p.shape, bool)
    for p in range(call.shape[0]):
        k = int(sizes[p])
        if k <= 0:
            continue
        k = min(k, call.shape[1])
        pv = np.where(np.isfinite(s.pred_p[p]), s.pred_p[p], np.inf)
        call[p, np.argsort(pv, kind="stable")[:k]] = True
    return call & s.adj


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
                nmae=s.nmae(), prec=float(np.mean(prec)))


def main() -> None:
    path = {}
    for tg in TAGS:
        for f in sorted(WORK.glob(f"{tg}.*.npz")):
            path[f.name[len(tg) + 1:-4]] = f
    arms = list(path)
    order = []
    for a in ("quasi", "rc_hybrid"):
        for m in (1, 2, 4, 8, 16, 32, 64):
            if f"{a}_m{m}" in arms:
                order.append(f"{a}_m{m}")
    arms = order + [a for a in arms if a not in order]
    meta = np.load(WORK / "ctrl_meta.npz", allow_pickle=True) if (WORK / "ctrl_meta.npz").exists() \
        else np.load(ROOT / "experiments/185-rank-controlled-emission/work/ctrl_meta.npz",
                     allow_pickle=True)
    cpm = meta["cpm"][meta["keep"]]
    diag = {"diag": {}, "failed": {}}
    for tg in TAGS:
        j = WORK / f"{tg}_diag.json"
        if j.exists():
            d = json.load(open(j))
            diag["diag"].update(d["diag"])
            diag["failed"].update(d.get("failed", {}))

    S = {a: Surface(path[a]) for a in arms}
    print(f"h1-test depth sweep, {len(S[arms[0]].perts)} perturbations, "
          f"{S[arms[0]].genes.size} gated genes; median n_real "
          f"{np.median(S[arms[0]].real_sig.sum(1)):.0f}\n")

    # The known-answer values are 185/191's, measured on the FULL 85-row cohort, so the gate
    # runs BEFORE the paired-cohort drop below.
    print("KNOWN-ANSWER ROWS (rule 10) -- full 85-row cohort, as 185/191 measured them")
    bad = []
    for a, want in KA.items():
        if a not in S:
            continue
        got = row(S[a])
        for k, w in want.items():
            g = got[k]
            tol = 0.5 if k == "n_pred" else 2e-5
            hit = abs(g - w) <= tol
            print(f"  {a:<14} {k:<7} {g:>10.5f}  expect {w:<10.5f} {'PASS' if hit else 'FAIL'}")
            if not hit:
                bad.append(f"{a}.{k}")
    if bad and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {bad}")
    print("  all known-answer rows PASS\n" if not bad else "")

    # PAIRED COHORT.  `rc_hybrid`'s greedy placement can fail to weld a row at depth; that
    # row's DE table is all-NaN.  A failed row must be dropped from EVERY arm, or the deep
    # rc rows are compared against a larger cohort than they were emitted on.  Detected from
    # the tables themselves, so it does not depend on the (possibly missing) diag json.
    T = len(S[arms[0]].perts)
    ok = np.ones(T, bool)
    for a in arms:
        ok &= np.isfinite(S[a].pred_p).any(1)
    if not ok.all():
        lost = [str(x) for x in S[arms[0]].perts[~ok]]
        print(f"PAIRED COHORT: dropping {int((~ok).sum())}/{T} perturbations that failed to "
              f"emit on at least one arm: {lost}\n")
        for a in arms:
            subset_rows(S[a], ok)

    res = {}
    for pool, mk in (("native", lambda a: S[a]),
                     ("wald89", lambda a: subset_rows(WaldSurface(path[a], cpm), ok))):
        print(f"=== POOL: {pool} ===")
        print(f"{'arm':<15}{'n_pred':>8}{'reach':>9}{'fid':>9}{'jac':>9}{'nmae':>9}"
              f"{'prec':>9}{'d overall':>11}")
        base = None
        for a in arms:
            s = mk(a)
            r = row(s)
            emitter = a.split("_m")[0]
            if base is None or a.endswith("_m1"):
                if a == f"{emitter}_m1":
                    res[(pool, emitter, "base")] = r
            b = res.get((pool, emitter, "base"), r)
            d = ((r["reach"] - b["reach"]) * INV["reach"] + (r["fid"] - b["fid"]) * INV["fid"]
                 + (r["jac"] - b["jac"]) * INV["jac"] - (r["nmae"] - b["nmae"]) * INV["nmae"]) / 6
            res[(pool, a, "native")] = dict(r, d_overall=d)
            print(f"{a:<15}{r['n_pred']:>8.0f}{r['reach']:>9.5f}{r['fid']:>9.5f}"
                  f"{r['jac']:>9.5f}{r['nmae']:>9.5f}{r['prec']:>9.5f}{d:>+11.5f}")
        print()

    print("=== POOL: native, CALL SET MATCHED to the m=1 arm's own per-row n_pred ===")
    print(f"{'arm':<15}{'n_pred':>8}{'reach':>9}{'fid':>9}{'jac':>9}{'nmae':>9}"
          f"{'prec':>9}{'d overall':>11}")
    for emitter in ("quasi", "rc_hybrid"):
        b1 = f"{emitter}_m1"
        if b1 not in S:
            continue
        sizes = (S[b1].pred_sig & S[b1].adj).sum(1)
        base = None
        for a in [x for x in arms if x.startswith(emitter + "_m")]:
            s = S[a]
            r = row(s, call=truncated_call(s, sizes))
            if base is None:
                base = r
            d = ((r["reach"] - base["reach"]) * INV["reach"]
                 + (r["fid"] - base["fid"]) * INV["fid"]
                 + (r["jac"] - base["jac"]) * INV["jac"]
                 - (r["nmae"] - base["nmae"]) * INV["nmae"]) / 6
            res[("matched", a, "native")] = dict(r, d_overall=d)
            print(f"{a:<15}{r['n_pred']:>8.0f}{r['reach']:>9.5f}{r['fid']:>9.5f}"
                  f"{r['jac']:>9.5f}{r['nmae']:>9.5f}{r['prec']:>9.5f}{d:>+11.5f}")
    print("\nEMISSION DIAGNOSTICS (median over perturbations)")
    print(f"{'arm':<15}{'emit L':>10}{'panel S':>13}{'nnz/block':>12}{'nnz ratio':>11}"
          f"{'max cell':>11}{'cell sd':>10}")
    base_nnz = {}
    for a in arms:
        d = diag["diag"].get(a)
        if not d:
            print(f"{a:<15}  (no diag on disk -- its sweep crashed after the tables were "
                  f"written; nnz for this emitter is measured in work/rcdepth.log)")
            continue
        d = [x for x in d if not x.get("failed")]
        med = {k: float(np.median([x[k] for x in d])) for k in
               ("emit_L", "panel_S", "nnz", "max_cell", "cell_total_sd")}
        e = a.split("_m")[0]
        base_nnz.setdefault(e, med["nnz"])
        print(f"{a:<15}{med['emit_L']:>10,.0f}{med['panel_S']:>13,.0f}{med['nnz']:>12,.0f}"
              f"{med['nnz']/base_nnz[e]:>11.4f}{med['max_cell']:>11,.0f}"
              f"{med['cell_total_sd']:>10,.1f}")
        res[("diag", a, "native")] = med
    json.dump({f"{k[0]}|{k[1]}": v for k, v in res.items()},
              open(WORK / "six196.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
