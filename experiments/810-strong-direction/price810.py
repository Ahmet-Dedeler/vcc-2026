"""810 -- pricing of PROPOSAL.md's candidates on one bulk fold through 801's engine (EC2 only; written before any of these
numbers exist; runs only after the judge's GO). STANDARD FORMAT out (stats800), read by check800 aggregate.

Arms (the incumbent's key w and per-target q are kept everywhere, so P2g terciles, the render and the ballast design differ
only through the rows' direction; every changed row keeps |fin_i|). PROPOSAL-2.md prices ONE candidate (sgu1); the rest are
comparators or known answers, not tries of the family:
  inc     October's incumbent on L folds (pct+P2g); on ABC folds the board's ALLSW bulk (cap6+P2g(abc)), as ka804.
  707     ka804's 707 row rebuilt with ka804's own functions (top z_src tertile over targets with a donor row, fin blended
          toward the fold sources' unit mean U at beta 0.25, non-panel genes). KA: == 804's work/run/bulk/<fold>/707.npz.
  sgu1    THE CANDIDATE: on the same top tertile, the composition's donor BASE (decomp810's part `base` = s W) is replaced
          by the unwhitened normrel'd donor mean `pre` at the base's norm, base' = |base| pre/|pre| (beta 1 of
          base' = |base| normalise((1-b) base/|base| + b pre/|pre|)); fin' = base' + cd4 + h1 + part (the residual parts
          untouched); fin' rescaled to |fin_i|.
  sgu1all the CLOSED global form (212 / 097's whitening relaxation at the base): the same swap on EVERY target with a
          donor row (U != 0). A known-answer comparator: expected net <= ~0 on BULK-ABC (PROPOSAL-2.md section 4).
  sgu0    beta 0 on the top tertile: a known answer (== inc), not saved.
Known answers (abort unless |d| < 1e-12): the parts (810parts/<fold>-parts.npz, decomp810) sum to F.fin; sgu0 equals inc
(x, a, pds, pdsl, K); inc and 707 equal 804's stored arrays on this fold.

    KA804_OUT=/work/ka804 P810_OUT=/work/p810 uv run python experiments/810-strong-direction/price810.py ABC/orion_hct116
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/804-known-answers"))
import numpy as np  # noqa: E402

import ka804 as K  # noqa: E402  (engine801, zsrc_fold, donor_unit_mean, p2g_for; guard801 via engine801)

E = K.E
OUT = Path(os.environ.get("P810_OUT", "/work/p810"))
PARTS = Path(os.environ.get("PARTS810", "/work/parts810"))
REF804 = ROOT / "experiments/804-known-answers/work/run/bulk"
TOL = 1e-12
TOP = 2


def tiers_have(z, have):
    tier = np.full(z.size, -1)
    r = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
    tier[np.flatnonzero(have)] = np.minimum((r * 3) // int(have.sum()), 2)
    return tier


def sgu_rows(P, fin, sel, beta, nonpanel):
    out = fin.copy()
    for i in np.flatnonzero(sel):
        b, p = P["base"][i, nonpanel], P["pre"][i, nonpanel]
        nb, npre = np.linalg.norm(b), np.linalg.norm(p)
        if nb == 0 or npre == 0:
            continue
        m = (1 - beta) * b / nb + beta * p / npre
        b2 = nb * m / np.linalg.norm(m)
        row = b2 + P["cd4"][i, nonpanel] + P["h1"][i, nonpanel] + P["part"][i, nonpanel]
        out[i, nonpanel] = row * (np.linalg.norm(fin[i]) / np.linalg.norm(row))
    return out


def main(fold_id: str):
    t0 = time.time()
    stem = fold_id.replace("/", "__")
    F = E.load(fold_id, check=True)
    abc = F.panel == "ABC"
    P = dict(np.load(PARTS / f"{stem}-parts.npz"))
    assert [str(t) for t in P["targets"]] == list(F.targets)
    fin = np.asarray(F.fin, np.float64)
    ka = {"parts sum == F.fin": float(np.max(np.abs(P["base"] + P["cd4"] + P["h1"] + P["part"] - fin)))}
    inc0 = E.named(F, "cap6+P2g(abc)" if abc else "pct+P2g")
    w, q = inc0.w, inc0.q
    z = K.zsrc_fold(F)
    U, have = K.donor_unit_mean(F)
    tier = tiers_have(z, have)
    sel = tier == TOP
    npn = F.nonpanel
    fin707 = fin.copy()
    for i in np.flatnonzero(sel):
        nb, nu = np.linalg.norm(fin[i, npn]), np.linalg.norm(U[i, npn])
        if nb > 0 and nu > 0:
            m = (1 - K.BETA) * fin[i, npn] / nb + K.BETA * U[i, npn] / nu
            fin707[i, npn] = nb * m / np.linalg.norm(m)
    arms = {"inc": E.Arm(w=w, q=q, name="inc"),
            "707": E.Arm(w=w, q=q, fin=fin707, name="707"),
            "sgu0": E.Arm(w=w, q=q, fin=sgu_rows(P, fin, sel, 0.0, npn), name="sgu0"),
            "sgu1": E.Arm(w=w, q=q, fin=sgu_rows(P, fin, sel, 1.0, npn), name="sgu1"),
            "sgu1all": E.Arm(w=w, q=q, fin=sgu_rows(P, fin, have, 1.0, npn), name="sgu1all")}
    res, arrs, lines = {}, {}, []
    for nm, arm in arms.items():
        s, arr = E.price(F, arm, base=None if nm == "inc" else arms["inc"])
        arr["targets"] = np.asarray(F.targets)
        arrs[nm] = arr
        res[nm] = dict(cos_agg=s["cos_agg"], pds=s["pds"], pds_lam05=s["pds_lam0.5"], K_mean=s["K_mean"], vs_inc=s.get("vs_base"))
        lines.append(f"[{fold_id}] {nm:6s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.5f}" +
              (f" r {100 * s['vs_base']['r']:+.3f}% dpds {s['vs_base']['dpds']:+.5f} dK {s['vs_base']['dK_rel']:+.4f}" if nm != "inc" else ""),
              )
    # known answers
    for k in ("x", "a", "pds", "pdsl", "K"):   # RULING-KA.md (08:55Z): x, a (unnormalised sums) read RELATIVE to max|ref|
        d = float(np.max(np.abs(np.asarray(arrs["sgu0"][k], float) - np.asarray(arrs["inc"][k], float))))
        if k in ("x", "a"):
            d /= max(float(np.max(np.abs(np.asarray(arrs["inc"][k], float)))), 1e-300)
        ka[f"sgu0 == inc: {k}"] = d
    ref = REF804 / stem
    for nm in ("inc", "707"):
        if (ref / f"{nm}.npz").exists():
            z0 = np.load(ref / f"{nm}.npz")
            for k in ("x", "a", "pds", "pdsl"):
                ka[f"{nm} == 804 {k}"] = float(np.max(np.abs(np.asarray(arrs[nm][k], float) - np.asarray(z0[k], float))))
    for k, d in ka.items():
        assert d < TOL, f"KA FAILED {fold_id} {k}: {d}"
    print("\n".join(lines), flush=True)
    d = OUT / "bulk" / stem
    d.mkdir(parents=True, exist_ok=True)
    for nm, arr in arrs.items():
        if nm == "sgu0":
            continue
        keep = {k: v for k, v in arr.items() if not k.startswith("base_")}
        np.savez_compressed(d / f"{nm}.npz", **keep)
    rec = dict(fold=fold_id, ka=ka, n=F.n, n_top=int(sel.sum()), n_have=int(have.sum()), arms=res, seconds=time.time() - t0)
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    (OUT / "meta" / f"{stem}.json").write_text(json.dumps(rec, indent=1, default=float))
    print(f"[{fold_id}] done, KA worst {max(ka.values()):.1e} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
