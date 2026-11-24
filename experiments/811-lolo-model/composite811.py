"""811 composite arms for LAB-DEF-RUNBOOK Amendment 1 (2026-10-17): 811's c0_l25 priced as an INCREMENT on top of what would ship.

On any engine801 fold, with every construction taken from its owner, unchanged:
  B0 = inc + 705's theta 0.25 rider           (R0's bulk: inc = pct+P2g; w705 above the floor; P2g re-cut on w705)
  B1 = B0 + 707's swap                         (beta 0.25 on the top z_src tertile, toward the sources' unit mean)
  C0 = B0 + 811's c0_l25 blend                 (lambda 0.25, m811's DAY bank list, centring off, the held group out)
  C1 = B1 + 811's c0_l25 blend                 (swap first, then 811; the blend acts on B1's composition rows)
Keys stay those of B0 (w705, P2g on w705). Every arm is priced by engine801.price vs inc and saved with ka804.save.

Where each piece comes from (nothing re-derived):
  z_src, the sources' unit mean, P2g cuts, FLOOR / THETA / BETA: ka804 (zsrc_fold, donor_unit_mean, p2g_for, constants);
  the 705 rider is ka804.main's own expression (asserted present, character for character, in ka804.py's source at run time);
  the 707 swap loop and 811's blend are base811.blend (bit-identical to ka804's loop: 811's KA on 22 folds, |d| 0.0);
  the consensus U is m811.fold_U(F, centre=False) (m811.py byte-identical, sha256 e00192e6...).
Two conventions (--mode):
  day   : inc = pct+P2g and P2g cut by p2g_oct on EVERY fold, ABCp included (R0's bulk; LAB-DEF-RUNBOOK Amendment 1);
  ka804 : ka804's own per-panel choice (cap6+P2g(abc) / p2g_abc on the 062 and ABCp folds), used ONLY to reproduce its arrays.
On L / LABC / kit folds the two modes coincide (ka804 picks pct+P2g / p2g_oct there).
Known answers (--ka; every failure aborts with rc 1):
  KA-705 : ka804-mode B0 == ka804's stored 705 arrays;  KA-707 : ka804-mode B1 with the rider OFF == ka804's stored 707 arrays
           (on L1/k562_gwps from 804's run, and on an ABCp fold from 801's ka-abcp run; byte copies in work/refcomp/);
  KA-L0  : C0 at lambda 0 == B0 and C1 at lambda 0 == B1, bit-exact (the composite path with lambda 0 returns its base rows;
           the blend arithmetic at lambda 0 is also run and its max |d| reported, not required to be 0);
  KA-811 : C0 with the rider OFF (inc's keys) == m811 stage 1's stored c0_l25 arrays on L1/k562_gwps (|d| <= 1e-12).

    uv run python experiments/811-lolo-model/composite811.py <fold-id> [--mode day] [--ka]     # writes $COMPOSITE811_OUT
env COMPOSITE811_BANKS: the compact-bank dir (default work/banks; the kit sets work/banks_ext after kit811ext).
"""
from __future__ import annotations

import argparse
import inspect
import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import m811 as M  # noqa: E402

B, E, K = M.B, M.E, M.K
OUT = Path(os.environ.get("COMPOSITE811_OUT", "/work/composite811"))
REFCOMP = HERE / "work" / "refcomp"
REF_C0 = HERE / "work" / "run-m811s1" / "bulk"
LAM = 0.25
TOL = 1e-12
RIDER_SRC = ("    at_floor = w0 <= FLOOR + 1e-12\n"
             "    w705 = np.where(at_floor, w0, np.maximum(FLOOR, w0 * np.exp(THETA * z)))\n")


def rider_705(w0: np.ndarray, z: np.ndarray) -> np.ndarray:
    """ka804.main's 705 expression (its two lines, asserted verbatim in ka804.py below)."""
    FLOOR, THETA = K.FLOOR, K.THETA
    at_floor = w0 <= FLOOR + 1e-12
    w705 = np.where(at_floor, w0, np.maximum(FLOOR, w0 * np.exp(THETA * z)))
    return w705


def top_tier(F, z, have):
    """ka804.main's tier rule for 707 (tertiles of z over the targets carrying U)."""
    tier = np.full(F.n, -1)
    r = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
    tier[np.flatnonzero(have)] = np.minimum((r * 3) // int(have.sum()), 2)
    return tier == 2


def compose(base_fin, U, lam, npn):
    """811's blend on a base composition; lambda 0 returns the base rows themselves (KA-L0)."""
    return base_fin if lam == 0 else B.blend(base_fin, U, lam, npn)


def build(F, mode: str, theta_on: bool = True, lam: float = LAM):
    abc_rules = mode == "ka804" and F.panel == "ABC"
    inc0 = E.named(F, "cap6+P2g(abc)" if abc_rules else "pct+P2g")
    inc = E.Arm(w=inc0.w, q=inc0.q, name="inc")
    p2g = (lambda w: E.p2g_abc(w)) if abc_rules else (lambda w: E.p2g_oct(w, F.carried))
    fin = np.asarray(F.fin, np.float64)
    npn = F.nonpanel
    z = K.zsrc_fold(F)
    w0 = np.asarray(inc.w, float)
    w = rider_705(w0, z) if theta_on else w0
    q = p2g(w) if theta_on else inc.q
    U7, have7 = K.donor_unit_mean(F)
    fin707 = B.blend(fin, U7, K.BETA, npn, rows=top_tier(F, z, have7))
    U, _ = M.fold_U(F, False)
    arms = {"inc": inc,
            "B0": E.Arm(w=w, q=q, name="B0"),
            "B1": E.Arm(w=w, q=q, fin=fin707, name="B1"),
            "C0": E.Arm(w=w, q=q, fin=compose(fin, U, lam, npn), name="C0"),
            "C1": E.Arm(w=w, q=q, fin=compose(fin707, U, lam, npn), name="C1")}
    info = dict(mode=mode, inc_arm="cap6+P2g(abc)" if abc_rules else "pct+P2g", theta_on=theta_on, lam=lam,
                n_floor=int((w0 <= K.FLOOR + 1e-12).sum()), n_w_changed=int(np.sum(w != w0)), n_top707=int(top_tier(F, z, have7).sum()),
                U_targets=int((np.linalg.norm(U, axis=1) > 0).sum()), shared_share_U=M.shared_share(U))
    return arms, info, (fin, fin707, U, npn)


def price_all(F, arms, names):
    out = {}
    for nm in names:
        b = None if nm == "inc" else "inc"
        s, arr = E.price(F, arms[nm], base=arms[b] if b else None)
        arr["targets"] = np.asarray(F.targets)
        out[nm] = (s, {k: v for k, v in arr.items() if not k.startswith("base_")})
    return out


def eq_arrays(a: dict, b: dict) -> float:
    worst = 0.0
    for k, v in b.items():
        if k == "targets":
            assert list(map(str, v)) == list(map(str, a["targets"])), "targets differ"
            continue
        assert k in a, f"missing array {k}"
        worst = max(worst, float(np.max(np.abs(np.asarray(a[k], float) - np.asarray(v, float)))) if np.size(v) else 0.0)
    return worst


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("fold")
    ap.add_argument("--mode", default="day", choices=("day", "ka804"))
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args(argv)
    t0 = time.time()
    log = lambda *x: print(f"[{a.fold}]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    assert RIDER_SRC in inspect.getsource(K), "ka804.py's 705 expression changed: composite811's copy is no longer verbatim"
    if os.environ.get("COMPOSITE811_BANKS"):
        B.BANKS = Path(os.environ["COMPOSITE811_BANKS"])
    F = E.load(a.fold, check=True)
    K.OUT = OUT
    fid = a.fold.replace("/", "__")
    rec = dict(fold=a.fold, n=F.n, held=F.held, sources=F.sources, ka_engine=F.ka, ka={},
               script_sha256=B.sha(Path(__file__)), m811_sha256=B.sha(HERE / "m811.py"), base811_sha256=B.sha(HERE / "base811.py"),
               ka804_sha256=B.sha(Path(inspect.getfile(K))), engine801_sha256=B.sha(Path(inspect.getfile(E))),
               banks=str(B.BANKS), banks_manifest_sha256=B.sha(B.BANKS / "manifest.json"))
    arms, info, (fin, fin707, U, npn) = build(F, a.mode)
    rec["info"] = info
    res = price_all(F, arms, ("inc", "B0", "B1", "C0", "C1"))
    for nm, (s, arr) in res.items():
        tag = nm if a.mode == "day" else f"{nm}_ka804"
        rec.setdefault("arms", {})[tag] = K.save(a.fold, tag, s, arr, None if nm == "inc" else "inc")
        if nm != "inc":
            v = s["vs_base"]
            log(f"{tag:9s} r {100 * v['r']:+.3f}% dpds {v['dpds']:+.5f} dpdsl {v['dpds_lam05']:+.5f} dK {100 * v['dK_rel']:+.2f}% "
                f"off-axis {v['share_unscored']:.2f}")
    if a.ka:
        fails = []

        def check(name, got, ref_path):
            try:
                ref = dict(np.load(ref_path, allow_pickle=True))
                d = eq_arrays(got, ref)
                ok = d <= TOL
            except Exception as e:                              # a missing reference is a failed KA
                d, ok = f"{type(e).__name__}: {e}", False
            rec["ka"][name] = dict(ref=str(ref_path.relative_to(HERE)), worst=d, ok=ok,
                                   ref_sha256=B.sha(ref_path) if ref_path.exists() else None)
            log(f"KA {name}: {'PASS' if ok else 'FAILED'} (worst |d| {d})")
            if not ok:
                fails.append(name)

        # KA-705 / KA-707 in ka804's own convention
        ka_arms, _, _ = build(F, "ka804")
        ka_off, _, _ = build(F, "ka804", theta_on=False)
        p = price_all(F, {"inc": ka_arms["inc"], "B0": ka_arms["B0"], "B1off": ka_off["B1"]}, ("inc", "B0", "B1off"))
        check("KA-705 (ka804-mode B0 == ka804 705)", p["B0"][1], REFCOMP / fid / "705.npz")
        check("KA-707 (ka804-mode B1, rider off == ka804 707)", p["B1off"][1], REFCOMP / fid / "707.npz")
        # KA-L0: lambda 0 identities (bit-exact), plus the blend arithmetic at lambda 0 reported
        l0, _, _ = build(F, a.mode, lam=0.0)
        for c, b_ in (("C0", "B0"), ("C1", "B1")):
            same = all(np.array_equal(getattr(l0[c], f) if getattr(l0[c], f) is not None else fin,
                                      getattr(arms[b_], f) if getattr(arms[b_], f) is not None else fin) for f in ("fin",)) \
                and np.array_equal(np.asarray(l0[c].w), np.asarray(arms[b_].w)) and np.array_equal(np.asarray(l0[c].q), np.asarray(arms[b_].q))
            pr = price_all(F, {"inc": l0["inc"], c: l0[c]}, ("inc", c))[c][1]
            d = eq_arrays(pr, res[b_][1])
            ok = bool(same) and d == 0.0
            rec["ka"][f"KA-L0 {c}(lam 0) == {b_}"] = dict(rows_equal=bool(same), priced_worst=d, ok=ok)
            log(f"KA-L0 {c}(lam 0) == {b_}: {'PASS' if ok else 'FAILED'} (rows equal {same}, priced worst |d| {d})")
            if not ok:
                fails.append(f"KA-L0 {c}")
            base_rows = fin if c == "C0" else fin707
            rec["ka"][f"info: blend arithmetic at lam 0 vs {b_} rows, max |d|"] = float(np.max(np.abs(B.blend(base_rows, U, 0.0, npn) - base_rows)))
        # KA-811: C0 with the rider off == m811's stored c0_l25 (L folds only: the stored arrays exist there)
        if F.panel.startswith("L") and (REF_C0 / fid / "c0_l25.npz").exists():
            off, _, _ = build(F, a.mode, theta_on=False)
            pc = price_all(F, {"inc": off["inc"], "C0": off["C0"]}, ("inc", "C0"))
            check("KA-811 (C0 rider off == m811 c0_l25)", pc["C0"][1], REF_C0 / fid / "c0_l25.npz")
        rec["ka_failed"] = fails
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    (OUT / "meta" / f"{fid}.json").write_text(json.dumps(rec, indent=1, default=float))
    assert not rec.get("ka_failed"), f"KA FAILED: {rec['ka_failed']}"
    log("done")


if __name__ == "__main__":
    main()
