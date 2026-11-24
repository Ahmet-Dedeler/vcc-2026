"""246 — the six-member route of 216 (`six216.py`) / 222 (`six222.py`) for the arms that
survived the sweep: 188's EXACT `nmae` identity on three pools, `key_precision` (the ordering
half of `reach`/`fid`), delivered amplitude, and `pds` on the a23 render -- 4 PLOLO folds,
challenge frame A, `kolf2` reported apart.

KNOWN-ANSWER ROWS (abort on failure, 5e-6):
  a15 f=0.10             nmae_wald89 / top89_precision / delivered_norm_ratio == 216 six216.json
  cd4 f=0.34             same, == 216
  STACK (h1union g=0.28) nmae_wald89 / top89_precision == 222 six222.json (0.98569 / 0.60312)

The `nmae` identity has mis-signed the board twice (222 §8, 216 §4.4): its delta is reported
as a MAGNITUDE BOUND, never as a sign.

    uv run python six246.py
"""
from __future__ import annotations

import gc
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
import lib246 as L  # noqa: E402
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from sweep246b import partial_row_matched_self  # noqa: E402,F401  (keeps the import graph honest)

plolo = L.plolo
NMAE_SPAN = 0.60
CANCER = L.CANCER


def arms_for(st):
    """(name, parts) -- everything through `L.compose`."""
    cd4, union = st["cd4"], st["union"]
    f010, _ = L.solve_fraction(st["whitened"], cd4, 0.10)
    out = [("a15 f=0.10 SHIPPED", [(cd4, f010)]),
           ("cd4 f=0.34 <- cd4f034", [(cd4, st["frac_cd4"])]),
           ("STACK h1union g=0.28", L.stack_parts(st))]
    for k in (2, 3, 5, 7, 10):
        d, _ = L.deflate(cd4, "svd", k)
        out.append((f"B2_cd4_svd{k}", [(d, st["frac_cd4"]), (union, st["frac_x"])]))
    for k in (3, 5):
        d, _ = L.deflate(cd4, "svd", k)
        fr, _ = L.solve_fraction(st["whitened"], d, L.CD4_F)
        out.append((f"B2_cd4_svd{k}_resolved", [(d, fr), (union, st["frac_x"])]))
    d, _ = L.deflate(cd4, "rand", 60, seed=L.RAND_SEED)
    out.append(("B2_cd4_rand60 (ctrl)", [(d, st["frac_cd4"]), (union, st["frac_x"])]))
    part = L.cd4_partial(st["targets"], equalise=False)
    if part.condition_details[0]["n_partial"] > 0:
        fp, _ = L.solve_fraction_carrying(st["whitened"], part, 0.34)
        out.append(("C_partial@0.34", L.stack_parts(st) + [(part, fp)]))
        for k in (3, 5):
            d, _ = L.deflate(cd4, "svd", k)
            out.append((f"E_cd4_svd{k}+partial@0.34",
                        [(d, st["frac_cd4"]), (union, st["frac_x"]), (part, fp)]))
    return out


def main() -> None:
    acc, prec, norms, pds = {}, {}, {}, {}
    names = None
    for held in L.FOLDS:
        st = L.stage(held)
        fold = st["fold"]
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
        pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
                 "lfc89": topn(np.abs(LR), POOL_N),
                 "h89": topn(np.abs(truth_h[:, sc]), POOL_N)}
        Ts = truth_h[:, sc]
        arms = arms_for(st)
        names = [n for n, _ in arms]
        for name, parts in arms:
            fin = L.compose(st, parts)
            H = L.plolo.deliver(fin, bA, gA, None)
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, mmask in pools.items():
                acc.setdefault((name, pk), []).append(umd(LP, LR, mmask))
            prec.setdefault(name, []).append(
                key_precision(np.abs(H[:, sc]), H[:, sc], Ts, n=POOL_N))
            norms.setdefault(name, []).append(
                float(np.median(np.linalg.norm(H[:, fold["wcols"]], axis=1))))
            P = L.plolo.deliver(fin, bA, gA, "minnorm")[:, sc]
            A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
            B = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-300)
            S = A @ B.T
            own = np.diag(S).copy()
            n = S.shape[0]
            pds.setdefault(name, []).append((held, float(np.mean(
                ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1))
                / (n - 1)))))
            del fin, H, LP, P, A, B, S
            gc.collect()
        print(f"[{held}] done ({len(arms)} arms)", flush=True)
        del st, truth_h, LR, Ts, pools
        gc.collect()

    base, inc = names[0], names[2]
    ref = {pk: np.concatenate([x[2] for x in acc[(base, pk)]])
           for pk in ("wald89", "lfc89", "h89")}
    b0 = {}
    for pk in ref:
        U = np.concatenate([x[0] for x in acc[(base, pk)]])
        M = np.concatenate([x[1] for x in acc[(base, pk)]])
        b0[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
    p0 = float(np.mean(prec[base]))
    a15n = float(np.mean(norms[base]))
    out = {}
    print("\n" + "=" * 128)
    print("188's EXACT nmae identity, ordering precision, delivered amplitude -- 4 PLOLO folds, "
          "frame A;  pds on the a23 render (4 folds / 3 cancer)")
    print("=" * 128)
    print(f"{'arm':28s}{'dlv/a15':>9}{'nmae w89':>10}{'d vs inc':>10}{'ovr nmae':>10}"
          f"{'e=dM/dU':>9}{'prec89':>9}{'d prec':>9}{'nmae lfc':>10}{'nmae h89':>10}"
          f"{'pds4f':>9}{'pds3f':>9}{'d pds3f':>10}")
    for name in names:
        vals = {}
        for pk in ref:
            U = np.concatenate([x[0] for x in acc[(name, pk)]])
            M = np.concatenate([x[1] for x in acc[(name, pk)]])
            vals[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
        nm, U, M = vals["wald89"]
        pr = float(np.mean(prec[name]))
        dn = float(np.mean(norms[name])) / a15n
        p4 = float(np.mean([v for _, v in pds[name]]))
        p3 = float(np.mean([v for f, v in pds[name] if f in CANCER]))
        out[name] = dict(nmae_wald89=nm, top89_precision=pr, delivered_norm_ratio=dn,
                         nmae_lfc89=vals["lfc89"][0], nmae_h89=vals["h89"][0],
                         pds_4folds=p4, pds_3cancer=p3, U=U.tolist() if False else None,
                         pds_per_fold={f: v for f, v in pds[name]})
        out[name].pop("U")
    ni, pi, di = (out[inc]["nmae_wald89"], out[inc]["top89_precision"],
                  out[inc]["delivered_norm_ratio"])
    Ui = np.concatenate([x[0] for x in acc[(inc, "wald89")]])
    Mi = np.concatenate([x[1] for x in acc[(inc, "wald89")]])
    for name in names:
        o = out[name]
        U = np.concatenate([x[0] for x in acc[(name, "wald89")]])
        M = np.concatenate([x[1] for x in acc[(name, "wald89")]])
        du = float(np.mean((U - Ui) / ref["wald89"]))
        dm = float(np.mean((M - Mi) / ref["wald89"]))
        o["dU_vs_inc"], o["dM_vs_inc"] = du, dm
        o["e_vs_inc"] = dm / du if abs(du) > 1e-12 else float("nan")
        o["d_nmae_vs_inc"] = o["nmae_wald89"] - ni
        o["overall_nmae_vs_inc"] = -(o["nmae_wald89"] - ni) / NMAE_SPAN / 6.0
        o["d_precision_vs_inc"] = o["top89_precision"] - pi
        o["d_delivered_vs_inc"] = o["delivered_norm_ratio"] - di
        o["d_pds_3cancer_vs_inc"] = o["pds_3cancer"] - out[inc]["pds_3cancer"]
        o["d_pds_4folds_vs_inc"] = o["pds_4folds"] - out[inc]["pds_4folds"]
        print(f"{name:28s}{o['delivered_norm_ratio']:>9.4f}{o['nmae_wald89']:>10.5f}"
              f"{o['d_nmae_vs_inc']:>+10.5f}{o['overall_nmae_vs_inc']:>+10.6f}"
              f"{o['e_vs_inc']:>9.3f}{o['top89_precision']:>9.5f}{o['d_precision_vs_inc']:>+9.5f}"
              f"{o['nmae_lfc89']:>10.5f}{o['nmae_h89']:>10.5f}{o['pds_4folds']:>9.5f}"
              f"{o['pds_3cancer']:>9.5f}{o['d_pds_3cancer_vs_inc']:>+10.6f}")

    # known-answer rows
    ka = []
    six216 = json.loads((ROOT / "experiments/216-cd4-hard/work/six216.json").read_text())
    six222 = json.loads((ROOT / "experiments/222-h1-union-build/work/six222.json").read_text())
    checks = [("a15 f=0.10 SHIPPED", six216["a15 f=0.10 SHIPPED"],
               ("nmae_wald89", "top89_precision", "delivered_norm_ratio")),
              ("cd4 f=0.34 <- cd4f034", six216["res f=0.34 <- INCUMBENT"],
               ("nmae_wald89", "top89_precision", "delivered_norm_ratio")),
              ("STACK h1union g=0.28", six222["arms"]["h1union g=0.28"],
               ("nmae_wald89", "top89_precision", "delivered_norm_ratio", "pds_4folds",
                "pds_3cancer"))]
    print("\nKNOWN-ANSWER rows:")
    fails = 0
    for mine, theirs, keys in checks:
        for key in keys:
            got, want = out[mine][key], theirs[key]
            ok = abs(got - want) <= 5e-6
            fails += (not ok)
            ka.append(dict(arm=mine, key=key, got=got, want=want, ok=bool(ok)))
            print(f"  {'PASS' if ok else 'FAIL'}  {mine:26s} {key:22s} got {got:.8f} "
                  f"want {want:.8f} |d| {abs(got-want):.2e}")
    (L.WORK / "six246.json").write_text(json.dumps(dict(arms=out, known_answer=ka),
                                                   indent=2))
    print(f"\nwrote {L.WORK/'six246.json'}")
    if fails:
        raise SystemExit(f"{fails} known-answer failures")


if __name__ == "__main__":
    main()
