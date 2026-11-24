"""248 -- the rider through 216's EXACT `nmae` identity (188's split), not through a rate card.

229's `six229.py` (= 222's `six222.py` = 216's `six216.py` with the ladder swapped), with the
arm construction routed through the PRODUCTION builder's own `apply_frozen_cd4_residual`: CD4
by its rule, union 1 (the H1 pair) at 0.28, union 2 at its dose, ONE merged add, ONE preserve.
Everything that computes a number -- 188's exact identity, the three pools, `key_precision`,
the delivered-amplitude readout, the `pds` on the a23 render -- is 216's, imported, not
re-derived.

KNOWN-ANSWER ROWS: `a15 f=0.10`, `cd4 f=0.34` against `216/work/six216.json`; `h1union g=0.28`
against `222/work/six222.json`; and `h1union`, `stack6`, `stack4`, `stack3` against 229's
`six229.json` (its stack arms re-expressed here as union 1 = the six / four / three at 0.35).
The builder's merged path must reproduce every one at 0.00e+00.

229 §6.4 / 222 §8: this instrument MIS-SIGNED its own family's board point twice (−0.000152
predicted vs +0.001762 realised on the CD4 dose step; −0.000108 vs the h1union step's realised
+0.00588 scaled `nmae`).  Read its output as a MAGNITUDE bound, reported separately.

    uv run python six248.py
"""
from __future__ import annotations

import gc
import importlib.util
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
os.environ["VCC_UNION2_SOURCES"] = "dld1_lowmoi"
os.environ["VCC_UNION2_G"] = "0.35"
for _k in ("VCC_CD4_SIZING", "VCC_CD4_TARGET_G", "VCC_UNION_SOURCES", "VCC_CD4_DEFLATE_K"):
    os.environ.pop(_k, None)

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
sys.path.insert(0, str(ROOT / "experiments/216-cd4-hard"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from plolo import (FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, deliver)

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
NMAE_SPAN = 0.60
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
H1_PAIR = ("gse295214_h1", "h1_2025")
STACK4 = ("dld1_lowmoi", "jiang_5line_noK562", "igvf_cardiomyocyte", "myllia_a549")
STACK6 = H1_PAIR + STACK4
STACK3 = H1_PAIR + ("dld1_lowmoi",)
# (name, cd4 residual f, union-1 g, union-1 sources, union-2 g, union-2 sources)
ARMS = [
    ("a15 f=0.10 SHIPPED", 0.10, None, None, None, None),
    ("cd4 f=0.34", 0.34, None, None, None, None),
    ("h1union g=0.28 <- INCUMBENT", 0.34, 0.28, H1_PAIR, None, None),
    ("CUM1 dld1@0.35 <- RIDER", 0.34, 0.28, H1_PAIR, 0.35, ("dld1_lowmoi",)),
    ("CUM4 four@0.35", 0.34, 0.28, H1_PAIR, 0.35, STACK4),
    ("CUM1 dld1@0.28", 0.34, 0.28, H1_PAIR, 0.28, ("dld1_lowmoi",)),
    ("stack6 g=0.35 (229 ARM6)", 0.34, 0.35, STACK6, None, None),
    ("stack4 g=0.35 (229 ARM4)", 0.34, 0.35, STACK4, None, None),
    ("stack3 g=0.35 (229 ARM3)", 0.34, 0.35, STACK3, None, None),
]


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("six248_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["six248_builder"] = module
    spec.loader.exec_module(module)
    return module


MOD = load_builder()


def sources_of(names):
    return tuple((n, MOD.UNION_SOURCE_REGISTRY[n]["path"]) for n in names)


def build_final(pre, op, fold, cd4, resf, u1, g1, u2, g2):
    """Production's own path: `apply_frozen_cd4_residual` with CD4 at `resf`, union 1 at g1,
    union 2 at g2 (each None = absent).  The module globals ARE the builder's knobs."""
    wh = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
    MOD.CD4_TARGET_F = resf
    MOD.H1_UNION_G = g1
    MOD.UNION2_G = g2
    audits = {}
    final, cd4_audit, u1_audit = MOD.apply_frozen_cd4_residual(
        wh, cd4, np.asarray(fold["wcols"], dtype=np.int64),
        np.asarray(fold["panel_positions"], dtype=np.int64),
        extra=u1, extra_audit=None if u1 is None else {},
        extra2=u2, extra2_audit=None if u2 is None else {}, audits_out=audits)
    final[:, fold["panel_positions"]] = 0.0
    fx1 = 0.0 if u1_audit is None else float(u1_audit["residual_frame"]["solved_fraction"])
    fx2 = 0.0 if "union2" not in audits else float(audits["union2"]["residual_frame"]["solved_fraction"])
    return final, float(cd4_audit["fraction"]), fx1, fx2


def main() -> None:
    acc: dict[tuple[str, str], list] = {}
    prec: dict[str, list] = {}
    norms: dict[str, list] = {}
    pds: dict[str, list] = {}
    fractions: dict[str, list] = {}

    for held in FOLDS:
        fold = build_fold(held)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        unions = {}
        for sel in {a[3] for a in ARMS if a[3] is not None} | {a[5] for a in ARMS if a[5] is not None}:
            unions[sel], _details = MOD.prepare_h1_union(
                fold["targets"], OFFICIAL, PANEL, max_rss_gib=40.0, sources=sources_of(sel))
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)

        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
        pools = {
            "wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
            "lfc89": topn(np.abs(LR), POOL_N),
            "h89": topn(np.abs(truth_h[:, sc]), POOL_N),
        }
        Ts = truth_h[:, sc]

        for name, resf, g1, sel1, g2, sel2 in ARMS:
            fin, frac, fx1, fx2 = build_final(
                pre, op, fold, cd4, resf, unions.get(sel1), g1, unions.get(sel2), g2)
            fractions.setdefault(name, []).append(dict(fold=held, cd4=frac, u1=fx1, u2=fx2))
            H = deliver(fin, bA, gA, None)
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, mmask in pools.items():
                acc.setdefault((name, pk), []).append(umd(LP, LR, mmask))
            prec.setdefault(name, []).append(
                key_precision(np.abs(H[:, sc]), H[:, sc], Ts, n=POOL_N))
            norms.setdefault(name, []).append(
                float(np.median(np.linalg.norm(H[:, fold["wcols"]], axis=1))))
            P = deliver(fin, bA, gA, "minnorm")[:, sc]
            A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
            B = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-300)
            S = A @ B.T
            own = np.diag(S).copy()
            n = S.shape[0]
            pds.setdefault(name, []).append(
                (held, float(np.mean(((S < own[:, None]).sum(1)
                                      + 0.5 * ((S == own[:, None]).sum(1) - 1))
                                     / (n - 1)))))
            del fin, H, LP, P, A, B, S
            gc.collect()
        print(f"[{held}] done", flush=True)
        del op, pre, cd4, unions, truth_h, fold, LR, Ts, pools
        gc.collect()

    base = ARMS[0][0]
    inc = ARMS[2][0]
    ref = {pk: np.concatenate([x[2] for x in acc[(base, pk)]])
           for pk in ("wald89", "lfc89", "h89")}
    out = {}
    b0 = {}
    for pk in ref:
        U = np.concatenate([x[0] for x in acc[(base, pk)]])
        M = np.concatenate([x[1] for x in acc[(base, pk)]])
        b0[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
    p0 = float(np.mean(prec[base]))
    a15n = float(np.mean(norms[base]))
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    P("=" * 124)
    P("188's EXACT `nmae` identity, ordering precision, delivered amplitude -- 4 PLOLO folds, "
      "challenge frame A; every arm through the builder's apply_frozen_cd4_residual")
    P("=" * 124)
    P(f"{'arm':<30}{'dlv/a15':>9}{'nmae wald89':>13}{'d vs a15':>11}{'ovr nmae':>11}{'top89 prec':>12}"
      f"{'d prec':>10}{'nmae lfc89':>12}{'nmae h89':>11}{'pds 4f':>9}{'pds 3f':>9}")
    for name, *_ in ARMS:
        vals = {}
        for pk in ref:
            U = np.concatenate([x[0] for x in acc[(name, pk)]])
            M = np.concatenate([x[1] for x in acc[(name, pk)]])
            vals[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
        nm = vals["wald89"][0]
        ov = -(nm - b0["wald89"][0]) / NMAE_SPAN / 6.0
        pr = float(np.mean(prec[name]))
        dn = float(np.mean(norms[name])) / a15n
        p4 = float(np.mean([v for _f, v in pds[name]]))
        p3 = float(np.mean([v for f, v in pds[name] if f in CANCER]))
        out[name] = dict(nmae_wald89=nm, d_nmae_vs_a15=nm - b0["wald89"][0],
                         overall_nmae_vs_a15=ov, top89_precision=pr,
                         d_precision_vs_a15=pr - p0, delivered_norm_ratio=dn,
                         nmae_lfc89=vals["lfc89"][0], nmae_h89=vals["h89"][0],
                         pds_4folds=p4, pds_3cancer=p3,
                         pds_per_fold={f: v for f, v in pds[name]},
                         fractions=fractions[name],
                         U_sum=float(vals["wald89"][1].sum()),
                         M_sum=float(vals["wald89"][2].sum()))
        P(f"{name:<30}{dn:>9.4f}{nm:>13.5f}{nm-b0['wald89'][0]:>+11.5f}{ov:>+11.6f}{pr:>12.5f}"
          f"{pr-p0:>+10.5f}{vals['lfc89'][0]:>12.5f}{vals['h89'][0]:>11.5f}{p4:>9.5f}{p3:>9.5f}")

    ni, pi, di = (out[inc]["nmae_wald89"], out[inc]["top89_precision"], out[inc]["delivered_norm_ratio"])
    P(f"\nsame table, deltas against the board incumbent `{inc}`:")
    P(f"{'arm':<30}{'d nmae':>10}{'ovr nmae':>11}{'d prec':>10}{'d dlv':>9}{'d pds 4f':>10}{'d pds 3f':>10}")
    for name, *_ in ARMS:
        o = out[name]
        o["d_nmae_vs_incumbent"] = o["nmae_wald89"] - ni
        o["overall_nmae_vs_incumbent"] = -(o["nmae_wald89"] - ni) / NMAE_SPAN / 6.0
        o["d_precision_vs_incumbent"] = o["top89_precision"] - pi
        o["d_delivered_vs_incumbent"] = o["delivered_norm_ratio"] - di
        o["d_pds_4folds_vs_incumbent"] = o["pds_4folds"] - out[inc]["pds_4folds"]
        o["d_pds_3cancer_vs_incumbent"] = o["pds_3cancer"] - out[inc]["pds_3cancer"]
        dU = o["U_sum"] - out[inc]["U_sum"]
        dM = o["M_sum"] - out[inc]["M_sum"]
        o["dU_vs_incumbent"] = dU
        o["dM_vs_incumbent"] = dM
        o["efficiency_dM_over_dU"] = (dM / dU) if dU != 0 else float("nan")
        P(f"{name:<30}{o['d_nmae_vs_incumbent']:>+10.5f}{o['overall_nmae_vs_incumbent']:>+11.6f}"
          f"{o['d_precision_vs_incumbent']:>+10.5f}{o['delivered_norm_ratio']-di:>+9.4f}"
          f"{o['d_pds_4folds_vs_incumbent']:>+10.6f}{o['d_pds_3cancer_vs_incumbent']:>+10.6f}")

    P(f"\n188's SPLIT -- does the arm move SIGNS or MAGNITUDES?  (vs `{inc}`; break-even e = 0.5; "
      f"224's arm died at 0.072; 229's ARM6 read -0.033)")
    P(f"{'arm':<30}{'dU':>14}{'dM':>14}{'e = dM/dU':>12}{'verdict':>28}")
    for name, *_ in ARMS:
        o = out[name]
        if name == inc:
            continue
        e = o["efficiency_dM_over_dU"]
        v = ("SIGNS (nmae improves)" if e > 0.5 else "MAGNITUDES (nmae worsens)" if e == e else "degenerate")
        P(f"{name:<30}{o['dU_vs_incumbent']:>+14.4f}{o['dM_vs_incumbent']:>+14.4f}{e:>12.4f}{v:>28}")

    # ---- KNOWN-ANSWER rows -------------------------------------------------------------------
    ka_report = []
    keys = ("nmae_wald89", "top89_precision", "delivered_norm_ratio")
    refs = [
        (ROOT / "experiments/216-cd4-hard/work/six216.json", None,
         [("a15 f=0.10 SHIPPED", "a15 f=0.10 SHIPPED"), ("cd4 f=0.34", "res f=0.34 <- INCUMBENT")]),
        (ROOT / "experiments/222-h1-union-build/work/six222.json", "arms",
         [("h1union g=0.28 <- INCUMBENT", "h1union g=0.28"), ("cd4 f=0.34", "cd4 f=0.34 <- INCUMBENT"),
          ("a15 f=0.10 SHIPPED", "a15 f=0.10 SHIPPED")]),
        (ROOT / "experiments/229-stack4-build/work/six229.json", "arms",
         [("h1union g=0.28 <- INCUMBENT", "h1union g=0.28 <- INCUMBENT"),
          ("stack6 g=0.35 (229 ARM6)", "stack6 g=0.35 <- ARM"), ("stack4 g=0.35 (229 ARM4)", "stack4 g=0.35"),
          ("stack3 g=0.35 (229 ARM3)", "stack3 g=0.35")]),
    ]
    for path, sub, pairs in refs:
        if not path.is_file():
            P(f"\n(no known-answer file {path})")
            continue
        k = json.loads(path.read_text())
        if sub:
            k = k.get(sub, k)
        P(f"\nKNOWN-ANSWER against {path.relative_to(ROOT)}:")
        for mine, theirs in pairs:
            if theirs not in k:
                P(f"  (no arm {theirs!r} there)")
                continue
            for key in keys + (("pds_3cancer",) if "pds_3cancer" in k[theirs] else ()):
                got, want = out[mine][key], k[theirs].get(key)
                if want is None:
                    continue
                ok = abs(got - want) <= 5e-6
                ka_report.append(dict(arm=mine, key=key, got=got, want=want, ok=bool(ok), source=path.name))
                P(f"  {'PASS' if ok else 'FAIL'}  {mine:30s} {key:22s} got {got:.8f}  ref {want:.8f}  |d| {abs(got-want):.2e}")
    n_ok = sum(r["ok"] for r in ka_report)
    P(f"\nknown-answer rows: {n_ok}/{len(ka_report)} pass")

    (WORK / "six248.json").write_text(json.dumps(dict(arms=out, known_answer=ka_report), indent=2))
    (WORK / "six248.log").write_text("\n".join(lines) + "\n")
    P(f"wrote {WORK/'six248.json'}")
    if n_ok != len(ka_report):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
