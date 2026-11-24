"""612 ds-c: the CALL-SET-AWARE fold model, exactly as specified in amendment A1 of PREREGISTRATION-draft-KGID.md
(committed 3efebfd7, 2026-10-02 12:52Z, before this file existed).

Production adds the hook's sign only on UNCALLED genes (rc_emit_k: use_s = ~head & ~blk & (asg != 0)). ds612b's fold
model applied it to every gene. Here every form takes the no-hook prediction on production's call set:
  head  = each target's top N_HEAD 400 genes by |donor consensus key| (cons0)
  block = the top N_CPM 2,200 of the rest by the surface's control expression (lib521's construction)
so called cells cancel in every contrast. Folds h1 / HCT116 / HEK293T / K562 at 521's pools f 1 / 0.5 / 0.35 (ds612b).
Forms: nohook, plain (= c4h), kgC4_05 (= production c4kg: CLEAN4 tau 0.05 on the weak half), kgID (= idkg), IDw_C4s.

Order of output (A1): (0) known answer: an EMPTY call set reproduces ds612b exactly; (1) the retrodiction verdict
from the h1 cell alone, criteria (a)-(d); (2) only then the held-out rows and A1's decision quantities.
Analysis only; < 3 GB; seconds. Output: work/ds612c.out, work/ds612c.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ds612 as D  # noqa: E402
import ds612b as B  # noqa: E402
from common480 import boot_delta  # noqa: E402

N_HEAD, N_CPM = 400, 2200
REGIMES = B.REGIMES
SURFS = B.SURFS
FORMS = ["plain", "kgC4_05", "kgID", "IDw_C4s"]
REAL = {"kgID-plain": {1.0: -0.00315, 0.5: -0.00019, 0.35: +0.00067},     # H1 real emitter, 3 seeds (score612.json)
        "kgC4_05-plain": {1.0: +0.00034, 0.5: +0.00126, 0.35: +0.00108}}  # 2 seeds


def callset(cons0, cnt, n_head, n_cpm):
    """Production's call set per target: head by |key|, then the CPM block by control expression among the rest."""
    C = np.zeros(cons0.shape, bool)
    if n_head == 0 and n_cpm == 0:
        return C
    c = np.broadcast_to(np.asarray(cnt, float).reshape(-1), (cons0.shape[1],))
    for i in range(cons0.shape[0]):
        head = np.zeros(cons0.shape[1], bool)
        head[np.argsort(-np.abs(cons0[i]), kind="stable")[:n_head]] = True
        blk = np.zeros(cons0.shape[1], bool)
        blk[np.argsort(-np.where(head, -np.inf, c), kind="stable")[:n_cpm]] = True
        C[i] = head | blk
    return C


def masked_forms(X, n_head, n_cpm):
    P, weak = X["P"], X["weak"]
    F = {"nohook": P["nohook"], "plain": P["plain"], "kgC4_05": np.where(weak, P["plain"], P["nohook"]),
         "kgID": P["kgID"], "IDw_C4s": P["IDw_C4s"]}
    C = callset(X["cons0"], X["S"]["cnt"], n_head, n_cpm)
    return {k: np.where(C, F["nohook"], v) for k, v in F.items()}, C


def table(name, X, p, n_head, n_cpm):
    Fm, C = masked_forms(X, n_head, n_cpm)
    res = {}
    for f in REGIMES:
        Tf = B.pool_T(X["T"], p, f)
        rows = {k: D.nmae_cells(v, Tf)[0] for k, v in Fm.items()}
        r = {"n": int(np.isfinite(rows["nohook"]).sum()),
             "vs_nohook": {k: boot_delta(rows[k], rows["nohook"]) for k in FORMS},       # + = the form is better
             # raw nmae differences, first - second (negative = the first form better), as the real emitter reports
             "kgID-plain": boot_delta(rows["plain"], rows["kgID"]),
             "kgC4_05-plain": boot_delta(rows["plain"], rows["kgC4_05"]),
             "IDw_C4s-kgID": boot_delta(rows["kgID"], rows["IDw_C4s"])}
        res[f] = r
    return res, float(C.sum(1).mean())


def main() -> None:
    BETA, prof = D.setup()
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    X, P_ = {}, {}
    for name in SURFS:
        X[name] = D.surface_forms(name, BETA, prof)
        P_[name] = B.pvals(name, X[name]["T"])
    # (0) known answer: the empty call set must reproduce ds612b's recorded rows exactly
    ref = json.loads((HERE / "work/ds612b.json").read_text())["surfaces"]
    ka_ok = True
    for name in SURFS:
        res0, _ = table(name, X[name], P_[name], 0, 0)
        for f in REGIMES:
            for k in ("plain", "kgID", "IDw_C4s"):
                got, want = res0[f]["vs_nohook"][k][0], ref[name][str(f)]["vs_nohook"][k][0]
                if got != want:
                    ka_ok = False
                    say(f"KA MISMATCH {name} f {f} {k}: {got!r} vs ds612b {want!r}")
    say(f"(0) KA: empty call set == ds612b on 4 folds x 3 regimes x 3 forms: {'PASS' if ka_ok else 'FAIL'}")
    if not ka_ok:
        raise SystemExit("KA failed: ds612c is void (A1)")
    full = {name: table(name, X[name], P_[name], N_HEAD, N_CPM) for name in SURFS}
    # (1) the retrodiction verdict, h1 cell only
    h = full["h1"][0]
    kp = {f: h[f]["kgID-plain"][0] for f in REGIMES}
    cp = {f: h[f]["kgC4_05-plain"][0] for f in REGIMES}
    crit = {"a: f1 kgID-plain in [-0.0063,-0.0016]": -0.0063 <= kp[1.0] <= -0.0016,
            "b: f0.35 kgID-plain > 0": kp[0.35] > 0,
            "c: f0.5 kgID-plain <= +0.0002 and monotone f1 < f0.5 < f0.35": (kp[0.5] <= 0.0002) and (kp[1.0] < kp[0.5] < kp[0.35]),
            "d: kgC4_05-plain > 0 at f0.5 and f0.35": (cp[0.5] > 0) and (cp[0.35] > 0)}
    say(f"\n(1) RETRODICTION on h1 (A1). Call set: mean {full['h1'][1]:.0f} called genes per target.")
    for f in REGIMES:
        say(f"  f {f:<4g} kgID - plain {kp[f]:+.5f} (real {REAL['kgID-plain'][f]:+.5f})   "
            f"kgC4_05 - plain {cp[f]:+.5f} (real {REAL['kgC4_05-plain'][f]:+.5f})")
    for k, v in crit.items():
        say(f"  {'PASS' if v else 'FAIL'}  {k}")
    verdict = all(crit.values())
    say(f"  VERDICT: ds612c {'PASSES' if verdict else 'FAILS'} A1's retrodiction criterion")
    # (2) the held-out rows and A1's decision quantities
    say("\n(2) Held-out rows (Δnmae vs no hook, + = better; contrasts are raw nmae first - second, negative = first better)")
    for name in SURFS:
        res, nc = full[name]
        say(f"== {name} (mean {nc:.0f} called genes per target)")
        for f in REGIMES:
            r = res[f]
            say(f"  f {f:<4g} n {r['n']:3d} | " + "  ".join(f"{k} {r['vs_nohook'][k][0]:+.5f}" for k in FORMS)
                + f" | kgID-plain {r['kgID-plain'][0]:+.5f} [{r['kgID-plain'][1]:+.5f},{r['kgID-plain'][2]:+.5f}]"
                + f"  IDw_C4s-kgID {r['IDw_C4s-kgID'][0]:+.5f}")
    reg = {"plain": [], "kgID": []}
    lead_all = True
    for name in SURFS:
        for f in (0.5, 0.35):
            a, b = full[name][0][f]["vs_nohook"]["plain"][0], full[name][0][f]["vs_nohook"]["kgID"][0]
            best = max(a, b)
            reg["plain"].append(best - a)
            reg["kgID"].append(best - b)
            lead_all &= b > a
    mr = {k: max(v) for k, v in reg.items()}
    say(f"\nA1 TIE rule (applies only if the verdict PASSES): max regret over 4 folds x f 0.5/0.35: plain {mr['plain']:.5f}, "
        f"KGID {mr['kgID']:.5f} -> TIE goes to {'KGID' if mr['kgID'] < mr['plain'] else 'plain'}")
    say(f"A1 BELOW rule: KGID ahead of plain on all 4 folds at both f 0.5 and f 0.35: {lead_all}")
    out = {"verdict": verdict, "criteria": crit, "ka": ka_ok, "max_regret": mr, "kgid_leads_all": bool(lead_all),
           "surfaces": {n: {str(f): full[n][0][f] for f in REGIMES} for n in SURFS},
           "called_per_target": {n: full[n][1] for n in SURFS}}
    (HERE / "work/ds612c.out").write_text("\n".join(lines) + "\n")
    (HERE / "work/ds612c.json").write_text(json.dumps(out, indent=1, default=float))
    say("wrote work/ds612c.out, work/ds612c.json")


if __name__ == "__main__":
    main()
