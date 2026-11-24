#!/usr/bin/env python3
"""658 reader tests: real captures plus synthetic entries built from them (base + exact member deltas).
usage: uv run python experiments/658-prereg-1007-builds/test658.py"""
from __future__ import annotations

import contextlib
import copy
import io
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import read658 as T  # noqa: E402

R = T.R
TD = HERE / "work/test"
TD.mkdir(parents=True, exist_ok=True)
NMAE_SCALE = None   # scaled-per-raw factor for nmae on A/B/C, measured from 446bj2z vs 446bj2


def save(e: dict, name: str) -> str:
    p = TD / f"{name}.json"
    p.write_text(json.dumps({k: v for k, v in e.items() if not k.startswith("_")}))
    return str(p)


def tweak(e: dict, member: str, raw: float = 0.0, scaled: float = 0.0, name: str = "t") -> dict:
    out = copy.deepcopy(e)
    out[R.RAW[member]] += raw
    out[R.SCALED[member]] += scaled
    out["score_avg"] += scaled / 6.0
    out["entry_id"] = name
    return out


def recipe(env: dict, name: str, q=("1.75", "0.5", "noprog")) -> str:
    p = TD / f"{name}.recipe.json"
    p.write_text(json.dumps({"env": env, "args": ["b", "w", "n", *q, "4"], "pending_treated_off": [], "problems": []}))
    return str(p)


def run(argv: list[str]) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = T.main(argv)
    return rc, buf.getvalue()


OK, BAD = [], []


def check(name: str, cond: bool, out: str = ""):
    (OK if cond else BAD).append(name)
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        print(out)


def main():
    L = {k: T.load(v) for k, v in dict(b446=T.E446B, bj2=T.E446BJ2, bj2z=T.E446BJ2Z, p479=T.E479P2G, c495=T.E495,
                                         p498=T.E498P2G, dr453=T.E453DR).items()}
    global NMAE_SCALE
    NMAE_SCALE = (L["bj2z"][R.SCALED["nmae"]] - L["bj2"][R.SCALED["nmae"]]) / (L["bj2z"][R.RAW["nmae"]] - L["bj2"][R.RAW["nmae"]])
    # a synthetic 496a: the hook moves nmae by -0.0025 raw on 446bj2 (inside 496a's predicted -0.0015 ... -0.0035)
    hook = tweak(L["bj2"], "nmae", raw=-0.0025, scaled=-0.0025 * NMAE_SCALE, name="syn496a")
    H = save(hook, "syn496a")

    # ---- allsw
    a1 = T.synth(L["p498"], (L["bj2z"], L["bj2"]), tag="allsw_additive")
    rc, out = run(["allsw", save(a1, "allsw_additive")])
    check("T1 allsw additive -> ADDITIVE, R1 PASS (v2), prediction 0.2754683, R1 window [-0.0008573, +0.0003427]",
          rc == 0 and "VERDICT: ADDITIVE" in out and "R1 (637 --validate, leak band v2: centre -0.0002, tol 0.0006): PASS" in out
          and "0.2754683" in out and "I in [-0.0008573, +0.0003427]" in out, out)
    a2 = tweak(a1, "nmae", raw=0.0018 / -NMAE_SCALE, scaled=0.0018, name="allsw_weak")      # +0.0003 overall interaction
    rc, out = run(["allsw", save(a2, "allsw_weak")])
    check("T2 allsw +0.0003 -> WEAK INTERACTION, R1 PASS under v2 (upper edge +0.0003427)", rc == 0 and "WEAK INTERACTION" in out
          and "): PASS" in out, out)
    a2c = tweak(a1, "nmae", raw=0.0024 / -NMAE_SCALE, scaled=0.0024, name="allsw_weak2")    # +0.0004 overall interaction
    rc, out = run(["allsw", save(a2c, "allsw_weak2")])
    check("T2c allsw +0.0004 -> WEAK INTERACTION, R1 FAIL (above +0.0003427)", rc == 0 and "WEAK INTERACTION" in out
          and "): FAIL" in out, out)
    a2b = tweak(a1, "nmae", raw=0.0036 / NMAE_SCALE, scaled=-0.0036, name="allsw_neg")       # -0.0006 overall
    rc, out = run(["allsw", save(a2b, "allsw_neg")])
    check("T2b allsw -0.0006 -> INTERACTION, R1 PASS (lower edge -0.0008573)", rc == 0 and "VERDICT: INTERACTION" in out
          and "): PASS" in out, out)
    a2d = tweak(a1, "nmae", raw=0.0054 / NMAE_SCALE, scaled=-0.0054, name="allsw_neg2")      # -0.0009 overall
    rc, out = run(["allsw", save(a2d, "allsw_neg2")])
    check("T2d allsw -0.0009 -> INTERACTION, R1 FAIL (below -0.0008573)", rc == 0 and "VERDICT: INTERACTION" in out
          and "): FAIL" in out, out)
    a3 = tweak(a1, "pds", raw=1e-5, name="allsw_f1")
    rc, out = run(["allsw", save(a3, "allsw_f1")])
    check("T3 allsw bulk moved 1e-5 raw -> F1 FALSIFIED (rc 3)", rc == 3 and "F1 pds" in out, out)
    a4 = tweak(a1, "reach", raw=1e-4, name="allsw_f2")
    rc, out = run(["allsw", save(a4, "allsw_f2")])
    check("T4 allsw zs-only, reach moved -> F2 FALSIFIED", rc == 3 and "F2 reach" in out, out)
    a5 = T.synth(L["p498"], (L["bj2z"], L["bj2"]), (hook, L["bj2"]), tag="allsw_hook")
    rc, out = run(["allsw", save(a5, "allsw_hook"), "--hook", H])
    check("T5 allsw with hook, additive -> ADDITIVE, R1 PASS", rc == 0 and "VERDICT: ADDITIVE" in out and "): PASS" in out, out)

    # ---- dr secondary
    rc, out = run(["dr", T.E453DR.replace(T.E453DR, "LrkhsIXWMjEsvRXSGToa"), "--recipe", recipe({"DEF_DAY_LODO": "1"}, "lodo_only")])
    check("T6 460dr vs 453dr, LODO only -> PASS (gap +0.0006356)", rc == 0 and "-> PASS" in out and "+0.0006356" in out, out)
    cur = {"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI_PCT": "96.64697522769474", "DEF_DAY_SHARPEN_PT": "tercile:2.5,2.0,1.75",
           "DEF_DAY_DE_MAG": "dn_k250_zs"}
    RC524 = recipe(cur, "r524")
    s524 = T.synth(L["dr453"], (L["bj2"], L["b446"]), (L["c495"], L["bj2"]), (L["p479"], L["bj2"]), (L["bj2z"], L["bj2"]),
                   tag="syn524dr")
    S524 = save(s524, "syn524dr")
    rc, out = run(["dr", S524, "--recipe", RC524])
    check("T7 synthetic 524dr = 453dr + Sigma -> PASS, gap 0", rc == 0 and "VERDICT: PASS" in out and "+0.0000000 -> PASS" in out, out)
    lag = tweak(s524, "mse", raw=0.012, scaled=-0.0126, name="syn524lag")
    rc, out = run(["dr", save(lag, "syn524lag"), "--recipe", RC524])
    check("T8 bulk shortfall -> LAG with the mse owners named", rc == 0 and "-> LAG" in out and "mse   d - Sigma" in out
          and "lodo_hi_pct" in out, out)
    fz = tweak(s524, "nmae", raw=0.0025, scaled=-0.0025 * NMAE_SCALE, name="syn524f")
    rc, out = run(["dr", save(fz, "syn524f"), "--recipe", RC524])
    check("T9 nmae beyond Sigma by 0.0025 raw -> FALSIFIED", rc == 3 and "DE member nmae" in out, out)
    rc, out = run(["dr", S524, "--recipe", recipe({**cur, "DEF_DAY_PROGRAMME_SOURCE": "x"}, "unknown")])
    check("T10 an unwired switch -> REFUSED (rc 2)", rc == 2 and "REFUSED" in out, out)
    rc, out = run(["dr", S524, "--recipe", recipe(cur, "q2", q=("2.0", "0.5", "noprog"))])
    check("T10b q differs from 524's -> REFUSED", rc == 2 and "REFUSED" in out, out)

    # ---- dr primary (524dr2 vs 524dr)
    hk = {**cur, "DEF_DAY_MDE_SIGN": "clean4:0.05"}
    tc = {**cur, "DEF_DAY_SHARPEN_PT": "tercilec:2.5,2.0,1.75"}
    both = {**tc, "DEF_DAY_MDE_SIGN": "clean4:0.05"}
    e1 = T.synth(s524, (hook, L["bj2"]), tag="syn524dr2_hook")
    rc, out = run(["dr", save(e1, "syn524dr2_hook"), "--recipe", recipe(hk, "r_hk"), "--prev", S524, "--prev-recipe", RC524, "--hook", H])
    check("T11 N1 hook new, delivered -> CARRIED, bulk pinned", rc == 0 and "branch N1" in out and "-> CARRIED" in out, out)
    rc, out = run(["dr", S524, "--recipe", recipe(hk, "r_hk2"), "--prev", S524, "--prev-recipe", RC524, "--hook", H])
    check("T12 N1 hook new, not delivered -> primary LAG; 524's nmae falsifier also fires (0.0025 > 0.002 raw) -> rc 3",
          rc == 3 and "branch N1" in out and "LAG (the day does not deliver" in out and "DE member nmae" in out, out)
    e3 = tweak(tweak(s524, "pds", raw=0.00003, scaled=0.0006, name="x"), "reach", raw=0.0008, scaled=0.0009, name="syn524dr2_tc")
    rc, out = run(["dr", save(e3, "syn524dr2_tc"), "--recipe", recipe(tc, "r_tc"), "--prev", S524, "--prev-recipe", RC524])
    check("T13 N2 tercilec new, bulk +0.0001, DE re-roll +0.00015 -> CARRIED", rc == 0 and "branch N2" in out
          and "CARRIED (tercilec)" in out and "VERDICT: CARRIED" in out, out)
    e4 = tweak(s524, "pds", raw=0.0002, scaled=0.003, name="syn524dr2_tcbad")
    rc, out = run(["dr", save(e4, "syn524dr2_tcbad"), "--recipe", recipe(tc, "r_tc2"), "--prev", S524, "--prev-recipe", RC524])
    check("T14 N2 tercilec new, bulk +0.0005 -> FALSIFIED", rc == 3 and "only tercilec new" in out, out)
    e5 = tweak(T.synth(s524, (hook, L["bj2"]), tag="y"), "reach", raw=0.0010, scaled=0.0012, name="syn524dr2_both")
    rc, out = run(["dr", save(e5, "syn524dr2_both"), "--recipe", recipe(both, "r_both"), "--prev", S524, "--prev-recipe", RC524,
                   "--hook", H])
    check("T15 N3 tercilec + hook -> CARRIED / CARRIED", rc == 0 and "branch N3" in out and "CARRIED (tercilec) / CARRIED" in out, out)
    rc, out = run(["dr", S524, "--recipe", RC524, "--prev", S524, "--prev-recipe", RC524])
    check("T16 N0 same recipe, identical members -> TREE-IDENTICAL", rc == 0 and "TREE-IDENTICAL" in out, out)
    rc, out = run(["dr", S524, "--recipe", recipe({**cur, "DEF_DAY_SHARPEN_PT": "tercile:3.0,2.0,1.75"}, "r_chg"), "--prev", S524,
                   "--prev-recipe", RC524])
    check("T17 a changed switch value -> REFUSED", rc == 2 and "changed ['sharpen_pt']" in out, out)
    rc, out = run(["dr", S524, "--recipe", recipe({**cur, "DEF_DAY_H1_UNION_G": "0.84"}, "r_h1"), "--alias", "DEF_DAY_H1_UNION_G=h1_dose",
                   "--h1", T.E453DR, "--prev", S524, "--prev-recipe", RC524])
    check("T18 an aliased H1-dose switch reads through --alias/--h1 (branch N2)", rc in (0, 3) and "branch N2" in out, out)

    # ---- bulkid658 on tiny submission-shaped h5ads
    import subprocess
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.default_rng(658)
    obs = pd.DataFrame({"context": pd.Categorical(["A"] * 6 + ["B"] * 6), "target_gene": pd.Categorical(["g1", "g1", "g1", "g2", "g2", "g2"] * 2)},
                       index=[f"c{i}" for i in range(12)])
    var = pd.DataFrame(index=[f"G{i}" for i in range(7)])
    X = rng.poisson(3, size=(12, 7)).astype(np.float32)
    X[X == 0] = 1
    def write(M, name):
        a = ad.AnnData(X=sp.csr_matrix(M), obs=obs.copy(), var=var.copy())
        path = TD / f"{name}.h5ad"
        a.write_h5ad(path)
        return str(path)
    pa = write(X, "bid_a")
    pb = write(X.copy(), "bid_b")
    M2 = X.copy(); M2[0, 2] += 1; M2[1, 2] -= 1          # a count moved between two cells of one (A, g1) block: bulk equal
    pc = write(M2, "bid_c")
    M3 = X.copy(); M3[4, 5] += 1                          # a count added: bulk differs
    pd_ = write(M3, "bid_d")
    bid = str(HERE / "bulkid658.py")
    r1 = subprocess.run([sys.executable, bid, pa, pb, "--exact"], capture_output=True, text=True)
    check("T19 bulkid: identical files -> exit 0, exact IDENTICAL", r1.returncode == 0 and "exact content: IDENTICAL" in r1.stdout, r1.stdout + r1.stderr)
    r2 = subprocess.run([sys.executable, bid, pa, pc, "--exact"], capture_output=True, text=True)
    check("T20 bulkid: a count moved within a block -> bulk identical, exact DIFFERENT (exit 1)",
          r2.returncode == 1 and "0 (block, gene) sums differ" in r2.stdout and "exact content: DIFFERENT" in r2.stdout, r2.stdout + r2.stderr)
    r3 = subprocess.run([sys.executable, bid, pa, pd_], capture_output=True, text=True)
    check("T21 bulkid: a count added -> bulk DIFFERENT (exit 1)", r3.returncode == 1 and "1 (block, gene) sums differ" in r3.stdout,
          r3.stdout + r3.stderr)

    # ---- branch658 (524dr2's recipe classifier)
    import branch658 as Bc
    base = {"env": dict(cur), "args": ["b", "w", "n", "1.75", "0.5", "noprog", "4"], "problems": []}
    def rr(env=None, args=None, problems=None):
        r = copy.deepcopy(base)
        if env is not None: r["env"] = env
        if args is not None: r["args"] = args
        if problems is not None: r["problems"] = problems
        return r
    cases = [
        ("T22 branch: identical -> N0", rr(), "N0 none"),
        ("T23 branch: + hook -> N1", rr(env=hk), "N1 mde_sign"),
        ("T24 branch: tercile -> tercilec -> N2", rr(env=tc), "N2 tercilec"),
        ("T25 branch: tercilec + hook -> N3", rr(env=both), "N3 tercilec,mde_sign"),
        ("T26 branch: + H1 key -> N2 h1_dose", rr(env={**cur, "DEF_DAY_H1_UNION_G": "0.84"}), "N2 h1_dose(DEF_DAY_H1_UNION_G)"),
    ]
    for name, new_r, want in cases:
        got = Bc.classify(new_r, base)
        check(name, got == want, got)
    stops = [
        ("T27 branch: dropped switch -> STOP", rr(env={k: v for k, v in cur.items() if k != "DEF_DAY_DE_MAG"})),
        ("T28 branch: other sharpen change -> STOP", rr(env={**cur, "DEF_DAY_SHARPEN_PT": "tercile:3.0,2.0,1.75"})),
        ("T29 branch: unwired new knob -> STOP", rr(env={**cur, "DEF_DAY_DE_UPRN": "b_k100_q20r3p10"})),
        ("T30 branch: hook with another tau -> STOP", rr(env={**cur, "DEF_DAY_MDE_SIGN": "clean4:0.10"})),
        ("T31 branch: def_day args changed -> STOP", rr(args=["b", "w", "n", "2.0", "0.5", "noprog", "4"])),
        ("T32 branch: recipe problems -> STOP", rr(problems=["x"])),
    ]
    for name, new_r in stops:
        got = Bc.classify(new_r, base)
        check(name, got.startswith("STOP"), got)

    # ---- driver DRY runs: the last log line is the single terminal line
    import os
    def last_line(path):
        return Path(path).read_text().rstrip("\n").splitlines()[-1]
    env = dict(os.environ, DRY="1")
    r = subprocess.run(["bash", str(HERE / "build_allsw658.sh")], env=env, capture_output=True, text=True)
    ll = last_line(HERE / "work/allsw/dry-build_allsw658.log")
    check("T33 build_allsw658 DRY: rc 0, terminal 'NOT READY DRY=1'", r.returncode == 0 and ll.endswith("NOT READY DRY=1 (plan only; nothing built)"), ll)
    live = TD / "live_prev.recipe.json"
    subprocess.run(["uv", "run", "python", "scripts/october_recipe.py", "command", "--bundle", "experiments/180-def-readiness/work/def-bundle",
                    "--work", str(HERE / "work/dr2/run-dry"), "--name", "524dr2", "--allow-pending", "--json"],
                   cwd=str(HERE.parents[1]), text=True, check=True, stdout=open(live, "w"), stderr=subprocess.DEVNULL)
    r = subprocess.run(["bash", str(HERE / "recipe_dr524dr2.sh")], env=dict(env, DR2_PREV_RECIPE=str(live)), capture_output=True, text=True)
    ll = last_line(HERE / "work/dr2/dry-recipe_dr524dr2.log")
    check("T34 recipe_dr524dr2 DRY vs the same recipe: branch N0, terminal NOT READY DRY", r.returncode == 0 and ll.endswith("(plan only; nothing built; branch N0)"), ll)
    lj = json.loads(live.read_text()); lj["env"].pop("DEF_DAY_DE_MAG", None); stale = TD / "stale_prev.recipe.json"; stale.write_text(json.dumps(lj))
    r = subprocess.run(["bash", str(HERE / "recipe_dr524dr2.sh")], env=dict(env, DR2_PREV_RECIPE=str(stale)), capture_output=True, text=True)
    ll = last_line(HERE / "work/dr2/dry-recipe_dr524dr2.log")
    check("T35 recipe_dr524dr2 DRY with an unwired change: terminal NOT READY ... amend 658b first", "NOT READY unwired recipe changes" in ll, ll)

    print(f"\n{len(OK)} passed, {len(BAD)} failed")
    (HERE / "work/test658.out").write_text(f"{len(OK)} passed, {len(BAD)} failed\n" + "\n".join("PASS " + x for x in OK)
                                           + "".join("\nFAIL " + x for x in BAD) + "\n")
    return 0 if not BAD else 1


if __name__ == "__main__":
    sys.exit(main())
