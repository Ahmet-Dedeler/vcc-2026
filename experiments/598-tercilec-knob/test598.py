#!/usr/bin/env python3
"""598 tests, on temp trees only (the live repo is never written).

The canonical POST-bundle tree is rebuilt the way post594 [1] does: git HEAD's copies + 538 -> 554 -> 590. 598 is applied
on a copy of it. Light: one 85 MB dump (~0.4 GB RSS), no builds, no scorer.

  T0  canonical tree: the builder is HEAD's (538/554/590 do not touch it), def_day.sh carries 554's mark
  T1  --check writes nothing; T2 apply, a second apply refused, a partial apply refused, a pre-bundle tree refused
  T3  py_compile on the five Python files, bash -n on def_day.sh
  T4  def_pin_panel_artifacts.rewrite (the day-of registry pin) commutes with the patch
  T5  synthetic arrays: tie blocks (floor at the lower cut), the cap tie (cap block at the upper cut), all-carried ==
      tercile (bit for bit, several seeds), efp-like 50 % non-carried (== 493h's w-space rule; MID holds every
      non-carried target, robust to 1e-14 noise), a zero-carrier target, the refusals; tercile / stercile of the new
      function == the old function bit for bit on every case
  T6  the real A/B/C dump (479p2g's 00ccfcba): the old function reproduces 479p2g's recorded audit exactly (instrument
      known answer); new tercile / stercile == old; dump_key_carried resolves lodo415_lo0.5_hi3 (8c82bef3) and the ratio
      holds to ~1e-14; wrong keys refuse; tercilec = 98 / 104 / 98 and moves exactly 4 targets, == 493h's w-space rule
      on the artifact's w; writes work/ka598c_expected.json (the expected audit for KA-c)
  T7  the env parse accepts tercilec (any case), refuses malformed qs and a missing ballast
  T8  def_day.sh: the case lines accept tercile / stercile / tercilec, refuse others; tercilec needs dump_key_carried
  T9  def_day_helpers.py meta_matches with a tercilec meta (the pre-598 helper rejects it)
  T10 def_pipeline.september_rows with tercilec (the pre-598 driver raises ValueError)
  T11 october_recipe check: tercilec accepted with the patched files, flagged with either unpatched file, malformed
      forms flagged; scripts/test_october_recipe.py passes on the canonical tree and on the patched one
  T12 gate524's tercilec row (538's fixture): PASS on a consistent meta, FAIL on a wrong key / maxdev / count
  T13 ka598.py (the coordinator's KA-c checker) on synthetic metas

  uv run python experiments/598-tercilec-knob/test598.py      (exit 0 = all PASS)
"""
from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import os
import py_compile
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/598-tercilec-knob"
W = X / "work"
PY = str(ROOT / ".venv/bin/python")
PATCH = X / "patch598.py"
DUMP = ROOT / "experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"
DUMP_SHA = "00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667"
REF479 = ROOT / "submissions/479p2g-cap15-partial50-clean4q80-pertscale3-ballast-sharpenp2g-stacklodo05load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json"
REF453 = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
CANON = ["scripts/def_day.sh", "scripts/def_day_helpers.py", "experiments/448-def-rehearsal-q175/relevel_iter.py",
         "scripts/def_sameline.py", "scripts/october_recipe.py", "docs/october-decisions.json",
         "experiments/524-october-recipe-dr/gate524.py", "scripts/def_postbuild.py", "scripts/def_pipeline.py",
         "scripts/def_preflight.py", "scripts/def_identity.py", "scripts/frozen_inputs.json", "scripts/make_submission_pm_cd4.py"]
# (script, apply args, check args): each step is skipped when its --check says the tree already carries it (after 10-05
# the bundle is on HEAD, and the canonical tree is then HEAD itself)
BUNDLE = [(ROOT / "experiments/538-cloud-bundle-extras/apply_october538.py", [], ["--check"]),
          (ROOT / "experiments/554-dose-guard/patch_guard554.py", [], ["--check"]),
          (ROOT / "experiments/590-identity-validation/patch_identity590.py", ["--apply"], ["--check"])]
BUNDLE_STATE: dict[str, str] = {}
SIX = ["scripts/make_submission_pm_cd4.py", "scripts/def_day.sh", "scripts/def_day_helpers.py", "scripts/def_pipeline.py",
       "scripts/october_recipe.py", "experiments/524-october-recipe-dr/gate524.py"]

RESULTS: list[tuple[str, bool, str]] = []


def rec(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'} {name}{': ' + detail if detail else ''}", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)


def canon_tree(dst: Path) -> None:
    for rel in CANON:
        (dst / rel).parent.mkdir(parents=True, exist_ok=True)
        (dst / rel).write_bytes(subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True, check=True).stdout)
    shutil.copymode(ROOT / "scripts/def_day.sh", dst / "scripts/def_day.sh")
    for script, apply_args, check_args in BUNDLE:
        c = run([PY, script, "--root", dst, *check_args])
        if c.returncode != 0 and "already" in (c.stdout + c.stderr):
            BUNDLE_STATE[script.name] = "on HEAD"
            continue
        q = run([PY, script, "--root", dst, *apply_args])
        if c.returncode != 0 or q.returncode != 0:
            raise SystemExit(f"bundle step {script.name} failed: {c.stdout[-300:]} {c.stderr[-300:]} {q.stdout[-300:]} {q.stderr[-300:]}")
        BUNDLE_STATE[script.name] = "applied here"


def copy_tree(src: Path, dst: Path) -> None:
    shutil.copytree(src, dst, symlinks=True)


# ------------------------------------------------------------------------------------------------ function extraction
def module_ns(builder_text: str) -> dict:
    """sharpen_q_by_tercile, dump_key_carried (if present), file_sha256 and PERT_SCALE_ARTIFACTS from a builder text,
    exec'd in a namespace whose ROOT is the live repo (the registry's artifacts are read from there)."""
    tree = ast.parse(builder_text)
    ns: dict = {"np": np, "Path": Path, "hashlib": hashlib, "ROOT": ROOT}
    want_f = {"sharpen_q_by_tercile", "dump_key_carried", "file_sha256"}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in want_f:
            exec(compile(ast.get_source_segment(builder_text, node), "<builder>", "exec"), ns)
        elif isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "PERT_SCALE_ARTIFACTS" for t in node.targets):
            exec(compile(ast.get_source_segment(builder_text, node), "<builder>", "exec"), ns)
    return ns


def parse_block(builder_text: str) -> str:
    """The module-level `if BULK_SHARPEN_Q is not None and ... in ("tercile", ...)` statement."""
    for node in ast.parse(builder_text).body:
        if isinstance(node, ast.If):
            src = ast.get_source_segment(builder_text, node)
            if src.startswith("if BULK_SHARPEN_Q is not None and BULK_SHARPEN_Q.split("):
                return src
    raise SystemExit("parse block not found")


# ----------------------------------------------------------------------------------------------------------- synthetic
def synth(w: np.ndarray, *, seed: int, c: float = 0.6, n_np: int = 40, noise: float = 0.0, zero: tuple = ()):
    """bulk = carriers * (c * w)[:, None] (exactly apply_pert_scale's row scale), with optional ~1e-14 noise on r.
    Genes: n_np non-panel genes + the panel's own genes (excluded from every norm)."""
    rng = np.random.default_rng(seed)
    n = len(w)
    perts = [f"T{i:03d}" for i in range(n)]
    genes = [f"G{j:03d}" for j in range(n_np)] + perts
    car = rng.normal(size=(n, len(genes))) * rng.uniform(0.5, 3.0, size=(n, 1))
    for i in zero:
        car[i, :n_np] = 0.0
    s = c * w * (1.0 + noise * rng.standard_normal(n))
    return car * s[:, None], car, genes, perts


def tiers(fn, bulk, car, genes, perts, mode, carried=None):
    kw = {} if carried is None else {"carried": carried}
    q, a = fn(bulk, car, genes, perts, {"mode": mode, "q": [0.0, 1.0, 2.0]}, **kw)
    return np.rint(q).astype(int), a


def rule493h(w: np.ndarray, carried: np.ndarray, floor_w: float) -> np.ndarray:
    """493h's w-space rule (comp493h.py): cut points over carried targets, non-carried -> MID, floor -> LOW."""
    floor = w <= floor_w + 1e-12
    ec = np.quantile(w[carried], [1 / 3, 2 / 3])
    t = np.where(~carried, 1, np.digitize(w, ec))
    return np.where(floor, 0, t)


def same_q(fa, fb, *args, mode, **kw) -> bool:
    qa, aa = fa(*args, {"mode": mode, "q": [2.5, 2.0, 1.75]}, **kw)
    qb, ab = fb(*args, {"mode": mode, "q": [2.5, 2.0, 1.75]}, **kw)
    return np.array_equal(qa, qb) and qa.tobytes() == qb.tobytes() and aa == ab


def t5(old: dict, new: dict) -> None:
    fo, fn = old["sharpen_q_by_tercile"], new["sharpen_q_by_tercile"]
    rng = np.random.default_rng(598)
    # S1 tie block: the floor (150 of 300) straddles the lower cut
    w1 = np.r_[np.full(150, 0.5), rng.uniform(0.6, 3.0, 150)]
    b, c, g, p = synth(w1, seed=1, noise=1e-14)
    t_pl, a_pl = tiers(fn, b, c, g, p, "tercile")
    t_c, a_c = tiers(fn, b, c, g, p, "tercilec", carried=np.ones(300, bool))
    rec("T5 S1 floor tie block -> LOW (tercile)", (t_pl[:150] == 0).all() and a_pl["floor_to_low"] == 150,
        f"counts {a_pl['counts']} tie_to_mid {a_pl['tie_to_mid']}")
    rec("T5 S1 all-carried tercilec == tercile", np.array_equal(t_pl, t_c) and a_c["edges"] == a_pl["edges"]
        and a_c["counts"] == a_pl["counts"] and a_c["moved_vs_tercile"] == {} and a_c["n_carried"] == 300)
    # S2 the cap tie: 120 targets at the cap (w = 3) straddle the upper cut
    w2 = np.r_[np.full(50, 0.5), rng.uniform(0.6, 2.9, 130), np.full(120, 3.0)]
    b, c, g, p = synth(w2, seed=2, noise=1e-14)
    t_pl, a_pl = tiers(fn, b, c, g, p, "tercile")
    t_c, a_c = tiers(fn, b, c, g, p, "tercilec", carried=np.ones(300, bool))
    rec("T5 S2 cap tie block -> HIGH (tercile)", (t_pl[180:] == 2).all() and a_pl["counts"][2] == 120,
        f"counts {a_pl['counts']}")
    rec("T5 S2 all-carried tercilec == tercile", np.array_equal(t_pl, t_c) and a_c["counts"] == a_pl["counts"])
    # cap tie among carried with non-carried present (tercilec's own tie handling)
    carried2 = np.ones(300, bool)
    carried2[rng.choice(np.arange(50, 180), 30, replace=False)] = False
    w2c = np.where(carried2, w2, 1.0)
    b, c, g, p = synth(w2c, seed=22, noise=1e-14)
    t_c, a_c = tiers(fn, b, c, g, p, "tercilec", carried=carried2)
    want = rule493h(w2c, carried2, 0.5)
    rec("T5 S2c cap tie with 30 non-carried: carried cap block HIGH, == 493h", (t_c[180:] == 2).all()
        and (t_c[~carried2] == 1).all() and np.array_equal(t_c, want), f"counts {a_c['counts']}")
    # S3 efp-like: 150 of 300 non-carried (w = 1 exactly, r = c up to 1e-14 noise), straddling the cuts
    carried3 = np.zeros(300, bool)
    carried3[rng.choice(300, 150, replace=False)] = True
    w3 = np.ones(300)
    w3[carried3] = np.r_[np.full(30, 0.5), rng.uniform(0.55, 3.0, 120)]
    b, c, g, p = synth(w3, seed=3, noise=1e-14)
    t_c, a_c = tiers(fn, b, c, g, p, "tercilec", carried=carried3)
    t_pl, a_pl = tiers(fn, b, c, g, p, "tercile")
    want = rule493h(w3, carried3, 0.5)
    car_counts = [int((t_c[carried3] == i).sum()) for i in range(3)]
    rec("T5 S3 efp-like: every non-carried target MID", (t_c[~carried3] == 1).all() and a_c["noncarried_to_mid"] == 150)
    rec("T5 S3 efp-like: carried terciles 50/50/50, floor LOW, == 493h's w-space rule",
        car_counts == [50, 50, 50] and (t_c[w3 <= 0.5] == 0).all() and np.array_equal(t_c, want),
        f"counts {a_c['counts']} carried {car_counts}; plain tercile counts {a_pl['counts']} (moves {len(a_c['moved_vs_tercile'])})")
    rec("T5 S3 moved_vs_tercile lists exactly the differing targets",
        set(a_c["moved_vs_tercile"]) == {p[i] for i in np.flatnonzero(t_c != t_pl)}
        and all(v == [int(t_pl[p.index(k)]), int(t_c[p.index(k)])] for k, v in a_c["moved_vs_tercile"].items()))
    # S3 noise robustness: a different noise draw on the non-carried block does not move anything
    b2, c2, g2, p2 = synth(w3, seed=3, noise=3e-14)
    t_c2, _ = tiers(fn, b2, c2, g2, p2, "tercilec", carried=carried3)
    rec("T5 S3 tercilec is insensitive to 1e-14 noise in the w = 1 block", np.array_equal(t_c, t_c2))
    # S4 all-carried == tercile on random keys with clip blocks and noise, several seeds (bit for bit, q and audit)
    eq = []
    for s in range(6):
        r = np.random.default_rng(100 + s)
        w = np.clip(np.exp(r.normal(0, 0.9, 300)), 0.3 + 0.2 * (s % 2), 3.0 + 3.0 * (s % 3 == 0))
        b, c, g, p = synth(w, seed=100 + s, noise=1e-14)
        qa, aa = fn(b, c, g, p, {"mode": "tercile", "q": [2.5, 2.0, 1.75]})
        qb, ab = fn(b, c, g, p, {"mode": "tercilec", "q": [2.5, 2.0, 1.75]}, carried=np.ones(300, bool))
        common = {k: ab[k] for k in aa if k != "mode"}
        eq.append(qa.tobytes() == qb.tobytes() and common == {k: aa[k] for k in aa if k != "mode"})
    rec("T5 S4 all-carried tercilec == tercile on 6 random keys (q bit for bit, audit)", all(eq), str(eq))
    # S5 a zero-carrier target (r non-finite): MID in every mode, carried or not
    w5 = np.r_[np.full(40, 0.5), rng.uniform(0.6, 3.0, 260)]
    b, c, g, p = synth(w5, seed=5, zero=(7, 200))
    carried5 = np.ones(300, bool)
    carried5[200] = False
    ok5 = []
    for mode, kw in (("tercile", {}), ("stercile", {}), ("tercilec", {"carried": carried5})):
        t, a = tiers(fn, b, c, g, p, mode, **kw)
        if mode == "stercile":   # strength = ||carriers|| = 0 is finite: the weakest rows, LOW (today's rule, unchanged)
            ok5.append(bool(t[7] == 0 and t[200] == 0 and a["nonfinite"] == 0))
        else:                    # the key ratio 0/0 is non-finite: MID
            ok5.append(bool(t[7] == 1 and t[200] == 1 and a["nonfinite"] == 2))
    rec("T5 S5 zero-carrier targets: MID in tercile / tercilec (r non-finite), LOW in stercile (strength 0)", all(ok5), str(ok5))
    # S6 refusals
    errs = []
    for kw in ({}, {"carried": np.ones(299, bool)}, {"carried": np.zeros(300, bool)}):
        try:
            fn(b, c, g, p, {"mode": "tercilec", "q": [2.5, 2.0, 1.75]}, **kw)
            errs.append(False)
        except SystemExit:
            errs.append(True)
    rec("T5 S6 tercilec refuses a missing / wrong-length / all-False carried", all(errs), str(errs))
    # S7 tercile / stercile: new function == old function bit for bit on every synthetic case above
    cases = [(w1, 1, ()), (w2, 2, ()), (w2c, 22, ()), (w3, 3, ()), (w5, 5, (7, 200))]
    ok7 = []
    for w, s, z in cases:
        b, c, g, p = synth(w, seed=s, noise=1e-14, zero=z)
        for mode in ("tercile", "stercile"):
            ok7.append(same_q(fo, fn, b, c, g, p, mode=mode))
    rec("T5 S7 tercile / stercile: new == old (q bytes and audit) on all synthetic cases", all(ok7), f"{sum(ok7)}/{len(ok7)}")


def t6(old: dict, new: dict) -> dict:
    fo, fn, dk = old["sharpen_q_by_tercile"], new["sharpen_q_by_tercile"], new["dump_key_carried"]
    rec("T6 dump sha", sha(DUMP) == DUMP_SHA, sha(DUMP)[:16])
    with np.load(DUMP, allow_pickle=False) as z:
        bulk = np.asarray(z["bulk"], np.float64)
        car = np.asarray(z["carriers"], np.float64)
        perts = [str(x) for x in z["perts"]]
        genes = [str(x) for x in z["genes"]]
        env = json.loads(str(z["env"]))
    pt = {"mode": "tercile", "q": [2.5, 2.0, 1.75]}
    q_old, a_old = fo(bulk, car, genes, perts, pt)
    ref = json.loads(REF479.read_text())["config"]["bulk_sharpen_pt"]
    rec("T6 instrument KA: the old function reproduces 479p2g's recorded audit exactly", a_old == ref, json.dumps(a_old))
    rec("T6 tercile new == old on the dump (q bytes, audit)", same_q(fo, fn, bulk, car, genes, perts, mode="tercile"))
    rec("T6 stercile new == old on the dump (q bytes, audit)", same_q(fo, fn, bulk, car, genes, perts, mode="stercile"))
    carried, key = dk(env, perts, bulk, car, genes)
    art = new["PERT_SCALE_ARTIFACTS"][key["carried_key"]]
    rec("T6 dump_key_carried resolves the dump's own key", key["carried_key"] == "lodo415_lo0.5_hi3"
        and key["carried_key_sha256"] == art[1] and int(carried.sum()) == 294 and key["key_ratio_maxdev"] <= 1e-12,
        f"{key['carried_key']} {key['carried_key_sha256'][:16]} c {key['key_ratio_c']:.13g} maxdev {key['key_ratio_maxdev']:.2e} "
        f"carried {int(carried.sum())}")
    refused = []
    for kid in ("wshare_lo0.3_hi3", "lodo415_lo0.5_hi6", "wshare_lo0.5_hi3", "nope", ""):
        try:
            dk({**env, "VCC_PERT_SCALE": kid}, perts, bulk, car, genes)
            refused.append((kid, False))
        except SystemExit as e:
            refused.append((kid, True))
            print(f"      refused {kid!r}: {str(e)[:150]}")
    rec("T6 a wrong / unknown / absent key refuses (ratio check or registry)", all(r for _, r in refused), str(refused))
    q_c, a_c = fn(bulk, car, genes, perts, {"mode": "tercilec", "q": [2.5, 2.0, 1.75]}, carried=carried)
    t_old = np.rint(fo(bulk, car, genes, perts, {"mode": "tercile", "q": [0.0, 1.0, 2.0]})[0]).astype(int)
    t_c = np.rint(fn(bulk, car, genes, perts, {"mode": "tercilec", "q": [0.0, 1.0, 2.0]}, carried=carried)[0]).astype(int)
    moved = a_c["moved_vs_tercile"]
    rec("T6 tercilec on A/B/C: 98 / 104 / 98, exactly 4 targets move (492/493h's prediction)",
        a_c["counts"] == [98, 104, 98] and len(moved) == 4 and a_c["n_carried"] == 294 and a_c["noncarried_to_mid"] == 6
        and int((t_c != t_old).sum()) == 4, f"counts {a_c['counts']} moved {moved}")
    with np.load(art[0], allow_pickle=False) as z:
        w = np.asarray(z["w"], np.float64)
    rec("T6 tercilec == 493h's w-space rule on the key artifact's w (all 300 targets)",
        np.array_equal(t_c, rule493h(w, carried, 0.5)))
    rec("T6 tercilec q follows its tiers", np.array_equal(q_c, np.asarray([2.5, 2.0, 1.75])[t_c]))
    expected = {**a_c, **key}
    W.mkdir(parents=True, exist_ok=True)
    (W / "ka598c_expected.json").write_text(json.dumps(expected, indent=1) + "\n")
    detail = {p: {"from": int(t_old[perts.index(p)]), "to": int(t_c[perts.index(p)]),
                  "w": float(w[perts.index(p)]), "r": float(np.linalg.norm(bulk[perts.index(p)][[g not in set(perts) for g in genes]])
                                                         / np.linalg.norm(car[perts.index(p)][[g not in set(perts) for g in genes]]))}
              for p in moved}
    return {"moved": detail, "edges_tercile": a_old["edges"], "edges_tercilec": a_c["edges"], "key": key,
            "counts_tercilec": a_c["counts"]}


def t7(old_text: str, new_text: str) -> None:
    blk_new = parse_block(new_text)
    out = []
    for val, ballast, want in (("tercilec:2.5,2.0,1.75", "1", ("tercilec", [2.5, 2.0, 1.75])),
                               ("TERCILEC:2.5,2,1.75", "1", ("tercilec", [2.5, 2.0, 1.75])),
                               ("tercile:2.5,2.0,1.75", "1", ("tercile", [2.5, 2.0, 1.75])),
                               ("stercile:2.5,2.0,1.75", "1", ("stercile", [2.5, 2.0, 1.75])),
                               ("tercilec:2.5,2.0", "1", None), ("tercilec:4,2,2", "1", None), ("tercilec:a,b,c", "1", None),
                               ("tercilec:2.5,2.0,1.75", None, None)):
        ns = {"BULK_SHARPEN_Q": val, "PERT_SCALE_BALLAST": ballast, "BULK_SHARPEN_PT": None}
        try:
            exec(blk_new, ns)
            got = (ns["BULK_SHARPEN_PT"]["mode"], ns["BULK_SHARPEN_PT"]["q"]) if ns["BULK_SHARPEN_PT"] else None
            ok = want is not None and got == want and ns["BULK_SHARPEN_Q"] is None
        except SystemExit:
            ok = want is None
        out.append(ok)
    rec("T7 env parse: tercilec (any case) parsed; short / out-of-range / non-numeric / no-ballast refused", all(out), str(out))
    ns = {"BULK_SHARPEN_Q": "tercilec:2.5,2.0,1.75", "PERT_SCALE_BALLAST": "1", "BULK_SHARPEN_PT": None}
    exec(parse_block(old_text), ns)
    rec("T7 (pre-598 parse leaves tercilec to the scalar path, which would float() it)", ns["BULK_SHARPEN_PT"] is None)


def t8(tree: Path, canon: Path) -> None:
    day = (tree / "scripts/def_day.sh").read_text().splitlines()
    i = next(k for k, ln in enumerate(day) if ln.startswith('  case "$DEF_DAY_SHARPEN_PT" in tercile:'))
    snippet = "\n".join(day[i:i + 4])
    rec("T8 the four lines are the case, the 479 grep, the 598 comment, the tercilec grep",
        "tercilec:*,*,*" in day[i] and "def sharpen_q_by_tercile" in day[i + 1] and day[i + 2].lstrip().startswith("# 598")
        and "def dump_key_carried" in day[i + 3])
    script = 'die(){ echo "DIE: $*"; exit 3; }\n' + snippet + '\necho OK\n'
    res = []
    for val, cwd, want in (("tercile:2.5,2.0,1.75", tree, "OK"), ("stercile:2.5,2.0,1.75", tree, "OK"),
                           ("tercilec:2.5,2.0,1.75", tree, "OK"), ("tercilec:2.5,2.0", tree, "DIE"),
                           ("terciled:2.5,2.0,1.75", tree, "DIE"), ("tercilec:2.5,2.0,1.75", canon, "DIE"),
                           ("tercile:2.5,2.0,1.75", canon, "OK")):
        p = subprocess.run(["/bin/bash", "-c", script], cwd=cwd, env={**os.environ, "DEF_DAY_SHARPEN_PT": val},
                           capture_output=True, text=True)
        res.append((val, cwd == tree, p.stdout.strip().startswith(want)))
    rec("T8 def_day.sh case lines (patched tree; the canonical builder lacks dump_key_carried -> tercilec dies)",
        all(r[2] for r in res), str([r for r in res if not r[2]] or "all as expected"))


def t9(tree: Path, canon: Path, tmp: Path) -> None:
    meta = tmp / "m_c.json"
    meta.write_text(json.dumps({"config": {"bulk_sharpen_pt": {"mode": "tercilec", "q": [2.5, 2.0, 1.75], "counts": [98, 104, 98]}}}))
    meta_t = tmp / "m_t.json"
    meta_t.write_text(json.dumps({"config": {"bulk_sharpen_pt": {"mode": "tercile", "q": [2.5, 2.0, 1.75]}}}))
    h, h0 = tree / "scripts/def_day_helpers.py", canon / "scripts/def_day_helpers.py"
    a = run([PY, h, "meta_matches", meta, "A=1,B=1,C=1", "tercilec:2.5,2.0,1.75"])
    b = run([PY, h, "meta_matches", meta_t, "A=1,B=1,C=1", "tercilec:2.5,2.0,1.75"])
    c = run([PY, h, "meta_matches", meta_t, "A=1,B=1,C=1", "tercile:2.5,2.0,1.75"])
    d = run([PY, h0, "meta_matches", meta, "A=1,B=1,C=1", "tercilec:2.5,2.0,1.75"])
    rec("T9 meta_matches: tercilec meta matches, a tercile meta does not, tercile unchanged; the pre-598 helper rejects",
        a.returncode == 0 and b.returncode == 1 and c.returncode == 0 and d.returncode != 0,
        f"rc {a.returncode}/{b.returncode}/{c.returncode}/pre-598 {d.returncode}")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.pop(0)
    return mod


def t10(tree: Path, canon: Path) -> None:
    new, old = _load(tree / "scripts/def_pipeline.py", "dp598"), _load(canon / "scripts/def_pipeline.py", "dp590")
    resolved = {"VCC_BULK_SHARPEN_Q": "tercilec:2.5,2.0,1.75"}
    cfg = {"bulk_sharpen_pt": {"mode": "tercilec", "q": [2.5, 2.0, 1.75], "counts": [98, 104, 98]}}
    rows = {r[0]: r for r in new.september_rows(cfg, resolved)}
    ok_new = rows.get("bulk_sharpen_pt", (None, None, None, False))[3]
    bad = {r[0]: r for r in new.september_rows({"bulk_sharpen_pt": {"mode": "tercile", "q": [2.5, 2.0, 1.75]}}, resolved)}
    try:
        old.september_rows(cfg, resolved)
        old_raises = False
    except ValueError:
        old_raises = True
    rec("T10 driver: tercilec row passes, a tercile meta fails it, the pre-598 driver raises ValueError",
        ok_new and not bad["bulk_sharpen_pt"][3] and old_raises, f"row {rows.get('bulk_sharpen_pt')}")


def t11(tree: Path, canon: Path, tmp: Path) -> None:
    d = json.loads((tree / "docs/october-decisions.json").read_text())

    def check(value, def_day, builder):
        dd = json.loads(json.dumps(d))
        dd["switches"]["sharpen_pt"].update({"status": "DECIDED", "value": value})
        f = tmp / "dec.json"
        f.write_text(json.dumps(dd))
        env = {**os.environ, "OCTOBER_DECISIONS": str(f), "OCTOBER_DEF_DAY": str(def_day), "OCTOBER_BUILDER": str(builder)}
        p = subprocess.run([PY, str(tree / "scripts/october_recipe.py"), "check", "--allow-pending"], env=env,
                           capture_output=True, text=True)
        return [ln for ln in p.stdout.splitlines() if ln.startswith("PROBLEM") and ("sharpen" in ln or "598" in ln or "tercilec" in ln)]

    T, C = tree / "scripts", canon / "scripts"
    cases = [("tercilec:2.5,2.0,1.75", T / "def_day.sh", T / "make_submission_pm_cd4.py", 0),
             ("tercile:2.5,2.0,1.75", T / "def_day.sh", T / "make_submission_pm_cd4.py", 0),
             ("stercile:2.5,2.0,1.75", T / "def_day.sh", C / "make_submission_pm_cd4.py", 0),
             ("tercilec:2.5,2.0,1.75", T / "def_day.sh", C / "make_submission_pm_cd4.py", 1),
             ("tercilec:2.5,2.0,1.75", C / "def_day.sh", T / "make_submission_pm_cd4.py", 1),
             ("tercilec:2.5,2.0", T / "def_day.sh", T / "make_submission_pm_cd4.py", 1),
             ("tercilec:4,2,2", T / "def_day.sh", T / "make_submission_pm_cd4.py", 1),
             ("tercilex:2.5,2.0,1.75", T / "def_day.sh", T / "make_submission_pm_cd4.py", 1)]
    got = [(v, len(check(v, dd, bu))) for v, dd, bu, _ in cases]
    rec("T11 october_recipe: tercilec clean with the patched pair, one problem with either unpatched file, malformed flagged",
        [g[1] for g in got] == [c[3] for c in cases], str(got))
    res = {}
    for name, root in (("canonical", canon), ("patched", tree)):
        t = root / "scripts/test_october_recipe.py"
        t.write_bytes((ROOT / "scripts/test_october_recipe.py").read_bytes())
        for rel in ("experiments/497-october-switches", "experiments/528-callprop-block", "src"):
            link = root / rel
            if not link.exists():
                link.symlink_to(ROOT / rel)
        p = subprocess.run([PY, str(t)], cwd=root, capture_output=True, text=True)
        fails = sorted(ln.split()[1] for ln in p.stderr.splitlines() if ln.startswith(("FAIL:", "ERROR:")))
        res[name] = (p.returncode, (p.stderr.strip().splitlines() or ["?"])[-1], fails)
        t.unlink()
    # The suite is stale on its own (1 failure on the live tree since 10-02's sharpen_pt decision, 5 on the post-bundle
    # tree from 538's new PENDING switches): 598 must change nothing in it.
    rec("T11 scripts/test_october_recipe.py: the patched tree fails exactly the tests the canonical tree fails",
        res["canonical"][2] == res["patched"][2] and res["canonical"][0] == res["patched"][0], str(res))


def t12(tree: Path, tmp: Path) -> None:
    gate = tree / "experiments/524-october-recipe-dr/gate524.py"
    (tree / "submissions").mkdir(exist_ok=True)
    shutil.copy(REF453, tree / "submissions" / REF453.name)
    base_pt = {"mode": "tercilec", "q": [2.5, 2.0, 1.75], "edges": [0.3, 0.9], "counts": [98, 104, 98], "floor_to_low": 97,
               "tie_to_mid": 0, "nonfinite": 0, "n_carried": 294, "noncarried_to_mid": 6, "moved_vs_tercile": {"X": [0, 1]},
               "carried_key": "lodo_lo0.5_pct_def", "carried_key_sha256": "ab" * 32, "key_ratio_c": 0.6,
               "key_ratio_maxdev": 2e-14}

    def go(env, pt, log_extra):
        work = Path(tempfile.mkdtemp(prefix="gate598-", dir=tmp))
        rd = work / "run"
        rd.mkdir()
        (work / "recipe524.json").write_text(json.dumps({"env": env, "args": ["B", "W", "N", "1.75", "auto", "pt", "4"],
                                                         "pending_treated_off": [], "problems": []}))
        (rd / "def_day.log").write_text("CONVERGED at D=1,E=1,F=1 after 0 levelled build(s)\n"
                                        "final: lvl0 factors D=1,E=1,F=1 -> out\n" + log_extra)
        (rd / "lvl0.ntclevel.json").write_text(json.dumps({"per_context": {c: {"rho": 0.5} for c in "DEF"}}))
        (rd / "postbuild.json").write_text(json.dumps({"checks": [{"check": "a", "pass": True}]}))
        (rd / "audit428.txt").write_text("ZERO block total (target excluded): 0 in 0 blocks\n")
        meta = json.loads(REF453.read_text())
        meta["config"].pop("bulk_sharpen_q", None)
        meta["config"]["bulk_sharpen_pt"] = pt
        stem = work / "pkg"
        Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(meta))
        blob = b"vcc-synthetic-598"
        Path(str(stem) + ".vcc").write_bytes(blob)
        (work / "chain.log").write_text(f"novelty rc=0\npackage pkg.vcc sha256 {hashlib.sha256(blob).hexdigest()}\n")
        p = subprocess.run([PY, str(gate), "--write", "--work", str(work), "--stem", str(stem), "--chainlog",
                            str(work / "chain.log")], capture_output=True, text=True)
        row = [ln for ln in p.stdout.splitlines() if "(598)" in ln]
        return p.returncode, row
    lodo = "2026-10-22T10:00:00Z LODO key ON for the stack dump: lodo_lo0.5_pct_def 0123456789abcdef; cap 6.0\n"
    env = {"DEF_DAY_SHARPEN_PT": "tercilec:2.5,2.0,1.75"}
    want_pass = go(env, base_pt, lodo)
    own = go(env, {**base_pt, "carried_key": json.loads(REF453.read_text())["config"].get("pert_scale")}, "")
    fails = [go(env, {**base_pt, "carried_key": "wshare_lo0.3_hi3_def"}, lodo),
             go(env, {**base_pt, "key_ratio_maxdev": 1e-3}, lodo),
             go(env, {**base_pt, "noncarried_to_mid": 5}, lodo),
             go(env, {k: v for k, v in base_pt.items() if k != "n_carried"}, lodo)]
    plain = go({"DEF_DAY_SHARPEN_PT": "tercile:2.5,2.0,1.75"},
               {k: base_pt[k] for k in ("q", "edges", "counts", "floor_to_low", "tie_to_mid", "nonfinite")} | {"mode": "tercile"}, lodo)
    rec("T12 gate524: tercilec row PASSES (LODO key from the log; the build's own key without one)",
        want_pass[0] == 0 and want_pass[1] and want_pass[1][0].strip().startswith("ok") and own[0] == 0,
        f"{want_pass[1]} | own-key rc {own[0]}")
    rec("T12 gate524: wrong key / maxdev / count / missing n_carried FAIL; the tercile mode has no 598 row",
        all(f[0] != 0 and f[1] and f[1][0].strip().startswith("FAIL") for f in fails) and plain[0] == 0 and not plain[1],
        str([f[0] for f in fails]) + f" plain rc {plain[0]}")


def t13(tmp: Path) -> None:
    ka = X / "ka598.py"
    ref = json.loads(REF479.read_text())
    exp = json.loads((W / "ka598c_expected.json").read_text())

    def go(edit):
        m = json.loads(json.dumps(ref))
        edit(m)
        f = tmp / "ka.json"
        f.write_text(json.dumps(m))
        return run([PY, ka, f]).returncode

    def good(m):
        m["config"]["bulk_sharpen_pt"] = dict(exp)
        m["config"]["de_sign_emitter_sha256"] = "0f56209d" + "0" * 56
        m["output"]["sha256"] = "1" * 64
    rc_good = go(good)
    rc_same = go(lambda m: (good(m), m["output"].__setitem__("sha256", ref["output"]["sha256"])))
    rc_moved = go(lambda m: (good(m), m["config"]["bulk_sharpen_pt"]["moved_vs_tercile"].pop("E2F4")))
    rc_extra = go(lambda m: (good(m), m["config"].__setitem__("bulk_sharpen_q", 1.75)))
    rc_ref = go(lambda m: None)
    rec("T13 ka598.py: PASS on the expected meta; FAIL on an identical artifact, a missing move, an extra config diff, 479p2g itself",
        rc_good == 0 and rc_same == 1 and rc_moved == 1 and rc_extra == 1 and rc_ref == 1,
        f"rc {rc_good}/{rc_same}/{rc_moved}/{rc_extra}/{rc_ref}")


def main() -> int:
    W.mkdir(parents=True, exist_ok=True)
    live_before = {rel: sha(ROOT / rel) for rel in SIX}
    with tempfile.TemporaryDirectory(prefix="t598-") as td:
        tmp = Path(td)
        canon = tmp / "canon"
        canon_tree(canon)
        b0 = sha(canon / "scripts/make_submission_pm_cd4.py")
        rec("T0 canonical tree: builder == HEAD's (the bundle does not touch it), def_day.sh carries 554's mark",
            b0 == sha(ROOT / "scripts/make_submission_pm_cd4.py") and "554: the K-class emitter" in (canon / "scripts/def_day.sh").read_text(),
            f"builder {b0[:16]}")
        pre = {rel: sha(canon / rel) for rel in SIX}
        for rel in SIX:
            print(f"      canonical {rel} {pre[rel][:16]}")
        # T1 / T2
        q = run([PY, PATCH, "--root", canon, "--check"])
        rec("T1 --check: clean, nothing written", q.returncode == 0 and {r: sha(canon / r) for r in SIX} == pre, q.stdout.strip().splitlines()[-1])
        tree = tmp / "tree"
        copy_tree(canon, tree)
        q = run([PY, PATCH, "--root", tree])
        post = {rel: sha(tree / rel) for rel in SIX}
        rec("T2 apply: rc 0, all six files changed", q.returncode == 0 and all(post[r] != pre[r] for r in SIX))
        for rel in SIX:
            print(f"      patched   {rel} {post[rel][:16]}")
        q2 = run([PY, PATCH, "--root", tree])
        rec("T2 a second apply is refused, nothing written", q2.returncode != 0 and "already patched" in (q2.stderr + q2.stdout)
            and {r: sha(tree / r) for r in SIX} == post)
        part = tmp / "part"
        copy_tree(tree, part)
        shutil.copy(canon / "scripts/def_pipeline.py", part / "scripts/def_pipeline.py")
        q3 = run([PY, PATCH, "--root", part])
        rec("T2 a partial apply is refused", q3.returncode != 0 and "partially patched" in (q3.stderr + q3.stdout))
        raw = tmp / "raw"                     # the canonical tree with a PRE-bundle def_day.sh (b425c8bc's, 0b53ba89)
        copy_tree(canon, raw)
        (raw / "scripts/def_day.sh").write_bytes(subprocess.run(["git", "show", "b425c8bc:scripts/def_day.sh"], cwd=ROOT,
                                                                capture_output=True, check=True).stdout)
        q4 = run([PY, PATCH, "--root", raw, "--check"])
        rec("T2 a pre-bundle tree (HEAD's def_day.sh, no 554 mark) is refused", q4.returncode != 0 and "554" in (q4.stderr + q4.stdout))
        # T3
        okc = True
        for rel in SIX:
            if rel.endswith(".py"):
                try:
                    py_compile.compile(str(tree / rel), cfile=str(tmp / "pyc" / (Path(rel).name + "c")), doraise=True)
                except py_compile.PyCompileError as e:
                    okc = False
                    print(e)
        bn = subprocess.run(["/bin/bash", "-n", str(tree / "scripts/def_day.sh")], capture_output=True, text=True)
        rec("T3 py_compile x5 and bash -n def_day.sh", okc and bn.returncode == 0, bn.stderr.strip())
        # T4
        pin = _load(ROOT / "scripts/def_pin_panel_artifacts.py", "pin598")
        new_text = (tree / "scripts/make_submission_pm_cd4.py").read_text()
        old_text = (canon / "scripts/make_submission_pm_cd4.py").read_text()
        args = ("union_pruned_def", "aa" * 32, "bb" * 32, {}, ("lodo_lo0.5_pct_def", "x/pert_scale_lodo_def_lo0.5_pct.npz", "cc" * 32))
        pinned_new = pin.rewrite(new_text, *args)
        mod = _load(PATCH, "patch598")
        pinned_old_then_patched = pinned_old = pin.rewrite(old_text, *args)
        for a, b in mod.edits()[mod.BUILDER]:
            pinned_old_then_patched = pinned_old_then_patched.replace(a, b)
        rec("T4 the day-of pin commutes with the patch (pin(patch(b)) == patch(pin(b)))",
            pinned_new == pinned_old_then_patched and pinned_new.count('"lodo_lo0.5_pct_def": (') == 1)
        old_ns, new_ns = module_ns(old_text), module_ns(new_text)
        rec("T4b the new builder defines dump_key_carried, the old one does not",
            "dump_key_carried" in new_ns and "dump_key_carried" not in old_ns)
        t5(old_ns, new_ns)
        info = t6(old_ns, new_ns)
        t7(old_text, new_text)
        t8(tree, canon)
        t9(tree, canon, tmp)
        t10(tree, canon)
        t11(tree, canon, tmp)
        t12(tree, tmp)
        t13(tmp)
        out = {"bundle": BUNDLE_STATE, "canonical": pre, "patched": post, "head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                                                           capture_output=True, text=True).stdout.strip(),
               "abc": info}
        (W / "test598.json").write_text(json.dumps(out, indent=1) + "\n")
        shutil.copy(tree / "scripts/make_submission_pm_cd4.py", W / "builder.post598.py")
    live_after = {rel: sha(ROOT / rel) for rel in SIX}
    rec("live tree untouched", live_before == live_after)
    n_ok = sum(ok for _, ok, _ in RESULTS)
    print(f"{n_ok}/{len(RESULTS)} PASS" if n_ok == len(RESULTS) else f"{n_ok}/{len(RESULTS)} PASS -- FAILED: "
          + "; ".join(n for n, ok, _ in RESULTS if not ok))
    return 0 if n_ok == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
