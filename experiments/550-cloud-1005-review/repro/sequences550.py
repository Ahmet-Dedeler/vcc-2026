#!/usr/bin/env python3
"""550: both 10-05 apply sequences, decide-survival, callblock equivalence, census preseed.

Temp copies only. Does not write the live nine files. Does not run vcc.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FILES = [
    Path("scripts/def_day.sh"),
    Path("scripts/def_day_helpers.py"),
    Path("experiments/448-def-rehearsal-q175/relevel_iter.py"),
    Path("scripts/def_sameline.py"),
    Path("scripts/october_recipe.py"),
    Path("docs/october-decisions.json"),
    Path("experiments/524-october-recipe-dr/gate524.py"),
    Path("scripts/def_postbuild.py"),
    Path("scripts/def_pipeline.py"),
]
APPLY = ROOT / "experiments/538-cloud-bundle-extras/apply_october538.py"
P530 = ROOT / "experiments/530-cloud-october-bundle/patch_bundle530.py"
P528 = ROOT / "experiments/528-callprop-block/patch_defday528.py"
PY = os.environ.get("SEQ550_PYTHON", sys.executable)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def live_hashes() -> dict[str, str]:
    return {str(rel): sha256(ROOT / rel) for rel in FILES}


def stage_nine(dest: Path, decisions_text: str | None = None) -> None:
    for rel in FILES:
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if rel.name == "october-decisions.json" and decisions_text is not None:
            target.write_text(decisions_text)
        else:
            shutil.copy(ROOT / rel, target)


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main() -> int:
    before = live_hashes()
    section("A. commits since abcb541d touching the nine files or their callees")
    r = run(["git", "log", "--since=2026-09-30 11:40:00 +0000", "--format=%h %ci %s", "--"] + [str(p) for p in FILES], cwd=ROOT)
    print("nine-file commits since 11:40Z:")
    print(r.stdout.strip() or "(none)")
    callees = [
        "scripts/def_pin_panel_artifacts.py",
        "scripts/def_lodo_key.py",
        "scripts/def_census.py",
        "scripts/make_submission_pm_cd4.py",
        "src/vcc2026/rc_emit_k.py",
        "experiments/528-callprop-block/callblock528.py",
        "experiments/528-callprop-block/feat528.py",
        "scripts/def_preflight.py",
        "scripts/build_submission.py",
    ]
    r = run(["git", "log", "--since=2026-09-30 11:40:00 +0000", "--format=%h %ci %s", "--"] + callees, cwd=ROOT)
    print("callee commits since 11:40Z:")
    print(r.stdout.strip() or "(none)")
    r = run(["git", "diff", "--stat", "abcb541d^..HEAD", "--"] + [str(p) for p in FILES] + callees, cwd=ROOT)
    print("diff stat abcb541d^..HEAD (nine + callees):")
    print(r.stdout.strip() or "(empty)")

    section("A. apply_october538.py --check on a temp copy of today's nine files")
    check_root = Path(tempfile.mkdtemp(prefix="550check-"))
    try:
        stage_nine(check_root)
        r = run([PY, str(APPLY), "--root", str(check_root), "--check"])
        print(r.stdout)
        print(r.stderr)
        print("check rc", r.returncode)
        if r.returncode != 0 or "patches apply cleanly (not written)" not in r.stdout:
            return 1
        # destination of --check is the temp root, which MUST change only inside the stage
        # the applier's --check must leave the temp root byte-identical
        after_check = {str(rel): sha256(check_root / rel) for rel in FILES}
        staged_before = {str(rel): sha256(ROOT / rel) for rel in FILES}
        if after_check != staged_before:
            print("FAIL: --check wrote the temp root")
            return 1
        print("temp root byte-identical after --check")
    finally:
        shutil.rmtree(check_root, ignore_errors=True)

    section("C. WRONG sequence: patch_bundle530.py then patch_defday528.py")
    wrong = Path(tempfile.mkdtemp(prefix="550wrong-"))
    try:
        stage_nine(wrong)
        day = wrong / "scripts/def_day.sh"
        pre = sha256(day)
        r530 = run([PY, str(P530), "--root", str(wrong)])
        print("530 rc", r530.returncode)
        print((r530.stdout or "")[-400:])
        print((r530.stderr or "")[-400:])
        if r530.returncode != 0:
            return 1
        mid = day.read_bytes()
        r528 = run([PY, str(P528), str(day)])
        text = (r528.stdout or "") + (r528.stderr or "")
        print("528 rc", r528.returncode)
        print(text.strip())
        print("def_day bytes unchanged by 528:", day.read_bytes() == mid)
        print("DEF_DAY_CALLSET_BLOCK in file:", "DEF_DAY_CALLSET_BLOCK" in day.read_text())
        print("anchor-miss text present:", "anchor found" in text)
        if r528.returncode == 0 or day.read_bytes() != mid or "DEF_DAY_CALLSET_BLOCK" in day.read_text():
            print("FAIL: wrong sequence did not miss cleanly")
            return 1
        if "anchor found" not in text:
            print("FAIL: expected the anchor-miss message")
            return 1
        print("pre-530 sha", pre)
        print("post-530 sha", sha256(day))
    finally:
        shutil.rmtree(wrong, ignore_errors=True)

    section("C. RIGHT sequence: apply_october538.py (530, 528 miss, 538)")
    right = Path(tempfile.mkdtemp(prefix="550right-"))
    mutated: Path | None = None
    try:
        stage_nine(right)
        r = run([PY, str(APPLY), "--root", str(right)])
        print(r.stdout)
        print(r.stderr)
        print("apply rc", r.returncode)
        if r.returncode != 0:
            return 1
        day_txt = (right / "scripts/def_day.sh").read_text()
        for needle in (
            "DEF_DAY_CALLSET_BLOCK",
            'callblock528.py --bundle "$B" --contexts "$CTXS" --out "$W/callblock_def.npz"',
            'export VCC_RC_CALLSET_BLOCK="$CSB"',
            'export VCC_DE_MDE_SIGN="$MDS"',
            "DONE. Next",
            "DEF_DAY_CENSUS_MAX_ZERO",
            'REPO="$(CDPATH=',
        ):
            print(f"contains {needle!r}:", needle in day_txt)
        if 'callblock528.py --bundle "$B" --contexts "$CTXS" --out "$W/callblock_def.npz"' not in day_txt:
            return 1
        # prop without a suffix must not append :swap
        if 'case "$DEF_DAY_CALLSET_BLOCK" in prop:*) CSB="$CSB:${DEF_DAY_CALLSET_BLOCK#prop:}"' not in day_txt:
            print("FAIL: missing the prop:* suffix append")
            return 1
        r2 = run([PY, str(APPLY), "--root", str(right)])
        print("second apply rc", r2.returncode)
        print((r2.stdout or "") + (r2.stderr or ""))
        if r2.returncode == 0 or "already patched" not in ((r2.stdout or "") + (r2.stderr or "")):
            print("FAIL: second apply was not refused")
            return 1
        # 528's own patch on the 538 tree must not rewrite
        day = right / "scripts/def_day.sh"
        before528 = day.read_bytes()
        r528 = run([PY, str(P528), str(day)])
        text = (r528.stdout or "") + (r528.stderr or "")
        print("528-after-538 rc", r528.returncode)
        print(text.strip()[:500])
        print("bytes unchanged:", day.read_bytes() == before528)
        if day.read_bytes() != before528:
            return 1
        bn = run(["bash", "-n", str(day)])
        print("bash -n rc", bn.returncode, bn.stderr.strip())
        if bn.returncode != 0:
            return 1

        section("C. morning decide of mde_sign survives the JSON patch")
        mutated = Path(tempfile.mkdtemp(prefix="550mut-"))
        raw = json.loads((ROOT / "docs/october-decisions.json").read_text())
        raw["switches"]["mde_sign"]["status"] = "DECIDED"
        raw["switches"]["mde_sign"]["value"] = "clean4:0.05"
        raw["switches"]["mde_sign"]["evidence"] = "550 synthetic morning decide"
        raw["switches"]["de_mag"]["evidence"] = "550 keep this evidence"
        stage_nine(mutated, decisions_text=json.dumps(raw, indent=1) + "\n")
        r = run([PY, str(APPLY), "--root", str(mutated)])
        print("mutated apply rc", r.returncode)
        if r.returncode != 0:
            print(r.stdout)
            print(r.stderr)
            return 1
        got = json.loads((mutated / "docs/october-decisions.json").read_text())
        sw = got["switches"]
        print("mde_sign", sw["mde_sign"]["status"], sw["mde_sign"]["value"], sw["mde_sign"].get("evidence"))
        print("de_mag evidence", sw["de_mag"].get("evidence"))
        print("lodo_hi", sw["lodo_hi"]["status"], sw["lodo_hi"]["value"])
        print("callset_block", sw["callset_block"]["status"], sw["callset_block"]["value"])
        print("cd4_stack", sw["cd4_stack"]["status"], sw["cd4_stack"]["value"])
        keys = list(sw)
        print("lodo_hi index", keys.index("lodo_hi"), "callset_block index", keys.index("callset_block"),
              "cd4_thresholds index", keys.index("cd4_thresholds"))
        if sw["mde_sign"]["status"] != "DECIDED" or sw["mde_sign"]["value"] != "clean4:0.05":
            print("FAIL: morning mde_sign was not preserved")
            return 1
        if sw["de_mag"].get("evidence") != "550 keep this evidence":
            print("FAIL: unrelated evidence rewritten")
            return 1
        if sw["lodo_hi"]["status"] != "PENDING" or sw["callset_block"]["status"] != "PENDING":
            return 1
        if not (keys.index("lodo_hi") < keys.index("cd4_thresholds") and keys.index("callset_block") < keys.index("cd4_thresholds")):
            print("FAIL: new switches not inserted before cd4_thresholds")
            return 1

        section("C. patched october_recipe.py check, then decide on the TEMP tree only")
        env = os.environ.copy()
        env.update({
            "OCTOBER_DECISIONS": str(right / "docs/october-decisions.json"),
            "OCTOBER_DEF_DAY": str(right / "scripts/def_day.sh"),
            "OCTOBER_BUILDER": str(ROOT / "scripts/make_submission_pm_cd4.py"),
            "OCTOBER_EMITTER": str(ROOT / "src/vcc2026/rc_emit_k.py"),
            "OCTOBER_CALLBLOCK": str(ROOT / "experiments/528-callprop-block/callblock528.py"),
            "OCTOBER_SAMELINE": str(right / "scripts/def_sameline.py"),
        })
        # unknown switch on the UNPATCHED resolver, writing only a temp decisions file
        unk_dec = Path(tempfile.mkdtemp(prefix="550unk-")) / "october-decisions.json"
        shutil.copy(ROOT / "docs/october-decisions.json", unk_dec)
        before_unk = unk_dec.read_bytes()
        uenv = os.environ.copy()
        uenv["OCTOBER_DECISIONS"] = str(unk_dec)
        u = run([PY, str(ROOT / "scripts/october_recipe.py"), "decide", "callset_block", "prop",
                 "--evidence", "550 must not write a new switch into the unpatched file"], env=uenv)
        print("unpatched decide callset_block rc", u.returncode)
        print((u.stdout or "") + (u.stderr or ""))
        print("unpatched decisions unchanged:", unk_dec.read_bytes() == before_unk)
        if u.returncode == 0 or unk_dec.read_bytes() != before_unk:
            print("FAIL: unpatched decide accepted callset_block or wrote the file")
            return 1

        chk = run([PY, str(right / "scripts/october_recipe.py"), "check", "--allow-pending"], env=env)
        print("patched check rc", chk.returncode)
        print(chk.stdout)
        print(chk.stderr)
        if chk.returncode != 0:
            return 1
        d1 = run([PY, str(right / "scripts/october_recipe.py"), "decide", "callset_block", "prop",
                  "--evidence", "550 temp only"], env=env)
        d2 = run([PY, str(right / "scripts/october_recipe.py"), "decide", "mde_sign", "clean4:0.05",
                  "--evidence", "550 temp only"], env=env)
        print(d1.stdout, d1.stderr, "rc", d1.returncode)
        print(d2.stdout, d2.stderr, "rc", d2.returncode)
        if d1.returncode or d2.returncode:
            return 1
        bare = run([PY, str(right / "scripts/october_recipe.py"), "command",
                    "--bundle", "data/challenge", "--work", "W", "--name", "rider", "--json"], env=env)
        print("command WITHOUT --allow-pending rc", bare.returncode)
        print(bare.stderr.strip())
        if bare.returncode != 2:
            print("FAIL: expected exit 2 while other switches stay PENDING")
            return 1
        cmd = run([PY, str(right / "scripts/october_recipe.py"), "command", "--allow-pending",
                   "--bundle", "data/challenge", "--work", "W", "--name", "rider", "--json"], env=env)
        print("command --allow-pending rc", cmd.returncode)
        print(cmd.stdout)
        print(cmd.stderr)
        if cmd.returncode != 0:
            return 1
        payload = json.loads(cmd.stdout)
        print("resolved env", payload["env"])
        if payload["env"].get("DEF_DAY_CALLSET_BLOCK") != "prop":
            return 1
        if payload["env"].get("DEF_DAY_MDE_SIGN") != "clean4:0.05":
            return 1
        # live decisions file must be untouched
        if sha256(ROOT / "docs/october-decisions.json") != before["docs/october-decisions.json"]:
            print("FAIL: live october-decisions.json changed")
            return 1

        section("F. census knob dies on a pre-seeded work dir")
        census = r"""
set -uo pipefail
W="$1"
die(){ echo "STOP: $*"; exit 1; }
if [ -n "${DEF_DAY_CENSUS_MAX_ZERO:-}" ]; then
  python3 -c 'import sys; n=sys.argv[1]; sys.exit(0 if n.isdigit() and not (len(n)>1 and n[0]=="0") and 6<=int(n)<=30 else 1)' "$DEF_DAY_CENSUS_MAX_ZERO" 2>/dev/null \
    || die "range"
  if [ -d "$W" ] && [ -n "$(ls -A "$W" 2>/dev/null)" ]; then
    die "fresh work dir ($W is not empty)"
  fi
fi
echo CENSUS_OK
"""
        preseed = Path(tempfile.mkdtemp(prefix="550pre-"))
        (preseed / "whitener_union_pruned_def.npz").write_bytes(b"not a real whitener")
        (preseed / "whitener_union_pruned_def.json").write_text("{}\n")
        empty = Path(tempfile.mkdtemp(prefix="550empty-"))
        c1 = run(["bash", "-c", census, "census", str(preseed)], env={**os.environ, "DEF_DAY_CENSUS_MAX_ZERO": "8"})
        print("preseed rc", c1.returncode, (c1.stdout or "") + (c1.stderr or ""))
        c2 = run(["bash", "-c", census, "census", str(empty)], env={**os.environ, "DEF_DAY_CENSUS_MAX_ZERO": "8"})
        print("empty rc", c2.returncode, (c2.stdout or c2.stderr or "").strip())
        c3 = run(["bash", "-c", census, "census", str(preseed)], env={k: v for k, v in os.environ.items() if k != "DEF_DAY_CENSUS_MAX_ZERO"})
        # parent env may not have the knob; the script checks unset
        env3 = os.environ.copy()
        env3.pop("DEF_DAY_CENSUS_MAX_ZERO", None)
        c3 = run(["bash", "-c", census, "census", str(preseed)], env=env3)
        print("unset-on-preseed rc", c3.returncode, (c3.stdout or "").strip())
        if c1.returncode == 0 or "not empty" not in (c1.stdout + c1.stderr):
            print("FAIL: preseed did not die")
            return 1
        if c2.returncode != 0 or c3.returncode != 0:
            print("FAIL: empty or unset path died")
            return 1
        shutil.rmtree(preseed, ignore_errors=True)
        shutil.rmtree(empty, ignore_errors=True)
    finally:
        shutil.rmtree(right, ignore_errors=True)
        if mutated is not None:
            shutil.rmtree(mutated, ignore_errors=True)

    section("D. callblock528.main twice on synthetic features: same npz bytes")
    import numpy as np
    sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
    import callblock528 as C

    g = 12
    rng = np.random.default_rng(0)
    mean = rng.random(g) * 80 + 1
    x = rng.random(g) * 0.02
    keep = mean > 5
    genes = np.array([f"G{i}" for i in range(g)])

    def fake(_path):
        return dict(genes=genes, mean_cpm=mean, x=x, keep=keep, n_ctrl=np.int64(40))

    C.features = fake
    bundle = Path(tempfile.mkdtemp(prefix="550bundle-"))
    (bundle / "gene_names.csv").write_text("gene_name\n" + "\n".join(genes.tolist()) + "\n")
    (bundle / "context_A.h5ad").write_bytes(b"synthetic-A")
    (bundle / "context_B.h5ad").write_bytes(b"synthetic-B")
    out1 = bundle / "k1.npz"
    out2 = bundle / "k2.npz"
    old = sys.argv
    try:
        sys.argv = ["callblock528.py", "--bundle", str(bundle), "--contexts", "A,B", "--out", str(out1)]
        C.main()
        sys.argv = ["callblock528.py", "--bundle", str(bundle), "--contexts", "A,B", "--out", str(out2)]
        C.main()
    finally:
        sys.argv = old
    b1, b2 = out1.read_bytes(), out2.read_bytes()
    print("npz bytes equal", b1 == b2, "len", len(b1), "sha", hashlib.sha256(b1).hexdigest())
    z1, z2 = np.load(out1), np.load(out2)
    print("key_A equal", np.array_equal(z1["key_A"], z2["key_A"]))
    # same function the artifact builder calls: key_from_features(features(...))
    F = fake(bundle / "context_A.h5ad")
    k, med = C.key_from_features(F, C.MEDIAN_X_H1, beta=C.BETA)
    print("key matches key_from_features", np.allclose(z1["key_A"], k), "median_x", med, "beta", C.BETA)
    if b1 != b2 or not np.array_equal(z1["key_A"], k):
        return 1

    section("F. feat528 empty library bins -> NaN mean")
    L = np.array([5000.0] * 20)
    nb = 10
    edges = np.quantile(L, np.linspace(0, 1, nb + 1))
    b = np.clip(np.searchsorted(edges, L, side="right") - 1, 0, nb - 1)
    Lb = np.array([L[b == k].mean() if np.any(b == k) else float("nan") for k in range(nb)])
    # the production line is L[b == k].mean() with no empty guard; numpy mean of empty is NaN and warns
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Lb_prod = np.array([L[b == k].mean() for k in range(nb)])
    print("edges", edges.tolist())
    print("occupied bins", sorted(set(b.tolist())), "prod Lb", Lb_prod.tolist())
    print("nonfinite bins", int(np.sum(~np.isfinite(Lb_prod))))
    if np.all(np.isfinite(Lb_prod)):
        print("FAIL: expected a non-finite bin mean on a degenerate library")
        return 1

    section("D. builder parse of prop:path:sha (no swap, no steer)")
    sha = "ab" * 32
    env_s = f"prop:/tmp/callblock_def.npz:{sha}"
    _cb = [s.strip() for s in env_s.split(":")]
    # path may contain no extra colons here; the production split is the same
    assert _cb[0] == "prop" and len(_cb) == 3 and len(_cb[2]) == 64
    swap = float(_cb[3]) if len(_cb) >= 4 else 1.0
    steer = _cb[4].lower() if len(_cb) == 5 else "kept"
    cpm = 2200
    block_keep = int(cpm - int(round(swap * cpm)))
    print("swap", swap, "steer", steer, "block_keep", block_keep)
    if swap != 1.0 or steer != "kept" or block_keep != 0:
        return 1
    # source lines that this duplicates
    src = (ROOT / "scripts/make_submission_pm_cd4.py").read_text().splitlines()
    window = "\n".join(src[2020:2030])
    print("builder window 2021-2030:\n" + window)
    if "else 1.0" not in window or 'else "kept"' not in window:
        return 1

    section("D. pipeline resolve_env fills MDE defaults under a parent MDE_SIGN export")
    code = r"""
import os, sys
sys.path.insert(0, "scripts")
os.environ.pop("VCC_DE_MDE", None)
os.environ["VCC_DE_MDE_SIGN"] = "clean4:0.05"
os.environ["VCC_RC_CALLSET_BLOCK"] = "prop:/tmp/callblock_def.npz:%s" % ("ab" * 32)
import def_pipeline as P
resolved, notes = P.resolve_env()
# Runner.run builds the child env as {**os.environ, **resolved} (def_pipeline.py:511-513).
# resolve_env only returns SHIPPED_ENV keys, so a parent export that is not a default survives
# the merge and a missing default is filled.
child = {**os.environ, **resolved}
for k in ("VCC_DE_MDE", "VCC_DE_ZERO_REPAIR", "VCC_DE_SIGN_PRIOR", "VCC_DE_MDE_SIGN", "VCC_RC_CALLSET_BLOCK"):
    print(f"resolved {k}={resolved.get(k)!r}")
    print(f"child {k}={child.get(k)}")
"""
    r = run([PY, "-c", code], cwd=ROOT)
    print(r.stdout)
    print(r.stderr)
    if r.returncode != 0:
        return 1
    if "VCC_DE_MDE=0.1,0.1,0.3,coarse" not in r.stdout:
        return 1
    if "VCC_DE_MDE_SIGN=clean4:0.05" not in r.stdout:
        return 1
    if "VCC_DE_SIGN_PRIOR=clean4q80" not in r.stdout:
        return 1

    section("E. EXPECT is a set (order does not matter)")
    snippet = r"""
want=set(x for x in "de_mde_sign,rc_callset_block".split(",") if x)
got=set(["rc_callset_block", "de_mde_sign"])
print("set equal", got == want)
print("string equal", "rc_callset_block,de_mde_sign" == "de_mde_sign,rc_callset_block")
"""
    # mirror run_arm_generic.sh lines 81-82
    r = run([PY, "-c", snippet])
    print(r.stdout)
    if "set equal True" not in r.stdout or "string equal False" not in r.stdout:
        return 1

    section("F. chain502 clobber-before-wait, on a temp tree, with the 10-05 clock")
    sim = Path(tempfile.mkdtemp(prefix="550chain-"))
    (sim / "scripts").mkdir()
    (sim / "work").mkdir()
    live_day = sim / "scripts/def_day.sh"
    live_day.write_text("PATCHED-OCTOBER\n")
    (sim / "scripts/def_day_helpers.py").write_text("PATCHED-HELPERS\n")
    (sim / "work/pre502.def_day.sh").write_text("PRE502-DAY\n")
    (sim / "work/pre502.def_day_helpers.py").write_text("PRE502-HELPERS\n")
    sim_sh = r"""
set -uo pipefail
DD="$1"; HP="$2"; PRED="$3"; PREH="$4"; NOW="$5"
DDSHA="$(sha256sum "$DD" | cut -c1-16)"; HPSHA="$(sha256sum "$HP" | cut -c1-16)"
case "$DDSHA $HPSHA" in
  "e29e793a7964df03 60b2a7771ef857d4"|"610db78277239f50 20142142e4574b7f") echo "kept";;
  *) cp "$PRED" "$DD"; cp "$PREH" "$HP"; echo "reverted (was $DDSHA $HPSHA)";;
esac
if [ "$NOW" -ge 202610031200 ]; then echo "cap hit, exit 1"; exit 1; fi
echo "would wait"
"""
    r = run(["bash", "-c", sim_sh, "sim", str(live_day), str(sim / "scripts/def_day_helpers.py"),
             str(sim / "work/pre502.def_day.sh"), str(sim / "work/pre502.def_day_helpers.py"), "202610050400"])
    print(r.stdout, r.stderr, "rc", r.returncode)
    print("day after", live_day.read_text().strip())
    if r.returncode == 0 or live_day.read_text().strip() != "PRE502-DAY":
        print("FAIL: clobber simulation did not revert then exit 1")
        return 1
    print("fix526 cap: 202610050400 >= 202610021200", 202610050400 >= 202610021200)

    section("F. partial extras refuses a second apply")
    partial = Path(tempfile.mkdtemp(prefix="550part-"))
    stage_nine(partial)
    # borrow one patched def_day by applying on a sibling and copying only def_day.sh
    donor = Path(tempfile.mkdtemp(prefix="550donor-"))
    stage_nine(donor)
    r = run([PY, str(APPLY), "--root", str(donor)])
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        return 1
    shutil.copy(donor / "scripts/def_day.sh", partial / "scripts/def_day.sh")
    r = run([PY, str(APPLY), "--root", str(partial), "--check"])
    text = (r.stdout or "") + (r.stderr or "")
    print("partial rc", r.returncode)
    print(text.strip()[:400])
    if r.returncode == 0 or "partial extras" not in text:
        print("FAIL: partial tree was not refused")
        return 1
    shutil.rmtree(partial, ignore_errors=True)
    shutil.rmtree(donor, ignore_errors=True)
    shutil.rmtree(sim, ignore_errors=True)
    shutil.rmtree(bundle, ignore_errors=True)

    section("live nine files unchanged")
    after = live_hashes()
    for rel, digest in before.items():
        same = after[rel] == digest
        print(rel, "unchanged" if same else "CHANGED", digest[:16])
        if not same:
            return 1
    print("ALL SEQUENCES PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
