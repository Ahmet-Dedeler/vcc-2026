"""644: run every command the DEF-RUNBOOK amendments cite, in --help / --check / --verify / dry-run form, against scratch
POST-BUNDLE trees (copies; the live repo is never written):
  T0 = 538 -> 554 -> 590              (the 10-05 bundle; 641's --check runs here)
  T1 = T0 + 641                       (step 6b; 642's --check runs here)
  T2 = T1 + 642                       (10-06; the day-of commands run here)
After the 10-05 / 10-06 applies, patches already on HEAD are skipped when the trees are built, and a patch --check row
for one reads "already on HEAD" (606's P1 rule; a PARTIALLY patched tree still fails).
Post-538, def_day.sh derives REPO from its own location (pre-538 it hard-coded the live repo), so in a scratch tree it acts
on that tree; here it is only syntax-checked and run without arguments (it exits at ${1:?bundle}). 644's skip patch is
checked on T2 and applied for T3 = T2 + skip; its own behaviour is test_skip644.py (538's stub harness). Writes work/verify644.json and work/verify644.md.

    uv run python experiments/644-runbook-amendments/verify644.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/554-dose-guard"))
import test_guard554 as G  # noqa: E402

PY = str(ROOT / ".venv/bin/python")


def build(dest: Path, upto: str) -> Path:
    for rel in G.REL538 + G.EXTRA + [Path("scripts/def_identity.py")]:
        t = dest / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, t)
    steps = [[PY, str(G.APPLIER538), "--root", str(dest)], [PY, str(G.PATCH), "--root", str(dest)],
             [PY, str(ROOT / "experiments/590-identity-validation/patch_identity590.py"), "--root", str(dest), "--apply"]]
    if upto in ("641", "642"):
        steps.append([PY, str(ROOT / "experiments/641-depth-row/patch_depth641.py"), "--root", str(dest)])
    for cmd in steps:
        if already(cmd):        # after the 10-05 / 10-06 applies the live files carry these patches (606's P1 rule)
            continue
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, (cmd, r.stdout[-1500:], r.stderr[-1500:])
    for d in ("data", "submissions"):
        os.symlink(ROOT / d, dest / d)
    (dest / "experiments").mkdir(exist_ok=True)
    for e in sorted((ROOT / "experiments").iterdir()):
        if not (dest / "experiments" / e.name).exists():
            os.symlink(e, dest / "experiments" / e.name)
    os.symlink(ROOT / ".venv", dest / ".venv")
    if upto == "642":
        cmd = [PY, str(ROOT / "experiments/642-hek293-contingency/patch_hek293_642.py"), "--root", str(dest)]
        if not already(cmd):
            r = subprocess.run(cmd, capture_output=True, text=True)
            assert r.returncode == 0, r.stdout + r.stderr
    return dest


def already(cmd: list[str]) -> bool:
    """606's P1 rule: a patch whose --check says it is ALREADY in the tree is skipped; a PARTIALLY patched tree is not."""
    c = subprocess.run([x for x in cmd if x != "--apply"] + ["--check"], capture_output=True, text=True)
    o = (c.stdout + c.stderr).lower()
    return c.returncode != 0 and "already" in o and "partially" not in o


def clean_env(**kw):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    env.update(kw)
    return env


def main():
    t0 = time.time()
    tmp = Path(tempfile.mkdtemp(prefix="verify644-"))
    T0, T1, T2 = build(tmp / "T0", "590"), build(tmp / "T1", "641"), build(tmp / "T2", "642")
    T3 = build(tmp / "T3", "642")
    _sk = [PY, str(HERE / "patch_skip644.py"), "--root", str(T3)]
    if not already(_sk):
        r = subprocess.run(_sk, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
    rows = []

    def check(cid, tree, cmd, expect_rc=0, must=(), mustnot=(), env=None, shell=False, note="", patch_check=False):
        p = subprocess.run(cmd, cwd=tree, env=env or clean_env(), capture_output=True, text=True, shell=shell)
        out = p.stdout + p.stderr
        ok = (expect_rc is None or p.returncode == expect_rc) and all(m in out for m in must) and not any(m in out for m in mustnot)
        if patch_check and not ok and p.returncode != 0 and "already" in out.lower() and "partially" not in out.lower():
            ok, note = True, "already on HEAD (post-apply): nothing left to check"
        rows.append(dict(id=cid, tree=(tree.name if tree != ROOT else "live"), cmd=cmd if isinstance(cmd, str) else " ".join(map(str, cmd)).replace(str(tmp), "<tmp>"),
                         rc=p.returncode, expect_rc=expect_rc, ok=ok, note=note,
                         tail=[ln for ln in out.strip().splitlines()[-3:]]))
        print(f"{'PASS' if ok else 'FAIL'} {cid} [{tree.name}] rc {p.returncode}", flush=True)

    py = PY
    # --- patch application (the 10-05 / 10-06 checklist) ---------------------------------------------------------
    check("A1 641 --check after 538+554+590", T0, [py, str(ROOT / "experiments/641-depth-row/patch_depth641.py"), "--root", ".", "--check"],
          must=("apply cleanly",), patch_check=True)
    check("A2 642 --check after 641", T1, [py, str(ROOT / "experiments/642-hek293-contingency/patch_hek293_642.py"), "--root", ".", "--check"],
          must=("apply cleanly",), patch_check=True)
    check("A3 642 refuses a second apply", T2, [py, str(ROOT / "experiments/642-hek293-contingency/patch_hek293_642.py"), "--root", "."],
          expect_rc=1, must=("refusing a second apply",))
    check("A4 644 skip --check after 642", T2, [py, str(HERE / "patch_skip644.py"), "--root", ".", "--check"], must=("apply cleanly",), patch_check=True)
    check("A5 644 skip: def_day.sh syntax after the patch", T3, ["bash", "-n", "scripts/def_day.sh"])
    check("A6 644 skip: the knob is logged at start", T3, ["grep", "-c", "-F", 'log "644: DEF_DAY_SAMELINE_SKIP=', "scripts/def_day.sh"], must=("1",))
    # --- def_day.sh -------------------------------------------------------------------------------------------------
    check("D1 def_day.sh syntax", T2, ["bash", "-n", "scripts/def_day.sh"])
    for k, pat in (("WIDE_LIBS", "log \"554: DEF_DAY_ALLOW_WIDE_LIBS"), ("LOW_DOSE", "log \"554: DEF_DAY_ALLOW_LOW_DOSE"),
                   ("DEPTH", "log \"641: the depth row"), ("NONCCLE", "log \"642: the HEK293 contingency")):
        check(f"D2 def_day.sh logs the {k} override at start", T2, ["grep", "-c", "-F", pat, "scripts/def_day.sh"], must=("1",))
    check("D3 def_day.sh without arguments stops at the bundle argument", T2, ["bash", "scripts/def_day.sh"], expect_rc=1, must=("bundle",))
    # --- preflight (554 spread + 641 depth) --------------------------------------------------------------------------
    check("P1 def_preflight --verify (A/B/C: both rows silent)", T2, [py, "scripts/def_preflight.py", "--verify"],
          must=("VERIFY PASS",), mustnot=("q75/q25 <=", "within 15%"))
    check("P2 def_preflight --help", T2, [py, "scripts/def_preflight.py", "--help"])
    # --- identity (590 + 642) ------------------------------------------------------------------------------------------
    check("I1 def_identity --verify (prints 590's min top-1)", T2, [py, "scripts/def_identity.py", "--verify"],
          must=("lowest top-1 correlation", "PASS"), mustnot=("[642]",))
    idj = tmp / "identity_hek.json"
    check("I2 def_identity with DEF_DAY_ALLOW_NONCCLE_REF=hek293 (642 stage 2, A/B/C not called)", T2,
          [py, "scripts/def_identity.py", "--bundle", "data/challenge", "--contexts", "A,B,C", "--out", str(idj)],
          env=clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293"), must=("[642] A:", "[642] B:", "[642] C:", "not called"), mustnot=("CALLED HEK293",))
    check("I3 def_identity refuses an unknown reference", T2, [py, "scripts/def_identity.py", "--bundle", "data/challenge", "--contexts", "A"],
          env=clean_env(DEF_DAY_ALLOW_NONCCLE_REF="rpe1"), expect_rc=1, must=("is not one of",))
    # --- same-line (590 floor + 642 gate) ------------------------------------------------------------------------------
    check("S1 def_sameline --help shows --min-top1", T2, [py, "scripts/def_sameline.py", "--help"], must=("--min-top1",))
    slo = tmp / "sameline_hek.npz"
    check("S2 def_sameline on the stage-2 identity file (no HEK293 call -> A on via Nadig Jurkat, B/C none)", T2,
          [py, "scripts/def_sameline.py", "--bundle", "data/challenge", "--identity", str(idj), "--out", str(slo)],
          must=("A: JURKAT", "590 gate passed", "B: HELA_CERVIX -> no registered bank"))
    # --- postbuild (554 dose guard) -------------------------------------------------------------------------------------
    check("B1 def_postbuild --help shows --rc-audit", T2, [py, "scripts/def_postbuild.py", "--help"], must=("--rc-audit",))
    inc = "submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
    m433 = "submissions/433a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse.h5ad.meta.json"
    d433 = "experiments/433-coarse-mde-arm/work/rc_audit_433a.jsonl"
    check("B2 def_postbuild --rc-audit on a healthy dump (433a: no 554 row)", T2,
          [py, "scripts/def_postbuild.py", "--meta", m433, "--bundle", "data/challenge", "--incumbent", inc, "--rc-audit", d433,
           "--out", str(tmp / "post433.json")], expect_rc=None, mustnot=("(554)",), note="rc is 1 for unrelated old-meta rows (554's note); only the 554 row is checked")
    # --- post-apply (10-06 APPLY DAY v2: 598 -> 612 -> 612dd -> 630 -> 630dd -> 642 -> 604 -> 644 v3): the knobs that were
    #     unverifiable before the apply. Run only when the live tree carries the apply (644's mark in def_day.sh).
    if "644: DEF_DAY_SAMELINE_SKIP" in (ROOT / "scripts/def_day.sh").read_text():
        rcw = tmp / "recipe-work"
        check("K1 october_recipe command emits the decided day (tercilec, LODO pct, zs, clean4:0.05)", T2,
              [py, "scripts/october_recipe.py", "command", "--bundle", "data/challenge", "--work", str(rcw), "--name", "verify644",
               "--allow-pending"],
              must=("DEF_DAY_LODO=1", "DEF_DAY_LODO_HI_PCT=96.64697522769474", "DEF_DAY_SHARPEN_PT=tercilec:2.5,2.0,1.75",
                    "DEF_DAY_DE_MAG=dn_k250_zs", "DEF_DAY_MDE_SIGN=clean4:0.05", "bash scripts/def_day.sh data/challenge"))
        for k, pat in (("598 tercilec pass-through", "tercilec"), ("612dd KEYGATE pass-through", "KEYGATE"),
                       ("630dd DE_UPRN pass-through", "DEF_DAY_DE_UPRN"), ("644 skip pass-through", '${SLSKIP:+--skip-contexts "$SLSKIP"}')):
            check(f"K2 def_day.sh carries the {k}", T2, ["grep", "-c", "-F", pat, "scripts/def_day.sh"], mustnot=("\n0\n",),
                  must=())
        check("K3 def_sameline --help shows --skip-contexts (644)", T2, [py, "scripts/def_sameline.py", "--help"], must=("--skip-contexts",))
        check("K4 def_sameline carries 604's identity programme source", T2, ["grep", "-c", "-F", 'a.programme_source == "identity"',
              "scripts/def_sameline.py"])
        for k in ("VCC_DE_MDE_SIGN_KEYGATE", "VCC_DE_UPRN", "tercilec"):
            check(f"K5 the builder carries {k}", T2, ["grep", "-c", "-F", k, "scripts/make_submission_pm_cd4.py"], mustnot=("^0$",))
        rows_k = [r for r in rows if r["id"].startswith(("K2", "K4", "K5"))]
        for r in rows_k:                       # grep -c prints the count: require a positive count
            r["ok"] = r["ok"] and r["rc"] == 0 and r["tail"] and r["tail"][-1].strip().isdigit() and int(r["tail"][-1]) > 0
    # --- recipe, capture, readers ------------------------------------------------------------------------------------------
    check("R1 october_recipe.py show", T2, [py, "scripts/october_recipe.py", "show"])
    check("R2 october_recipe.py check --allow-pending", T2, [py, "scripts/october_recipe.py", "check", "--allow-pending"])
    check("C1 capture_board_entry.py --help (live repo: no bundle patch touches it)", ROOT, [py, "scripts/capture_board_entry.py", "--help"])
    check("F1 531 readdef.py --help", T2, [py, "experiments/531-cloud-final-phase-tooling/readdef.py", "--help"])
    check("F2 637 readdef_factorial637.py (no args prints usage)", T2, [py, "experiments/637-final-phase-factorial/readdef_factorial637.py"],
          expect_rc=None, must=("--factorial",))
    check("F3 637 tests", T2, [py, "experiments/637-final-phase-factorial/test_factorial637.py"], must=("factorial checks passed", "PASS 531's suite still passes"))
    res = dict(seconds=round(time.time() - t0, 1), trees={"T0": "538+554+590", "T1": "+641", "T2": "+642", "T3": "+644 skip", "note": "post-apply: every tree == live (patches skipped as already on HEAD)"}, rows=rows,
               all_ok=all(r["ok"] for r in rows),
               not_verifiable=[
                   "612 KEYGATE / 630 UPRN / 604 identity prior VALUES: the code paths are verified (K2-K5) but none is in the "
                   "recipe yet (their reads are 10-08+); a decided value must be re-checked through K1's command",
                   "the D/E/F bundle itself: `vcc datasets download controls` after Oct 22 (and the CLI upgrade)",
                   "a real def_day.sh run on a bundle (~2 h; the coordinator's rehearsals cover it; 644's skip is stub-tested)",
                   "readdef637.diff dispatch into readdef.main (UNAPPLIED until the final-phase tooling freeze)"])
    (HERE / "work/verify644.json").write_text(json.dumps(res, indent=1))
    md = ["| id | tree | command | rc | ok |", "|---|---|---|---:|---|"]
    md += [f"| {r['id']} | {r['tree']} | `{r['cmd'][:150]}` | {r['rc']} | {'PASS' if r['ok'] else 'FAIL'} |" for r in rows]
    md += ["", "Not verifiable yet:"] + [f"- {x}" for x in res["not_verifiable"]]
    (HERE / "work/verify644.md").write_text("\n".join(md) + "\n")
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"all_ok {res['all_ok']}  ({res['seconds']} s)")


if __name__ == "__main__":
    main()
