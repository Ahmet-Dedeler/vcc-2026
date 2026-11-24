#!/usr/bin/env python3
"""630dd tests (council 600, [4d0207] for [d4cd5f]). Light: no data, no build. Each tree is test598's canonical
post-bundle tree with 641 -> 598 -> 612 -> 630 applied, like 606's chain.
  T1  --check, apply, and a second apply refused; refusals without 598 / without 630 / with another module pin
  T2  612dd and 630dd commute: both orders give the same def_day.sh bytes
  T3  bash -n on every patched def_day.sh
  T4  behaviour of the six insertions, run in a bash harness with stubbed die/log (the inserted text itself, cut from the
      patched file):
      - unset knob: no export, no meta check, no final.env line, no ALLOW change (the KA for "today's command unchanged");
      - validation accepts the builder grammar with K in [0.5, 3] and refuses everything else, and refuses without
        DEF_DAY_DE_MAG;
      - the meta check passes only with the requested variant AND the reviewed module sha.
usage: uv run python experiments/630-rank-neutral-up/test630dd.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/598-tercilec-knob"))
sys.path.insert(0, str(HERE))
import test598 as T598  # noqa: E402
import patch630dd as P  # noqa: E402

PY = sys.executable
PRE = ["experiments/641-depth-row/patch_depth641.py", "experiments/598-tercilec-knob/patch598.py",
       "experiments/612-key-gated-hook/patch612.py", "experiments/630-rank-neutral-up/patch_builder630.py"]
DD612 = "experiments/612-key-gated-hook/patch612dd.py"
DD630 = "experiments/630-rank-neutral-up/patch630dd.py"
RESULTS: list[tuple[bool, str]] = []


def rec(ok: bool, name: str) -> None:
    RESULTS.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(script: str, root: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([PY, str(ROOT / script), "--root", str(root), *extra], capture_output=True, text=True)


def tree(td: Path, name: str, pre=PRE) -> Path:
    t = td / name
    T598.canon_tree(t)
    (t / "data/identity_refs").mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", t / "data/identity_refs/hek293_nonflex_ref642.npz")
    (t / "src/vcc2026").mkdir(parents=True, exist_ok=True)
    for f in (ROOT / "src/vcc2026").glob("*.py"):
        shutil.copy2(f, t / "src/vcc2026" / f.name)
    for s in pre:
        c = run(s, t, "--check")
        co = (c.stdout + c.stderr).lower()
        if c.returncode != 0 and "already" in co and "partially" not in co:
            continue   # B2: already on HEAD (after the 10-05 / 10-06 applies)
        r = run(s, t)
        if r.returncode != 0:
            raise SystemExit(f"setup {s} failed: {(r.stdout + r.stderr)[-400:]}")
    return t


def inserted(new: str, old: str) -> str:
    """The text an edit appends after (or puts before) its anchor."""
    return new[len(old):] if new.startswith(old) else new[: len(new) - len(old)]


def harness(snippet: str, env: dict[str, str], builder_has_mark: bool, meta: dict | None, td: Path,
            pre: str = "", post: str = "", audit: str | None = None, keep: list | None = None) -> subprocess.CompletedProcess:
    w = Path(tempfile.mkdtemp(dir=td))
    if keep is not None:
        keep.append(w)
    if audit is not None:
        (w / "final.rc_audit.jsonl").write_text(audit)
    (w / "scripts").mkdir()
    (w / "scripts/make_submission_pm_cd4.py").write_text("# 630: rank-neutral UP\n" if builder_has_mark else "# none\n")
    if meta is not None:
        (w / "final.h5ad.meta.json").write_text(json.dumps({"config": meta}))
    script = ('set -u\ndie(){ echo "DIE: $*"; exit 7; }\nlog(){ echo "LOG: $*"; }\n'
              f'W="{w}"; FINAL="{w}/final"; ALLOW="base"\n{pre}\n{snippet}\n{post}\n'
              'echo "ALLOW=$ALLOW"; echo "EXPORTED=${VCC_DE_UPRN-<unset>}"; '
              '[ -f "$W/def_day.final.env" ] && cat "$W/def_day.final.env"; echo END\n')
    full_env = {"PATH": "/usr/bin:/bin:/usr/sbin:/sbin", "HOME": str(Path.home()), **env}
    return subprocess.run(["bash", "-c", script], cwd=w, capture_output=True, text=True, env=full_env)


def main() -> int:
    live = sha(ROOT / "scripts/def_day.sh")
    with tempfile.TemporaryDirectory(prefix="test630dd-") as tds:
        td = Path(tds)
        # T1
        t = tree(td, "t1")
        c = run(DD630, t, "--check")
        rec(c.returncode == 0 and "630dd check OK" in c.stdout, f"T1 --check: {c.stdout.strip()[:120]}")
        before = (t / "scripts/def_day.sh").read_text()
        a = run(DD630, t)
        after = (t / "scripts/def_day.sh").read_text()
        rec(a.returncode == 0 and after.count(P.MARK) == 7 and after != before, f"T1 apply: {a.stdout.strip()}; marks {after.count(P.MARK)}")
        b = run(DD630, t)
        rec(b.returncode != 0 and "second apply" in (b.stdout + b.stderr), "T1 second apply refused")
        t0 = tree(td, "t1b", pre=PRE[:2])          # 641 + 598 only: no 630 in the builder
        r = run(DD630, t0, "--check")
        rec(r.returncode != 0 and "apply 630 first" in (r.stdout + r.stderr), "T1 refused without 630")
        t00 = tree(td, "t1c", pre=PRE[:1])          # 641 only: def_day.sh lacks 598's mark
        r = run(DD630, t00, "--check")
        rec(r.returncode != 0 and "apply 598 first" in (r.stdout + r.stderr), "T1 refused without 598")
        t3 = tree(td, "t1d")
        bp = t3 / "scripts/make_submission_pm_cd4.py"
        bp.write_text(bp.read_text().replace(P.UPRN_MODULE_SHA256, "0" * 64))
        r = run(DD630, t3, "--check")
        rec(r.returncode != 0 and "reviewed" in (r.stdout + r.stderr), "T1 refused when the builder pins another rank_up module")
        # T2 + T3
        ta, tb = tree(td, "ord_a"), tree(td, "ord_b")
        ra = [run(DD612, ta), run(DD630, ta)]
        rb = [run(DD630, tb), run(DD612, tb)]
        ok = all(x.returncode == 0 for x in ra + rb)
        sa, sb = sha(ta / "scripts/def_day.sh"), sha(tb / "scripts/def_day.sh")
        rec(ok and sa == sb, f"T2 612dd/630dd commute: {sa[:12]} == {sb[:12]}")
        for name, x in (("630dd only", t / "scripts/def_day.sh"), ("612dd+630dd", ta / "scripts/def_day.sh")):
            bn = subprocess.run(["bash", "-n", str(x)], capture_output=True, text=True)
            rec(bn.returncode == 0, f"T3 bash -n ({name}) {bn.stderr.strip()[:200]}")
        # T4: cut the inserted text from the 630dd-only tree's file
        day = after
        ins = {name: inserted(new, old) for name, old, new in P.EDITS}
        for name in ins:
            rec(day.count(ins[name]) == 1, f"T4 inserted {name} present once")
        v, e, cc, f, al = ins["V validation"], ins["E export"], ins["C meta check"], ins["F final env"], ins["A allow"]
        body = "\n".join([v, e, cc, f, al])
        good_meta = {"de_uprn": {"variant": "b_k100_q20r3p10", "module_sha256": P.UPRN_MODULE_SHA256}}
        # unset knob: nothing happens, even with a meta that lacks de_uprn
        h = harness(body, {}, True, {"other": 1}, td)
        out = h.stdout
        rec(h.returncode == 0 and "ALLOW=base\n" in out and "EXPORTED=<unset>" in out and "VCC_DE_UPRN=" not in out
            and "DIE" not in out, "T4 unset: no export, no check, no final.env line, ALLOW unchanged")
        # set + good meta: export, check passes, final.env line, ALLOW gains de_uprn
        env = {"DEF_DAY_DE_UPRN": "b_k100_q20r3p10", "DEF_DAY_DE_MAG": "dn_k250_zs"}
        h = harness(v + "\n" + e + "\n" + cc + "\n" + f + "\n" + al, env, True, good_meta, td)
        out = h.stdout
        rec(h.returncode == 0 and "ALLOW=base,de_uprn" in out and "EXPORTED=b_k100_q20r3p10" in out
            and "VCC_DE_UPRN=b_k100_q20r3p10" in out, "T4 set: export, meta check passes, final.env line, ALLOW += de_uprn")
        for label, meta in (("wrong variant", {"de_uprn": {"variant": "b_k075_q20r3p10", "module_sha256": P.UPRN_MODULE_SHA256}}),
                            ("wrong module sha", {"de_uprn": {"variant": "b_k100_q20r3p10", "module_sha256": "8e1f05a6" + "0" * 56}}),
                            ("no de_uprn in meta", {"de_mag": {}})):
            h = harness(e + "\n" + cc, env, True, meta, td)
            rec(h.returncode == 7 and "DIE: the final meta does not carry VCC_DE_UPRN" in h.stdout, f"T4 meta check refuses: {label}")
        for variant in ("b_k100_q20r3p10", "s_k100", "b_k075_q20r3p10", "b_k50", "b_k300", "s_k125"):
            h = harness(v, {"DEF_DAY_DE_UPRN": variant, "DEF_DAY_DE_MAG": "dn_k250_zs"}, True, None, td)
            rec(h.returncode == 0 and "DIE" not in h.stdout, f"T4 validation accepts {variant}")
        for variant in ("b_k20", "b_k400", "b_k301", "b_k49", "x_k100", "B_k100", "b_k100_q20", "b_k1000", "b_k100 ", "b-k100",
                        "b_k100_q20r0p10", "b_k100_q0r3p10", "b_k100_q20r3p0", "b_k100_q020r3p10", "b_k100_q20r03p10"):
            h = harness(v, {"DEF_DAY_DE_UPRN": variant, "DEF_DAY_DE_MAG": "dn_k250_zs"}, True, None, td)
            rec(h.returncode == 7 and "must be <s|b>_k" in h.stdout, f"T4 validation refuses {variant!r}")
        h = harness(v, {"DEF_DAY_DE_UPRN": "b_k100_q20r3p10"}, True, None, td)
        rec(h.returncode == 7 and "needs DEF_DAY_DE_MAG" in h.stdout, "T4 validation refuses without DEF_DAY_DE_MAG")
        # T5: the non-fatal per-context report from the final rc_audit dump
        def row(c, t, **kw):
            base = dict(context=c, target=t, ur_n_sel=10.0, ur_frac_got=0.3, ur_n_moved=900.0, ur_guard_failed=0.0,
                        ur_invariant_failed=0.0, ur_n_no_pos_donor=0.0)
            return json.dumps({**base, **kw})
        audit = "\n".join([row("D", "T1"), row("D", "T2", ur_guard_failed=1.0, ur_frac_got=0.0, ur_n_moved=0.0),
                            row("E", "T3", ur_invariant_failed=1.0, ur_frac_got=0.0, ur_n_moved=0.0),
                            row("F", "T4", ur_n_no_pos_donor=3.0), json.dumps(dict(context="F", target="T5", zs_x=1.0))]) + "\n"
        keep: list = []
        h = harness(e + "\n" + cc, env, True, good_meta, td, audit=audit, keep=keep)
        rep = json.loads((keep[-1] / "uprn_audit.json").read_text()) if (keep[-1] / "uprn_audit.json").exists() else {}
        rec(h.returncode == 0 and "LOG: 630dd: uprn D:" in h.stdout and rep.get("D", {}).get("guard_failed") == 1
            and rep.get("D", {}).get("targets") == 2 and rep.get("E", {}).get("invariant_failed") == 1
            and rep.get("F", {}).get("no_pos_donor_genes") == 3 and rep.get("F", {}).get("targets") == 1
            and rep.get("D", {}).get("capacity_floor_met") is False and rep.get("F", {}).get("capacity_floor_met") is True,
            f"T5 report: per-context counts logged and written, not fatal ({json.dumps(rep)[:160]})")
        rec("LOG: 630dd: uprn FAILSAFE_TOTAL 2" in h.stdout, "T5 report: the greppable FAILSAFE_TOTAL line counts guard + invariant failures (2)")
        offf = keep[-1] / "uprn_effectively_off.txt"
        rec("LOG: 630dd: 630 EFFECTIVELY OFF in context D (1/2)" in h.stdout and "EFFECTIVELY OFF in context E (1/1)" in h.stdout
            and "EFFECTIVELY OFF in context F" not in h.stdout and offf.exists() and offf.read_text() == "D (1/2)\nE (1/1)\n",
            "T5 report: contexts with > 10 % fail-safe targets are flagged EFFECTIVELY OFF (logged + file)")
        sins = inserted(P.S_NEW, P.S_OLD)
        hs = harness(sins, {}, True, None, td, pre=f'cp "{offf}" "$W/uprn_effectively_off.txt"')
        rec(hs.returncode == 0 and "LOG: 630dd: WARNING 630 EFFECTIVELY OFF" in hs.stdout and "D (1/2) E (1/1)" in hs.stdout,
            "T6 summary: the EFFECTIVELY OFF warning is repeated before DONE")
        hs = harness(sins, {}, True, None, td)
        rec(hs.returncode == 0 and "WARNING" not in hs.stdout, "T6 summary: silent when no context is off")
        h = harness(e + "\n" + cc, env, True, good_meta, td, audit="{not json\n")
        rec(h.returncode == 0 and "END" in h.stdout and "LOG: 630dd: uprn FAILSAFE_TOTAL -1" in h.stdout,
            "T5 report: a malformed dump is logged as FAILSAFE_TOTAL -1, not fatal")
        h = harness(e + "\n" + cc, env, True, good_meta, td)
        rec(h.returncode == 0 and "END" in h.stdout and "LOG: 630dd: uprn FAILSAFE_TOTAL -1" in h.stdout,
            "T5 report: a missing dump is logged as FAILSAFE_TOTAL -1, not fatal")
        h = harness(e + "\n" + cc, env, True, good_meta, td, audit=json.dumps(dict(context="D", target="T1", zs_x=1.0)) + "\n")
        rec(h.returncode == 0 and "WARNING no ur_ diagnostics" in h.stdout and "LOG: 630dd: uprn FAILSAFE_TOTAL -1" in h.stdout,
            "T5 report: a dump with no ur_ rows warns and reports FAILSAFE_TOTAL -1")
        clean = "\n".join([row("D", "T1"), row("E", "T2")]) + "\n"
        keep2: list = []
        h = harness(e + "\n" + cc, env, True, good_meta, td, audit=clean, keep=keep2)
        rec(h.returncode == 0 and "LOG: 630dd: uprn FAILSAFE_TOTAL 0" in h.stdout and "EFFECTIVELY OFF" not in h.stdout
            and not (keep2[-1] / "uprn_effectively_off.txt").exists(), "T5 report: FAILSAFE_TOTAL 0 and no OFF flag on a clean dump")
        h = harness(v, env, False, None, td)
        rec(h.returncode == 7 and "no VCC_DE_UPRN (630)" in h.stdout, "T4 validation refuses a builder without 630")
    rec(sha(ROOT / "scripts/def_day.sh") == live, "live scripts/def_day.sh untouched")
    n_ok = sum(ok for ok, _ in RESULTS)
    print(f"\n{n_ok}/{len(RESULTS)} PASS" if n_ok == len(RESULTS) else f"\n{len(RESULTS) - n_ok} FAILED of {len(RESULTS)}")
    return 0 if n_ok == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
