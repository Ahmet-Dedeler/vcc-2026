"""551's instrument on R7a (552's reserve: panels T1-T4, fresh in targets): the R7a access driver.

It lives in 551 because 551's write guard (guard551) confines every output of the instrument to this directory: the T
panels' specs, CD4 banks, operators, fold caches and prices go to work/T1 .. work/T4, this driver's state and logs to
work/r7a/. Nothing here changes a byte of 551's own code: lib551 / full551 / comp551 / price551 run unchanged, exactly as
the efp known answers validated them (commit e7be21ca); `plan` proves that by git blob.

  plan     light, writes only work/r7a/plan.json: 552's four bundles == their draw record (draws.json: the three control
           files, pert_counts, gene_names, manifest by sha256); every module the instrument imports is byte-identical to
           its blob at e7be21ca (the efp KA); no T output exists yet; disk; the pinned external scripts' shas.
  ka-cd4   light (~5 s, < 0.5 GB): the CD4 step below, run on efp's bundle into work/r7a/ka-cd4, must equal 238's efp
           cd4_compact file for file, byte for byte (the same extractor, source and flags as 238 and 534).
  record   (the owner, after the claim is written and its sha256 is in the ledger): work/r7a/claim.json = the claim file's
           path and sha256 (asserted present in docs/CONFIRMATION-LEDGER.md), the bulk configuration F (a designs551.CONFIGS
           key), the coordinator's window and message verbatim. Refuses a second record.
  access   THE R7a ACCESS. Per panel, T1 .. T4 in order: CD4 (gates on, then relaxed, as 534) -> spec.json -> the four
           operator fits, FRESH (534's commands without --expect-sha256; each sha256 recorded) -> axis check (each
           operator's axis is its own panel's; not A/B/C's, not a tuned-on panel's, not R6's) -> full551 T --config F ->
           comp551 T <3 folds> --level --config F (no --repro: R7a has no C10a values) -> operators deleted (sha256 and
           bytes kept). Then aggregate_r7a551.py (light). Heavy steps wait for 551's gate (memory, builders, chains, vcc
           preps, disk) at C10b's floors (22 GiB + pressure 55 % + swap 400 MiB) inside the recorded window only.

    uv run python experiments/551-r6-c10b/r7a551.py plan | ka-cd4 | access
    uv run python experiments/551-r6-c10b/r7a551.py record --claim experiments/552-r7-reserve/CLAIM-R7.md --config M6+P2g \
        --window-start 2026-10-08T06:00Z --window-end 2026-10-08T18:00Z --received <UTC> --message-file <file>
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import gate551 as G  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REL = HERE.relative_to(ROOT)
WORK = HERE / "work"
R7A = WORK / "r7a"
STATE = R7A / "state.json"
CLAIM = R7A / "claim.json"
E552 = ROOT / "experiments/552-r7-reserve"
DRAWS = E552 / "work/draws.json"
LEDGER = ROOT / "docs/CONFIRMATION-LEDGER.md"
PANELS = ("T1", "T2", "T3", "T4")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FITS = {"held_k562_gwps": ("orion_hct116", "orion_hek293t"), "held_orion_hct116": ("k562_gwps", "orion_hek293t"),
        "held_orion_hek293t": ("k562_gwps", "orion_hct116"), "full3": ("k562_gwps", "orion_hct116", "orion_hek293t")}
PY = ["uv", "run", "python"]
KA_COMMIT = "e7be21ca"                   # 551's efp known answers PASS (|d| 0.0) on this code
CD4_SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
CD4_EXTRACT = ROOT / "experiments/047-cd4-direct/extract_cd4.py"
CD4_CONDITIONS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
WHITENER = ROOT / "scripts/build_union_whitener.py"
R238 = ROOT / "experiments/238-essential-free-panel/work"
W534 = ROOT / "experiments/534-generalization-audit/work/partB"


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def say(*a):
    print(f"[r7a551 {utc()}]", *a, flush=True)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


# ------------------------------------------------------------------------------------------------------------- plan
def instrument_files() -> list[str]:
    """Every repo file the instrument's processes import (lib551 + full551's and comp551's imports + price551), found by
    importing them in a child process (which writes nothing) and listing sys.modules under the repo."""
    code = ("import sys, json; sys.path.insert(0, %r); sys.argv = ['x']\n"
            "import lib551, designs551, price551, gate551\n"
            "import importlib.util as u\n"
            "for m in ('full551', 'comp551', 'aggregate551'):\n"
            "    s = u.spec_from_file_location(m + '_probe', %r + '/' + m + '.py'); mod = u.module_from_spec(s)\n"
            "    src = open(s.origin).read(); i = src.find('\\ndef '); exec(compile(src[:i], s.origin, 'exec'), mod.__dict__)\n"
            "root = %r\n"
            "print(json.dumps(sorted({m.__file__ for m in list(sys.modules.values())\n"
            "      if getattr(m, '__file__', None) and m.__file__.startswith(root + '/') and '/.venv/' not in m.__file__})))\n"
            % (str(HERE), str(HERE), str(ROOT)))
    r = subprocess.run([*PY, "-c", code], cwd=ROOT, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"instrument import probe failed:\n{r.stderr[-2000:]}")
    files = json.loads(r.stdout.strip().splitlines()[-1]) + [str(HERE / "full551.py"), str(HERE / "comp551.py")]
    return sorted({str(Path(f).relative_to(ROOT)) for f in files})


def blob_at(commit: str, rel: str) -> str | None:
    r = subprocess.run(["git", "rev-parse", f"{commit}:{rel}"], cwd=ROOT, capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def bundle_check() -> dict:
    d = json.loads(DRAWS.read_text())
    out = {}
    for p in PANELS:
        rec = d["draws"][p]["bundle"]
        b = ROOT / rec["path"]
        got = {}
        for f, want in {**rec["control_files_identical_to_199"], **rec["files"]}.items():
            h = sha256(b / f)
            assert h == want, f"{p}/{f}: sha256 {h[:12]} != draws.json {want[:12]}"
            got[f] = h
        targets = [ln.split(",")[0] for ln in (b / "pert_counts.csv").read_text().splitlines()[1:] if ln.strip()]
        assert sorted(targets) == sorted(d["draws"][p]["targets"]) and len(targets) == 300, f"{p}: target list differs"
        out[p] = dict(bundle=rec["path"], files=got)
    return out


def plan(write: bool = True) -> dict:
    t0 = time.time()
    rec = dict(at=utc(), ka_commit=KA_COMMIT, head=git("rev-parse", "HEAD"))
    rec["bundles"] = bundle_check()
    say("552's four bundles == draws.json (controls, pert_counts, gene_names, manifest by sha256; 300 targets each)")
    files = instrument_files()
    bad, pins = [], {}
    for rel in files:
        now = git("hash-object", rel)
        then = blob_at(KA_COMMIT, rel)
        pins[rel] = now
        if now != then:
            bad.append((rel, now[:12], (then or "absent")[:12]))
    rec["instrument"] = dict(n_files=len(files), blobs=pins, differs_from_ka=bad)
    assert not bad, f"instrument files differ from the efp-KA'd code at {KA_COMMIT}: {bad}"
    say(f"instrument: {len(files)} imported repo files, every one byte-identical to its blob at {KA_COMMIT} (the efp KA)")
    rec["external"] = {str(p.relative_to(ROOT)): git("hash-object", str(p.relative_to(ROOT)))
                       for p in (WHITENER, CD4_EXTRACT, HERE / "r7a551.py", HERE / "aggregate_r7a551.py") if p.exists()}
    rec["cd4_source"] = dict(path=str(CD4_SOURCE.relative_to(ROOT)), bytes=CD4_SOURCE.stat().st_size)
    exist = [str((WORK / p).relative_to(ROOT)) for p in PANELS if (WORK / p).exists()]
    rec["t_outputs_exist"] = exist
    assert not exist or STATE.exists(), f"T outputs exist before any access state: {exist}"
    rec["disk_free_gib"] = round(shutil.disk_usage(ROOT).free / 2**30, 1)
    rec["armed"] = CLAIM.exists()
    rec["seconds"] = round(time.time() - t0, 1)
    if write:
        R7A.mkdir(parents=True, exist_ok=True)
        (R7A / "plan.json").write_text(json.dumps(rec, indent=1) + "\n")
    say(f"plan OK ({rec['seconds']} s): T outputs {exist or 'none'}; disk {rec['disk_free_gib']} GiB; "
        f"{'ARMED (claim recorded)' if rec['armed'] else 'not armed (no claim recorded)'}")
    return rec


# -------------------------------------------------------------------------------------------------------------- CD4
def cd4_run(bundle: Path, out: Path, relax: bool, log: Path) -> dict:
    """534's cd4_534.run with the paths as arguments (the same extractor, source and flags)."""
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [*PY, str(CD4_EXTRACT), "--input", str(CD4_SOURCE), "--output", str(out),
           "--panel-csv", str(bundle / "pert_counts.csv"), "--genes-csv", str(bundle / "gene_names.csv"),
           "--max-rss-gib", "4.0", "--overwrite"]
    if relax:
        cmd += ["--min-targets", "0", "--min-union-targets", "0"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(r.stdout + "\n----- stderr -----\n" + r.stderr)
    return dict(relaxed=relax, rc=r.returncode, seconds=round(time.time() - t0, 1), log=str(log.relative_to(ROOT)),
                err=(r.stderr.strip().splitlines() or [""])[-2:])


def cd4_coverage(root: Path) -> dict:
    per, union = {}, set()
    for cond in CD4_CONDITIONS:
        perts = [x.strip() for x in (root / cond / "perts.txt").read_text().splitlines() if x.strip()]
        per[cond] = len(perts)
        union |= set(perts)
    genes = len([x for x in (root / CD4_CONDITIONS[0] / "genes.txt").read_text().splitlines() if x.strip()])
    return dict(per_condition=per, union=len(union), genes=genes)


def cd4_panel(bundle: Path, out: Path, logdir: Path, tag: str) -> dict:
    res = {"gates_on": cd4_run(bundle, out, False, logdir / f"cd4-{tag}.log")}
    res["relaxed"] = cd4_run(bundle, out, True, logdir / f"cd4-{tag}-relaxed.log")
    assert res["relaxed"]["rc"] == 0, f"cd4 extraction failed: {res['relaxed']}"
    res["coverage"] = cd4_coverage(out)
    return res


def tree_shas(root: Path) -> dict:
    return {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*")) if p.is_file()}


def _strip_rss(x):
    """A json record without its run-specific `peak_rss_gib` measurements (the extractor writes them into its metas)."""
    if isinstance(x, dict):
        return {k: _strip_rss(v) for k, v in x.items() if k != "peak_rss_gib"}
    if isinstance(x, list):
        return [_strip_rss(v) for v in x]
    return x


def ka_cd4() -> None:
    """Arrays and tables byte for byte; the json records equal once each `peak_rss_gib` measurement is removed."""
    out = R7A / "ka-cd4" / "cd4_compact"
    ref_root = R238 / "cd4/efp/cd4_compact"
    res = cd4_panel(R238 / "ef-bundle", out, R7A / "ka-cd4", "efp")
    mine, ref = tree_shas(out), tree_shas(ref_root)
    diff = sorted(set(mine) ^ set(ref))
    rss_only = []
    for k in sorted(set(mine) & set(ref)):
        if mine[k] == ref[k]:
            continue
        if k.endswith(".json") and _strip_rss(json.loads((out / k).read_text())) == _strip_rss(json.loads((ref_root / k).read_text())):
            rss_only.append(k)
        else:
            diff.append(k)
    res.update(files=len(mine), reference_files=len(ref), byte_identical=len(mine) - len(rss_only) - len(diff),
               json_equal_but_peak_rss=rss_only, differing=diff, passed=not diff)
    (R7A / "ka-cd4.json").write_text(json.dumps(res, indent=1) + "\n")
    say(f"KA-CD4 (efp): {len(mine)} files vs 238's {len(ref)}: " + (f"PASS ({res['byte_identical']} byte-identical, "
        f"{len(rss_only)} json records equal but for peak_rss_gib)" if not diff else "FAIL " + str(diff[:5])))
    if diff:
        raise SystemExit(1)
    for f in sorted(out.rglob("*"), key=lambda x: len(x.parts), reverse=True):   # the copy is not kept (== 238's);
        f.unlink() if f.is_file() else f.rmdir()     # absolute paths: shutil.rmtree's dir_fd calls trip guard551
    out.rmdir()


# --------------------------------------------------------------------------------------------------------- record
def record(a) -> None:
    import designs551 as D
    assert not CLAIM.exists(), f"{CLAIM} exists: the R7a claim is recorded once"
    assert a.config in D.CONFIGS, f"--config must be one of {sorted(D.CONFIGS)}"
    claim = (ROOT / a.claim).resolve()
    h = sha256(claim)
    assert h in LEDGER.read_text(), f"the claim's sha256 {h} is not in {LEDGER.relative_to(ROOT)}: record it there first"
    ws, we = (dt.datetime.strptime(x, "%Y-%m-%dT%H:%MZ").replace(tzinfo=dt.timezone.utc) for x in (a.window_start, a.window_end))
    assert we > ws, "window end before start"
    rec = dict(claim=str(claim.relative_to(ROOT)), claim_sha256=h, config=a.config, base=D.CONFIGS[a.config][0],
               step=D.CONFIGS[a.config][1], window=dict(start=a.window_start, end=a.window_end), received=a.received,
               message=Path(a.message_file).read_text(), recorded_at=utc(), by='session "Explore contribution areas"',
               plan=plan(write=False))
    R7A.mkdir(parents=True, exist_ok=True)
    CLAIM.write_text(json.dumps(rec, indent=1) + "\n")
    say(f"recorded: claim {rec['claim']} sha256 {h[:12]} (in the ledger), F {a.config}, window {a.window_start} - {a.window_end}")


# --------------------------------------------------------------------------------------------------------- access
def state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def done(step: str, rec: dict) -> None:
    s = state()
    s[step] = dict(rec, done=utc())
    STATE.write_text(json.dumps(s, indent=1))


def arm_gate(claim: dict) -> None:
    """551's gate with R7a's window and C10b's floors; its logs go to work/r7a/ (module globals, read at call time)."""
    ws, we = (dt.datetime.strptime(claim["window"][k], "%Y-%m-%dT%H:%MZ").replace(tzinfo=dt.timezone.utc)
              for k in ("start", "end"))

    def window_ok(now=None):
        now = now or dt.datetime.now(dt.timezone.utc)
        if not (ws <= now < we):
            return False, f"R7a's window is {claim['window']['start']} - {claim['window']['end']} (the coordinator's)"
        return True, "R7a window open"
    G.window_ok = window_ok
    G.thresholds = lambda now=None: (22.0, 55.0, 400.0, "R7a access: 22 GiB + pressure 55 % + swap 400 MiB (C10b's floors)")
    G.GATE_LOG, G.HEAVY_LOG = R7A / "gate.log", R7A / "heavy.log"


def heavy(step: str, cmd: list[str], log: Path) -> dict:
    if step in state():
        say(f"skip {step} (done)")
        return state()[step]
    say(f"{step}: waiting for the gate")
    r = G.heavy(cmd, step, log)
    say(f"{step}: rc {r['rc']}, {r['wall_s']:.0f} s, peak RSS {r['peak_rss_gib']:.2f} GiB")
    if r["rc"]:
        raise SystemExit(f"{step} failed (rc {r['rc']}); see {log}")
    done(step, r)
    return r


def axis_check(p: str, paths: dict) -> dict:
    """axes534's rule in a child process (534's module, which writes nothing): each operator's json axis == the axis
    recomputed from ITS panel's bundle and donor set; not A/B/C's axis, not a tuned-on panel's, not any R6 panel's."""
    bundle = ROOT / json.loads(DRAWS.read_text())["draws"][p]["bundle"]["path"]
    code = ("import sys, json; sys.path.insert(0, %r)\n"
            "import axes534 as X\nfrom pathlib import Path\n"
            "ops = json.loads(sys.argv[1]); bundle = Path(sys.argv[2]); tuned = X.tuned_on_axes(); out = {}\n"
            "r6 = {X.axis_sha(Path(%r) / s / 'bundle', X.DONOR_SETS[t])[0]: s + '/' + t for s in ('S1','S2','S3','S4') for t in X.DONOR_SETS}\n"
            "for tag, path in ops.items():\n"
            "    meta = json.loads(Path(path).with_suffix('.json').read_text())\n"
            "    sha, n = X.axis_sha(bundle, X.DONOR_SETS[tag])\n"
            "    assert meta['axis_sha256'] == sha and int(meta['n_cols']) == n and list(meta['donors']) == list(X.DONOR_SETS[tag]), (tag, meta['axis_sha256'][:12], sha[:12])\n"
            "    assert sha not in X.ABC_AXES and sha not in tuned and sha not in r6, (tag, 'axis equals an existing panel axis')\n"
            "    out[tag] = dict(axis_sha256=sha, n_cols=n, op_sha256=meta['sha256'])\n"
            "print(json.dumps(out))\n" % (str(W534), str(W534)))
    r = subprocess.run([*PY, "-c", code, json.dumps({k: str(v) for k, v in paths.items()}), str(bundle)],
                       cwd=W534, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"{p}: axis check failed:\n{r.stderr[-1500:]}")
    return json.loads(r.stdout.strip().splitlines()[-1])


def panel(p: str, F: str) -> None:
    pw = WORK / p
    ops = pw / "ops"
    bundle = ROOT / json.loads(DRAWS.read_text())["draws"][p]["bundle"]["path"]
    if f"{p}:cd4" not in state():
        res = cd4_panel(bundle, pw / "cd4" / "cd4_compact", R7A / "logs", p)
        (pw / "cd4_r7a.json").write_text(json.dumps(res, indent=1) + "\n")
        done(f"{p}:cd4", dict(coverage=res["coverage"], gates_on_rc=res["gates_on"]["rc"]))
        say(f"{p}: CD4 banks (relaxed; gates-on rc {res['gates_on']['rc']}); coverage {res['coverage']}")
    ops.mkdir(parents=True, exist_ok=True)
    spec = dict(bundle=str(bundle), rules="def", cd4_root=str(pw / "cd4/cd4_compact"), full_op=str(ops / f"{p}__full3.npz"),
                folds=list(CANCER), synthetic=True, meta=None, fold_ops={h: str(ops / f"{p}__held_{h}.npz") for h in CANCER})
    (pw / "spec.json").write_text(json.dumps(spec, indent=1))
    paths = {}
    for tag, donors in FITS.items():
        out = ops / f"{p}__{tag}.npz"
        paths[tag] = out
        step = f"{p}:fit:{tag}"
        if step in state():
            continue
        free = shutil.disk_usage(ROOT).free / 2**30
        assert free > 10.0, f"disk free {free:.1f} GiB"
        cmd = [*PY, str(WHITENER.relative_to(ROOT)), "--bundle", str(bundle), "--donors", ",".join(donors),
               "--basis-prune", "essential", "--out", str(out), "--scratch", str(ops / ".scratch"), "--max-rss-gib", "14"]
        heavy(step, cmd, R7A / f"logs/fit-{p}-{tag}.log")
        got = json.loads(out.with_suffix(".json").read_text())
        s = state()
        s[step].update(op_sha256=got["sha256"], axis_sha256=got["axis_sha256"], command=" ".join(cmd))
        STATE.write_text(json.dumps(s, indent=1))
        say(f"{p} {tag}: operator fitted fresh, sha256 {got['sha256'][:12]}, axis {got['axis_sha256'][:12]}")
    if f"{p}:axes" not in state():
        rec = axis_check(p, paths)
        (pw / "axes_r7a.json").write_text(json.dumps(rec, indent=1) + "\n")
        done(f"{p}:axes", rec)
        say(f"{p}: four axes are the panel's own (not A/B/C's, not tuned-on, not R6's)")
    lv = str(REL / f"work/{p}/level551.json")
    heavy(f"{p}:full551", [*PY, str(REL / "full551.py"), p, "--config", F], R7A / f"logs/full551-{p}.log")
    heavy(f"{p}:comp551", [*PY, str(REL / "comp551.py"), p, ",".join(CANCER), "--level", lv, "--config", F],
          R7A / f"logs/comp551-{p}.log")
    for f in CANCER:
        c = json.loads((pw / f"comp551-{f}.json").read_text())
        assert "peak_rss_gib" in c and c["op_axis_checked"], (p, f)
    if f"{p}:ops_deleted" not in state():
        rec = {}
        for tag, path in paths.items():
            rec[tag] = dict(sha256=json.loads(path.with_suffix(".json").read_text())["sha256"], bytes=path.stat().st_size)
            path.unlink()
        done(f"{p}:ops_deleted", rec)
        say(f"{p}: priced; operators deleted ({sum(r['bytes'] for r in rec.values()) / 2**30:.1f} GiB freed; jsons kept)")


def access() -> None:
    assert CLAIM.exists(), "not armed: record the claim first (`r7a551.py record ...`, after its sha256 is in the ledger)"
    claim = json.loads(CLAIM.read_text())
    assert sha256(ROOT / claim["claim"]) == claim["claim_sha256"], "the claim file changed after it was recorded"
    assert claim["claim_sha256"] in LEDGER.read_text(), "the claim's sha256 is no longer in the ledger"
    p0 = plan(write=False)
    assert p0["instrument"]["blobs"] == claim["plan"]["instrument"]["blobs"], "the instrument changed since the record"
    assert p0["external"] == claim["plan"]["external"], "an external script changed since the record"
    arm_gate(claim)
    say(f"R7a ACCESS: F {claim['config']}, claim {claim['claim']} ({claim['claim_sha256'][:12]})")
    for p in PANELS:
        panel(p, claim["config"])
    if "aggregate" not in state():
        log = R7A / "logs/aggregate_r7a551.log"
        with open(log, "w") as fh:
            rc = subprocess.run([*PY, str(REL / "aggregate_r7a551.py")], cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
        if rc:
            raise SystemExit(f"aggregate_r7a551 failed (rc {rc}); see {log}")
        done("aggregate", dict(rc=rc, log=str(log.relative_to(ROOT))))
    say("R7a access: complete")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("plan", "ka-cd4", "record", "access"))
    ap.add_argument("--claim")
    ap.add_argument("--config")
    ap.add_argument("--window-start")
    ap.add_argument("--window-end")
    ap.add_argument("--received")
    ap.add_argument("--message-file")
    a = ap.parse_args()
    R7A.mkdir(parents=True, exist_ok=True)
    try:
        {"plan": lambda: plan(), "ka-cd4": ka_cd4, "record": lambda: record(a), "access": access}[a.what]()
    except G.WindowClosed as e:
        say(f"WINDOW CLOSED ({e}); re-run `uv run python {REL}/r7a551.py access` when it reopens (done steps are skipped)")
        raise SystemExit(3)


if __name__ == "__main__":
    main()
