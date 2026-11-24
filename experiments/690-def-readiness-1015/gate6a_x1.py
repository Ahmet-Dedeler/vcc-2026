#!/usr/bin/env python3
"""690 gate6a_x1 (readiness [401d9b], 2026-10-22): OCT22-CHECKLIST §C 6a's lines for X1, judged on the FINISHED X1 build against R0.

X1 = R0's day recipe + DEF_DAY_H1_STACK=x3 + DEF_DAY_CD4_STACK=x1.5 + DEF_DAY_DUMP_HOOK=dumphook707.sh (day708's FINAL, the
coordinator's GO, SIXA d98c1c48). R0 = experiments/692-def-r0/work/r0 (entry k2OKr0zPvJq48BCtdUx8). Checks, each fail-closed:
  (a) the content-aware config diff vs R0 (scripts/def_metadiff.py) is within {bulk_load, bulk_load_sha256, bulk_ctx_scale,
      bulk_sharpen_pt} and contains bulk_load + bulk_load_sha256 (a dose and a hook change the dump by design);
  (a') bulk_sharpen_pt differs from R0's only in its DATA-DERIVED fields {edges, key_ratio_c, key_ratio_maxdev} (regate673's
      rule, the W1 lesson): mode, q, counts, carried key + sha, per-target moves, n_carried and every tally are equal;
  (a'') the switches really are in the build: the stack dump's env has VCC_H1_UNION_G 0.84 and VCC_CD4_TARGET_G
      0.5165482917974462, def_day.log shows the 707 hook ran and replaced the dump, and the final meta's bulk_load is the hook's;
  (b) R0's gates again: the ready log says READY with the package's sha, the run log has "[<stem>] PASS", "[<stem>] CLIFF OK"
      and "novelty rc=0"; the final ntclevel is in band (def_day_helpers inband: every rho <= 1 and within 1 % of rho*);
  (e) the same tree as R0: hashes.implementation equal, config.de_sign_emitter_sha256 equal, and def_day.log's builder
      restored to e99bbd4a;
  (d) bookkeeping only: the package is novel (the run log's novelty line) and not in the submitted ledger.
Writes <out>/gate6a_<candidate>.json and prints "6A PASS" or "6A FAIL: <reasons>"; exit 0 / 1.
  uv run python experiments/690-def-readiness-1015/gate6a_x1.py --ready-dir experiments/692-def-r0/work/x1-ready [--stem x1def]

--candidate x2 (added 10-22 07:2xZ; SIXA 1ebd3a1f: X2 = C0 = R0's day recipe + DEF_DAY_DUMP_HOOK=dumphook811c0.sh, no weights move,
no swap). Its lines: (a) the diff vs R0 within {bulk_ctx_scale, bulk_load, bulk_load_sha256} needing bulk_load_sha256 (699's on
A/B/C); bulk_sharpen_pt is tolerated ONLY under (a')'s data-fields rule and is then NOTED for the coordinator; (a'') the stack dump
keeps R0's doses (VCC_H1_UNION_G 0.56, VCC_CD4_TARGET_G 0.3443655278649641) and equals R0's stack dump ada1baa2 BY CONTENT (bulk,
carriers, perts, genes), the 811 hook ran and replaced it, and the used dump equals the hook applied to R0's stack dump BY CONTENT
(--hook-ref: that output, made before the build); (b), (e), (d) as X1.
  uv run python experiments/690-def-readiness-1015/gate6a_x1.py --candidate x2 --work experiments/692-def-r0/work/x2 --stem x2def \
      --ready-dir experiments/692-def-r0/work/x2-ready --hook-ref <hookdir>/dump_h1x2_cd4x1_s0.25.npz
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R0W = ROOT / "experiments/692-def-r0/work/r0"
R0M = R0W / "r0def_final.h5ad.meta.json"
ALLOW = ["bulk_load", "bulk_load_sha256", "bulk_ctx_scale", "bulk_sharpen_pt"]
NEED = ["bulk_load", "bulk_load_sha256"]
DATA_FIELDS = {"edges", "key_ratio_c", "key_ratio_maxdev"}
H1G, CD4G = "0.84", "0.5165482917974462"
CANDS = {
    "x1": dict(allow=ALLOW, need=NEED, h1g="0.84", cd4g="0.5165482917974462", sharpen_free=True,
               hook="experiments/707-strong-tier-direction/dumphook707.sh"),
    "x2": dict(allow=["bulk_load", "bulk_load_sha256", "bulk_ctx_scale"], need=["bulk_load_sha256"], h1g="0.56",
               cd4g="0.3443655278649641", sharpen_free=False, hook="experiments/811-lolo-model/dumphook811c0.sh"),
}
R0_DUMP = R0W / "dump_h1x2_cd4x1_s0.25.npz"


def same_content(p: Path, q: Path) -> bool:
    a, b = np.load(p, allow_pickle=True), np.load(q, allow_pickle=True)
    keys = ("bulk", "carriers", "perts", "genes")
    return all(k in a.files and k in b.files and a[k].shape == b[k].shape and np.array_equal(a[k], b[k]) for k in keys)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, default=ROOT / "experiments/692-def-r0/work/x1")
    ap.add_argument("--stem", default="x1def")
    ap.add_argument("--ready-dir", type=Path, required=True, help="the dir holding <stem>.ready.log / <stem>.run.log (x1-ready or a re-gate's)")
    ap.add_argument("--pkg", type=Path, default=None)
    ap.add_argument("--candidate", choices=sorted(CANDS), default="x1")
    ap.add_argument("--hook-ref", type=Path, default=None, help="x2: the hook's output on R0's stack dump (content reference)")
    a = ap.parse_args()
    K = CANDS[a.candidate]
    # 10-22 10:12Z: resolve every path against the repo (the chain passes relative paths; hook dumps are compared by relative_to(ROOT))
    for f in ("work", "ready_dir", "pkg", "hook_ref"):
        v = getattr(a, f)
        if v is not None and not v.is_absolute():
            setattr(a, f, ROOT / v)
    W, N = a.work, a.stem
    meta_p = W / f"{N}_final.h5ad.meta.json"
    pkg = a.pkg or ROOT / f"submissions/{N}.vcc"
    fails, notes = [], []
    meta, r0 = json.loads(meta_p.read_text()), json.loads(R0M.read_text())
    c, c0 = meta["config"], r0["config"]
    dl = (W / "def_day.log").read_text()

    # (a) content-aware config diff vs R0
    md_json = a.ready_dir / "gate6a_metadiff.json"
    p = subprocess.run([sys.executable, str(ROOT / "scripts/def_metadiff.py"), str(meta_p), str(W), str(R0M), str(R0W),
                        "--allow", ",".join(K["allow"] + (["bulk_sharpen_pt"] if not K["sharpen_free"] else [])),
                        "--need", ",".join(K["need"]), "--json", str(md_json)],
                       capture_output=True, text=True, cwd=ROOT)
    last = (p.stdout.strip().splitlines() or [""])[-1]
    notes.append(f"(a) {last}")
    if p.returncode != 0 or not last.startswith("[metadiff] PASS"):
        fails.append(f"(a) config diff vs R0: rc {p.returncode}: {last[:300]}")

    # (a') bulk_sharpen_pt: only the data-derived fields may move
    sp, sp0 = c.get("bulk_sharpen_pt") or {}, c0.get("bulk_sharpen_pt") or {}
    moved = sorted(k for k in set(sp) | set(sp0) if sp.get(k) != sp0.get(k))
    behav = [k for k in moved if k not in DATA_FIELDS]
    notes.append(f"(a') bulk_sharpen_pt moved {moved}")
    if behav:
        fails.append(f"(a') bulk_sharpen_pt behavioural fields moved: {behav}")
    elif moved and not K["sharpen_free"]:
        notes.append(f"(a') NOTE for the coordinator: bulk_sharpen_pt's data-derived fields moved {moved} (not in the X2 spec; the "
                     "W1 rule: behavioural fields identical)")

    # (a'') the switches are in the build
    dump_name = "dump_h1x2_cd4x1_s0.25.npz"
    if a.candidate == "x1":
        m = re.search(r"dump ([^ ]+\.npz) VCC_H1_UNION_G=([0-9.]+)", dl)
        if not m or m.group(2) != K["h1g"]:
            fails.append(f"(a'') no 'H1 stack x3 ... VCC_H1_UNION_G={K['h1g']}' line in def_day.log")
            dump_name = None
        else:
            dump_name = m.group(1)
        if f"VCC_CD4_TARGET_G={K['cd4g']}" not in dl:
            fails.append(f"(a'') no 'CD4 stack x1.5 ... VCC_CD4_TARGET_G={K['cd4g']}' line in def_day.log")
    else:
        if "H1 stack" in dl or "CD4 stack" in dl:
            fails.append("(a'') def_day.log shows an H1 / CD4 stack move (X2 keeps R0's doses)")
    if dump_name and (W / dump_name).exists():
        env = json.loads(str(np.load(W / dump_name, allow_pickle=True)["env"]))
        if env.get("VCC_H1_UNION_G") != K["h1g"] or env.get("VCC_CD4_TARGET_G") != K["cd4g"]:
            fails.append(f"(a'') the stack dump's env has H1 {env.get('VCC_H1_UNION_G')} CD4 {env.get('VCC_CD4_TARGET_G')}")
        else:
            notes.append(f"(a'') stack dump {dump_name}: VCC_H1_UNION_G {K['h1g']}, VCC_CD4_TARGET_G {K['cd4g']}")
        if a.candidate == "x2":
            if same_content(W / dump_name, R0_DUMP):
                notes.append("(a'') X2's stack dump == R0's stack dump ada1baa2 by content")
            else:
                fails.append("(a'') X2's stack dump differs from R0's by content")
    elif dump_name:
        fails.append(f"(a'') stack dump {dump_name} missing")
    if f"690g: dump hook {K['hook']}" not in dl or "690g: dump replaced by the hook's" not in dl:
        fails.append(f"(a'') def_day.log does not show {K['hook']} running and replacing the dump")
    hk = W / "dumphook" / (dump_name or "?")
    if dump_name and hk.exists():
        # the meta records bulk_load as the file NAME (the hook writes under the stack dump's name), so the sha decides
        bl = str(c.get("bulk_load", ""))
        if Path(bl).name != dump_name:
            fails.append(f"(a'') the final meta's bulk_load {bl} is not named {dump_name}")
        elif c.get("bulk_load_sha256") != sha256(hk):
            fails.append("(a'') bulk_load_sha256 != the hook dump's sha")
        elif c.get("bulk_load_sha256") == sha256(W / dump_name):
            fails.append("(a'') bulk_load_sha256 is the UN-hooked stack dump's sha")
        else:
            notes.append(f"(a'') bulk_load_sha256 == the hook's dump {sha256(hk)[:16]} (!= the stack dump {sha256(W / dump_name)[:16]})")
        if a.candidate == "x2":
            if not a.hook_ref or not a.hook_ref.exists():
                fails.append("(a'') x2 needs --hook-ref (the hook's output on R0's stack dump)")
            elif same_content(hk, a.hook_ref):
                notes.append(f"(a'') the used dump == the 811 hook on R0's stack dump, by content ({a.hook_ref})")
            else:
                fails.append("(a'') the used dump differs from the 811 hook applied to R0's stack dump")
    elif dump_name:
        fails.append(f"(a'') no hook dump {hk}")

    # (b) R0's gates again
    rdy = (a.ready_dir / f"{N}.ready.log").read_text() if (a.ready_dir / f"{N}.ready.log").exists() else ""
    rl = (a.ready_dir / f"{N}.run.log").read_text() if (a.ready_dir / f"{N}.run.log").exists() else ""
    mm = re.findall(rf"READY {N} ([0-9a-f]{{64}})", rdy)
    if "NOT READY" in rdy or not mm:
        fails.append(f"(b) the ready log has no clean READY line: {rdy.strip()[-200:]}")
    elif not pkg.exists() or sha256(pkg) != mm[-1]:
        fails.append("(b) the package on disk is not the READY sha")
    else:
        notes.append(f"(b) READY {mm[-1][:16]} == package on disk")
    for pat in (rf"^\[{N}\] PASS", rf"^\[{N}\] CLIFF OK", r"novelty rc=0"):
        if not re.search(pat, rl, re.M):
            fails.append(f"(b) run log lacks {pat}")
    ft = re.findall(r" final: ([^ ]+) factors ([^ ]+) -> ", dl)
    if not ft:
        fails.append("(b) no final line in def_day.log")
    else:
        nj = W / f"{ft[-1][0]}.ntclevel.json"
        q = subprocess.run([sys.executable, str(ROOT / "scripts/def_day_helpers.py"), "inband", str(nj)], cwd=ROOT)
        rho = subprocess.run([sys.executable, str(ROOT / "scripts/def_day_helpers.py"), "rho", str(nj)], cwd=ROOT,
                             capture_output=True, text=True).stdout.strip()
        notes.append(f"(b) final {ft[-1][0]} at {ft[-1][1]}: {rho} (inband rc {q.returncode})")
        if q.returncode != 0:
            fails.append(f"(b) the final build is not in band (rho {rho})")

    # (e) the same tree. 10-22 10:1xZ (found on X1): the PINNED builder's sha can never equal R0's across work dirs: the day's pin
    # embeds the work-dir paths and the sha of each cd4_compact/<cond>/meta.json, which records a measured peak_rss_gib. So the
    # builder entry is judged by content: each pinned builder must be the meta's builder (provenance), the two sources must be
    # identical after normalising the work dir and the cd4 meta.json pins, and every cd4 meta.json may differ from R0's only in
    # VOLATILE_CD4 keys. Every other implementation entry must be byte-equal.
    BK = "scripts/make_submission_pm_cd4.py"
    VOLATILE_CD4 = {"peak_rss_gib"}
    im, im0 = meta["hashes"]["implementation"], r0["hashes"]["implementation"]
    other = sorted(k for k in set(im) | set(im0) if k != BK and im.get(k) != im0.get(k))
    if other:
        fails.append(f"(e) hashes.implementation differ from R0's: {other}")
    if im.get(BK) != im0.get(BK):
        pb, pb0 = W / "builder.pinned.py", R0W / "builder.pinned.py"
        if not (pb.exists() and pb0.exists()) or sha256(pb) != im.get(BK) or sha256(pb0) != im0.get(BK):
            fails.append("(e) a pinned builder is missing or is not its meta's builder")
        else:
            def norm(t: str, w: Path) -> str:
                t = t.replace(str(w.relative_to(ROOT)), "WORK")
                return re.sub(r'("cd4/[^"]+/meta\.json": \(\s*")[0-9a-f]{64}"', r'\1META"', t)
            if norm(pb.read_text(), W) != norm(pb0.read_text(), R0W):
                fails.append("(e) the pinned builders differ beyond the work dir and the cd4 meta.json pins")
            else:
                bad = []
                for mj in sorted((R0W / "cd4_compact").glob("*/meta.json")):
                    mx = W / "cd4_compact" / mj.parent.name / "meta.json"
                    j0, j1 = json.loads(mj.read_text()), json.loads(mx.read_text()) if mx.exists() else None
                    if j1 is None:
                        bad.append(f"{mj.parent.name}: missing")
                        continue
                    ks = {k for k in set(j0) | set(j1) if j0.get(k) != j1.get(k)}
                    if ks - VOLATILE_CD4:
                        bad.append(f"{mj.parent.name}: {sorted(ks - VOLATILE_CD4)}")
                if bad:
                    fails.append(f"(e) cd4 meta.json differ beyond {sorted(VOLATILE_CD4)}: {bad}")
                else:
                    notes.append(f"(e) pinned builder {im[BK][:8]} vs R0's {im0[BK][:8]}: identical after normalising the work dir and "
                                 f"the cd4 meta.json pins; those metas differ only in {sorted(VOLATILE_CD4)}")
    if c.get("de_sign_emitter_sha256") != c0.get("de_sign_emitter_sha256"):
        fails.append("(e) de_sign_emitter_sha256 differs from R0's")
    rb = re.findall(r"builder after restore: ([0-9a-f]{16}) \(want ([0-9a-f]{16})\)", dl)
    if not rb or rb[-1][0] != rb[-1][1] or rb[-1][0] != "e99bbd4a63ad834c":
        fails.append(f"(e) builder not restored to e99bbd4a: {rb[-1:] }")
    else:
        notes.append("(e) implementation + emitter == R0's; builder restored e99bbd4a")

    # (d) bookkeeping
    led = json.loads((ROOT / "submissions/SUBMITTED_LEDGER.json").read_text())
    if mm and mm[-1] in json.dumps(led):
        fails.append("(d) the package sha is already in the submitted ledger")

    out = {"candidate": a.candidate, "stem": N, "work": str(W), "ready_dir": str(a.ready_dir), "pass": not fails, "fails": fails,
           "notes": notes}
    (a.ready_dir / f"gate6a_{a.candidate}.json").write_text(json.dumps(out, indent=1))
    for n in notes:
        print("  " + n)
    print("6A PASS" if not fails else "6A FAIL: " + " | ".join(fails))
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
