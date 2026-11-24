"""432: known-answer rows for the PROPOSED `def_pin_panel_artifacts.py --pert-scale` (never applied; nothing written).

  R1  dry run of the proposed script's main() (ROOT/BUILDER pointed at the repo) with the day's operator, CD4 banks and
      key: the diff adds exactly ONE PERT_SCALE_ARTIFACTS entry + ONE WHITENER_REGISTRY entry and replaces only CD4
      hash lines (<= 15: bank files identical across the two panels keep their hash)
  R2  the rewritten source compiles
  R3  rewrite() is idempotent on its own output (a second pin of the same key is a no-op)
  R4  a re-pin with a DIFFERENT key sha replaces the entry in place (still exactly one entry)
  R5  a hand-written entry under the same id is refused
  R6  the builder on disk is byte-identical before and after; R7 the scratch HEAD copy too (dry run)
  (the rewrite runs on git HEAD's builder blob, copied to work/pintest/, never on the live file)
usage: uv run python proposals/test_pin432.py [--set r|ef]
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
W = HERE.parent / "work"
HEADB = "2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f"
SETS = {
    "ef": dict(op=ROOT / "experiments/244-october-rehearsal/work/run1/whitener_union_pruned_def.npz",
               cd4=ROOT / "experiments/244-october-rehearsal/work/run1/cd4_compact",
               key=W / "key/pert_scale_ef432_lo0.3_hi3.npz"),
    "r": dict(op=W / "base/whitener_union_pruned_def.npz", cd4=W / "base/cd4_compact",
              key=W / "key/pert_scale_def432_lo0.3_hi3.npz"),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="ef", choices=sorted(SETS))
    a = ap.parse_args()
    s = SETS[a.set]
    live = ROOT / "scripts/make_submission_pm_cd4.py"
    b0 = hashlib.sha256(live.read_bytes()).hexdigest()
    # rewrite HEAD's builder, never the file on disk (another run's pin may be live there): a scratch copy of the blob
    import subprocess
    blob = subprocess.run(["git", "show", "HEAD:scripts/make_submission_pm_cd4.py"], cwd=ROOT, capture_output=True,
                          check=True).stdout
    assert hashlib.sha256(blob).hexdigest() == HEADB, "git HEAD's builder moved; update HEADB"
    builder = W / "pintest/make_submission_pm_cd4.head.py"
    builder.parent.mkdir(parents=True, exist_ok=True)
    builder.write_bytes(blob)
    spec = importlib.util.spec_from_file_location("pin432", HERE / "def_pin_panel_artifacts.proposed.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.ROOT, mod.BUILDER = ROOT, builder
    rows = []

    def row(name, ok, note=""):
        rows.append(dict(row=name, result="PASS" if ok else "FAIL", note=note))
        print(f"{'PASS' if ok else 'FAIL'}  {name}  {note}", flush=True)

    argv = ["pin432", "--whitener", str(s["op"]), "--whitener-metadata", str(s["op"].with_suffix(".json")),
            "--whitener-id", "union_pruned_def", "--cd4-root", str(s["cd4"]), "--pert-scale", str(s["key"])]
    buf = io.StringIO()
    old = sys.argv
    sys.argv = argv
    try:
        with contextlib.redirect_stdout(buf):
            mod.main()
    finally:
        sys.argv = old
    out = buf.getvalue()
    (W / f"test_pin432_{a.set}.diff").write_text(out)
    plus = [ln for ln in out.splitlines() if ln.startswith("+") and not ln.startswith("+++")]
    minus = [ln for ln in out.splitlines() if ln.startswith("-") and not ln.startswith("---")]
    n_ps = sum('"wshare_lo0.3_hi3_def": (' in ln for ln in plus)
    n_wr = sum('"union_pruned_def": (' in ln for ln in plus)
    key_sha = hashlib.sha256(s["key"].read_bytes()).hexdigest()
    row("R1 dry-run diff: one key entry, one operator entry, only CD4 hash lines replaced (<= 15; files equal across panels keep theirs)",
        n_ps == 1 and n_wr == 1 and 0 < len(minus) <= 15 and all(re.fullmatch(r'-\s+"[0-9a-f]{64}"', ln) for ln in minus)
        and any(key_sha in ln for ln in plus) and "dry run: nothing written" in out,
        f"+{len(plus)} -{len(minus)} key_entries {n_ps} operator_entries {n_wr}")
    # rebuild the rewrite() inputs exactly as main() did, then exercise rewrite() directly
    cd4_sha = {f"cd4/{b}/{f}": mod.sha256(s["cd4"] / b / f) for b in mod.CD4_NAMES for f in mod.CD4_FILES}
    npz_sha, json_sha = mod.sha256(s["op"]), mod.sha256(s["op"].with_suffix(".json"))
    ps = ("wshare_lo0.3_hi3_def", str(s["key"].resolve().relative_to(ROOT)), key_sha)
    src = builder.read_text()
    new1 = mod.rewrite(src, "union_pruned_def", npz_sha, json_sha, cd4_sha, ps)
    try:
        compile(new1, str(builder), "exec")
        row("R2 rewritten source compiles", True)
    except SyntaxError as e:
        row("R2 rewritten source compiles", False, str(e))
    new2 = mod.rewrite(new1, "union_pruned_def", npz_sha, json_sha, cd4_sha, ps)
    row("R3 idempotent", new2 == new1)
    new3 = mod.rewrite(new1, "union_pruned_def", npz_sha, json_sha, cd4_sha, (ps[0], ps[1], "0" * 64))
    row("R4 re-pin replaces in place", new3.count('"wshare_lo0.3_hi3_def": (') == 1 and "0" * 64 in new3 and key_sha not in new3)
    hand = src.replace("PERT_SCALE_ARTIFACTS = {\n",
                       'PERT_SCALE_ARTIFACTS = {\n    "wshare_lo0.3_hi3_def": (\n        ROOT / "x.npz",\n        "' + "1" * 64 + '",\n    ),\n', 1)
    try:
        mod.rewrite(hand, "union_pruned_def", npz_sha, json_sha, cd4_sha, ps)
        row("R5 hand-written entry refused", False)
    except SystemExit as e:
        row("R5 hand-written entry refused", "hand-written" in str(e), str(e)[:90])
    b1 = hashlib.sha256(live.read_bytes()).hexdigest()
    row("R6 the builder on disk is untouched by the dry run", b0 == b1, f"{b1[:16]} ({'HEAD' if b1 == HEADB else 'a live pin'})")
    row("R7 the scratch copy is untouched too (dry run writes nothing)", hashlib.sha256(builder.read_bytes()).hexdigest() == HEADB)
    (W / f"test_pin432_{a.set}.json").write_text(json.dumps(rows, indent=1))
    bad = [r for r in rows if r["result"] != "PASS"]
    print(f"{len(rows) - len(bad)}/{len(rows)} rows PASS")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
