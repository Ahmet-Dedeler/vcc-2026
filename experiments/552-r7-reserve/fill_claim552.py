#!/usr/bin/env python3
"""552 — freeze the R7 claim from CLAIM-R7-DRAFT.md, mechanically, by the draft's own rules (§2 R7A_F, §3 the pairs rule).

Inputs are the coordinator's decisions (10-07 ~04:30Z message): the bulk configuration F (a designs551.CONFIGS key) and the
decided DE switches (block B, hook H, call-set size K). Everything else follows from the draft's rule table; nothing is a
choice made here. The script reads no R7 truth: it checks only the prediction side (each arm's stats on r7bd5000 at seeds
0-2, the veto-eligible arms on the 10 control halves at seed 0) and git (553's blobs committed, unchanged).

    uv run python experiments/552-r7-reserve/fill_claim552.py --bulk-config M6+P2g --block on --hook on --k default \\
        --decisions "<the coordinator's decided-switch lines, verbatim>" [--out experiments/552-r7-reserve/CLAIM-R7.md]

It refuses to overwrite an existing claim, refuses a recipe outside the 10-07 tree, and refuses if any FILL survives. It
prints the claim's sha256 (-> the ledger's `## Reserve R7` section) and the 10-07 commands with every variable expanded.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E553 = ROOT / "experiments/553-r7b-harness"
DRAFT = HERE / "CLAIM-R7-DRAFT.md"
CONFIGS = ("M5", "M6", "M5+P2g", "M6+P2g", "M5+S_P2", "M6+S_P2")        # designs551.CONFIGS (asserted below)
HALF_CORE = ("cpm", "lcap-kept", "cpm+c4h", "lcap-kept+c4h", "cpm-n1600", "lcap-kept-n1600+c4h")
TREE = {(True, False, "default"), (False, True, "default"), (True, True, "default"), (True, True, "2000"),
        (False, False, "2000"), (False, False, "1600"), (False, False, "default")}
BLOBS = ("access553.py", "lib553.py", "band553.py", "run553.py", "prep553.py")


def mode(b: bool, h: bool, k: str) -> str:
    """The draft's MODE(B, H, K)."""
    return ("lcap-kept" if b else "cpm") + {"default": "", "2000": "-n1600", "1600": "-n1200"}[k] + ("+c4h" if h else "")


def plan_arms(b: bool, h: bool, k: str) -> dict:
    """The draft's pairs rule (§3): veto-eligible pairs on core = MODE(B, H, default); K's pair and the sign control report-only."""
    if not (b or h):
        return dict(access=False, reason="no veto-eligible switch (block / hook) decided: R7b is not accessed")
    final, core = mode(b, h, k), mode(b, h, "default")
    veto = []
    if b:
        veto.append(f"{core}:{mode(False, h, 'default')}")
    if h:
        veto.append(f"{core}:{mode(b, False, 'default')}")
    report = []
    if k != "default":
        report.append(f"{final}:{core}")
    if b and h:
        report.append("lcap-kept+c4h-s:lcap-kept+c4h")
    pairs = veto + report
    arms = []
    for x in [final, core, "lcap-kept"] + [s for p in pairs for s in p.split(":")]:
        if x != "cpm" and x not in arms:
            arms.append(x)
    half_arms = []
    for x in [s for p in veto for s in p.split(":")] + ([final, core] if k == "2000" else []):
        if x != "cpm" and x not in half_arms:
            half_arms.append(x)
    half_pairs = veto + ([f"{final}:{core}"] if k == "2000" else [])       # K on CPH is on the halves: banded, still report-only
    assert all(x in HALF_CORE for x in half_arms), half_arms
    return dict(access=True, final=final, core=core, veto=veto, report=report, arms=arms, pairs=pairs,
                half_arms=half_arms, half_pairs=half_pairs)


def check_emitted(p: dict) -> list[str]:
    """Prediction side only: every arm at seeds 0-2 on r7bd5000, every half arm at seed 0 on the 10 halves (zs form)."""
    missing = []
    for arm in ["cpm"] + p["arms"]:
        for s in (0, 1, 2):
            if not list((E553 / "work" / f"r7bd5000_{arm}-xBST2").glob(f"stats-base_s{s}_*_dnk250zs.npz")):
                missing.append(f"r7bd5000 {arm} s{s}")
    for arm in ["cpm"] + p["half_arms"]:
        for k in range(5):
            for hh in "AB":
                if not list((E553 / "work" / f"r7bd5000c{k}{hh}_{arm}-xBST2").glob("stats-base_s0_*_dnk250zs.npz")):
                    missing.append(f"r7bd5000c{k}{hh} {arm} s0")
    return missing


def rel(p: Path) -> str:
    """the path relative to the repo (the claim lives there); absolute for a test output elsewhere"""
    try:
        return str(p.resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def git(*a) -> str:
    return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk-config", required=True, choices=CONFIGS)
    ap.add_argument("--block", required=True, choices=("on", "off"))
    ap.add_argument("--hook", required=True, choices=("on", "off"))
    ap.add_argument("--k", required=True, choices=("default", "2000", "1600"))
    ap.add_argument("--decisions", required=True, help="the coordinator's decided-switch lines, verbatim (quoted in §0)")
    ap.add_argument("--out", default=str(HERE / "CLAIM-R7.md"))
    a = ap.parse_args()
    sys.path.insert(0, str(ROOT / "experiments/551-r6-c10b"))
    import designs551 as D
    assert tuple(sorted(D.CONFIGS)) == tuple(sorted(CONFIGS)), "designs551.CONFIGS changed"
    out = Path(a.out)
    assert not out.exists(), f"{out} exists: the claim is frozen once"
    b, h, k = a.block == "on", a.hook == "on", a.k
    assert (b, h, k) in TREE, f"block {a.block} / hook {a.hook} / K {k} is outside the 10-07 tree"
    p = plan_arms(b, h, k)
    if p["access"]:
        miss = check_emitted(p)
        assert not miss, f"not emitted (prediction side): {miss[:8]}"
    # 553's code: committed and unchanged (the blobs go into the claim)
    dirty = git("status", "--porcelain", "--", *[f"experiments/553-r7b-harness/{f}" for f in BLOBS])
    assert not dirty, f"553 code has uncommitted changes:\n{dirty}"
    blobs = {f: git("rev-parse", f"HEAD:experiments/553-r7b-harness/{f}") for f in BLOBS}
    t = DRAFT.read_text()
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    t = t.replace("# R7 claim — DRAFT (not the claim: nothing here is recorded in the ledger, no R7 number exists)",
                  f"# R7 claim — FROZEN {stamp} from CLAIM-R7-DRAFT.md by fill_claim552.py (its sha256 goes into the ledger "
                  f"before any R7 step)", 1)
    t = t.replace("\n## 1. What R7 confirms", "\n## 0. The decisions this claim is filled from (verbatim)\n\n"
                  + "\n".join("> " + ln for ln in a.decisions.strip().splitlines())
                  + f"\n\nFilled: bulk F `{a.bulk_config}`; DE: block {a.block}, hook {a.hook}, call-set size {k}.\n\n## 1. What R7 confirms", 1)
    t = re.sub(r"(?m)^R7A_F = FILL$", f"R7A_F = {a.bulk_config}", t, count=1)
    if p["access"]:
        rows = {"  - ARMS = FILL\n": f"  - ARMS = `{','.join(p['arms'])}` (final `{p['final']}`; `cpm` is always included)\n",
                "  - PAIRS = FILL (every switch-on-recipe pair)\n":
                    f"  - PAIRS = `{','.join(p['pairs'])}` (veto-eligible: `{','.join(p['veto'])}`; report-only: "
                    f"`{','.join(p['report']) or 'none'}`)\n",
                "  - HALF_ARMS / HALF_PAIRS = FILL:": f"  - HALF_ARMS = `{','.join(p['half_arms'])}`; HALF_PAIRS = "
                                                       f"`{','.join(p['half_pairs'])}`:"}
    else:
        rows = {"  - ARMS = FILL\n": f"  - ARMS: none. {p['reason']}.\n",
                "  - PAIRS = FILL (every switch-on-recipe pair)\n": "  - PAIRS: none.\n",
                "  - HALF_ARMS / HALF_PAIRS = FILL:": "  - HALF_ARMS / HALF_PAIRS: none (R7b not accessed):"}
    for old, new in rows.items():
        assert t.count(old) == 1, old
        t = t.replace(old, new)
    old = "plus the 555 fixes (commit 90ec5d3d, 2026-10-02)"
    assert t.count(old) == 1
    t = t.replace(old, old + "; the blobs at freeze (HEAD " + git("rev-parse", "--short", "HEAD") + "): "
                  + ", ".join(f"{f} `{v[:8]}`" for f, v in blobs.items()))
    t = t.replace("**Arms (FILL, in 553's MODE grammar):**", "**Arms (in 553's MODE grammar):**", 1)
    left = [ln for ln in t.splitlines() if re.search(r"= FILL\b|FILL:|\(FILL", ln)]   # placeholders; prose may say FILL
    assert not left, f"a placeholder is left in the claim: {left}"
    out.write_text(t)
    h = hashlib.sha256(out.read_bytes()).hexdigest()
    print(f"wrote {rel(out)}  sha256 {h}")
    print(f"ledger line (## Reserve R7): - **R7 claim recorded {stamp}**: `{rel(out)}` sha256 `{h}` "
          f"(F `{a.bulk_config}`; DE block {a.block}, hook {a.hook}, K {k}).")
    if p["access"]:
        common = (f"--exclude-targets BST2 --arms {','.join(p['arms'])} --seeds 0,1,2 --pairs {','.join(p['pairs'])} "
                  f"--half-arms {','.join(p['half_arms'])} --half-pairs {','.join(p['half_pairs'])}")
        print("R7b dry run:  uv run python experiments/553-r7b-harness/access553.py r7bd5000 --dry-run " + common
              + " --halves 0,1,2,3,4")
        print(f"R7b access:   /usr/bin/time -l uv run python experiments/553-r7b-harness/access553.py r7bd5000 --owner-go "
              f"{rel(out)} " + common + " --ctrl-splits 5 --retro work/access553_rehearsal5 --halves 0,1,2,3,4")
    else:
        print("R7b: not accessed (" + p["reason"] + ")")
    print(f"R7a record:   uv run python experiments/551-r6-c10b/r7a551.py record --claim {rel(out)} "
          f"--config {a.bulk_config} --window-start 2026-10-07T06:00Z --window-end 2026-10-07T12:00Z --received <UTC> "
          f"--message-file <message.txt>")


if __name__ == "__main__":
    main()
