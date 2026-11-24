"""438 DAY-OF: re-price the bulk stack on the day's D/E/F bundle and print the recommended dump recipe.

    uv run python experiments/438-october-stack-repricing/dayof438.py --bundle <D/E/F bundle> --op <day's 3-donor
        whitener .npz> --cd4-root <day's cd4_compact> [--tag dayof] [--fold-ops prod|held --held-ops <dir>]

What it runs (strictly serial, each step behind the memory gate, every process <= 8.8 GB RSS; rehearsal: 27.5 min):
  1. full438   the full panel (300 targets, 3 donors, the day's operator): D/E/F composition rules, the 410 key and
               LODO-V1 re-fitted on the panel, every grid design's realised render energy, the anchored rho per context
               (claims from the bundle's own NTC cells, ntclevel417's f = 1 rule), and each design's re-level factor.
  2. comp438   PLOLO on the three cancer folds (targets = panel n held coverage): every design + its re-levelled twin.
               --fold-ops prod (default) uses the day's 3-donor operator on every fold -- no extra operator fits; its
               error vs held-out operators is measured in RESULTS.md (section 6). --fold-ops held needs the three
               held-out fits in --held-ops (build_union_whitener.py --donors <the other two>, ~5 min / 13 GB each).
  3. the table the runbook reads: halves, stack, neighbours, q, keys, the H1 x CD4 surface (levelled = the day's
     default, since the re-level is mandatory), and the recommendation with its dump recipe.

Decision rule (fixed before Oct 22; recommend438.py / RESULTS.md section 5): incumbent H1 x2, CD4 x1 (the lowest-regret
point on the rehearsal panels); move only to a grid point that beats it on the LEVELLED k = 1 price by >= 0.0010, on 3/3
folds, and still after leave-top-3-out; q is the board's, vetoed only if the calibrated q 1.5 step is <= 0.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
CD4_G = 0.3443655278649641
H1_G = 0.28


def gate():
    return subprocess.run(["bash", str(HERE / "gate438.sh")], capture_output=True, text=True).returncode == 0


def run(cmd, log):
    for _ in range(240):
        if gate():
            break
        time.sleep(30)
    else:
        raise SystemExit("memory gate never opened")
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.run(["/usr/bin/time", "-l", "uv", "run", "python", *cmd], cwd=HERE, stdout=fh, stderr=subprocess.STDOUT,
                            env={**__import__("os").environ, "PYTHONDONTWRITEBYTECODE": "1"}).returncode
    txt = Path(log).read_text()
    rss = [ln for ln in txt.splitlines() if "maximum resident set size" in ln]
    print(f"  {' '.join(cmd[:3])}: rc {rc}, {time.time() - t0:.0f} s, {rss[0].split()[0] if rss else '?'} bytes max RSS", flush=True)
    if rc:
        raise SystemExit(f"{cmd} failed; see {log}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--op", required=True, type=Path, help="the day's 3-donor whitener (driver stage `whitener`)")
    ap.add_argument("--cd4-root", required=True, type=Path, help="the day's cd4_compact (driver stage `cd4-extract`)")
    ap.add_argument("--tag", default="dayof")
    ap.add_argument("--fold-ops", default="prod", choices=("prod", "held"))
    ap.add_argument("--held-ops", type=Path, default=None, help="dir with held_<fold>.npz (+ .json) for --fold-ops held")
    ap.add_argument("--synthetic", action="store_true", help="the bundle's controls are A/B/C relabelled (rehearsal)")
    ap.add_argument("--skip-runs", action="store_true", help="only re-print the table from existing outputs")
    ap.add_argument("--board-q", default="1.25", choices=("1.25", "1.5"), help="the q A/B/C's board settled (433k: 1.5)")
    a = ap.parse_args()
    tag = a.tag
    (WORK / tag).mkdir(parents=True, exist_ok=True)
    spec = dict(bundle=str(a.bundle.resolve()), rules="def", cd4_root=str(a.cd4_root.resolve()),
                full_op=str(a.op.resolve()), folds=["k562_gwps", "orion_hct116", "orion_hek293t"], synthetic=a.synthetic,
                meta=None, fold_ops={})
    if a.fold_ops == "held":
        spec["fold_ops"] = {h: str((a.held_ops / f"held_{h}.npz").resolve()) for h in spec["folds"]}
    (WORK / tag / "spec.json").write_text(json.dumps(spec, indent=1))
    t0 = time.time()
    if not a.skip_runs:
        print(f"[dayof438] {tag}: full panel", flush=True)
        run(["full438.py", tag, "--tag", tag], WORK / tag / "full438.log")
        lv = str(WORK / tag / "level438.json")
        for f in spec["folds"]:
            print(f"[dayof438] {tag}: fold {f} ({a.fold_ops} operators)", flush=True)
            run(["comp438.py", tag, f, "--fold-ops", a.fold_ops, "--level", lv], WORK / tag / f"comp438-{f}.log")
    sys.path.insert(0, str(HERE))
    import recommend438
    recommend438.main(tag, a.fold_ops, a.board_q)
    print(f"[dayof438] total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
