"""The Oct-22 path as ONE command, stage by stage, timed, fail-closed. Experiment 199.

    preflight -> identity -> cd4-extract -> census -> whitener -> repin-frozen -> pin-artifacts
             -> repin-panel-support -> build -> postbuild -> prep

Every stage is a subprocess wrapped in `/usr/bin/time -l`, so its wall clock and peak RSS are
captured into `<work>/timings.json` as it goes; the run stops at the first non-zero exit.
Nothing here submits: `vcc submit` is never invoked and the driver refuses to run it.

The two heavy stages (`repin-panel-support`, `build`) need ~17-20 GiB (AGENTS.md: the three-donor
model load is the floor; 180 §11). The driver refuses to start them while another builder is
running or while free+inactive memory is under `--min-free-gib`.

Rehearsal (synthetic bundle, throwaway frozen table, throwaway whitener id, builder restored):

    uv run python scripts/def_pipeline.py \
        --bundle experiments/199-def-runbook/work/def-bundle \
        --work experiments/199-def-runbook/work/run1 \
        --frozen-inputs experiments/199-def-runbook/work/run1/frozen_inputs.json \
        --whitener-id union_pruned_def199 --name def199-rehearsal \
        --cd4-thresholds relaxed --restore-builder-after

The round (real bundle, the production table, the builder edit COMMITTED afterwards):

    uv run python scripts/def_pipeline.py --bundle data/challenge \
        --work experiments/<round-dir>/work --name <opaque-name> --whitener-id union_pruned_def

Stage subsets: `--stages preflight,identity,census` or `--from build`. The exact command of
every stage is printed and logged (`<work>/logs/<stage>.log`), so the runbook and the run can
never disagree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGES = ["preflight", "identity", "cd4-extract", "census", "whitener", "repin-frozen",
          "pin-artifacts", "repin-panel-support", "build", "postbuild", "prep"]
# `vcc prep` is heavy too: measured 21.7 GiB max RSS on a 13.2 GB artifact (the CLI reads the
# whole X into memory), so it gets the same "alone, and enough free memory" guard as the builds.
HEAVY = {"repin-panel-support", "build", "prep"}
CD4_SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
CD4_SOURCE_SHA256 = "c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62"
CD4_SOURCE_BYTES = 16_786_240_107
ABC_OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
# The incumbent's build flags, byte-for-byte from experiments/191-rc-hybrid-build/build191.sh
# (reused unchanged by 201 and 202). The env pair is exported for the effects run too.
INCUMBENT_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.1"}   # the ONLY env the incumbent sets
INCUMBENT_FLAGS = ["--min-cpm", "10", "--tail-spread", "1.0", "--keep-top", "8000",
                   "--keep-spread", "0.25", "--overdispersion", "4.0"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def free_gib() -> float:
    out = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    pages = {}
    for line in out.splitlines():
        m = re.match(r"Pages (free|inactive|speculative|purgeable):\s+(\d+)\.", line)
        if m:
            pages[m.group(1)] = int(m.group(2))
    return sum(pages.values()) * 16384 / 2**30


def other_builders() -> list[str]:
    """Python processes actually executing the builder or the operator fit -- not shells whose
    command line merely mentions them (the first rehearsal tripped on its own launcher)."""
    out = subprocess.run(["ps", "-axo", "pid=,args="], capture_output=True, text=True).stdout
    hits = []
    for ln in out.splitlines():
        parts = ln.split()
        if len(parts) < 3 or not re.search(r"python[0-9.]*$", parts[1]):
            continue
        if any(a.endswith(("scripts/make_submission_pm_cd4.py", "scripts/build_union_whitener.py"))
               for a in parts[2:4]):
            hits.append(ln.strip())
    return hits


class Runner:
    def __init__(self, work: Path, dry: bool):
        self.work, self.dry = work, dry
        self.logs = work / "logs"
        self.logs.mkdir(parents=True, exist_ok=True)
        self.timings_path = work / "timings.json"
        self.timings = json.loads(self.timings_path.read_text()) if self.timings_path.exists() else {}

    def run(self, stage: str, cmd: list[str], env: dict[str, str] | None = None) -> None:
        assert "submit" not in cmd, "this driver never submits"
        full = ["/usr/bin/time", "-l", *cmd]
        shown = " ".join(f"{k}={v}" for k, v in (env or {}).items()) + (" " if env else "") + " ".join(map(str, cmd))
        print(f"\n=== [{stage}] {shown}", flush=True)
        if self.dry:
            return
        log = self.logs / f"{stage}.log"
        t0 = time.time()
        with open(log, "w") as fh:
            fh.write(f"$ {shown}\n\n")
            fh.flush()
            proc = subprocess.run(full, cwd=ROOT, env={**os.environ, **(env or {})},
                                  stdout=fh, stderr=subprocess.STDOUT)
        wall = time.time() - t0
        text = log.read_text()
        m = re.search(r"(\d+)\s+maximum resident set size", text)
        rss = int(m.group(1)) / 2**30 if m else None
        self.timings[stage] = {"seconds": round(wall, 1), "peak_rss_gib": round(rss, 2) if rss else None,
                               "returncode": proc.returncode, "log": str(log.relative_to(ROOT)),
                               "command": shown}
        self.timings_path.write_text(json.dumps(self.timings, indent=2) + "\n")
        tail = "\n".join(ln for ln in text.splitlines()[-25:]
                         if not re.match(r"^\s*\d+\s+(average|page|swaps|block|messages|signals|voluntary|involuntary|instructions|cycles|peak memory)", ln))
        print(tail, flush=True)
        print(f"--- [{stage}] {wall:.1f}s, peak RSS {rss:.2f} GiB, rc={proc.returncode}" if rss
              else f"--- [{stage}] {wall:.1f}s rc={proc.returncode}", flush=True)
        if proc.returncode != 0:
            raise SystemExit(f"stage {stage} failed (rc={proc.returncode}); see {log}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, default=ROOT / "data/challenge")
    ap.add_argument("--work", type=Path, required=True, help="everything this run writes goes here")
    ap.add_argument("--name", required=True, help="artifact stem; the .h5ad is <work>/<name>.h5ad")
    ap.add_argument("--frozen-inputs", type=Path, default=ROOT / "scripts/frozen_inputs.json",
                    help="the frozen table to re-pin and build against (rehearsal: a throwaway copy)")
    ap.add_argument("--whitener-id", default="union_pruned_def")
    ap.add_argument("--whitener-out", type=Path, default=None,
                    help="default <work>/whitener_<id>.npz; pass an existing fit to skip the stage")
    ap.add_argument("--cd4-source", type=Path, default=CD4_SOURCE)
    ap.add_argument("--cd4-out", type=Path, default=None, help="default <work>/cd4_compact")
    ap.add_argument("--cd4-thresholds", choices=["default", "relaxed"], default="default",
                    help="047's gates (275 per condition, 290 union) or none. Relax ONLY for a synthetic panel")
    ap.add_argument("--emitter", default="rc_callset",
                    help="the incumbent's emitter: rc_callset (202, board 0.1918638, entry 31oB9yWy8CSNBDEkY4iW); "
                         "rc_hybrid is 191's (0.1834190)")
    ap.add_argument("--incumbent-meta", type=Path, default=ROOT / "submissions/202-rc_callset.h5ad.meta.json",
                    help="the meta whose config the D/E/F build must differ from in exactly contexts + whitener_id")
    ap.add_argument("--max-rss-gib", type=float, default=24.0, help="builder guard, as build191.sh")
    ap.add_argument("--min-free-gib", type=float, default=26.0)
    ap.add_argument("--stages", default=None, help="comma-separated subset, in pipeline order")
    ap.add_argument("--from", dest="from_stage", default=None)
    ap.add_argument("--skip-source-hash", action="store_true", help="skip the 30 s sha256 of the CD4 source")
    ap.add_argument("--restore-builder-after", action="store_true",
                    help="REHEARSAL ONLY: byte-restore the builder from git HEAD after `prep`")
    ap.add_argument("--dry-run", action="store_true", help="print every command, run nothing")
    args = ap.parse_args()

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    bundle = args.bundle.resolve()
    manifest = json.loads((bundle / "manifest.json").read_text())
    contexts = list(manifest["contexts"])
    ctx = ",".join(contexts)
    is_abc = contexts == ["A", "B", "C"]
    frozen = args.frozen_inputs.resolve()
    wout = (args.whitener_out or work / f"whitener_{args.whitener_id}.npz").resolve()
    wmeta = wout.with_suffix(".json")
    cd4_out = (args.cd4_out or work / "cd4_compact").resolve()
    out_h5ad = work / f"{args.name}.h5ad"
    out_meta = Path(str(out_h5ad) + ".meta.json")
    py = ["uv", "run", "python"]

    stages = STAGES if not args.stages else [s.strip() for s in args.stages.split(",")]
    for s in stages:
        if s not in STAGES:
            ap.error(f"unknown stage {s}; order is {STAGES}")
    if args.from_stage:
        stages = [s for s in STAGES if STAGES.index(s) >= STAGES.index(args.from_stage)]
    r = Runner(work, args.dry_run)
    print(f"bundle {bundle}  panel {manifest.get('panel_id')}  contexts {contexts}\n"
          f"work {work}\nfrozen table {frozen}\nwhitener {wout} (id {args.whitener_id})\n"
          f"cd4 banks {cd4_out}\nartifact {out_h5ad}\nstages {stages}", flush=True)
    t_all = time.time()

    for stage in stages:
        if stage == "preflight":
            r.run(stage, [*py, "scripts/def_preflight.py", "--bundle", str(bundle),
                          *(["--verify"] if is_abc else []), "--out", str(work / "preflight.json")])
        elif stage == "identity":
            r.run(stage, [*py, "scripts/def_identity.py", "--bundle", str(bundle),
                          *(["--verify"] if is_abc else ["--contexts", ctx]),
                          "--out", str(work / "identity.json")])
        elif stage == "census":
            cd4_root = cd4_out if (cd4_out / "cd4_rest" / "perts.txt").exists() else None
            r.run(stage, [*py, "scripts/def_census.py", "--bundle", str(bundle),
                          *(["--verify"] if is_abc else []),
                          *(["--cd4-root", str(cd4_root)] if cd4_root else []),
                          "--out", str(work / "census.json")])
        elif stage == "cd4-extract":
            src = args.cd4_source.resolve()
            if not src.is_file() or src.stat().st_size != CD4_SOURCE_BYTES:
                raise SystemExit(f"CD4 source missing or wrong size: {src} (expect {CD4_SOURCE_BYTES} bytes)")
            if not args.skip_source_hash and not args.dry_run:
                t0 = time.time()
                got = sha256(src)
                if got != CD4_SOURCE_SHA256:
                    raise SystemExit(f"CD4 source sha256 {got} != pinned {CD4_SOURCE_SHA256}")
                print(f"CD4 source sha256 verified in {time.time() - t0:.0f}s", flush=True)
            thr = [] if args.cd4_thresholds == "default" else ["--min-targets", "0", "--min-union-targets", "0"]
            r.run(stage, [*py, "experiments/047-cd4-direct/extract_cd4.py", "--input", str(src),
                          "--output", str(cd4_out), "--panel-csv", str(bundle / "pert_counts.csv"),
                          "--genes-csv", str(bundle / "gene_names.csv"), "--overwrite", *thr])
        elif stage == "whitener":
            if wout.exists() and wmeta.exists():
                print(f"\n=== [whitener] exists, skipping: {wout}", flush=True)
                continue
            r.run(stage, [*py, "scripts/build_union_whitener.py", "--bundle", str(bundle),
                          "--basis-prune", "essential", "--out", str(wout),
                          "--scratch", str(work / ".whitener_scratch"), "--max-rss-gib", "14"])
            if not args.dry_run and not is_abc and sha256(wout) == ABC_OPERATOR_SHA:
                raise SystemExit("the re-fit operator is byte-identical to the A/B/C one: the panel did not change")
        elif stage == "repin-frozen":
            if not frozen.exists():
                shutil.copyfile(ROOT / "scripts/frozen_inputs.json", frozen)
                print(f"seeded throwaway table {frozen} from scripts/frozen_inputs.json", flush=True)
            r.run(stage, [*py, "scripts/make_submission_pm_cd4.py", "--repin-frozen-inputs",
                          "--bundle", str(bundle), "--frozen-inputs", str(frozen)])
        elif stage == "pin-artifacts":
            r.run(stage, [*py, "scripts/def_pin_panel_artifacts.py", "--whitener", str(wout),
                          "--whitener-metadata", str(wmeta), "--whitener-id", args.whitener_id,
                          "--cd4-root", str(cd4_out), "--apply", "--out", str(work / "pin-artifacts.json")])
        elif stage in ("repin-panel-support", "build"):
            if not args.dry_run:
                others = other_builders()
                if others:
                    raise SystemExit("another builder is running; refusing to start a ~19 GiB stage:\n  " + "\n  ".join(others))
                fg = free_gib()
                if fg < args.min_free_gib:
                    raise SystemExit(f"free+inactive memory {fg:.1f} GiB < {args.min_free_gib} GiB; refusing")
                print(f"memory check: {fg:.1f} GiB free+inactive, no other builder", flush=True)
            common = [*py, "scripts/make_submission_pm_cd4.py", "--bundle", str(bundle),
                      "--frozen-inputs", str(frozen), "--expect-contexts", ctx,
                      "--whitener", str(wout), "--whitener-metadata", str(wmeta),
                      "--whitener-id", args.whitener_id, "--cd4-root", str(cd4_out),
                      "--max-rss-gib", str(args.max_rss_gib)]
            if stage == "repin-panel-support":
                r.run(stage, [*common, "--repin-panel-support", "--min-cpm", "10"], INCUMBENT_ENV)
            else:
                if out_h5ad.exists():
                    raise SystemExit(f"{out_h5ad} exists; the builder refuses to overwrite. Pick a new --name")
                r.run(stage, [*common, *INCUMBENT_FLAGS, "--emitter", args.emitter, "-o", str(out_h5ad)],
                      INCUMBENT_ENV)
        elif stage == "postbuild":
            r.run(stage, [*py, "scripts/def_postbuild.py", "--meta", str(out_meta), "--bundle", str(bundle),
                          "--incumbent", str(args.incumbent_meta.resolve()), "--cd4-root", str(cd4_out),
                          "--allow-config-diff", "whitener_id", "--out", str(work / "postbuild.json")])
        elif stage == "prep":
            if not args.dry_run:
                others = other_builders()
                if others:
                    raise SystemExit("another builder is running; refusing to start vcc prep (~22 GiB):\n  " + "\n  ".join(others))
                fg = free_gib()
                if fg < args.min_free_gib:
                    raise SystemExit(f"free+inactive memory {fg:.1f} GiB < {args.min_free_gib} GiB; refusing")
            r.run(stage, ["vcc", "prep", str(out_h5ad), "-g", str(bundle / "gene_names.csv"),
                          "--perts", str(bundle / "pert_counts.csv"), "--contexts", ctx,
                          "-o", str(work / f"{args.name}.vcc"), "--dry-run"])
            if args.restore_builder_after:
                r.run("restore-builder", [*py, "scripts/def_pin_panel_artifacts.py", "--restore",
                                          "--record", str(work / "pin-artifacts.json")])

    total = time.time() - t_all
    if not args.dry_run:
        r.timings["_total"] = {"seconds": round(total, 1), "stages": stages}
        r.timings_path.write_text(json.dumps(r.timings, indent=2) + "\n")
        print("\nstage                   seconds   peak RSS GiB")
        for s, t in r.timings.items():
            if s.startswith("_"):
                continue
            print(f"{s:<22} {t['seconds']:>9.1f}   {t.get('peak_rss_gib') or 0:>6.2f}")
        print(f"{'TOTAL (this run)':<22} {total:>9.1f}")
        print(f"timings -> {r.timings_path}")


if __name__ == "__main__":
    main()
