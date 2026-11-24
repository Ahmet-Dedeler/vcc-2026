#!/usr/bin/env python3
"""Capture a scored entry's COMPLETE `vcc --json status` record to disk before it disappears.

    uv run python scripts/capture_board_entry.py <entry-id> [--out PATH] [--package PKG --sha256 HEX]
                                                 [--wait] [--attempts N] [--interval S] [--note TEXT]

Why (2026-09-05, written ~35 min before the 237 slot): the board API resolves ONLY the team's
current entry.  Every earlier id returns `not_found` (experiments/253-*/work/vcc-status.log: six of
six).  `vcc submit --wait` prints the six SCALED members but not the raw ones (`pds_cosine`,
`expr_mse_unbiased_capped_norm`, ...), and `scripts/arm_submit_after.sh` used to fetch the full
record only to grep one bit out of it.  So the first slot's raw `mse` -- the one number 243's
class-model reading needs -- would have been lost the moment the second slot published.

What it does: runs `vcc --json status <id>`, validates that the record is a TERMINAL entry for that
id (and, when `published`, that all twelve score fields and the panel/anchor stamp are present),
and writes an atomic JSON wrapper:

    {"captured_at": ..., "vcc_version": ..., "package": ..., "sha256": ..., "note": ..., "entry": {...}}

Default output: data/board-entries/<entry-id>.json (tracked by git; `submissions/*` is not).  An
existing VALID capture is never overwritten by a worse one.  Exit 0 only when a valid terminal
record is on disk; anything else exits 1 so a caller can fail closed.

Final phase (531 finding 8, 2026-09-29): a FINAL-phase entry (`is_final` true, or a partition other than `val`) that is
terminal with EVERY score field null is the "members hidden" state of the dark final board, not a transient: it is
captured with `members_hidden: true` in the wrapper and exits 0, so callers do not retry it as a failure. A validation
(A/B/C) entry with null scores still fails closed.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_DIR = REPO / "data" / "board-entries"

TERMINAL = {"published", "failed", "hidden_admin", "superseded"}
RAW = ("pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
       "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
       "de_wilcoxon_sig_jaccard")
SCALED = ("score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac")
STAMP = ("panel_id", "anchor_version", "partition")


def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def fetch(entry_id: str) -> tuple[dict | None, str]:
    """Return (entry-or-None, diagnostic).  Never raises on a CLI failure."""
    try:
        p = subprocess.run(["vcc", "--json", "status", entry_id], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, f"vcc status could not run: {exc}"
    out = (p.stdout or "").strip()
    err = (p.stderr or "").strip()
    if not out:
        return None, f"no JSON on stdout (rc={p.returncode}); stderr: {err[:300]}"
    try:
        entry = json.loads(out)
    except json.JSONDecodeError as exc:
        return None, f"stdout is not JSON ({exc}); head: {out[:200]}"
    if not isinstance(entry, dict):
        return None, "stdout JSON is not an object"
    if "error" in entry and "status" not in entry:
        return None, f"API error: {entry.get('code')}: {entry.get('error')}"
    return entry, f"rc={p.returncode}"


def members_hidden(entry: dict) -> bool:
    """A terminal FINAL-phase record whose every score field is null: the dark final board (531 finding 8)."""
    final = entry.get("is_final") is True or (entry.get("partition") not in (None, "", "val"))
    return (final and entry.get("status") in ("published", "superseded")
            and entry.get("score_avg") is None and all(entry.get(k) is None for k in RAW + SCALED))


def validate(entry: dict, entry_id: str) -> list[str]:
    """Problems that make this record unusable as THE record of the entry.  Empty = valid."""
    problems: list[str] = []
    if entry.get("entry_id") != entry_id:
        problems.append(f"entry_id mismatch: {entry.get('entry_id')!r} != {entry_id!r}")
    status = entry.get("status")
    if status not in TERMINAL:
        problems.append(f"status {status!r} is not terminal ({sorted(TERMINAL)})")
    if members_hidden(entry):
        return problems                      # the dark final board: nothing more to require
    if status in ("published", "superseded"):
        if not _num(entry.get("score_avg")):
            problems.append("score_avg missing")
        for k in RAW + SCALED:
            if not _num(entry.get(k)):
                problems.append(f"{k} missing or non-numeric")
        for k in STAMP:
            if not entry.get(k):
                problems.append(f"{k} missing (scores are not comparable across panels)")
        if status == "published" and not _num(entry.get("rank")):
            problems.append("rank missing")
    return problems


def existing_valid(path: Path, entry_id: str) -> bool:
    try:
        old = json.loads(path.read_text())
        return not validate(old["entry"], entry_id)
    except Exception:
        return False


def write_atomic(out: Path, wrapper: dict) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f"{out.name}.{os.getpid()}.tmp")  # pid-unique: two capturers must not share a tmp
    tmp.write_text(json.dumps(wrapper, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, out)


def fill_metadata(out: Path, a) -> None:
    """Backfill package / sha256 / note on an existing valid capture when this caller knows them and
    the file does not (the second-slot handoff captures the prior entry without them). Never
    overwrites a differing sha256: that is reported, not silently replaced. The entry is untouched."""
    try:
        wrapper = json.loads(out.read_text())
    except Exception:
        return
    changed = False
    for key in ("package", "sha256", "note"):
        new = getattr(a, key, None)
        old = wrapper.get(key)
        if new and not old:
            wrapper[key] = new
            changed = True
        elif new and old and key == "sha256" and new != old:
            print(f"WARNING: capture {out} records sha256 {old}, caller says {new}; left as recorded", file=sys.stderr)
    if changed:
        try:
            write_atomic(out, wrapper)
            print(f"backfilled metadata on {out}")
        except Exception as exc:
            print(f"WARNING: could not backfill metadata on {out}: {exc}", file=sys.stderr)


def vcc_version() -> str:
    try:
        return subprocess.run(["vcc", "--version"], capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return "unknown"


def summary(entry: dict) -> str:
    if members_hidden(entry):
        return f"status={entry.get('status')} MEMBERS HIDDEN (final phase) partition={entry.get('partition')} is_final={entry.get('is_final')}"
    if entry.get("score_avg") is None:
        return f"status={entry.get('status')} rank={entry.get('rank')} error_info={entry.get('error_info')}"
    s = [f"status={entry['status']} rank={entry.get('rank')} overall={entry['score_avg']:.7f}"]
    for r, c in zip(RAW, SCALED):
        s.append(f"{c[6:]}={entry[r]:.5f}/{entry[c]:+.5f}")
    s.append(f"panel={entry.get('panel_id')} anchors={entry.get('anchor_version')}")
    return " ".join(s)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entry_id")
    ap.add_argument("--out", type=Path, default=None, help="default data/board-entries/<entry-id>.json")
    ap.add_argument("--package", default=None, help="package path, recorded in the wrapper")
    ap.add_argument("--sha256", default=None, help="package sha256, recorded in the wrapper")
    ap.add_argument("--note", default="", help="free text recorded in the wrapper")
    ap.add_argument("--wait", action="store_true", help="keep polling while the entry is not terminal")
    ap.add_argument("--attempts", type=int, default=20, help="tries on API/validation failure (default 20)")
    ap.add_argument("--interval", type=float, default=30.0, help="seconds between tries (default 30)")
    a = ap.parse_args()

    out = a.out or (DEFAULT_DIR / f"{a.entry_id}.json")
    if existing_valid(out, a.entry_id):
        fill_metadata(out, a)
        print(f"already captured (valid, terminal): {out}")
        return 0

    last = "never fetched"
    for attempt in range(1, max(1, a.attempts) + 1):
        entry, diag = fetch(a.entry_id)
        if entry is None:
            last = diag
            print(f"attempt {attempt}/{a.attempts}: {diag}", file=sys.stderr)
        else:
            problems = validate(entry, a.entry_id)
            if not problems:
                wrapper = {
                    "captured_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                    "vcc_version": vcc_version(),
                    "package": a.package,
                    "sha256": a.sha256,
                    "note": a.note,
                    "members_hidden": members_hidden(entry),
                    "entry": entry,
                }
                try:
                    write_atomic(out, wrapper)
                except Exception as exc:
                    # Two capturers (the second-slot handoff and the coordinator's backup watcher) may race
                    # on the same target. The tmp name is pid-unique, but whatever happens, the file on
                    # disk decides: a valid capture written by the other process is a success here too.
                    if existing_valid(out, a.entry_id):
                        fill_metadata(out, a)
                        print(f"captured by a concurrent writer ({exc.__class__.__name__}); valid on disk: {out}")
                        return 0
                    print(f"attempt {attempt}/{a.attempts}: write failed: {exc}", file=sys.stderr)
                    last = f"write failed: {exc}"
                    if attempt < a.attempts:
                        time.sleep(a.interval)
                    continue
                # Re-read what we wrote: the file on disk is the evidence, not the variable.
                if not existing_valid(out, a.entry_id):
                    print(f"ERROR: wrote {out} but it does not re-validate", file=sys.stderr)
                    return 1
                print(f"captured {out}: {summary(entry)}")
                return 0
            last = "; ".join(problems)
            not_terminal = entry.get("status") not in TERMINAL
            print(f"attempt {attempt}/{a.attempts}: {summary(entry)} -- {last}", file=sys.stderr)
            if not_terminal and not a.wait:
                print("not terminal and --wait not given; giving up", file=sys.stderr)
                return 1
        if attempt < a.attempts:
            time.sleep(a.interval)
    if existing_valid(out, a.entry_id):  # the other capturer may have succeeded meanwhile
        fill_metadata(out, a)
        print(f"captured by a concurrent writer; valid on disk: {out}")
        return 0
    print(f"ERROR: no valid terminal record for {a.entry_id} after {a.attempts} attempts: {last}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
