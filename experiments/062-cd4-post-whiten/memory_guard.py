#!/usr/bin/env python3
"""Linux/macOS process-tree RSS guard with a machine-readable report."""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

GIB = 1024**3


def process_tree_rss(root_pid: int) -> tuple[int, int]:
    output = subprocess.check_output(
        ["ps", "-axo", "pid=,ppid=,rss="], text=True
    )
    rows: dict[int, tuple[int, int]] = {}
    children: dict[int, list[int]] = {}
    for line in output.splitlines():
        try:
            pid_text, parent_text, rss_text = line.split()
            pid, parent, rss_kib = (
                int(pid_text),
                int(parent_text),
                int(rss_text),
            )
        except ValueError:
            continue
        rows[pid] = (parent, rss_kib)
        children.setdefault(parent, []).append(pid)
    stack = [root_pid]
    seen: set[int] = set()
    total_kib = 0
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        if pid in rows:
            total_kib += rows[pid][1]
        stack.extend(children.get(pid, ()))
    return total_kib * 1024, len(seen)


def stop_group(process: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rss-gib", type=float, required=True)
    parser.add_argument("--poll-seconds", type=float, default=2.0)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("provide a command after --")
    process = subprocess.Popen(command, start_new_session=True)
    peak = 0
    killed = False
    while process.poll() is None:
        rss, count = process_tree_rss(process.pid)
        peak = max(peak, rss)
        print(
            f"[memory-guard] rss={rss / GIB:.3f} GiB "
            f"peak={peak / GIB:.3f} GiB processes={count}",
            file=sys.stderr,
            flush=True,
        )
        if rss > args.max_rss_gib * GIB:
            killed = True
            stop_group(process)
            break
        time.sleep(args.poll_seconds)
    return_code = process.wait()
    report = {
        "max_rss_gib": args.max_rss_gib,
        "peak_process_tree_rss_gib": peak / GIB,
        "return_code": 137 if killed else return_code,
        "killed_for_memory": killed,
        "command": command,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True), file=sys.stderr, flush=True)
    return int(report["return_code"])


if __name__ == "__main__":
    raise SystemExit(main())

