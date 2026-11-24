#!/usr/bin/env python3
"""Run one command under a process-tree RSS ceiling and write an audit report."""

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


def process_tree_rss_bytes(root_pid: int) -> tuple[int, int]:
    result = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss="],
        check=True,
        capture_output=True,
        text=True,
    )
    rows: dict[int, tuple[int, int]] = {}
    children: dict[int, list[int]] = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        pid, parent, rss_kib = map(int, fields)
        rows[pid] = (parent, rss_kib)
        children.setdefault(parent, []).append(pid)
    pending = [root_pid]
    seen: set[int] = set()
    total_kib = 0
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        if pid in rows:
            total_kib += rows[pid][1]
        pending.extend(children.get(pid, ()))
    return total_kib * 1024, len(seen)


def terminate_tree(process: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=10)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rss-gib", type=float, required=True)
    parser.add_argument("--poll-seconds", type=float, default=0.25)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        raise SystemExit("missing command")
    if args.max_rss_gib <= 0:
        raise SystemExit("--max-rss-gib must be positive")

    process = subprocess.Popen(command, start_new_session=True)
    peak_bytes = 0
    killed = False
    while process.poll() is None:
        try:
            rss_bytes, process_count = process_tree_rss_bytes(process.pid)
        except subprocess.CalledProcessError:
            rss_bytes, process_count = 0, 0
        peak_bytes = max(peak_bytes, rss_bytes)
        print(
            f"[memory-guard] rss={rss_bytes / GIB:.3f} GiB "
            f"peak={peak_bytes / GIB:.3f} GiB processes={process_count}",
            file=sys.stderr,
            flush=True,
        )
        if rss_bytes >= args.max_rss_gib * GIB:
            killed = True
            terminate_tree(process)
            break
        time.sleep(args.poll_seconds)
    return_code = process.wait()
    report = {
        "command": command,
        "max_rss_gib": args.max_rss_gib,
        "peak_process_tree_rss_gib": peak_bytes / GIB,
        "killed_for_memory": killed,
        "return_code": 137 if killed else return_code,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    return int(report["return_code"] or 0)


if __name__ == "__main__":
    raise SystemExit(main())
