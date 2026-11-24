"""Run one command while enforcing a process-tree RSS ceiling."""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import time
from pathlib import Path


def process_tree_rss_bytes(root_pid: int) -> int:
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
        pending.extend(children.get(pid, []))
    return total_kib * 1024


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
    parser.add_argument("--poll-seconds", type=float, default=2.0)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        raise SystemExit("missing command")
    limit_bytes = int(args.max_rss_gib * 1024**3)
    process = subprocess.Popen(command, start_new_session=True)
    peak_bytes = 0
    killed = False
    while process.poll() is None:
        try:
            rss_bytes = process_tree_rss_bytes(process.pid)
        except subprocess.CalledProcessError:
            rss_bytes = 0
        peak_bytes = max(peak_bytes, rss_bytes)
        if rss_bytes >= limit_bytes:
            killed = True
            terminate_tree(process)
            break
        time.sleep(args.poll_seconds)
    return_code = process.wait()
    report = {
        "command": command,
        "max_rss_gib": args.max_rss_gib,
        "peak_process_tree_rss_gib": peak_bytes / 1024**3,
        "killed_for_memory": killed,
        "return_code": return_code,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    if killed:
        return 137
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
