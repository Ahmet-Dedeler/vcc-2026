#!/usr/bin/env python3
"""Run one command and terminate its process group before it exhausts the Mac.

The VCC scorer can spawn children, so watching only the parent PID is unsafe.
This wrapper sums resident memory over the full descendant tree and also checks
macOS's system-wide free-memory percentage.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


GIB = 1024**3


def process_table() -> tuple[dict[int, tuple[int, int]], dict[int, list[int]]]:
    out = subprocess.check_output(
        ["ps", "-axo", "pid=,ppid=,rss="], text=True
    )
    rows: dict[int, tuple[int, int]] = {}
    children: dict[int, list[int]] = {}
    for line in out.splitlines():
        try:
            pid_s, ppid_s, rss_s = line.split()
            pid, ppid, rss_kib = int(pid_s), int(ppid_s), int(rss_s)
        except ValueError:
            continue
        rows[pid] = (ppid, rss_kib)
        children.setdefault(ppid, []).append(pid)
    return rows, children


def descendant_pids(root_pid: int) -> list[int]:
    rows, children = process_table()
    stack = [root_pid]
    seen: set[int] = set()
    ordered: list[int] = []
    while stack:
        pid = stack.pop()
        if pid in seen or pid not in rows:
            continue
        seen.add(pid)
        ordered.append(pid)
        stack.extend(children.get(pid, ()))
    return ordered


def process_rss(root_pid: int) -> tuple[int, int]:
    """Return (RSS bytes, process count) for root_pid and all descendants."""
    rows, children = process_table()

    stack = [root_pid]
    seen: set[int] = set()
    total_kib = 0
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        row = rows.get(pid)
        if row is None:
            continue
        total_kib += row[1]
        stack.extend(children.get(pid, ()))
    return total_kib * 1024, sum(pid in rows for pid in seen)


def free_percent() -> int | None:
    try:
        out = subprocess.check_output(
            ["memory_pressure", "-Q"],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"free percentage:\s*(\d+)%", out)
    return int(match.group(1)) if match else None


def terminate_group(proc: subprocess.Popen[bytes], reason: str) -> None:
    print(f"\n[memory-guard] stopping command: {reason}", file=sys.stderr, flush=True)
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=15)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def terminate_process_tree(root_pid: int, reason: str) -> None:
    print(f"\n[memory-guard] stopping process tree: {reason}", file=sys.stderr)
    monitor_pid = os.getpid()
    targets = [
        pid
        for pid in reversed(descendant_pids(root_pid))
        if pid != monitor_pid
    ]
    for sig in (signal.SIGTERM, signal.SIGKILL):
        for pid in targets:
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                continue
        if sig == signal.SIGTERM:
            time.sleep(2.0)


def monitor_existing_process(args: argparse.Namespace) -> int:
    if args.monitor_pid <= 1:
        raise ValueError("--monitor-pid must name a non-system process")
    if args.ready_file is None or args.stop_file is None or args.report_file is None:
        raise ValueError(
            "--ready-file, --stop-file, and --report-file are required "
            "with --monitor-pid"
        )
    ready_path = Path(args.ready_file)
    stop_path = Path(args.stop_file)
    report_path = Path(args.report_file)
    started = time.monotonic()
    peak = 0
    samples = 0
    minimum_free: int | None = None

    ready = {
        "status": "active",
        "root_pid": args.monitor_pid,
        "monitor_pid": os.getpid(),
        "max_rss_gib": args.max_rss_gib,
        "min_free_percent": args.min_free_percent,
        "poll_seconds": args.poll_seconds,
    }
    atomic_write_json(ready_path, ready)
    try:
        while True:
            if stop_path.exists():
                rss, count = process_rss(args.monitor_pid)
                peak = max(peak, rss)
                free = free_percent()
                if free is not None:
                    minimum_free = (
                        free if minimum_free is None else min(minimum_free, free)
                    )
                atomic_write_json(
                    report_path,
                    {
                        **ready,
                        "status": "complete",
                        "peak_rss_gib": peak / GIB,
                        "samples": samples + 1,
                        "last_process_count": count,
                        "minimum_free_percent": minimum_free,
                        "elapsed_seconds": time.monotonic() - started,
                    },
                )
                return 0

            rss, count = process_rss(args.monitor_pid)
            if count == 0:
                atomic_write_json(
                    report_path,
                    {
                        **ready,
                        "status": "root_exited_without_stop",
                        "peak_rss_gib": peak / GIB,
                        "samples": samples,
                        "minimum_free_percent": minimum_free,
                        "elapsed_seconds": time.monotonic() - started,
                    },
                )
                return 2
            samples += 1
            peak = max(peak, rss)
            free = free_percent()
            if free is not None:
                minimum_free = (
                    free if minimum_free is None else min(minimum_free, free)
                )
            if not args.quiet:
                print(
                    f"[memory-guard] rss={rss / GIB:.2f} GiB "
                    f"peak={peak / GIB:.2f} GiB procs={count} "
                    f"system_free={free if free is not None else '?'}%",
                    file=sys.stderr,
                    flush=True,
                )
            violation: str | None = None
            if rss > args.max_rss_gib * GIB:
                violation = (
                    f"process tree exceeded {args.max_rss_gib:.1f} GiB RSS"
                )
            elif free is not None and free < args.min_free_percent:
                violation = (
                    "system free memory fell below "
                    f"{args.min_free_percent}%"
                )
            if violation is not None:
                atomic_write_json(
                    report_path,
                    {
                        **ready,
                        "status": "limit_exceeded",
                        "reason": violation,
                        "peak_rss_gib": peak / GIB,
                        "samples": samples,
                        "last_process_count": count,
                        "minimum_free_percent": minimum_free,
                        "elapsed_seconds": time.monotonic() - started,
                    },
                )
                terminate_process_tree(args.monitor_pid, violation)
                return 137
            time.sleep(args.poll_seconds)
    except BaseException as error:
        atomic_write_json(
            report_path,
            {
                **ready,
                "status": "monitor_error",
                "error": repr(error),
                "peak_rss_gib": peak / GIB,
                "samples": samples,
                "minimum_free_percent": minimum_free,
                "elapsed_seconds": time.monotonic() - started,
            },
        )
        terminate_process_tree(args.monitor_pid, f"monitor failure: {error}")
        return 2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rss-gib", type=float, default=14.0)
    parser.add_argument("--min-free-percent", type=int, default=15)
    parser.add_argument("--poll-seconds", type=float, default=3.0)
    parser.add_argument("--monitor-pid", type=int)
    parser.add_argument("--ready-file", type=Path)
    parser.add_argument("--stop-file", type=Path)
    parser.add_argument("--report-file", type=Path)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.max_rss_gib <= 0:
        parser.error("--max-rss-gib must be positive")
    if not 0 <= args.min_free_percent <= 100:
        parser.error("--min-free-percent must lie in 0..100")
    if args.poll_seconds <= 0:
        parser.error("--poll-seconds must be positive")
    if args.monitor_pid is not None:
        if args.command:
            parser.error("monitor mode does not accept a command")
        try:
            return monitor_existing_process(args)
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2

    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("provide a command after --")

    proc = subprocess.Popen(command, start_new_session=True)
    peak = 0
    try:
        while proc.poll() is None:
            rss, count = process_rss(proc.pid)
            peak = max(peak, rss)
            free = free_percent()
            print(
                f"[memory-guard] rss={rss / GIB:.2f} GiB "
                f"peak={peak / GIB:.2f} GiB procs={count} "
                f"system_free={free if free is not None else '?'}%",
                file=sys.stderr,
                flush=True,
            )
            if rss > args.max_rss_gib * GIB:
                terminate_group(
                    proc,
                    f"process tree exceeded {args.max_rss_gib:.1f} GiB RSS",
                )
                return 137
            if free is not None and free < args.min_free_percent:
                terminate_group(
                    proc,
                    f"system free memory fell below {args.min_free_percent}%",
                )
                return 137
            time.sleep(args.poll_seconds)
    except KeyboardInterrupt:
        terminate_group(proc, "interrupted")
        return 130

    print(
        f"[memory-guard] command exited {proc.returncode}; "
        f"peak RSS={peak / GIB:.2f} GiB",
        file=sys.stderr,
        flush=True,
    )
    return int(proc.returncode or 0)


if __name__ == "__main__":
    raise SystemExit(main())
