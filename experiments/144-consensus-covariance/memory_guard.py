#!/usr/bin/env python3
"""Fail-closed external memory guard for Experiment 144.

The guarded command starts in a new process group. Monitoring failures, an
unavailable free-memory measurement, limit violations, interrupts, and a root
process that exits while descendants remain all terminate that process group.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import secrets
import select
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Sequence


GIB = 1024**3
ACTIVE_ENV = "VCC144_MEMORY_GUARD_ACTIVE"
PATH_ENV = "VCC144_MEMORY_GUARD_PATH"
MAX_RSS_ENV = "VCC144_MEMORY_GUARD_MAX_RSS_GIB"
MIN_FREE_ENV = "VCC144_MEMORY_GUARD_MIN_FREE_PERCENT"
REQUEST_FD_ENV = "VCC144_MEMORY_GUARD_REQUEST_FD"
RESPONSE_FD_ENV = "VCC144_MEMORY_GUARD_RESPONSE_FD"
TOKEN_ENV = "VCC144_MEMORY_GUARD_TOKEN"


def process_table() -> list[tuple[int, int, int, int]]:
    output = subprocess.check_output(
        ["ps", "-axo", "pid=,ppid=,pgid=,rss="],
        text=True,
    )
    rows: list[tuple[int, int, int, int]] = []
    for line in output.splitlines():
        try:
            pid_s, ppid_s, pgid_s, rss_s = line.split()
            rows.append((int(pid_s), int(ppid_s), int(pgid_s), int(rss_s)))
        except ValueError:
            continue
    return rows


def process_rss(root_pid: int) -> tuple[int, int]:
    """Return RSS bytes and count for a process and all descendants."""
    rows = process_table()
    by_pid = {pid: (ppid, rss_kib) for pid, ppid, _, rss_kib in rows}
    children: dict[int, list[int]] = {}
    for pid, ppid, _, _ in rows:
        children.setdefault(ppid, []).append(pid)

    stack = [root_pid]
    seen: set[int] = set()
    total_kib = 0
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        row = by_pid.get(pid)
        if row is None:
            continue
        total_kib += row[1]
        stack.extend(children.get(pid, ()))
    return total_kib * 1024, sum(pid in by_pid for pid in seen)


def process_group_members(process_group: int) -> list[tuple[int, int]]:
    """Return ``(pid, RSS KiB)`` for the isolated guarded process group."""
    return [
        (pid, rss_kib)
        for pid, _, pgid, rss_kib in process_table()
        if pgid == process_group
    ]


def process_group_rss(process_group: int) -> tuple[int, int]:
    """Return RSS bytes and count for the isolated guarded process group."""
    members = process_group_members(process_group)
    return sum(rss_kib for _, rss_kib in members) * 1024, len(members)


def free_percent() -> int | None:
    try:
        output = subprocess.check_output(
            ["memory_pressure", "-Q"],
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"free percentage:\s*(\d+)%", output)
    return int(match.group(1)) if match else None


def measured_free_percent() -> int:
    value = free_percent()
    if value is None:
        raise RuntimeError("cannot measure system free-memory percentage")
    return value


def terminate_group(
    process: subprocess.Popen[bytes],
    reason: str,
) -> None:
    print(
        f"\n[memory-guard] stopping command: {reason}",
        file=sys.stderr,
        flush=True,
    )
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.poll()
        return
    if wait_for_group_exit(process, 15.0):
        return

    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        process.poll()
        return
    if wait_for_group_exit(process, 5.0):
        return
    raise RuntimeError(f"could not terminate guarded process group {process.pid}")


def process_group_exists(process_group: int) -> bool:
    try:
        os.killpg(process_group, 0)
    except ProcessLookupError:
        return False
    return True


def wait_for_group_exit(
    process: subprocess.Popen[bytes],
    timeout_seconds: float,
) -> bool:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        process.poll()
        if not process_group_exists(process.pid):
            return True
        time.sleep(0.05)
    process.poll()
    return not process_group_exists(process.pid)


def service_handshake(
    request_fd: int,
    response_fd: int,
    expected_token: str,
) -> bool:
    readable, _, _ = select.select([request_fd], [], [], 0.0)
    if not readable:
        return False
    request = os.read(request_fd, len(expected_token) + 2)
    if request != f"{expected_token}\n".encode():
        raise RuntimeError("invalid memory-guard handshake request")
    response = f"{expected_token}\n".encode()
    written = os.write(response_fd, response)
    if written != len(response):
        raise RuntimeError("incomplete memory-guard handshake response")
    return True


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


def write_report(
    path: Path | None,
    *,
    status: str,
    started: float,
    peak_rss: int,
    minimum_free: int | None,
    samples: int,
    returncode: int | None = None,
    reason: str | None = None,
) -> None:
    if path is None:
        return
    payload: dict[str, Any] = {
        "schema_version": 1,
        "status": status,
        "peak_rss_gib": peak_rss / GIB,
        "minimum_free_percent": minimum_free,
        "samples": samples,
        "elapsed_seconds": time.monotonic() - started,
    }
    if returncode is not None:
        payload["returncode"] = returncode
    if reason is not None:
        payload["reason"] = reason
    atomic_write_json(path, payload)


def monitor_process(
    process: subprocess.Popen[bytes],
    *,
    max_rss_gib: float,
    min_free_percent: float,
    poll_seconds: float,
    report_file: Path | None,
    quiet: bool,
    request_fd: int | None = None,
    response_fd: int | None = None,
    handshake_token: str | None = None,
) -> int:
    started = time.monotonic()
    peak = 0
    samples = 0
    minimum_free: int | None = None
    handshake_complete = request_fd is None
    require_handshake = request_fd is not None
    if require_handshake:
        if response_fd is None or handshake_token is None:
            raise ValueError("incomplete handshake configuration")
    try:
        while True:
            if not handshake_complete:
                handshake_complete = service_handshake(
                    request_fd,
                    response_fd,
                    handshake_token,
                )
            members = process_group_members(process.pid)
            rss = sum(rss_kib for _, rss_kib in members) * 1024
            count = len(members)
            root_returncode = process.poll()
            if root_returncode is not None:
                if require_handshake and not handshake_complete:
                    reason = "guarded command exited before live guard handshake"
                    terminate_group(process, reason)
                    write_report(
                        report_file,
                        status="monitor_error",
                        started=started,
                        peak_rss=peak,
                        minimum_free=minimum_free,
                        samples=samples,
                        returncode=2,
                        reason=reason,
                    )
                    return 2
                descendants = [pid for pid, _ in members if pid != process.pid]
                if descendants:
                    reason = (
                        "guarded root exited while process-group descendants remained"
                    )
                    terminate_group(process, reason)
                    write_report(
                        report_file,
                        status="monitor_error",
                        started=started,
                        peak_rss=peak,
                        minimum_free=minimum_free,
                        samples=samples,
                        returncode=2,
                        reason=reason,
                    )
                    return 2
                write_report(
                    report_file,
                    status="complete",
                    started=started,
                    peak_rss=peak,
                    minimum_free=minimum_free,
                    samples=samples,
                    returncode=int(root_returncode),
                )
                return int(root_returncode)

            if count == 0:
                raise RuntimeError("guard lost the live process group")

            free = measured_free_percent()
            samples += 1
            peak = max(peak, rss)
            minimum_free = free if minimum_free is None else min(minimum_free, free)
            if not quiet:
                print(
                    f"[memory-guard] rss={rss / GIB:.2f} GiB "
                    f"peak={peak / GIB:.2f} GiB procs={count} "
                    f"system_free={free}%",
                    file=sys.stderr,
                    flush=True,
                )

            violation: str | None = None
            if rss >= max_rss_gib * GIB:
                violation = f"process group reached {max_rss_gib:.1f} GiB RSS"
            elif free < min_free_percent:
                violation = f"system free memory fell below {min_free_percent:.1f}%"
            if violation is not None:
                terminate_group(process, violation)
                write_report(
                    report_file,
                    status="limit_exceeded",
                    started=started,
                    peak_rss=peak,
                    minimum_free=minimum_free,
                    samples=samples,
                    returncode=137,
                    reason=violation,
                )
                return 137
            time.sleep(poll_seconds)
    except KeyboardInterrupt:
        terminate_group(process, "interrupted")
        write_report(
            report_file,
            status="interrupted",
            started=started,
            peak_rss=peak,
            minimum_free=minimum_free,
            samples=samples,
            returncode=130,
        )
        return 130
    except BaseException as error:
        reason = f"monitor failure: {error!r}"
        terminate_group(process, reason)
        try:
            write_report(
                report_file,
                status="monitor_error",
                started=started,
                peak_rss=peak,
                minimum_free=minimum_free,
                samples=samples,
                returncode=2,
                reason=reason,
            )
        except BaseException as report_error:
            print(
                f"[memory-guard] could not write failure report: {report_error!r}",
                file=sys.stderr,
                flush=True,
            )
        return 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rss-gib", type=float, required=True)
    parser.add_argument("--min-free-percent", type=float, required=True)
    parser.add_argument("--poll-seconds", type=float, default=0.25)
    parser.add_argument("--report-file", type=Path)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.max_rss_gib <= 0.0:
        raise SystemExit("--max-rss-gib must be positive")
    if not 0.0 <= args.min_free_percent < 100.0:
        raise SystemExit("--min-free-percent must lie in [0,100)")
    if args.poll_seconds <= 0.0:
        raise SystemExit("--poll-seconds must be positive")
    if args.report_file is not None and args.report_file.exists():
        raise SystemExit("--report-file already exists")
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        raise SystemExit("provide a command after --")

    environment = os.environ.copy()
    environment.update(
        {
            ACTIVE_ENV: "1",
            PATH_ENV: str(Path(__file__).resolve()),
            MAX_RSS_ENV: repr(float(args.max_rss_gib)),
            MIN_FREE_ENV: repr(float(args.min_free_percent)),
        }
    )
    request_read, request_write = os.pipe()
    response_read, response_write = os.pipe()
    handshake_token = secrets.token_hex(32)
    environment.update(
        {
            REQUEST_FD_ENV: str(request_write),
            RESPONSE_FD_ENV: str(response_read),
            TOKEN_ENV: handshake_token,
        }
    )
    process: subprocess.Popen[bytes] | None = None
    try:
        process = subprocess.Popen(
            command,
            start_new_session=True,
            env=environment,
            pass_fds=(request_write, response_read),
        )
        os.close(request_write)
        request_write = -1
        os.close(response_read)
        response_read = -1
        return monitor_process(
            process,
            max_rss_gib=float(args.max_rss_gib),
            min_free_percent=float(args.min_free_percent),
            poll_seconds=float(args.poll_seconds),
            report_file=args.report_file,
            quiet=bool(args.quiet),
            request_fd=request_read,
            response_fd=response_write,
            handshake_token=handshake_token,
        )
    except BaseException:
        if process is not None and process.poll() is None:
            terminate_group(process, "wrapper failure")
        raise
    finally:
        for descriptor in (
            request_read,
            request_write,
            response_read,
            response_write,
        ):
            if descriptor >= 0:
                try:
                    os.close(descriptor)
                except OSError:
                    pass


if __name__ == "__main__":
    raise SystemExit(main())
