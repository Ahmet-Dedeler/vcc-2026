#!/usr/bin/env python3
"""Read-only status for the three recovered research workers.

A printed PID, an early RESULTS.md, or process disappearance is not completion.
Completion requires a successful terminal trace event AND the final artifact.
This script never launches, retries, stops, submits, or reads expression data.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = (
    "303-second-overexpression-library",
    "304-f15-audit",
    "305-t7-expression-test",
)


def read_trace(path: Path) -> dict[str, Any]:
    """Ignore an incomplete final JSON line while a worker is writing."""
    summary: dict[str, Any] = {
        "events": 0,
        "tool_calls": 0,
        "terminal": None,
        "latest_message": "",
    }
    if not path.is_file():
        return summary
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            try:
                event = json.loads(line)
            except (ValueError, UnicodeError):
                continue
            if not isinstance(event, dict):
                continue
            summary["events"] += 1
            if event.get("type") == "result":
                summary["terminal"] = {
                    "subtype": event.get("subtype"),
                    "is_error": event.get("is_error"),
                }
            if event.get("type") == "assistant":
                for item in event.get("message", {}).get("content", []):
                    if not isinstance(item, dict):
                        continue
                    if item.get("type") == "tool_use":
                        summary["tool_calls"] += 1
                    elif item.get("type") == "text":
                        summary["latest_message"] = str(item.get("text", ""))[-900:]
    return summary


def classify(
    *,
    launched: bool,
    alive: bool,
    events: int,
    tools: int,
    terminal: dict[str, Any] | None,
    result_exists: bool,
    preaccess_only: bool,
) -> str:
    if terminal is not None:
        if terminal.get("is_error") is not False or terminal.get("subtype") != "success":
            return "TERMINAL_ERROR"
        if not result_exists:
            return "TERMINAL_MISSING_ARTIFACT"
        if preaccess_only:
            return "PREACCESS_REPORT_READY"
        return "RESULT_READY_FOR_AUDIT"
    if alive:
        return "RUNNING" if tools else "STARTING"
    if not launched and not events:
        return "NOT_LAUNCHED"
    return "INTERRUPTED_NO_TERMINAL" if events else "LAUNCH_FAILED_NO_EVENTS"


def inspect(directory: Path) -> dict[str, Any]:
    pidfile = directory / "claude-run.pid"
    pid: int | None = None
    alive = False
    try:
        pid = int(pidfile.read_text().strip())
        if pid > 0:
            os.kill(pid, 0)
            alive = True
    except (FileNotFoundError, ValueError, ProcessLookupError):
        pass
    except PermissionError:
        # The process exists, but no signal or modification is attempted.
        alive = True

    trace = read_trace(directory / "claude-run.jsonl")
    report = directory / "RESULTS.md"
    report_text = report.read_text(encoding="utf-8") if report.is_file() else ""
    preaccess = (
        "PREACCESS ONLY" in report_text.upper()
        or "PRE-ACCESS ONLY" in report_text.upper()
    )
    status = classify(
        launched=pidfile.exists(),
        alive=alive,
        events=trace["events"],
        tools=trace["tool_calls"],
        terminal=trace["terminal"],
        result_exists=report.is_file(),
        preaccess_only=preaccess,
    )
    return {
        "experiment": directory.name,
        "status": status,
        "pid": pid,
        "pid_exists": alive,
        "events": trace["events"],
        "tool_calls": trace["tool_calls"],
        "terminal": trace["terminal"],
        "report_exists": report.is_file(),
        "preaccess_report": (directory / "PREACCESS.md").is_file(),
        "latest_message": trace["latest_message"],
        "caution": "PID existence alone does not verify worker identity or progress.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    rows = [inspect(ROOT / "experiments" / name) for name in EXPERIMENTS]
    if args.compact:
        for row in rows:
            print(
                f"{row['experiment']}: {row['status']}; "
                f"tools={row['tool_calls']}; terminal={row['terminal']}; "
                f"report={row['report_exists']}"
            )
    else:
        print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
