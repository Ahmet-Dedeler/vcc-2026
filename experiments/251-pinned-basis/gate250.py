"""250 — the memory gate, checked before every fold load (mandate; AGENTS.md §11.1b).

Proceed when `memory_pressure` >= 40% free AND `vm_stat` free+inactive >= 20 GiB AND swap
free >= 400 MiB.  233's `whiten233.memory_ok`, verbatim in substance.

    uv run python gate250.py            # print the reading, exit 0/1
    uv run python gate250.py --wait     # poll every 60 s for up to 4 h, exit 0 when clear
"""
from __future__ import annotations

import subprocess
import sys
import time


def memory_ok(min_free_pct: float = 40.0, min_free_gib: float = 20.0,
              min_swap_free_mb: float = 400.0) -> tuple[bool, str]:
    mp = subprocess.run(["memory_pressure"], capture_output=True, text=True).stdout
    free_pct = 100.0
    for ln in mp.splitlines():
        if "free percentage" in ln:
            free_pct = float(ln.rsplit(":", 1)[1].strip().rstrip("%"))
    sw = subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout
    free_sw = float(sw.split("free =")[1].split("M")[0].strip())
    vm = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    page = 16384

    def pages(tag):
        for ln in vm.splitlines():
            if ln.startswith(tag):
                return int(ln.split(":")[1].strip().rstrip("."))
        return 0
    free_gib = (pages("Pages free") + pages("Pages inactive")
                + pages("Pages speculative") + pages("Pages purgeable")) * page / 2**30
    msg = (f"memory_pressure free {free_pct:.0f}%, free+inactive {free_gib:.1f} GiB, "
           f"swap_free {free_sw:.0f} MiB")
    return (free_pct >= min_free_pct and free_gib >= min_free_gib
            and free_sw >= min_swap_free_mb), msg


def wait_for_memory(max_minutes: int = 240, tag: str = "") -> str:
    for i in range(max_minutes):
        ok, msg = memory_ok()
        if ok:
            return msg
        print(f"[gate{tag}] HOLD ({msg})", flush=True)
        time.sleep(60)
    raise SystemExit(f"memory gate never cleared in {max_minutes} min")


if __name__ == "__main__":
    if "--wait" in sys.argv:
        print(f"[gate] CLEAR ({wait_for_memory()})")
    else:
        ok, msg = memory_ok()
        print(f"[gate] {'CLEAR' if ok else 'HOLD'} ({msg})")
        sys.exit(0 if ok else 1)
