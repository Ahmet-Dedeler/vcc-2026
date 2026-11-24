#!/usr/bin/env python3
"""Load or unpack the lossless 313-snapshot VCC 2026 validation leaderboard archive.

Usage:
    python scripts/reconstruct_leaderboard.py --summary
    python scripts/reconstruct_leaderboard.py --unpack out/leaderboard-json/
    python scripts/reconstruct_leaderboard.py --team "Ahmet Dedeler"
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

DEFAULT_DIR = Path(__file__).resolve().parents[1] / "data" / "leaderboard"


def _snap_sort_key(snap: dict) -> str:
    digits = "".join(re.findall(r"\d+", snap["file"]))
    return digits or snap["file"]


def load_archive(lb_dir: Path = DEFAULT_DIR) -> tuple[dict[str, dict], list[dict]]:
    entries: dict[str, dict] = {}
    with gzip.open(lb_dir / "entries.jsonl.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            obj = json.loads(line)
            h = obj.pop("_h")
            entries[h] = obj

    snapshots: list[dict] = []
    with gzip.open(lb_dir / "snapshots.jsonl.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            snapshots.append(json.loads(line))
    snapshots.sort(key=_snap_sort_key)
    return entries, snapshots


def iter_snapshots(lb_dir: Path = DEFAULT_DIR):
    entries, snapshots = load_archive(lb_dir)
    for snap in snapshots:
        yield snap["file"], [entries[h] for h in snap["entry_hashes"]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--summary", action="store_true", help="Print archive statistics")
    ap.add_argument("--unpack", type=Path, help="Unpack all 313 JSON files into this directory")
    ap.add_argument("--team", type=str, help="Print score trajectory for a team name (substring match)")
    args = ap.parse_args()

    entries, snapshots = load_archive(args.dir)

    if args.summary or (not args.unpack and not args.team):
        teams = {e.get("team_name") for e in entries.values() if e.get("team_name")}
        print(f"Snapshots:       {len(snapshots)}")
        print(f"Unique versions: {len(entries)}")
        print(f"Distinct teams:  {len(teams)}")
        print(f"Files:           {snapshots[0]['file']} .. {snapshots[-1]['file']}")

    if args.team:
        q = args.team.lower()
        seen = set()
        for snap in snapshots:
            for h in snap["entry_hashes"]:
                e = entries[h]
                name = str(e.get("team_name") or "")
                if q in name.lower():
                    key = (e.get("id"), e.get("score_avg"))
                    if key not in seen:
                        seen.add(key)
                        print(
                            f"{snap['file']:28s}  id={e.get('id')}  "
                            f"avg={e.get('score_avg'):.6f}  "
                            f"pds={e.get('score_pds')}  mse={e.get('score_mse')}  "
                            f"model={e.get('model_name')!r}"
                        )

    if args.unpack:
        args.unpack.mkdir(parents=True, exist_ok=True)
        for snap in snapshots:
            out = args.unpack / snap["file"]
            payload = {
                "entries": [entries[h] for h in snap["entry_hashes"]],
                "num_entries": snap["num_entries"],
            }
            out.write_text(json.dumps(payload), encoding="utf-8")
        print(f"Unpacked {len(snapshots)} snapshots to {args.unpack}")


if __name__ == "__main__":
    main()
