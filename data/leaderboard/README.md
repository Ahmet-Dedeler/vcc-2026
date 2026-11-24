# VCC 2026 Validation Leaderboard Snapshots

Lossless deduplicated archive of **313 public validation-phase leaderboard snapshots** captured from `virtualcellchallenge.org` between **2026-08-21** and **2026-10-22** (24,189 unique entry versions across 2,077 teams, compressed from 415 MB of raw JSON to 7.3 MB).

## Files

- `entries.jsonl.gz` (6.9 MB): one JSON object per unique entry version, keyed by `_h` (first 16 hex chars of SHA-1 of the canonical entry JSON).
- `snapshots.jsonl.gz` (0.37 MB): one JSON object per snapshot (`file`, `num_entries`, `entry_hashes`).

## Usage

```bash
# Print summary stats
python scripts/reconstruct_leaderboard.py --summary

# Inspect a team's validation-phase score history
python scripts/reconstruct_leaderboard.py --team "Ahmet Dedeler"

# Unpack all 313 original JSON files
python scripts/reconstruct_leaderboard.py --unpack /tmp/vcc2026-leaderboard/
```

These snapshots contain only public leaderboard records published on the competition website during the validation phase (before the leaderboard went dark on Oct 22, 2026).
