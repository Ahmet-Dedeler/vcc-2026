# bulkid658 smoke test on a real package (2026-10-04 07:20Z, after C10b's driver exited)

- **Pair:** 497dr's built h5ad (`submissions/497dr-october-path-lodopct-spt-abc.h5ad`, 11.9 GB) vs its own package, restored
  with `scripts/restore_prediction_h5ad.py` (10 s; streamed zstd).
- **Result: BULK IDENTICAL and exact content IDENTICAL** (rc 0; `bulkid_smoke.json`).
  - The comparison read 900 (context, target) blocks and 1,488,991,072 stored entries; the meta's per-context sum is
    1,488,991,072. Total counts are 14,330,628,083.
  - So `vcc prep` keeps row order and the X arrays bit for bit: `--exact` is reliable between a built h5ad and a restored
    package, not only between two restored packages.
- **Cost:** 28 s at 20,000 rows per read, peak 5.4 GB. The read block is now 4,000 rows: 21.6 s for two passes, peak
  **1.16 GB** (`/usr/bin/time -l`).
- **Disk:** a restore needs 12 GB free while it exists; both drivers delete it right after the check.
