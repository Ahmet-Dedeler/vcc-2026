# 455 — the October driver's SHIPPED q follows the board: 1.25 → 1.75, row 1 re-proved on 446b (PASS, 2026-09-26 03:02Z)

- **Patch** (`patch455.py`): `scripts/def_pipeline.py` SHIPPED_ENV `VCC_BULK_SHARPEN_Q` 1.25 → **1.75**. The
  ARM_BOARD_DELTA text now lists the three levelled q steps and their board gains (417a +0.0069612, 433k +0.0046475,
  446b +0.0022834). Driver `5463fe28` → **`62b5110d`**, as precomputed on a scratch copy before the patch.
- **Row 1 (the day's path builds what we ship).** The driver's defaults, with only the two day-of values exported
  (446b's stack dump `bb19f0e5…` and its re-level A 1.1004 / B 1.1517 / C 1.1720), rebuild **446b's artifact byte for
  byte**: `9559923a774eb248…`. The `[arms]` guard read 34/34 checks passing. The builder was restored to `4feae27d`.
- `run_arm_generic.sh` accepts the driver before or after this patch. An arm exports its own q, so SHIPPED q never
  changes an arm build.
- `def_day.sh` passes q explicitly (arg 4), so the October day never depended on this default. This keeps the
  driver's defaults equal to the shipped arm (432's T4), so a manual driver run cannot silently fall back to q 1.25.
