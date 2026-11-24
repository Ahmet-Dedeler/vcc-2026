# 453dr BOARD RESULT — 0.2663613 (−0.0036816 vs 446b): the October path is VALIDATED end to end

Entry `q3VuTImh4FRGKXH2waCp`, submitted 2026-09-28 12:11Z by `456/arm0928.sh` slot 1 ("the lab notebook ran a marathon",
description empty), scored and captured 12:37Z, rank 7. Package `submissions/453dr-october-path-abc-q1.75.vcc` sha256
`a76283bd…`. Counterpart 446b `UDgjOKL9gI1udAexq6ue` 0.2700428 (flags kept `noprog`, so 446b per the preregistration).
Bins: `PREREGISTRATION.md` (written 09-26 00:25Z). Reader: `456-zs-riders/read0928.py dress`.

The slot fired 12 h late: the 09-27 ~20:23Z reboot killed its waiter before the 09-28 reset; a new coordinator session
re-armed it with the same arguments (NEXT-SLOT-DECISION, 09-28 12:20Z). Nothing about the package or the rule changed.

453dr = `scripts/def_day.sh` (unmodified, b517b6d2) on 180's relabelled bundle (A/B/C as D/E/F), relabelled back.

| member | Δ raw vs 446b | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | −0.001345 | −0.002714 | −0.000452 |
| mse | +0.015033 | −0.015805 | −0.002634 |
| nmae | −0.000188 | +0.000355 | +0.000059 |
| fid | −0.000036 | −0.000116 | −0.000019 |
| reach | −0.003436 | −0.003797 | −0.000633 |
| jac | −0.000004 | −0.000011 | −0.000002 |
| **overall** | | | **−0.0036816** |

Scaled members: pds 0.81201, mse 0.34939, nmae 0.17167, fid 0.00546, reach 0.23127, jac 0.02836.

## Against the bins
- **Bulk channel −0.00309: VALIDATED** (bin −0.0045 … −0.0005). Predicted −0.0024 (CD4 ×1 instead of ×1.5 at −0.0017 …
  −0.0022, the NTC-proxy re-level ≈ −0.0002, re-fitted operator/key/CD4 banks ≈ 0). The realised cost sits inside the band,
  ~0.0007 below the point prediction; it is on mse (+0.015 raw), as the CD4 half predicted (pds −0.0013 raw, mse +0.003 …
  +0.010 raw expected; mse landed just past the top).
- **DE channel −0.00059: pinned** (|d| ≤ 0.009). The DEF carrier path equals the tuned one within emission noise
  (reach −0.0034 raw is inside the 0.0045 ordering sd).
- **Hard failures: none** (no ρ above 1, overall far above −0.02).

## Consequences
- **The October pipeline reproduces the tuned A/B/C arm at its priced cost.** Every day-of component (identity, CD4
  extraction, operator refit, key refit, NTC-proxy re-level, secant loop, relabel) is now board-validated in one entry.
- For `rule1004.py` / `rule1005.py`: **O = OK** (not BAD: bulk −0.0031 ≥ −0.0060, |DE| 0.0006 ≤ 0.009, overall ≥ −0.02).
- The −0.0031 bulk cost is the CD4 ×1 October default (438: ×1.5 loses on October-like panels). 457dr (09-29 slot 1, the
  CD4 ×1.5 twin of 453dr) reads that choice directly on the board tomorrow; `cd4_stack` stays ×1 unless 457dr reads ≥ +0.0040.
- The October path's A/B/C board number (0.2664) is the reference for every later "dr" arm (460dr, 461dr, 469dr, 497dr).
