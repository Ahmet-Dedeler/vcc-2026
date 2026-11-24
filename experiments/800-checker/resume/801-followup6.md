FROM THE LAB LEAD (2026-10-16 12:01Z by the clock, queued while follow-up 5 ran): follow-up 6, small. Budget ≤ $1, EC2 only.
Deadline 10-18 06:00Z (science's deadline is 10-18 18:00Z).

**Ask (science, 713 part 2, registered at experiments/713-labc-dose-partners/REGISTRATION-H1.json @ ee791df7).** Rebuild the
three LABC dumps with dump801's extra-design path, `--designs h2.5c1,h3c1` (H1 ×2.5 and ×3 at CD4 ×1, beside fin = h2c1).
Keep the same spec, the same operators (ABCp's day fits) and the same CD4 root as follow-up 4.
- **KA (abort unless it holds):** every non-design array of the new dump equals the current LABC dump byte-for-byte or to
  |d| < 1e-12, and fin equals fin_h2c1. Current shas:
  - dump708-LABC-k562_gwps.npz `eba095c1…`
  - dump708-LABC-orion_hct116.npz `e97f1fbd…`
  - dump708-LABC-orion_hek293t.npz `34ecd80e…`
- **Write under a NEW tag** so the originals stay untouched: `experiments/801-lab-folds/work/dumps/dump708-LABCd-<held>.npz`
  + `.json` (dump708 format, `designs: ["h2.5c1","h3c1"]`). Mirror them to S3 under the same dumps prefix.
- Record the three new sha256 in RESULTS.md §16 and stop. The lab lead sends them to science, which adds them to its
  registration before reading anything. Do not price the design arms yourself.
