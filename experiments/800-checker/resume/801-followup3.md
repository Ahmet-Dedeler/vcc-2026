FROM THE LAB LEAD (2026-10-16 09:28Z): follow-up 3, small. Deadline: done by 10-17 20:00Z, before 710's board read
(10-18 s1). Budget ≤ $2, EC2 only.

**Why.** Row 710 (H1 ×2.5 vs ×2 at CD4 ×1.5) on your ABCp folds failed its KA, as it should have. Run
`20261016T092059Z-ka804-710p` used `experiments/804-known-answers/price710p.py`. 062's A/B/C design `fin_h2c1.5` differs from
the ABCp fold's `fin` by 1.20 / 3.75 / 4.03 (k562 / hct116 / hek293t) because the compositions live in the operator's whitened
space. Logs are under `experiments/800-checker/cloud/runs/20261016T092059Z-ka804-710p/out/jobs/`.

**Task.** Produce `experiments/801-lab-folds/work/designs/ABCp__<held>.npz` (+ .json) for the three cancer held lines. Use the
designs801.py ABC path (`h2c1.5` and `h2.5c1.5`, A/B/C rules), with the held-out operator swapped for the day operator the same
way dump518p.py does for ABCp. Same keys and targets as the ABC__ files.
- KA: `h2c1.5` == the ABCp fold's `fin` (`engine801.load("ABCp/<held>").fin`), |d| < 1e-12.
- KA: with 062's operator, the same code reproduces the ABC__ file's `h2c1.5` row sha256.
- Record npz sha256, row sha256 and the operator sha256 in the .json.
- Mirror the files to the S3 input path the other designs use.

Do not run price710 / price710p. The lab lead runs `price710p.py` once the files exist; it prefers `ABCp__<held>.npz`
automatically. Your current window/hold rules apply: `window800.sh` may print "WINDOW hold (local jobs only)" until R7a ends.
That allows EC2 but no local job. Report in RESULTS.md §13 with the runs and $, then stop.
