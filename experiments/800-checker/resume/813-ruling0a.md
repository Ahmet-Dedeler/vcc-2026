FROM THE LAB LEAD (10:00Z): read `experiments/813-public-data/RULING-0a.md`. It is the decision on your stage 0a.
In short:
- Run ONE corrected 0a as `freeze813-v2.json` (≤ $0.15, EC2 via labcloud.sh). The code diff is limited to uncapped strata and
  the matched-incidence twin, exactly as PROPOSAL-2 specifies. Write the diff and the new sha256 in RESULTS.md before the
  run. The bars are unchanged.
- KA-SG takes the G5 fallback: status "uninformative (G5)", no bounds, and the caveat written down. No new statistic.
- If v2 fails a bar, record it and STOP; there is no third 0a. If it passes, record freeze813-v2.json's sha256 in RESULTS.md
  and STOP. I commit the sha before 0b.
- Record the v1 capped-strata leak (0.200) as a negative finding about coarsened-k estimators.
- The local-job HOLD is still in force: EC2 only, no local compute.
