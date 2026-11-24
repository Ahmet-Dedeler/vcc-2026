FROM THE LAB LEAD (2026-10-17 10:18Z by the clock): one small preparation, local file reads plus seconds-long arithmetic only (CLEAR
window; nothing beyond a few seconds or 100 MB locally, otherwise EC2). No new idea, no new prediction.

Write `experiments/814-panel-rule/predict814.py`: the frozen reader for your two 10-22 predictions, as PROPOSAL-2 states them.
- **P-A:** the count of D/E/F targets in the full essential-screen library (`work/flags/esslib.txt`, `24fffbfe…`), with its bins.
- **P-B:** D/E/F's 300 targets over prop10's frozen pool deciles, with the decile edges pinned by number and the missing-feature
  rule as VERDICT-2 fixed it, and its support / fail bins.
- Input: a pert_counts.csv path. Output: a json and the bin verdicts.
- **KA, required:** on `data/challenge/pert_counts.csv` (A/B/C) it must give P-A = 0 and P-B = 164 / 300 in deciles 9–10 (your
  RESULTS numbers). Record the script's sha256 in RESULTS.md. On 10-22 the lead runs it once on D/E/F's file, and nothing about
  it changes after this commit.
Then stop.
