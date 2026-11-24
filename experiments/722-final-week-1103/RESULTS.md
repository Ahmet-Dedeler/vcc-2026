# 722 — Final-week state check (2026-11-03, ~01:35–01:45Z)

Coordinator session restarted by Ahmet after an 11-day gap (last commit 10-22 23:49Z; the Mac rebooted ~10-27,
uptime 6 d). Nothing ran in between, and nothing needed to: the counting final entry was already in.

## 1. State
- **Counting final entry: X2 = R0 + 811 C0, `S5RxOrm4ZbqtnPRIX8NX`.** `vcc status` (vcc 0.3.0) at 01:36Z: published, mode
  FINAL, partition test, panel vcc2026-test-1, anchors testD/E/F-r1. Package `submissions/x2def.vcc` sha256 e39d2ba6… ==
  SUBMITTED_LEDGER. Members hidden until the Nov 19/20 ceremony.
- Arc's 11-02 email: final deadline **Nov 5 23:59 UTC**; 2 final submissions a day, only the last counts; an upload that
  starts before the deadline and finishes processing after it still counts.
- No live waiter, chain, worker, cron or launch agent that could submit (ps, crontab, LaunchAgents checked). AC power.
- AWS: no lab or production instance running (the one running instance is `peak-gateway`, t4g.small since 09-04, not
  VCC). The 10-02 snapshots (5, kept on purpose for the credits) are unchanged.
- Public A/B/C board, frozen since 10-22 (`work/leaderboard-20261103.json`, 1,635 rows): our displayed row is W1 0.2746
  (rank 90); our best entry was C0 0.2781. The top is Illumina AI 0.4720, then 0.3527 / 0.3447 / 0.3441. Versus #2–#12 the
  gap is mse (~0.50 vs our 0.41 scaled), nmae (~0.25–0.30 vs 0.185) and fid (0.05–0.09 vs 0.005).

## 2. Independent semantic check of X2's prediction (`sanity_x2.py`, sha256 048966d5…)
Input: `experiments/692-def-r0/work/x2/x2def_final.h5ad`, sha256 **94ae91b6…** == the ledger's h5ad_sha256 for
e39d2ba6, i.e. exactly the submitted content. Controls: `data/challenge-def`. Peak RSS 2.4 GB. Output
`work/sanity_x2.json`. This checks what `def_postbuild` does not (postbuild checks config and provenance).

| check | D | E | F |
|---|---:|---:|---:|
| pseudobulk log-CPM r with its OWN controls | 0.99989 | 0.99995 | 0.99992 |
| r with the other contexts' controls | 0.736 / 0.747 | 0.736 / 0.916 | 0.747 / 0.916 |
| median library, prediction / controls | 19,953 / 20,236 | 19,210 / 19,681 | 19,553 / 19,845 |
| own-target CPM ratio vs controls, median | 0.774 | 0.812 | 0.806 |
| targets with own-gene ratio < 0.5 / > 1 | 25.0 % / 0.7 % | 22.7 % / 0.3 % | 24.3 % / 1.3 % |

- 360,000 rows; each context has exactly the 300 `pert_counts.csv` targets × 400 cells; var order == gene_names.csv.
- Counts: 0 non-integer, 0 negative entries; max cell total 996,556 (< 1,000,000); max entry 674,810 (ballast cells).
- So contexts are not swapped, target labels carry their own knockdown (not shuffled), and depths match.

## 3. Decision: X2 stands; no further final submission
- No candidate exists. Every lab lane after 811 is closed or WEAK (805 FAIL, 810 no dossier, 813 no donor candidate,
  811's other configurations WEAK on ABCp), and nothing was queued after 10-22.
- The dose lever on 811 does not qualify: λ 0.50 (all-lines form) read HARM on 062-ABC (−0.0074, 0/3 folds). Picking an
  interior λ from A/B/C reads is the A/B/C-specific tuning §11.1c.3 forbids.
- With the board dark there is no instrument to read an X3 with except the lab's folds, and the sealed reserve LZ is
  spent. A replacement would carry build and gate risk on the only entry that counts. On the A/B/C board a few
  thousandths moves rank by a few places and is nowhere near the 0.07 gap to the prize places.
- This follows the 10-22 rule (NEXT-SLOT-DECISION): no further final submission unless a NEW decision, fixed in advance and
  not based on any D/E/F score, justifies one. None does.

## 4. After the deadline
- Revert the sleep override Ahmet set on 10-15: `sudo pmset -a disablesleep 0` (OCT22-CHECKLIST §G).
- Nov 19/20: capture R0's and X2's members if they are revealed (`scripts/capture_board_entry.py`, team endpoint with
  `get_final=true`). That gives the only D/E/F read of the recipe and of 811's increment, the calibration record for the
  write-up.
- Storage (65 GiB free, 97 %): after Nov 5, the X1 / X2 build trees can be cleaned with `git clean -xdf` per AGENTS §10.
  The .vcc packages restore their h5ads.

One honest sentence: nothing changed the model; this check confirms the submitted final is the intended, well-formed
prediction and records why no replacement is made.
