# 723 — Post-deadline disk cleanup (prepared 2026-11-06)

The final deadline passed with X2 as the last final entry, so the bulk arrays no longer serve a decision. All records were
committed first (def79420, a11002cd). `cleanup723.py` deletes only git-ignored files, from the reviewed list
`work/delete-list.txt.gz` (91,286 entries, ~831 GiB by unique inode). Ahmet runs the delete himself.

| deleted | GiB |
|---|---:|
| experiments/ ignored arrays and vendored packages | 486.7 |
| data/external/ public datasets | 327.4 |
| submissions/ old .vcc packages | 28.3 |
| data/def-download, data/challenge-def.unzipped (duplicates) | 1.2 |

Kept: everything git tracks; `submissions/x2def.vcc`, `r0def.vcc`; 692's build inputs (dump, whitener, scales; not its h5ads,
which `scripts/restore_prediction_h5ad.py` restores from the .vcc); `811/work/banks/`; small text under 5 MB; data/donors, harness,
challenge, challenge-def, leaderboard, board-entries; .venv.

**What this gives up:** rebuilding an arm or re-running an experiment now needs the public data re-downloaded and the intermediates
recomputed. Scoring X2 / R0 against D/E/F truth, if Arc releases it, needs only what is kept. This overrides AGENTS §10's competition-time
"keep data/external and 199's GWCD4i.DE_stats.h5ad" rule. Copies on S3 were not audited and are untouched.

Tested on a throwaway repo: refuses with no list, a dirty tree, unpushed commits, a tampered package and a changed list; the real delete
removed exactly the 9 listed entries (including a nested git repo), left all kept files, and git status stayed clean.

```bash
uv run --no-project python experiments/723-post-final-cleanup/cleanup723.py --delete
```

## Result (Ahmet ran `--delete`; verified 2026-11-06 12:10Z)
- 91,237 of 91,286 entries removed, 758.8 GiB freed (below the ~831 GiB estimate: APFS clones). Free space 829 GiB, data volume 55 %. Repo now 47 GB.
- 49 files (0.06 GiB) stayed: they sit in read-only preregistration folders (mode 555, experiments 107, 117, 118, 122, 125, 128, 135). Left in place.
- Verified after the delete: x2def.vcc e39d2ba6 and r0def.vcc 5c82ea9e match the ledger; X2 is still published and FINAL; no tracked file deleted; git clean and
  level with origin; data/donors, harness, challenge, challenge-def, leaderboard, board-entries, .venv, 811 banks and 692's build inputs all present.
