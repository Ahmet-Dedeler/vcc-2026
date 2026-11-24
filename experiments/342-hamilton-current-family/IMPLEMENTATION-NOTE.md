# 342 implementation note (Claude, September 9, 2026; nothing executed)

**Coordinator pre-freeze review:** the initial eleven synthetic tests passed.
Codex found that the arm-level integrity check incorrectly required floor to
restore rounded mass, defeating the intended unchanged baseline. That check
now uses the appropriate mass for each rule; an additional pure regression
test covers it. The time budget is also checked before every target. Neither
change modifies rounding, emission, scientific gates or any measured result.
The historical worker note below describes its unexecuted initial delivery.

Files in this directory, all new, nothing outside it touched:

| file | role |
|---|---|
| `PREREGISTRATION.md` | prospective registration; freeze hashes it |
| `rounding342.py` | pure rules (floor / Hamilton / random_up), residual count, per-block `BlockRounder` adapter |
| `seam342.py` | AST seam on production `emit_counts`: one statement in the quasi branch; attestation |
| `run342.py` | `freeze`, `verify`, `run`, `attest` |
| `report342.py` | paired 10,000-target bootstrap, six members, gates, `work/results.json` |
| `audit342.py` | optional independent replay: sufficient-statistic rescoring, dense MSE bootstrap, 0/1 totals contract |
| `test342.py` | pure tests; `verify` runs them; also `uv run python test342.py` or pytest |

Suggested sequence for Codex (each step is a separate decision):

```
uv run python experiments/342-hamilton-current-family/test342.py
uv run python experiments/342-hamilton-current-family/run342.py attest
uv run python experiments/342-hamilton-current-family/run342.py freeze
git add experiments/342-hamilton-current-family && git commit -m "Freeze Hamilton current-family remeasurement 342"
uv run python experiments/342-hamilton-current-family/run342.py verify
uv run python scripts/run_with_memory_guard.py --max-rss-gib 16 \
    --report-file experiments/342-hamilton-current-family/work/guard.json \
    uv run python experiments/342-hamilton-current-family/run342.py run
uv run python experiments/342-hamilton-current-family/audit342.py
```

Expected cost, from 323's cap15 record (142 s per 85-block arm, 5.04 GiB peak for
six arms in one process): three arms plus reference load, verify hashing and the
report in roughly 10 to 15 minutes, well inside the 2,400 s in-process budget.

What `run` does, in order: verify (hashes, committed and clean own files, seam
attestation, scorer files, pure tests, 323's own verify); load 323's reference
and cap15 effects (hash-checked); floor arm; exact known-answer against 323's
cap15 outputs (stops on any inequality, receipt preserved); hamilton arm;
random_up arm; report; `work/execution.json`.

Known limits stated plainly:

- The emission stream is shared across the 85 targets within an arm, exactly as
  in production. Hamilton's totals change the scatter draw shapes, so from the
  second target on, the emission draws differ from floor's even where totals
  agree. The random_up control shares this property.
- The floor known answer requires exact equality with 323's recorded cap15
  outputs, including per-target rows. If any check fails, the seam is not
  accepted; do not relax it.
- Nothing here was run. The end-to-end `rc_hybrid_block` fixture in the tests is
  synthetic and small; if its per-cell weld fails for fixture reasons, that is a
  test-fixture defect to report, not evidence about the seam.
