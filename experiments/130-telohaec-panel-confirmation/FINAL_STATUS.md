# Experiment 130 final status

Status: **CLOSED — NO ELIGIBLE CANDIDATE**

Experiment 130 remains a metadata-only, fail-closed preparation. Its original
pre-value implementation freeze is preserved byte-for-byte.

Experiment 125 has completed with an immutable **NO-GO** result:

- macro signed-credit delta: `-0.015049899010796045`;
- macro precision delta: `+0.012174572820953597`;
- signed-credit delta was negative in every family:
  - Replogle K562: `-0.01214610823191394`;
  - Orion: `-0.026082415627716488`;
  - CD4: `-0.006921173172757709`;
- precision failed both complete-procedure null gates:
  - observed: `+0.012174572820953597`;
  - target-label null p95: `+0.015527271023240793`;
  - response-label null p95: `+0.013516444470243838`.

Experiment 130 required an immutable Experiment 125 `GO` with every promotion
gate passing. That prerequisite is false. Therefore:

- no candidate satisfies Experiment 130's unlock conditions;
- `CANDIDATE_BINDING_AMENDMENT.json` must remain absent;
- no TeloHAEC expression value may be opened under Experiment 130;
- no TeloHAEC confirmation, null, bootstrap, or transport evaluation is
  authorized;
- Experiment 130 is terminally closed unless a materially new hypothesis is
  assigned a new experiment number.

No TeloHAEC `X` or layer value, compact TeloHAEC effect/control value, HCASMC,
DLD1, Pan, KOLF, H1/Orion held outcome, challenge truth, or challenge-context
expression value was accessed while closing this experiment. No `cell-eval2`,
count generation, challenge build/submission, commit, or push occurred.

The machine-readable evidence and hashes are recorded in
`TERMINAL_CLOSURE.v1.json`.
