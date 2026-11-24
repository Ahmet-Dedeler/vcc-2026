# Experiment 116 source-only conversion/QC protocol v1

This protocol is frozen before either LFC matrix payload is parsed. It creates
one canonical-axis source artifact while keeping Low1 and Low2 distinct.

The archival HDF5 preserves every finite author LFC value that maps to an exact
canonical response-gene symbol, in float64, for all 6,451 original columns.
NTC and unsupported target columns remain present there. Missingness is stored
as a separate exact Boolean mask. The accepted gzip files preserve every
original byte, including source response rows outside the VCC canonical axis.

Target labels use one rule fixed from metadata alone: `NTC_*` is a control;
otherwise the exact substring before the first underscore (or the complete
label if there is no underscore) must equal a canonical gene symbol. No alias,
case, transcript, punctuation, or fuzzy repair is allowed.

All promoter columns remain individually preserved. For a canonical target
with multiple promoters, the derived target-level value is the coordinatewise
median of finite promoter values, separately in Low1 and Low2. Observed
effects never select or delete a promoter.

The explicit transfer arrays are response-gene by target. Every one of the 300
VCC panel response rows is forced to exact zero and `measured=false`. The
unmasked archival HDF5 is labeled non-transferable.

Reliability thresholds create masks and diagnostics only. They do not alter
the preserved source values or remove rows. The predeclared recommendation
gate in `CONVERSION_PROTOCOL.v1.json` decides only whether a new, separately
preregistered held-context experiment is justified.

After acquisition:

```bash
uv run python experiments/116-dld1-lowmoi/convert_compact_v1.py build \
  --confirm I_CONFIRM_EXPERIMENT_116_SOURCE_ONLY_CONVERSION_V1
```

Verification is source-only:

```bash
uv run python experiments/116-dld1-lowmoi/convert_compact_v1.py verify
```
