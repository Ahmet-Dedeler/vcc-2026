# 334 recovery — fill empty primary descriptions, not vectors or labels

September 9, 2026, before any response-model fit or accuracy result.
The original NCBI-only coverage failure is preserved in
`COVERAGE-FAILURE.md`, with all original code and hashes unchanged.

Only one change: for an unavailable/empty usable NCBI description,
look up the SAME exact primary key in
`NCBI_UniProt_summary_of_genes.json`, in the already-verified archive,
SHA `319897d7e2f8a0306bccdf2d022982b3c8321dbec8273739094021d57e414c31`.
Require the same exact `Gene Symbol KEY` prefix and token boundary and
remove that header. Do not change a nonempty valid NCBI body, infer an
alias, use a pretrained vector, generate text, select genes again,
lower a gate, change a model or inspect an outcome.

All original target/coordinate identities, models, controls, random
seeds, transforms, splits, time/memory limits and advancement gates
remain fixed. Remaining missing rows are retained with zero features.
The feature-source receipt counts primary, fallback and missing bodies.
This is a different text artifact from the original, even though it
uses the same released archive and dimensions.

Run in `recovery-work/`, never overwrite original outputs.
`run_recovery334.py` explicitly replaces only the feature-preparation
function and verification wrapper in the frozen runner. Its input
manifest binds the original manifest, this policy, the failure record,
recovery implementation and tests. The original runner, model and
analysis remain hash-verified. The wrapper adds a recovery attestation
with terminal result and execution hashes after the original analysis.
