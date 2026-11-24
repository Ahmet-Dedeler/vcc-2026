# 333 — original schema gate failed, no biological fit

September 9, 2026. Filename recovery 001 was frozen at **f78057ce**.
The guarded conversion used only the approved NumPy unpickler globals and
reached `validate_mapping`, where it rejected the declared requirement that
every gene key be a nonempty uppercase string. The original gate is
**NO-GO as written**; no `embeddings.npz` or terminal coverage result was
published. No unrestricted pickle fallback was used.

Exact error:
`ValueError: symbols must be unique nonempty uppercase strings`.

Outer wall 3.37 seconds; macOS maximum RSS 4,248,109,056 bytes. The guard's
two-second sample saw 2.77 GiB and missed the brief maximum. Exit 1.
All activity was local. Embedding objects were deserialized, but no
biological response, feature-accuracy result or model fit was accessed.

Next safe diagnostic is the source's key/type inventory: count non-string,
empty, mixed-case and whitespace-bearing keys and case-fold collisions,
with a few identifier examples. Do not rewrite identifiers, discard rows,
change feature-coverage gates or infer a biological failure from this error.
Any necessary identifier-policy correction must be disclosed and frozen
before feature conversion resumes.
