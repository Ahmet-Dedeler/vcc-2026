# Experiment 137 result — terminal NO-GO

Date: **2026-08-30**

## Decision

The high-MOI Jurkat tangent-cone candidate failed the frozen development
promotion gates. The untouched 79-target test remains sealed and must not be
opened. Do not retune the cap grid, split, masks, fallback, source scale, or
geometry against these development outcomes.

## Frozen development result

- final cap selected on all five development folds: `6.0°`;
- nested selected caps: `6.0°, 1.5°, 6.0°, 6.0°, 0.0°`;
- nested macro raw-PDS delta: `-0.000921658986175089`;
- held-fold raw-PDS deltas:
  `+0.0002560164, -0.0002560164, -0.0056323605, +0.0010240655, 0.0`;
- nested macro signed-credit delta: `+0.0004293694830849792`;
- held-fold signed-credit deltas:
  `+0.0010772782, +0.0006804831, +0.0005421732, -0.0001530870, 0.0`;
- 20,000-draw cluster-bootstrap 95%:
  `[-0.0043010753, +0.0059395801]`;
- 99-permutation maximum: `+0.0002560164`;
- promotion: `false`.

The exact-identity, exact-fallback, and maximum-angle invariants passed. The
positive gain was not stable enough to define a finite concentration fraction,
because the nested macro PDS delta was nonpositive.

## Runtime amendment

The first authorized development process opened only the 315 development rows
and then terminated before printing a score or writing a result artifact. The
failure was a bytewise identity assertion attempting a `uint8` view of a
non-contiguous array. The repair wrapped both operands in
`numpy.ascontiguousarray`; it changed no candidate value, score, selection,
gate, seed, split, or outcome access.

The amendment was bound in `RUNTIME_AMENDMENT.v1.json`, independently audited
twice, and resumed under the exact original development lock. The 79 test rows
were never accessed.

## Artifact bindings

- development summary SHA-256:
  `8e7a898367c2a9930a5c72ab361cb4ba86823853c7e14d6ed454c52ee23515d2`;
- development arrays SHA-256:
  `bbfb9dda83b190e14933ba4448d8a8e2f0aa72f75d0e53c601d537f840f011b3`;
- source candidate NPZ SHA-256:
  `417efe32b11c45b5c8a517c47a6a9c741de374fe2f46df2a96aa75ff248ec3b6`;
- shared operator NPZ SHA-256:
  `52374bcb27dbd0feafe2ff09a2135bb0a57b3f14e7ffc93252ed4d78baaf0cc4`;
- runtime amendment SHA-256:
  `2551a15548cabee6a5156494d851e140ee22603531a082f58eff95382d5615df`.

No candidate freeze, test lock, test result, count build, `cell-eval2` run,
challenge submission, commit, or push was created.
