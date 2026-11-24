# 328 — Numerical continuation passes the source gates

September 9, 2026. **EARNS_H1_COUNT_STAGE_REGISTRATION**.
This does not clear a package, production change or submission.
V1's `SOURCE_PILOT_NO_GO` and its files remain unchanged.

The solver-only amendment was frozen at **ddd3c3cd** before continuation.
The joint fit converged at **588 total iterations**, with final mean
log-likelihood increment **9.9892485e-7**, below the unchanged 1e-6 threshold.
Its entire v1 trace is an exact prefix. Both already-converged controls,
including all their weights and metadata, remain exactly unchanged.

| comparison | target-macro log-density delta, nats/triple | paired target 95% CI | fraction favoring joint |
|---|---:|---|---:|
| joint minus diagonal | +.01337746 | [+.01034952, +.01685980] | .704 |
| joint minus shuffled joint | +.02232457 | [+.01422700, +.03556943] | .740 |

All four registered gates pass. This is repeated source-development evidence
after a disclosed post-result numerical amendment, not untouched confirmation.
Predictive density of noisy source coefficients does not establish latent
effect accuracy or recipient-context performance.

Verification:

- V1 training endpoints reproduced exactly; per-target held densities within
  4.45e-16. The recovery script checked these before continuation.
- A separate read-only audit verified unchanged controls, the trace prefix and
  convergence, then recomputed the 10,000 paired target intervals from the
  saved CSV; maximum discrepancy below 1.8e-17.
- Ten `unittest` tests pass, including resumed versus uninterrupted EM.
  `pytest` is absent; no dependency was installed.
- Recovery process exited zero: 33.16 s inside the runner, 39.19 s under the
  outer timer; sampled tree peak 2.30 GiB, macOS maximum RSS 2,546,352,128 bytes.
  No guard report file is claimed: command-mode guard output is terminal-only.
- Recovery result SHA-256:
  `dc4bb7eb502a1fc2ee4ae3edc30a2ce2029e1e165260fca810a16234768461b0`.
- Joint fitted-model SHA-256:
  `109a1502305b7624cc3715afce3e7f252239658b1e39e89c7c5323bc85ec2285`.

Next: separately freeze experiment 329 before constructing or scoring its
current-family H1 count predictions. No H1 outcomes were read by this recovery.
