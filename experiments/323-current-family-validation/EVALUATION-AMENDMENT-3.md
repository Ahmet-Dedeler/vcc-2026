# Numerical parity repair — September 9, 2026

The historical emitted panel-total digest reproduced, as did PDS, MSE and
fidelity. NMAE, reach and Jaccard failed their frozen tolerances. No new arm
was scored.

Diagnosis on the known-answer target MLLT1 found a normalization-order tie
difference. Scanpy computes `counts / (library / 1e6)`, whereas the production
emitter's fast DE helper computes `counts * (1e6 / library)`. Floating-point
differences change some ties. FGFR1's adjusted p-value was 0.04999425835
through Scanpy versus 0.05000068924 through the first evaluator: one reference
call crosses 0.05, changing three aggregate metrics.

On this full 18,400-control/400-perturbation-cell target, feeding the fast
rank-sum kernel Scanpy-normalized values matches the official p-values to
3.33e-16 and adjusted p-values to 8.88e-16, with zero changed calls.
Evidence: `work/de-parity-diagnosis.json` and `work/official-de-MLLT1.parquet`.

Add `exact_de323.py` for evaluation only. Keep **two separate control blocks**:
the unchanged production one for emitting predictions, and the
Scanpy-arithmetic-order one for evaluating them. Compute arithmetic LFC
directly from normalized counts, as cell-eval2 does. A synthetic fixture
compares all p, adjusted p and LFC values against actual `compute_de(scanpy)`.

Recompute reference statistics and rerun the same known-answer gate. Preserve
the failed metrics JSON and log. No model, emission, target cohort, seed,
tolerance or scientific decision threshold changes.
