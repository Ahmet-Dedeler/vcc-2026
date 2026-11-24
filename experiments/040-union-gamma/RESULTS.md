# 040 — Clean HCT116 union-gamma result

Date: 2026-08-23

## Verdict

The production-shaped union-axis `gamma=0.65` baseline completed successfully
on the full 132-target Orion HCT116 count harness. This is the clean replacement
for the quarantined H1 comparison described in `README.md`.

The result is now the fixed comparator for experiments 053, 061, and 065.
It does not itself authorize a challenge submission.

## Six-metric result

Official-anchor scaling uses the midpoint anchors already frozen by experiments
031/040.

| metric | raw | scaled |
|---|---:|---:|
| PDS cosine | `0.7978834143` | `0.6583058879` |
| MSE normalized | `2.0537529303` | `0.0000000000` |
| direction fidelity | `0.5281190914` | `0.0487303047` |
| direction reach | `0.0150248423` | `-0.0635883457` |
| DE Jaccard | `0.0337623585` | `0.0128712392` |
| LFC NMAE | `0.9570910100` | `0.0735223515` |
| **six-metric average** | — | **`0.1216402396`** |

The negative reach score confirms that the HCT116 harness remains a difficult
and imperfect proxy for the challenge contexts. Use this result for paired
candidate deltas, not as an estimate of the live leaderboard score.

## Execution and preservation

- Worker: `i-08d90ae88c7b210ed`, on-demand `r7i.2xlarge`.
- Scorer: `cell-eval2==0.16.0`, `vcc2026`, Scanpy backend.
- Wall time: `3:27:46`.
- Maximum resident set: `17,835,020 KiB` (about `17.0 GiB`).
- Exit code: `0`.
- Durable S3 prefix:
  `s3://<your-s3-bucket>/output-clean-hct/g0.65_q1.0_a1.08_spread0.25/`.

Locally preserved compact artifacts:

| file | SHA-256 |
|---|---|
| `agg_results.csv` | `f1011cd12786d2196949daa6ce3bf328d52c0afc9d46966e8d3e882f6200cccf` |
| `exit-code.txt` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `score.log` | `cb77f984edd25ca31127439a508839da2bc81503a54b3c1df2258bd5dece510c` |

The first experiment-053 watcher invocation failed before transformation
because its systemd environment did not include `uv` on `PATH`. The watcher
was patched and relaunched; the full transformation then completed under
`0.5 GiB` RSS and its scorer started. Experiment 061 also started scoring
after a temporary narrow S3 read permission was added for this baseline prefix
and immediately removed once the aggregate was copied.
