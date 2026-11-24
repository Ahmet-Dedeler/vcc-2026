# PLAN.md — 353-recipient-conditioning-scout

## Hypothesis
The held-out-line sign deficit (343: 0.60 |lfc|-weighted sign accuracy at design
stage; up-regulated reference genes a coin flip at 0.516) is addressable by
conditioning the *transferred response* on the recipient line's own basal state —
information our pipeline currently uses only for means, a CPM gate and a noise
proxy. Published zero-shot models (STATE, TxPert, CPA/scLAMBDA, X-Cell, scDEBART)
all consume the target cells' unperturbed state; the literature may contain a
recipient-conditioning mechanism that (a) is fittable on our existing banks,
(b) predicts per-gene *direction* on a held-out line, and (c) is not already
closed by 325–342.

## Exact inputs (read-only; nothing written outside experiments/353-*/)
- Local: `docs/RESEARCH-PROGRAMME-2026-09.md`, `docs/SPRINT-HANDOFF-2026-09-09.md`
  (ledger), `experiments/343-*/RESULTS.md`, `experiments/344-*/RESULTS.md`,
  `experiments/246-*/RESULTS.md` §0–§1, `experiments/315-*/RESULTS.md` §0,
  `experiments/314-*/RESULTS.md` §4–§5, `experiments/341-*/core341.py`,
  `experiments/334-*/selection.json`,
  `experiments/323-current-family-validation/work/{metrics-cap15.json,per-target-cap15.csv,stats-cap15.npz,reference_de.parquet}`.
- Web (cite URL + access date; ≤15 words quoted per source): STATE (Arc 2025),
  TxPert (Valence 2025), GEARS, CPA, scLAMBDA, scGPT/scFoundation benchmarks,
  Tahoe-100M / Tx1 / TX1-CD, X-Cell/Pisces 2026, scDEBART 2026, AROMA 2026,
  TwinCell 2026, LGR-Tahoe 2026, IRIS 2026, PerturbNet, GPerturb, PS-score
  (Nat Cell Biol 2025), basal-expression-predicts-response literature.

## Estimator / method
No fitted estimator. Literature + repo synthesis into a ranked proposal memo
(RESULTS.md) with ≤3 experiments. Each specifies: mechanism paragraph,
concrete estimator + exact bank paths, primary endpoint (|lfc|-weighted sign
accuracy on held-out line's reference-significant genes, retrieval as guard),
cheap falsifier + kill number, and why 325–341 do not already answer it.

## Controls / falsifiers
- Known-answer reproduction (own code, tiny): recompute 343's delivered
  |lfc|-weighted sign accuracy (0.588) and mean NMAE (0.9508) from 323's saved
  `stats-cap15.npz` + `reference_de.parquet`. If mismatch >1e-6, stop and report.
- For each proposed experiment: a named cheap falsifier (<1 h) with a numeric
  kill threshold, plus wrong-recipient / shuffled-label controls per programme
  rule 4.
- CLOSED-list check: every candidate mechanism grepped against the protocol's
  CLOSED list and the 325–342 ledger before inclusion.

## What kills the idea (proposal-level)
- If the literature shows recipient conditioning only helps with ≥10 diverse
  training lines (TxPert's own ablation: no basal encoder best at 4 lines) and
  our 6-line frame cannot support it → rank such proposals last or reject.
- If a candidate mechanism reduces to a CLOSED construction (text features
  334/336, TabICL 332, marginal maps 339, neighborhoods 340, own-line labels
  341, basal gating beyond 10 CPM per 246-H3, sign-consistent consensus per
  246-H4) → excluded with citation.
- If known-answer reproduction fails → no proposal numbers are trusted; report
  the failure instead.

## Budget
≤3 h total. No run beyond the tiny check (<5 min, <2 GiB RSS).
