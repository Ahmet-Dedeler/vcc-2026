# 190 — The regulatory-prior ordering key: take a competitor's board-measured +0.0141 without its `nmae` bill

You are a research subagent. Work ONLY in `experiments/190-regulon-ordering/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory (network downloads
of public network resources go to `data/external/regulon/` — that is the ONE exception).
Peak RSS ≤ 10 GB (four other agents are running). `uv run python`.
**Post nothing anywhere and put none of our method into any query or public field** (`AGENTS.md` §14).

## READ FIRST
`AGENTS.md` §4, §11.1 rules 1–13 with the audit-181 amendments. Then:
- `experiments/186-intel-0903/RESULTS.md` §M2 — the board evidence this mandate is built on.
- `experiments/108-signed-regulon/RESULTS.md` — **our own DoRothEA attempt, NO-GO.** Know exactly
  what it did; the whole point of this experiment is that it used the prior for a different purpose.
- `experiments/172-sign-confidence-key/RESULTS.md`, `171-reach-prefix/RESULTS.md`,
  `154-confidence-key-screen/RESULTS.md` — the ordering-key record.
- `docs/BOARD-LEDGER.md` — our `reach` trajectory and the shrunk-cohort control.

## The evidence

`tydymy` (Complex Disorders Alliance) published their config verbatim on the public leaderboard:

> "reference store from X-Atlas/Orion HCT116 and HEK293T plus Replogle K562 genome-wide,
> decomposition transfer, template_weight 0, conserved_alpha 1.30, shrinkage 1.0, denoiser off,
> seed 0, **plus a CollecTRI signed TF to target correction at 25 percent of backbone L2**."

The single arm that introduced it (`xfer-orion-genefill-k0-v4` → `regprop-collectri-v1`, 2026-08-26)
moved, **on the board**:

| member | before | after | Δ raw | Δ overall |
|---|---:|---:|---:|---:|
| `reach` | 0.2196 | 0.2549 | **+0.0353** | **+0.0066** |
| `fid` | 0.4789 | 0.4922 | **+0.0133** | **+0.0075** |
| `jac` | | | +0.0030 | +0.0014 |
| `pds` | | | −0.0013 | −0.0005 |
| `nmae` | 0.9792 | 1.0213 | **−0.0421** | **−0.0116** |
| `log mse` | | | +0.108 | 0 (clamped both ends) |
| | | | | **net +0.0021** |

**The prior bought the ORDER and the CALLS and paid for it in the MAGNITUDES.** That is the exact
signature of an ordering lever, from a team with no connection to us.

Independent same-week confirmation: **D²R² (arXiv 2608.15288, 2026-08-15)** reformulates perturbation
prediction as regulation-guided gene-wise *ordered* generation, initialises the ordering policy from
**a GRN inferred from control cells**, and reports controlled ablations — generator and budget fixed
— in which biological-prior ordering beats random ordering, beats uncertainty-based heuristics, and
**degrades every metric when reversed**. Evaluated on Norman19 and **VCC-H1**. Two properties make it
ours to use: the prior comes from **control cells only**, which is all we have for A/B/C and all we
will have for D/E/F; and its own ablation says the *ordering* is where it pays.

Weaker third corroboration: `raft.bio` ships a per-(target, gene) permutation-null weight ranked
against 18.3k genome-wide knockdowns and sits at `reach` 0.2300.

## The hypothesis, and why our own 108 does not close it

**108 blended the regulon into the DIRECTION** — regulon-only PDS 0.585/0.584 against 0.685/0.647,
sign-preserving reweight −0.00035 macro, every promotion gate failed. That is the same use that costs
`tydymy` their `nmae`, and it is not this hypothesis.

**This experiment uses the prior ONLY as a confidence ORDER over genes.** `reach` reads nothing but
the order of the reference's DE pool under our confidence and the sign we assign; `fid` reads which
genes we call. Neither reads the profile. So if the order can be set independently of the emitted
profile — which is `experiments/185-rank-controlled-emission`'s question, and its INBOX 3 shows our
half of the DE table is exactly computable on the real contexts — then **the `reach` and `fid` gain
is collectable with `pds`, `mse` and `nmae` provably unchanged**, and `tydymy`'s −0.0116 `nmae` bill
is not incurred at all. Their +0.0141 gross becomes the whole of it.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- `|H|` key h1-test `reach` **0.23533**; production **0.18283**; H-oracle **0.62119**; oracle
  **0.64553** (172/171). Reproduce with `experiments/153-de-bench/debench.py`'s `Surface.reach`.
- Donors on h1-test: **`k562_gwps` + `orion_hek293t` only**. `data/donors/h1_2025` covers 85/85 of
  h1-test's targets and is pure oracle leakage.
- Pools scaled **proportionally** per perturbation; never `Surface(n_real_target=N)` — it deletes the
  small-`n_real` stratum that carries 98% of our reach.

### Phase 1 — build the prior, two independent versions
1. **CollecTRI** signed TF→target regulons (`saezlab/CollecTRI`, via OmniPath or `decoupleR`):
   45,856 signed interactions, 1,183 TFs, `source/target/mor`. **Check the licence before use** — it
   aggregates 12 upstream resources and OmniPath separates academic from commercial. Record what you
   find; if it is not clean for a prize competition, say so and rely on (2).
2. **A control-derived co-expression GRN**, D²R²-shaped and licence-clean, built from the challenge's
   own 18,400 NTC cells per context (`data/challenge/context_{A,B,C}.h5ad`) and/or from h1-test's
   controls. This one is available for D/E/F on day one, which the external network also is, but this
   one is *context-specific*, which is strictly more information.

From either, form a **per-(target, gene) regulon-agreement score** `R(t, g)`: does the sign our model
predicts for `g` under knockdown of `t` agree with what the network says a loss of `t` should do to
`g`, and with what confidence. Define it precisely and write the definition down before scoring.

### Phase 2 — price it as a KEY, which is the cheapest possible read
`Surface.reach(key = f(|H|, R))` on h1-test against the `|H|` incumbent **0.23533**. No emitter, no
build. Sweep the mixing weight from 0 (identity, the known-answer row) upward. Report also top-1% and
top-5% sign accuracy over the pool — 172's calibration is that **reach reads only the front, so
precision at the top is the objective, not rank correlation**.

Then `fid` as a function of the call set chosen by the same score, remembering
`F_p = k / max(n_pred, n_real)` and that **`n_real` is a distribution** (audit 181: board median
90–110, heavy tail). Report against assumed `n_real`, not one value.

### Phase 3 — the falsifiers, all mandatory
- **Shuffled-regulon control** — permute the network's targets within each TF, and separately permute
  the signs. This is the control `tydymy` did not report and it is the one that decides the result.
- **Coverage honesty:** how many of the 300 panel targets are TFs with a regulon at all? CollecTRI
  has 1,183 TFs; if only ~20 of our 300 targets are covered, the arm's reach is capped at that
  fraction and you must say so in overall units, not as a per-covered-target effect.
- **Does it beat `|H|` × expression?** 189 measured that control expression alone predicts per-gene
  cross-line sign reliability with a **+0.0505** held-out decile spread. A regulon prior must beat
  that cheap baseline, not just beat random.
- **Retrodict `tydymy`.** Their arm is a board point of this exact family. Your model of it should
  predict `reach` up and `nmae` down; if you cannot reproduce the *shape*, say so (audit-181 rule 1).

### Phase 4 — the direction blend, only if Phase 2 is positive
Separately price a small-weight additive regulon correction to the direction on **PLOLO** `pds`
(162's instrument, the only one validated against a board `pds` point) and on h1-test `nmae`. This is
108's use and it is expected to lose; measure it so the two uses are cleanly separated in the record.

## Deliverable
`RESULTS.md`: known-answer rows, the prior's construction and coverage, the key sweep against
`|H|`, the `fid` curve, all four falsifiers, and a verdict in overall units. If the honest answer is
"the regulon adds nothing over `|H|` × expression", that is a full result and it closes M2.

## Territory
184 owns Tahoe. 185 owns the emission law and how an order is delivered — **you supply a key, it
supplies the mechanism; read its RESULTS.md if it lands first.** 188 owns per-gene magnitude.
189 owns the per-gene prior fitted from donor cross-line agreement — **you own the per-(target, gene)
REGULATORY prior.** Do not duplicate 189's expression baseline; use it as your comparator.
