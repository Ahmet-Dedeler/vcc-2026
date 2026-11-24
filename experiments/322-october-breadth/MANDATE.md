# 322 — October robustness: does a BROADER source stack degrade better on a disjoint panel, even when it loses today? Measure it for specific sources, on 233's own instrument.

You are a research subagent. Work in `experiments/322-october-breadth/`. **Edit nothing outside it.**
**NEVER run `vcc submit`; never build into `submissions/`; do not edit `scripts/make_submission_pm_cd4.py` or
`scripts/def_pipeline.py`** (they carry three shipped arms and a package is staged for tomorrow's slot).
`uv run python`. **One `cell-eval2` instance at a time, ever.** Cloud is allowed via `scripts/cloud_score.sh`
**capped at $12** — but read §"instrument" first, because 233's own harness may answer this for free.
**Time box: 5 hours.** A well-supported NO-GO closes the question.

## Why this is the highest-value open question
**The board we are ranked on is not the board that pays.** Only the last entry after **Oct 22** is prize-eligible,
and Oct 22 releases a **new panel and new cell lines** (D/E/F). **233 measured what a disjoint panel costs us:**
93 % of aggregate cosine retained but only **64 % of `pds`** and **60 % of ordering**, projecting 0.2118 →
**0.137 / 0.119** — *and the loss is carried by the TARGETS, not by coverage*. That is an explicit **lower bound**
(a random draw's transferability ceiling is 1.3–2.0× below Arc's, and Arc's protocol applies again).

**233 §7 item 3 says: "If you want insurance, buy breadth, not a new mechanism. Coverage explains 78 % of the
cosine loss."** It then names the mechanism — *"the component whose contribution collapses is the one whose
coverage collapses (CD4, 24.7 % → 10.7 %), and the one that is coverage-robust by construction (the H1 union,
sized on carrying rows) does not move"* — **but it never measured whether ADDING a broad source would help
October.** That is this mandate.

**The candidate that makes this concrete:** `kolf2` covers **273 of 300** panel targets (278's inventory) — far
broader than `dld1_lowmoi` (67) or `gse295214_h1` (64). **264 already measured it as a base donor on today's
cancer folds and it LOSES** (equal mean −0.0195; whitened-amplitude-matched −0.0107, 3/3 folds). **The question is
whether a source that loses today degrades better tomorrow** — which is exactly the trade 233's coverage argument
predicts and nobody has priced.

## The question, precisely
**Under 233's disjoint-panel protocol, does adding a broad-coverage source to the stack improve RETENTION, and is
the retention gain larger than the level loss it costs today?** Report both halves. A source that costs −0.0195
today and retains 5 % better is still a loss; say so.

## Instrument
**Read `experiments/233-disjoint-panel/RESULTS.md` first and reuse its protocol and its draws** — it built random
disjoint 300-target panels and re-scored; its per-component attribution table (CD4 24.72 % → 10.71 %, H1 union
flat) is the object you are extending. **Prefer its free local harness over cloud.** If a six-member reading is
genuinely needed, 317 §1 documents the cloud path ($1.56/score, workers self-terminate) — but justify it, and
respect the $12 cap. Carry 233's own known-answer row before trusting any number you produce.

## Arms and controls
- **Baseline:** the shipped stack, on today's panel and on 233's disjoint draws.
- **+`kolf2`** as a fourth base donor (264's construction — reuse it, do not re-derive), same two surfaces.
- **At least one narrow source as the contrast** (`dld1_lowmoi`, 67 targets) — 233's argument predicts a
  broad source retains better *because* of coverage, so a narrow one must not show the same gain.
- **Controls:** a coverage-matched *random* source assignment (does the gain track coverage, or the identity of
  the source?); and the same measurement with coverage held fixed, to separate "more targets covered" from
  "different biology".
- Report **retention** (disjoint / today) per component, not only levels — that is 233's unit and it is what
  October pays.

## Rules
Programme rules 1–9; name the unit (the **panel draw** is the resampled unit here, and 233's draw count bounds
every interval). **This is a measurement, not a build** — nothing competition-facing is selected, no arm is
proposed for a slot. **Do not quote a number measured on a synthetic fixture as a delivered value** (Correction 5
in `docs/LOG.md`), and **name the surface and stage of every magnitude you report.**

## Deliverable
`RESULTS.md`: 233's KA row reproduced; the retention table (component × surface × arm); whether breadth buys
October retention and at what cost today; the controls; and one line — **322 (broad-source October insurance):
WORTH IT / NOT WORTH IT / UNRESOLVABLE (retention Δ …, today's cost …)**. Report any spend. Then stop.
