# 488 — Does our emitter's model of the scorer match cell-eval2 exactly? (Cursor Cloud, Grok 4.7 extra-high)

Your experiment directory is `experiments/488-cloud-scorer-conformance/`.

## Why

Our DE channel is designed around a precise model of Arc's scorer (`cell-eval2`, preset `vcc2026`): which genes pass the
control filter, how the Wilcoxon p-values and BH correction are computed, how ties and zeros are handled, how log2 fold
changes use arithmetic means of per-cell CPM with ε 1e-9, how `reach` orders genes (significant first, then adjusted p,
p, |lfc|, gene name), what `fid` counts, how `nmae` handles missing genes, and how `pds` / `mse` pseudobulk and jackknife
work. Several emitter components claim exactness against this model (for example an "exact-U" pass that moves
arithmetic means while keeping every Wilcoxon U statistic unchanged: `src/vcc2026/rank_mag.py`). A 2026-09-26 review
already found one mismatch: two different "CPM > 5" gates (pseudobulk CPM in `scripts/make_submission_pm_cd4.py`
`emission_rank_axis` vs mean per-cell CPM in `src/vcc2026/rc_emit.py`). Any assumption that differs from the real
scorer can silently cost score.

## Method

1. Read `docs/reference/vcc2026-metrics-BRIEF.md` and `vcc2026-metrics-FULL.md`, then the INSTALLED `cell-eval2` source
   in the uv environment (find it with `uv run python -c "import cell_eval2, os; print(os.path.dirname(cell_eval2.__file__))"`
   or similar). Record the version.
2. List every assumption about the scorer that the emitter code makes (grep `src/vcc2026/` and the builder for Wilcoxon,
   BH, CPM, epsilon, ties, ordering, gate, target exclusion, jackknife, ρ). For each: the code location, the claim, and
   the scorer's actual behaviour (source location).
3. Test each with small synthetic AnnData inputs run through the real `cell-eval2` functions (not re-implementations),
   and through the emitter's own helpers, and compare numbers. Emphasise: the control filter definition, tie handling in
   the rank-sum and its normal approximation (continuity correction? tie correction?), BH scope, the lfc formula,
   `reach`'s ordering and tie-breaks, `fid`'s treatment of calls without direction, `nmae`'s gate and missing genes, and
   whether the "exact-U" moves in `rank_mag.py` really leave the scorer's p-values unchanged.
4. Keep synthetic inputs small (the full scorer needs ~25 GB on real data; your VM may be smaller).

## Output

`RESULTS.md`: a table of assumptions with CONFORMS / DIFFERS / UNTESTED, and for every DIFFERS a reproduction, the
likely direction and size of the score effect, and a proposed fix (patch file in your directory).

## Ground rules for every cloud teammate

- You run on a Cursor Cloud VM with this repo and its `uv` environment (`uv run python ...`; `cell-eval2` is installed).
  The gitignored data (`data/challenge`, `data/donors`, `data/harness`, `data/external`, `submissions/`) is NOT on your
  VM; do not assume it exists. Production runs on the lead's macOS machine (bash 3.2, BSD grep, APFS); you are on Linux,
  so reason explicitly about any shell or tool behaviour that differs (AGENTS.md §6.2 lists the traps).
- Read `AGENTS.md` and `docs/FINDINGS.md` §0–§2 first. The record is large: grep `docs/` and `experiments/*/RESULTS.md`
  before claiming something is new.
- Write ONLY inside your experiment directory. Commit only that directory to your branch. Do not open a PR, do not merge,
  do not edit production files (`scripts/`, `src/`, `docs/`, other experiments); put proposed patches as separate files
  in your directory.
- Nothing public: no posting, issues, PRs, accounts, emails or messages. Never run `vcc`.
- Evidence standard: every claim carries a reproduction (a test, a script, a command and its output) that the lead can
  re-run. Label each finding CONFIRMED (reproduced) or SUSPECTED (reasoned). A well-supported negative is a result.
- Commit `RESULTS.md` early and update it as you go; finish with a complete `RESULTS.md` whose first lines give the
  answer. Budget: up to ~8 hours.
