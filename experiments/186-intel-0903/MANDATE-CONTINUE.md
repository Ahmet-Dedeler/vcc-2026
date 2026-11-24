# 186 (continuation) — finish the intel sweep from artifacts already on disk

A previous run of this mandate gathered the raw material and then died mid-run on a model-side
classifier before writing anything up. **The data is already collected.** Your job is the analysis
and the write-up. Do not re-download what is already here.

**Read `experiments/186-intel-0903/MANDATE.md` first — it is your full mandate and the questions in
it are still the questions.** Then note what the previous run already saved in your directory:

```
board/leaderboard_2026-09-03.json     the live 612-team board
board/descriptions_all_2026-09-03.md  every model_name + description on the board
board/trajectories.md                 per-team submission trajectories
web/pertmind_v1.pdf, pertmind_v1.txt  the PertMind paper (arXiv 2608.16419)
web/arxiv_recent.xml                  recent arXiv listing
web/gh-repos-vcc.txt, gh-repos-more.txt   GitHub repo search results
web/state-commits.txt, cell-eval2-commits.txt   upstream commit logs
parse_board.py, clusters.py, diff_desc.py, traj.py   its scripts
```

Work ONLY in `experiments/186-intel-0903/`. Read-only elsewhere. No `vcc`, no builds, no
`cell-eval2`. **Post nothing anywhere, contact nobody, create no accounts, and put none of our
method into any query or public field** (`AGENTS.md` §14). Additional network research is allowed
and encouraged where the saved artifacts leave a question open. `uv run python`.

Deliver `RESULTS.md` exactly as the original mandate specifies, with the ranked list of concrete
testable mechanisms for **"what gives raw `reach` ≥ 0.25 at raw `pds` ≥ 0.84?"** as the centrepiece,
plus an explicit list of what you could not determine.
