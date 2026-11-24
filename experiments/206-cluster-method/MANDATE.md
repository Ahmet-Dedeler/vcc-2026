# 206 — Four accounts share one prediction with 1.66× our directional accuracy, and it is SPREADING. Find it.

You are an intel subagent. Work ONLY in `experiments/206-cluster-method/`.
Read-only on the repo. No `vcc submit`, no builds, no `cell-eval2`.
**Post nothing anywhere, contact nobody, create no accounts, join no Discord, and put none of our
method into any query, form, issue or public field** (`AGENTS.md` §14). Read-only research only.
Peak RSS ≤ 8 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §8, §8.2, §14. Then `experiments/186-intel-0903/RESULTS.md` (the previous sweep — **do
not repeat it**; it already identified this triple as *"almost certainly one shared pipeline"* and
identified `lingshu` as Lingshu-Cell, arXiv 2603.25240, Alibaba DAMO), and
`experiments/197-cosine-inversion/RESULTS.md` (the inversion and its identifiability bound).

## The observation

Live board, 2026-09-03 16:10 UTC, 640 teams. **Four accounts now report a six-tuple that is
identical to three or four significant figures on every member:**

| team | model | overall | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` | subs |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Hongjian Zhou | `pm-r14e` | 0.21064 | **0.8331** | 0.7533 | 0.9289 | 0.5087 | 0.2524 | 0.0223 | 18 |
| jw_13288 | `juniper_46` | 0.20949 | **0.8333** | 0.7515 | 0.9287 | 0.5085 | 0.2441 | 0.0223 | 15 |
| Andrew Liu | `AL-1466` | 0.20793 | **0.8330** | 0.7599 | 0.9296 | 0.5084 | 0.2469 | 0.0222 | 12 |
| MT Zeng | `Kestrel2` | 0.20744 | **0.8330** | 0.7591 | 0.9308 | 0.5087 | 0.2442 | 0.0224 | 12 |
| **us** | | 0.18342 | **0.8360** | **0.8738** | **0.9442** | 0.5037 | 0.2299 | 0.0302 | 26 |

`pds` agreeing to **0.8330–0.8333** across four independent accounts is the same prediction.
Hongjian Zhou's organisation reads **"Oxford / Prim…"**. **And it is spreading:** `MT Zeng` gained
**+0.046** today to join it and `jw_13288` **+0.0285**, both almost entirely on `mse`.

## The anomaly, stated precisely — this is what you are explaining

**Their `pds` is LOWER than ours (0.8330 against 0.8360) while their `mse` is far better (0.7515
against 0.8738) and their `nmae` is also better (0.9287 against 0.9442).**

197 resolves the apparent contradiction and it is the reason this matters: inverting 096's `mse`
decomposition, their raw `mse` implies a transfer cosine **≥ 0.1312 against our 0.079249 —
1.66× ours** — while `pds` is a *rank* statistic on which we extract **26% more retrieval per unit
of cosine** than the median shrunk team, because of whitening. **So they have a materially better
direction and we are compensating on `pds` with geometry, which cannot compensate on `mse` or
`nmae`.**

Per `AGENTS.md` §8.2, **direction is the entire remaining gap** (`mse` + `pds` + `nmae` account for
110% of it; we already lead on `fid` and `jac`), it is worth **+0.00365 overall per +0.01 of
cosine**, and every internal route to it has closed: 166 (no fifth panel-covering line), 184
(adaptation, on 50 lines), 193 (the operator), 194 (eleven estimators), 195 (learned transforms, at
the cell-line boundary), 203 (hierarchical shrinkage). **A shared, spreading method that holds
1.66× our cosine is the single most valuable thing on this board, and unlike our own experiments it
is external information we can look up.**

## What to find, in priority order

1. **Whether the cluster's method is public.** Model-name series are your best lever: `pm-r14a…e`
   (Hongjian Zhou), `AL-*` (Andrew Liu — note `AL-2574` on 09-02 became `AL-1466` on 09-03, so the
   number is not a version), `juniper_*`, `Kestrel*`, and cqawesome's `m0902a/b`, `m0903c`. Mine
   every `model_name` and `description` across all 640 entries and diff against 186's Aug-28 and
   Sep-03 snapshots (`experiments/186-intel-0903/board/`).
2. **The joining diff.** `MT Zeng` and `jw_13288` moved onto the cluster's six-tuple **today**.
   Their own submission trajectories before and after are the closest thing to a controlled
   experiment on this method that exists. Reconstruct both trajectories member by member and say
   **which member moved and by how much at the moment they joined** — that identifies the component,
   even if the code is never found.
3. **Attribution.** `gh search repos/code`, HF models and datasets, arXiv and bioRxiv since
   2026-08-20, Papers-with-Code, and any public write-up. `Oxford / Prim…` is a lead; resolve the
   organisation string. **Do not contact anyone.**
4. **Lingshu-Cell's actual role.** 186 measured that three teams swapped architectures onto it and
   every member moved < 0.006, all landing at `reach` ≈ 0.208 — so Lingshu-Cell is **not** the
   cluster's mechanism, and `huhuhu` (`(lingshu+fofj21)`, 0.19185) and `Leo` (`(CPA+i0bhu2)`,
   0.19120) sit well below it. Confirm or refute that separation; do not re-litigate 186's finding.
5. **The `cqawesome` question, still open from 186.** They lead at 0.22082 with `pds` **0.8703** —
   0.037 above the cluster — and remain untraceable. Their `m0903c` landed today (+0.0029).
   Report anything new, and if there is nothing, say so plainly.

## The discriminating question you should try hardest to answer

**Is the cluster's cosine advantage DATA or METHOD?** The tweet in this project's brief argues the
challenge *"becomes luck whether or not you have in-house data that's closer or farther from the
test cell lines"*. If the cluster's edge is a proprietary or newly-released atlas, we may be able to
get the same atlas (experiment 205 is chasing X-Atlas/Pisces and X-Cell for exactly this reason —
read its RESULTS if it lands first). If it is a method on public data, it is reproducible.

**Evidence that would settle it:** a named dataset in any description or repo; a cosine that is too
high for four public lines under 159's `1/k` law (197 puts the cluster at **~10.9 effective donor
lines** against our 3.09, and only four panel-covering lines exist publicly — **that arithmetic is
itself the strongest hint that the edge is data, and it is the thing to test**); or a public method
whose reported numbers match.

## Deliverable
`RESULTS.md` with: the refreshed cluster table and the two joining diffs member by member; whatever
attribution you can establish, with URLs; a ranked list of **concrete, testable hypotheses** for the
cosine advantage, each with how we would test it on PLOLO or h1-test and what it would be worth in
overall units at 197's +0.00365 per +0.01 cosine; and an explicit list of what you could not
determine. **"The method is not public and the arithmetic says it is data" is a full result** — it
would tell us to stop spending on internal direction work and concentrate on data acquisition.
