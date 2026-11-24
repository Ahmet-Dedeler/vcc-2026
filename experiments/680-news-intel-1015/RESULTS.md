# 680 — News and intel sweep, 2026-10-14 22:10Z → 10-15 ~00:40Z

**Lane:** news + intel in the five-session council of 10-14/15 (coordinator [400b6e], readiness [17adab],
science [bf003b], this session [020377]). **Deliverable:** `docs/NEWS-2026-10-15.md` (the full write-up); this
file records what was done, the method, and the one-line findings.

## What changed, one line each
- Scorer unchanged: cell-eval2 0.18.0 (10-01) moves no scored number; board `anchor_version` r4 since 09-02.
  Watch: a private Arc branch `fix/expr-mse-signal-fraction-floor` would change the mse rule; not shipped.
- Rules revised 2026-09-16 (not in the repo until now): "solely by ML models" dropped; blended classical
  components must be described, not hidden. No "not able to win" clause exists for State/Stack on the live FAQ.
- Arc released PIE (10-05): zero-shot per-gene `p_de` / `lfc_pred` / `delta_p_pred` for an unseen context from
  its control mean; public checkpoint trained on Tahoe, Jiang, Orion, VCC-2025 H1. Handed to [bf003b] (exp 700).
- Board: 19 teams ≥ 0.30 (7 on 10-08); 61 above our best; we are #66. The five "clones" are one operator's five
  accounts since 09-16 (lock-step scores, constant submission-count offsets, bulk identical to 12 digits); six
  more multi-account or shared-code families found. The honest movers gained mostly bulk (pds/mse). Illumina
  0.4729 gained mse/nmae/fid/reach at flat pds. No public recipe, notebook or post explains any of it.
- D/E/F release is 07:00Z Oct 22 by the site's own timeline JS (midnight Pacific); "Submissions Close" marker is
  Nov 5 08:00Z vs the rules' 23:59Z. `GET /api/leaderboard?get_final=true` is gated ("not available yet").
- Discord READ (01:00–01:30Z) through Ahmet's logged-in Aside browser (`aside repl`, read-only scraper
  `li[id^="chat-messages-"]`, raw JSON in `work/discord_*.json`): the 10-06 queue stall was a platform-wide "resource
  exhaustion" incident, resolved 10-07; the final-phase logistics announcement has not been posted yet; Arc has answered
  none of the open forum questions (NTC resampling, test genes, cell counts); multi-accounting is being enforced
  (consolidate by Oct 16). No recipe or notebook shared.
- Inbox READ (`gwsa <account> gmail +triage --query …`): the 09-16 rules revision was emailed on 09-17; Arc posts
  competitor mini-profiles (incl. Ahmet's) on 10-15; a 10-03 model-trade offer from another team sits unanswered.

## Method (so it can be re-run on Oct 21/22)
1. Board: `curl https://virtualcellchallenge.org/api/leaderboard` → `work/leaderboard-<ts>.json`; the analysis
   scripts are inline in this session's log; outputs `work/board-analysis-2214Z.txt`, `work/top25-movement.txt`
   (daily series from `data/leaderboard/2026-10-*.json`), `work/multiaccount-groups.txt` (raw-member ties).
2. Scorer: `gh api repos/ArcInstitute/cell-eval2/{releases,commits?since=…,contents/CHANGELOG.md}`; compare
   `anchor_version` across snapshots.
3. CLI: PyPI `vcc-cli` JSON; source under `~/.local/share/uv/tools/vcc-cli/lib/python3.14/site-packages/vcc/`.
4. Site: it is a Next.js app; curl returns "Loading...". Read /app/rules, /app/faq (click each question;
   `find` → `computer left_click ref`), /app/datasets in the built-in browser. The timeline constants live in
   `_next/static/chunks/app/(protected)/app/page-*.js` (saved in `work/site/chunks/`).
5. Repos: `gh search repos … --sort updated`, `gh search code "vcc2026"|"context_A.h5ad"|…`; HF via
   `hf datasets ls --search` and the `/api/datasets?author=arcinstitute` endpoint.
6. Papers: Arc news index via fetch; PIE preprint PDF downloaded to `work/papers/` and text-extracted with pypdf
   (run with `python -I` in a separate directory).
7. Discord precheck (added 10-16): `curl -s -m 10 -o /dev/null -w '%{http_code}' https://discord.com/` must print 200 before a
   read; `000` means the Aside client sits on its loading splash (HTTPS to discord.com stalls on this Mac behind WARP), and
   `discord_scrape.js` now throws "discord client not connected" rather than returning []. Reads 10-15 21:21Z → 10-16 21:50Z are void.
8. Board churn (10-17): a team's row is ABSENT from `/api/leaderboard` while its latest entry is unpublished (in flight,
   processing or failed) and returns with the new score; median 2 rows vanish per 3 h step, the leader did it 4 times before 10-17.
   Never read a vanished row as a removal without the absence-episode history (script in the 10-17 12:47Z heartbeat).

## Not done / caveats
- X unread (login wall in the built-in browser; the Aside browser has X logged in for next time).
- Daily routine through Nov 5 (owned by this lane; coordinator's ask): inbox query, Discord #announcements /
  #challenge-status / #changelog via Aside, Arc blog (the 10-15 profiles post), em2.bio, cell-eval2 changelog and board
  `anchor_version`; one dated line in docs/NEWS-2026-10-15.md per day.
- The Molina & Zhang (July 2026) summary is second-hand (a participant README); not read first-hand.
- Snapshots stop at 10-08 15:50Z, so the 10-08 → 10-14 movement has no daily resolution.
- Nothing was posted, starred, forked, joined or downloaded from untrusted sources; the only download is the
  bioRxiv PDF (2.5 MB) and public site JS (0.8 MB). Disk after: 200 GiB free.

## Peak memory
Negligible (< 1 GB; curl, python3 JSON parsing, gh).
