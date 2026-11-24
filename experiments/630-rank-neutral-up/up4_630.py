"""630 v4 capacity candidates (driver-only; the pinned module src_rank_up.py, 7e6b3e52, is NOT touched).

`up_pass4` is src_rank_up.up_pass2 with three independent flags, made by source transformation so it tracks the module
exactly. All flags off == up_pass2 bit for bit (KA in rnup630: variant `<v3>_f` vs `<v3>`). From [7772dc]'s review (c)
and the oracle (full delivery at K 1.0 is worth ~2x today's pass on H1):
  a  above_budget  a gene without a recipient takes an above-control-max carrier only if it still has sub-gate budget
                   (v3 picks the max-cw one even when drained, then moves nothing);
  z  free_room     in mode b, a donor's zero-cost in-interval room stays usable after the crossing budget is spent
                   (v3 skips every donor once budget <= 0);
  r  rerecip       a gene whose recipient has no sub-gate budget left takes a new one (an above-cmax carrier with budget,
                   else a fresh promotion) instead of stalling for the rest of the rounds.
Variant grammar (rnup630): <v3 variant>_f<flags>, flags a subset of "azr" in that order, e.g. b_k100_q20r3p10_fazr.
"""
from __future__ import annotations

import inspect

import src_rank_up as M


def _build():
    src = inspect.getsource(M.up_pass2)
    edits = [
        ("def up_pass2(", "def up_pass4("),
        ("quantum=None, rounds=1, rec_pool=1):",
         "quantum=None, rounds=1, rec_pool=1, above_budget=False, free_room=False, rerecip=False):"),
        ("        if j in rec_of:\n            a = rec_of[j]\n",
         "        if above_budget or rerecip:\n"
         "            above = above[sub_left[above] > 0]\n"
         "        if j in rec_of and (not rerecip or sub_left[rec_of[j]] > 0):\n            a = rec_of[j]\n"),
        ("                if budget <= 0 or xb[bi] <= 0 or not elig[bi]:\n",
         "                if (budget <= 0 and not free_room) or xb[bi] <= 0 or not elig[bi]:\n"),
        ("                kb = _max_drop_k(int(xb[bi]), float(s[bal[bi]]), c, budget)\n",
         "                kb = _max_drop_k(int(xb[bi]), float(s[bal[bi]]), c, max(budget, 0))\n"),
    ]
    for old, new in edits:
        if src.count(old) != 1:
            raise AssertionError(f"up4: anchor not found exactly once: {old[:60]!r}")
        src = src.replace(old, new)
    ns = dict(M.__dict__)
    exec(compile(src, "<up_pass4>", "exec"), ns)
    return ns["up_pass4"]


up_pass4 = _build()
FLAGS = {"a": "above_budget", "z": "free_room", "r": "rerecip"}


def flag_kwargs(flags: str) -> dict:
    if any(f not in FLAGS for f in flags) or list(flags) != sorted(set(flags), key="azr".index):
        raise ValueError(f"v4 flags must be a subset of 'azr' in that order: {flags!r}")
    return {FLAGS[f]: True for f in flags}
