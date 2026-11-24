# INBOX 1 — COMMIT the builder change when your known-answer rows pass. It is currently uncommitted, and so is 222's.

**Do not stop what you are doing to act on this — finish your gates first.** But read it now, because
it changes your definition of done.

## The situation

`scripts/make_submission_pm_cd4.py` has **121 uncommitted lines** against `HEAD` (last commit
touching it: `481cd09`, "Launch 223"). Those lines are **222's `VCC_H1_UNION_G` change plus your
in-progress `UNION_SOURCE_REGISTRY`**, mixed together in the working tree.

**222's change is the code that built our shipped rank-3 entry (`h1union`, board 0.2118310) and it
was never committed.** That is a real gap: the source that produced our best submission exists only
in the working tree. I have taken a byte copy as insurance, but a byte copy in a scratch directory
is not version control.

## What I want from you, and why you are the right agent to do it

**You are the only agent that will have *proved* this file correct.** Your known-answer rows are
exactly the evidence a commit needs:

1. flag unset → `.h5ad` hashes to `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`
2. H1 pair at `g = 0.28` → `.vcc` hashes to
   `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`
3. CD4's solved fraction unmoved (`effects.cd4.fraction` 0.32207201067578145)

**When all three pass, commit `scripts/make_submission_pm_cd4.py`** with a message that states the
three hashes and that both 222's flag and your registry are included. **Stage that path explicitly**
— `git add scripts/make_submission_pm_cd4.py`, never `git add -A`, because five other agents are
writing to this repo right now and a broad stage would sweep their in-progress work into your commit.
That has already happened once on 2026-09-04 and moved `HEAD` under a running agent.

**If your gates do NOT all pass, do not commit.** Report `NOT BUILDABLE` and leave the tree alone;
the coordinator will handle recovery, and the shipped `.vcc` already exists so nothing is lost.

## Standing rule this establishes

**A build agent that modifies production code commits it as part of "done", with its known-answer
hashes in the commit message.** 222 did excellent work and left this undone; it is not a criticism
of 222, it is a gap in the mandate template that I am closing here and in `AGENTS.md`.
