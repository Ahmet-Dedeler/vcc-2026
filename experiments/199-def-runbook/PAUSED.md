# 199 PAUSED 2026-09-03 12:05 UTC — memory contention, not a scientific decision

Stopped mid-run to protect two concurrent ~20 GB submission builds (191 measured a build peak of
**19.83 GiB**). 199's own rehearsal stage does a full D/E/F build, and two of those plus eight
research agents does not fit in 64 GB.

**Nothing about the mandate changed and nothing it found is invalidated.** It had reached "building
the pre-staged CCLE matcher". Relaunch verbatim once 201's build is submitted and the 192–200 wave
has drained:

```
CLAUDE_SUBAGENT_MODEL=opus bash scripts/claude_subagent.sh 199-def-runbook experiments/199-def-runbook/MANDATE.md
```

It is the least time-critical of the ten (Oct 22 is 49 days out) and the only one that would have
run a second concurrent build, which is why it was the one stopped.
