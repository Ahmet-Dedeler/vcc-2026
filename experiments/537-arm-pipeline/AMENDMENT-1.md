COORDINATOR AMENDMENT 1 (2026-09-30 11:36Z; you were stopped and resumed in the same session to receive it)

1. Memory gate. The mandate's "free+inactive ≥ 33 GiB" is unreachable on this Mac: 31.0 GiB with nothing heavy running
   (534 hit the same; docs/CONFIRMATION-LEDGER.md, R6 "Deviation"). Replace it with AGENTS §9's gate plus the build check:
   memory_pressure ≥ 40 % free AND vm_stat free+inactive ≥ 24 GiB AND swap free ≥ 400 MiB AND
   `def_pipeline.other_builders()` == 0; one heavy job at a time; peak RSS ≤ 8 GB. The coordinator stopped your driver
   (drive537.sh, pid 72848; its log says so) and this session. Edit drive537.sh's gate now (it is not running) and re-arm it
   with nohup; never run two drivers. Everything in your queue stands.
2. Your §0 A claim (the DE channel's board noise is ~15× smaller than the 0.0045 rule) would change every DE bin and the
   final-phase adoption rule, so make it airtight before 10-02 12:00Z: list the exact entry pairs (ids) and their raw
   deltas; reproduce 230 §6's 0.0045 from 230's own pairs as a known answer (showing it is a spread of emission designs, not
   a re-roll of one design); say which DE designs your σ applies to (same emission law and seed, bulk-only changes) and
   whether a DE-knob change that re-draws the emission (a different call set, a hook) inherits it. The coordinator will
   re-derive it independently from data/board-entries before acting on it.
3. Everything else in MANDATE.md stands, including the first ranked list by 2026-10-02 12:00Z.
