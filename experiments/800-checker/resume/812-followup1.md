FROM THE LAB LEAD (2026-10-16 08:15Z). Your context is intact. Your NO-GO on the four emission hypotheses is accepted and
logged (docs/LAB-LOG.md 08:12Z). Your reading, that the DE gap is sign skill on responding genes, is the most useful DE
result the lab has. Turn it into an instrument, which is the next step:

1. Finish RESULTS.md (the four screens, their KAs, tries, and rates812's class table).
2. Build `skill812.py`, a DE sign-skill scorer for ANY model's predictions. Input: per-target predicted lfc or sign rows on H1's
   targets (h1-test 85 + h1-trainval 126) from a model that never saw H1 (VCC25 H1 out of every bank it uses). Output: per-class
   sign accuracy (a0 called non-DE, a1 called DE, head a1, u1 uncalled DE, by expression bin), with paired target-bootstrap
   CIs against CLEAN4, and rates812's overall-equivalent price at H1's per-cell rates. Report it in checker units with the
   caveat that the rates come from H1 (one line).
   - KA: CLEAN4 through skill812 reproduces its recorded class accuracies (a0 0.517 on called non-DE, etc.; cite the source
     row), and the production emitter's realized signs reproduce rates812's own numbers.
3. Run it on two existing sign sources, as information only (no proposal yet): (a) the per-target mean of the non-H1 public
   lines' rows (811's line banks, read-only, H1 out), plain and per-line-centred; (b) the production bulk delta's own sign.
   Say where, if anywhere, either beats CLEAN4 by more than its CI.
4. Stop with a short message. 811 will use skill812 as its DE objective. If (a) beats CLEAN4 in a class that pays, write a
   PROPOSAL for a DE-H1 price (the hook/prior swap through 803's de803.py) and the lab lead will judge it.
Budget: $10 of EC2 (H1 truths and 803's members are on S3); local only sub-second < 100 MB aggregations, else gate602 or EC2;
window800 rules.
