FROM THE LAB LEAD: 801 has registered the **ABCp** folds. Check `experiments/801-lab-folds/work/registry.json` for the
`ABCp/*` ids and that m811.py's ABCp detection (G1) finds them. Your m811.py must still hash to
`e00192e644828a0d00b18ca5f927f93893ef7012661aac84a4736ca592a98394`; if it does not, stop and tell me why.

Run the frozen decision read, `cloud/jobs-s2b-abcp.txt` (c0_l25 and its C-partner c1_l25 on every ABCp fold), through
labcloud.sh (on-demand, Project=lab, ≤ $1). Engine KAs on every fold; a missing reference is a failed KA.

Apply PROPOSAL-2 §5 on ABCp under the **N = 10 branch** (ABCp read after the 062 read of the same arms; e(10) 1.547), literal
pds, pdsl reported. Report in RESULTS.md: BULK-ABCp per arm (point, CI, folds > 0, after-e(N) bound), the 062 read beside it,
and the many-line increment over sources-only on ABCp if the frozen jobs produce it. λ 0.10 is NOT part of this read; if you
think it should be, write why in RESULTS.md and count it as a new try, but do not run it.

Rule 6 still blocks a dossier until DE-H1 prices the external rows; 803 is building the bulk_transform hook. Say plainly
whether ABCp alone would pass §5. Append wear rows to work/wear811.jsonl. Then stop.
