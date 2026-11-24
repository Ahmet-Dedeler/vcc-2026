FROM THE LAB LEAD: 801 has registered the **ABCp** folds (A/B/C held-out operators refitted with the day recipe). Check
`experiments/801-lab-folds/work/registry.json` for the `ABCp/*` fold ids and run the engine KAs (`E.load(check=True)`) on each;
a failed or missing KA stops the read.

Run the funded information read from `LEAD-DECISION-1.md`: the same five arms as the GO run, on every ABCp fold, through
labcloud.sh (≤ $0.5, on-demand, Project=lab). Before any number is read, record price810.py's and decide810.py's sha256 and
confirm both match the GO run's (no re-tune). ABCp folds keep `F.panel == "ABC"` (801 RESULTS §12.3), so price810's A/B/C
branch applies exactly as on 062's folds (same incumbent, same 707 comparison); record which incumbent it used.

Report in RESULTS.md as a new section "ABCp information read":
- per arm: BULK-ABCp point, CI, folds > 0, pds share of the mse gain; and the same arms on 062's BULK-ABC beside it;
- decide810 applied as frozen to ABCp;
- the one sentence that answers the question this read buys: does the W/ABC split shrink when A/B/C's operator is the day
  operator (i.e. was part of the −0.0020 the 062 whitener)?
- wear rows: one per arm read on ABCp, with N counted as BULK-ABC reads of these arms (062 then ABCp; the N = 6 branch plus
  this read).
This is a read, not a reopening. A DOSSIER still needs every G1 condition on ABCp. Then stop.
