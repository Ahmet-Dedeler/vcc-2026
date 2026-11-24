"""644: write the DEF-RUNBOOK amendment as a full copy (DEF-RUNBOOK.644.md) and a unified diff against the live
docs/DEF-RUNBOOK.md (DEF-RUNBOOK.644.diff). The shared doc is never written; the coordinator merges.

    uv run python experiments/644-runbook-amendments/make_amendment644.py
"""
from __future__ import annotations

import difflib
import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "docs/DEF-RUNBOOK.md"

TOP_ANCHOR = "## THE DAY IN ONE PAGE (2026-09-25; supersedes §2's command block where they differ)\n\n"
TOP = r"""> **2026-10-02 (644: the guards and switches added since 09-30; these override anything older on this page).**
> - **Patch order:** 538 → 554 → 590 on 10-05, 641 at checklist step 6b, then 598 (tercilec) → 612 → 642 on 10-06.
> - **Every new guard is silent on PASS.** An override is legitimate only in the case its row names.
> - **To override:** export it in the shell before `def_day.sh`, in a FRESH work dir (a tripped guard stops the day and
>   nothing resumes in place), and write the reason into the round's RESULTS.md before the build starts. `def_day.sh`
>   logs each override at start and records it in `<work>/def_day.final.env`.
>
> | guard (source) | the row that trips | meaning | override, and when it is legitimate | record |
> |---|---|---|---|---|
> | preflight library spread (554) | `[X] control library q75/q25 <= 3.0 (554: …)` | The K-class emitter's library classes are the control quartiles. Above q75/q25 ~3.5 its steering dose (CLEAN4 block prior, MDE, hook) drops to 0 with no error (553); by ~7 the placement crashes. A/B/C sit at 1.57–1.95 | `DEF_DAY_ALLOW_WIDE_LIBS=1`, only if every context reads ≤ 3.5; the 554 dose row below must then PASS. Above 3.5 in any context: STOP (no tested fallback exists for the real day) | the three spreads (preflight.json `library_spread_q75_q25`) |
> | preflight depth (641) | `[X] control library median within 15% of A/B/C's 20,030 (641: …)` | D/E/F were released at another depth than A/B/C's ~20k (Arc's full data is 29,482 UMI/cell). Every per-gene quantity is refitted from the controls, but the DE call-set size was tuned at 20k | All three contexts off by the same factor: price that ratio on H1's regime instrument (553's h1r7 how-to), then `DEF_DAY_ALLOW_DEPTH=1`. One context off alone: STOP (a damaged or mislabelled file) | the three medians and `library_median_ratio`, plus the price |
> | identity floor (590) | `[sameline] X: … min top-1 correlation 0.6xx < 0.65; context unchanged (590)` | A line absent from CCLE is matched, unanimously, to its nearest CCLE line. Correct calls sit at ≥ 0.68: A/B/C 0.772 / 0.719 / 0.733; real-Flex K562 0.693, HAP1 0.747. Absent lines on Flex: HEK293 0.647 (split), H1 0.653–0.665 (split) | none. The abstention is the design; a non-CCLE line switches on only through 642 | every context's call, n/8, `min_top1_score` and consensus top 5 (identity.json) |
> | postbuild dose (554) | `emitter steering dose per context: median prior_lambda >= 0.9, mde_lambda >= 0.9, call_set_missed <= 100 (554)` | A context's DE design was delivered at low or zero dose (553; 458's programme builds ran at λ 0.24–0.44). Healthy: prior 1.0, MDE 0.951–1.0 (A/D lowest), missed 4–21.5 | `DEF_DAY_ALLOW_LOW_DOSE=1`, only when the sole trip is an MDE median of 0.85–0.90 with prior ≥ 0.9 and missed ≤ 100 (≈ −0.0003 overall). Never for a prior below 0.9 or any dose near 0. Postbuild must run on the pinned builder, and the failed day's exit trap has already restored it: rerun the day in a fresh dir with the override exported (~2 h) | postbuild.json `dose_guard` (per-context values and margins) |
> | HEK293 contingency (642, OFF) | `[642] X: with the hek293 reference row: … -> CALLED HEK293` | Stage 2 named the context HEK293: 8/8, min_top1 ≥ 0.75, and the reference row beats every CCLE line | see "The HEK293 path" below; the coordinator's decision only | identity.json `nonccle_ref`, plus the decision text |
>
> **The Oct 22 identity reading rule (590).**
> - Before the build, copy into the round's RESULTS.md, for every context: the call, n/8 methods, `min_top1_score`, and
>   the consensus top 5 with lineages. `def_identity.py` prints all of them.
> - Read any unanimous call with min_top1 in 0.60–0.70 by hand. The top 5 should be one lineage, show no sister-line or
>   duplicate pattern, and the call's tissue should be plausible.
> - If in doubt, do not switch same-line rows on. Abstaining costs at most the same-line gain; a false switch-on puts
>   another line's rows into a context.
> - **How to switch a context off:** rerun in a fresh dir with `DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>` (644; the last 10-06
>   step).
>   - `def_sameline --skip-contexts` marks those contexts unmatched whatever their call, so each takes exactly an
>     unmatched context's path: no row blend, no `pt` signs.
>   - With a programme source (a bank, or 604's identity prior) a skipped context still gets the source's programme at
>     w 0, like any unmatched context.
>   - Unknown names stop the day.
>   - Do not use sameline w `0` for this: at w 0 the builder skips the blend, but `pt` still applies the bank's per-target
>     signs.
> - Never `--force` an identity on the day without a written reason in RESULTS.md: `--force` bypasses unanimity and the
>   floor.
>
> **The HEK293 path (642; OFF by default, applied 10-06).** It applies only to a context that abstained (not registered, or
> below the 590 floor) and whose top 5 or lineage suggests HEK293 (kidney or neural, embryonic-like).
> 1. **Diagnostic, read-only:**
>    `DEF_DAY_ALLOW_NONCCLE_REF=hek293 uv run python scripts/def_identity.py --bundle data/challenge --contexts D,E,F --out <work>/identity-hek293.json`
> 2. **Read** `nonccle_ref.per_context[X]`: `accepted` (8/8 and min_top1 ≥ 0.75), `row_top1_methods`, and the top 5.
>    Shallow or odd chemistries reached 6/8 at 0.67–0.69 without being HEK293 (642), so the bar is real, not a formality.
> 3. **The coordinator decides, in writing.** The open question stays open: Orion HEK293T's rows have never been priced on
>    an Arc-Flex HEK293 context. 450's same-line value was measured on H1 only.
> 4. **The day:** prefix the usual `def_day.sh` call with
>    `DEF_DAY_ALLOW_NONCCLE_REF=hek293 DEF_DAY_NONCCLE_DECISION="<who, when, which context, why>"`, in a fresh work dir.
>    - `def_sameline` then reads `orion_hek293t` for that context (294/300 on the A/B/C panel, w 2).
>    - A stage-2 HEK293 call vetoes that context's CCLE call, whatever it was.
>    - Without the decision text, the context is unchanged.
> 5. **The reference is pinned:** `data/identity_refs/hek293_nonflex_ref642.npz`, sha256 `1b829de0…`, built from non-Flex
>    controls only. Never rebuild it on the day, and never from Flex data: a Flex-built row pulled Flex HAP1 to 8/8.
>
> **10-06 knobs.** `DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>` (644, last): a doubtful context takes the unmatched same-line path;
> logged and recorded. 598's `DEF_DAY_SHARPEN_PT=tercilec:2.5,2.0,1.75` replaces the tercile form once its KAs pass and the
> recipe decides it; `october_recipe.py command` then emits it. Re-run 644's `verify644.py` after 10-06, and again on
> Oct 21.
>
> **PENDING switches** (612, 630, 639/604) are not part of the day until their reads decide them. See the PENDING
> section at the end.

"""

FB_ANCHOR = "   - If nothing comes back, the dark-board plan above stands unchanged.\n"
FB = r"""   - **Factorial reads (637, 2026-10-02; reader accepted, validation PENDING).** If members come back, each entry after
     R0 may carry ONE bulk switch AND ONE prior-only DE switch:
     - the bulk switch is a row 2/3/4, read from pds + mse;
     - the DE switch is a row 6/7, read from the DE members after a −0.0003 ± 0.0002 bulk-leak allowance.

     Command: `uv run python experiments/637-final-phase-factorial/readdef_factorial637.py --factorial BULK_ROW DE_ROW ENTRY BASE`.
     - A verdict counts only if it is the same at every leak in the band; otherwise it reads `inconclusive_leak_band`
       and nothing is adopted.
     - Refused: row 5 (moves both channels), rows 8–9 (ballast-moving), rows 10–12 (combinations).
     - Use it only after its A/B/C validation arm reads PASS (`--validate FACT X_ALONE Y_ALONE BASE`, exit 0; the
       all-switches arm, 10-08 … 10-10). Until then, one switch per entry.
"""

TBL_ANCHOR = "| source download needed | it is not: the source is pinned locally (§1). If the file is gone: `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad`, 16.8 GB, then verify `c355f535…` | |\n"
TBL = r"""| `FAIL [X] control library q75/q25 <= 3.0 (554: …)` | wide control libraries: the emitter's steering dose would collapse | the top block's 554 row: ≤ 3.5 everywhere → `DEF_DAY_ALLOW_WIDE_LIBS=1` in a fresh dir; above → STOP |
| `FAIL [X] control library median within 15% of A/B/C's 20,030 (641: …)` | D/E/F released at another depth | the top block's 641 row |
| `[sameline] X: … min top-1 correlation … < 0.65; context unchanged (590)` | the line is probably absent from CCLE | expected; read the top 5. Only 642 can switch a non-CCLE line on |
| `postbuild FAILED …` with a `(554)` row in postbuild.json | a context's DE design was delivered at low dose | the top block's 554 row |
| `STOP: 554: the final build wrote no rc_audit dump` | the build ran without `VCC_RC_AUDIT_DUMP` (a pre-554 `def_day.sh`) | check `def_day.sh`'s sha against the post-bundle table (594), then rerun in a fresh dir |
| `642: … is not the pinned …` or `642: reference … lists Flex (or no) sources` | the HEK293 reference file changed | `git checkout -- data/identity_refs/hek293_nonflex_ref642.npz`; never rebuild it on the day |
| `[sameline] X: the non-CCLE reference called HEK293; …; context unchanged (642)` | the knob or the decision text is missing | by design; see the HEK293 path |
| `642: DEF_DAY_ALLOW_NONCCLE_REF='…' is not one of ['hek293']` | a typo in the knob | fix it |
| `STOP: 644: DEF_DAY_SAMELINE_SKIP names X, which is not a manifest context (…)` | a typo, or a label not in the manifest | use the manifest's context names |
"""

FILES_ANCHOR = "- `experiments/199-def-runbook/work/run-r/` — the realistic bundle's cheap stages.\n"
FILES = r"""- 2026-10-02 additions (644):
  - `experiments/554-dose-guard/` (spread and dose guard);
  - `experiments/590-identity-validation/` (identity floor);
  - `experiments/641-depth-row/`;
  - `experiments/642-hek293-contingency/` and `data/identity_refs/` (the HEK293 contingency);
  - `experiments/637-final-phase-factorial/` (factorial reader);
  - `experiments/644-runbook-amendments/`: `verify644.py` re-runs every command these rows cite against scratch
    post-bundle trees (24/24 on 10-02), and `DECISION-CARD.md` is the one-page day-of card.
"""

PENDING = r"""
## PENDING switches: placeholders (644, 2026-10-02). Not part of the day until a read decides them.
If adopted, each one needs:
- `october_recipe.py decide <switch> <value> --evidence …`;
- a `def_day.sh` pass-through;
- an inert KA and a rehearsal.

Only then does its row move into the top block.

- **612, the key-gated uncalled-gene hook** (council [21dd4c]).
  - Builder knob `VCC_DE_MDE_SIGN_KEYGATE=<q>`; needs `VCC_DE_MDE_SIGN`.
  - A robustness lever for unseen lines: worst-case regret 0.0020 against 0.0115 for the plain hook, ~0 EV on
    CLEAN4-friendly lines.
  - Read: a 10-08+ arm. No DEF_DAY pass-through yet.
- **630, rank-neutral UP magnitude** (council [d4cd5f]).
  - Builder knob `VCC_DE_UPRN=<s|b>_k<K×100>`, applied after zs.
  - H1: nmae −0.0038 … −0.0066 raw, with the other five members exact and 0 call flips.
  - Read: 10-08+. No pass-through yet.
- **639 / 604, the identity-axis sign prior** (council; [4d0207] does the plumbing).
  - Delivered through the same-line path (`VCC_SAMELINE=<npz>:<sha>:0`, w 0, no rows) and computed on the day from
    D/E/F's own controls. Ship point λ 0.5.
  - Hard gate before any slot: per-context median prior_lambda ≥ 0.9 (554's dump).
  - The day-of `def_day.sh` patch and its knobs: TBD.
"""


DL_OLD = ("**Get Ahmet's go-ahead for the downloads before Oct 22** (an agent's rules require explicit permission, with filename,\n"
          "source and size, for any download). The day needs two: `vcc datasets download controls` (the D/E/F bundle,\n"
          "~1.2 GB like A/B/C's) and a `vcc` CLI upgrade. Ask on Oct 21 so the day does not wait on a reply.\n")
DL_NEW = ("**Downloads need no go-ahead (Ahmet, 2026-10-02: \"never ask my permission again for such a basic thing\").** The day\n"
          "needs two: `vcc datasets download controls` (the D/E/F bundle, ~1.2 GB like A/B/C's) and a `vcc` CLI upgrade. Run\n"
          "them, and record each file, its source, size and sha256 in the round's RESULTS.md.\n")


def insert_after(text: str, anchor: str, block: str) -> str:
    assert text.count(anchor) == 1, anchor[:80]
    return text.replace(anchor, anchor + block)


def main():
    src = SRC.read_text()
    out = insert_after(src, TOP_ANCHOR, TOP)   # line 77 (downloads) and §6 item 8 (H-631) are the coordinator's, already done
    out = insert_after(out, FB_ANCHOR, FB)
    out = insert_after(out, TBL_ANCHOR, TBL)
    out = insert_after(out, FILES_ANCHOR, FILES)
    out = out.rstrip("\n") + "\n" + PENDING
    (HERE / "DEF-RUNBOOK.644.md").write_text(out)
    diff = difflib.unified_diff(src.splitlines(keepends=True), out.splitlines(keepends=True),
                                fromfile="a/docs/DEF-RUNBOOK.md", tofile="b/docs/DEF-RUNBOOK.md")
    (HERE / "DEF-RUNBOOK.644.diff").write_text("".join(diff))
    print("base docs/DEF-RUNBOOK.md sha256", hashlib.sha256(src.encode()).hexdigest())
    print("amended sha256", hashlib.sha256(out.encode()).hexdigest(), f"(+{len(out.splitlines()) - len(src.splitlines())} lines)")


if __name__ == "__main__":
    main()
