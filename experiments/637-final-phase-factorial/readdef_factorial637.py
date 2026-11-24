"""637 — factorial reads of a D/E/F entry: one bulk switch AND one DE switch in the same entry (H-630).

  uv run python experiments/637-final-phase-factorial/readdef_factorial637.py --factorial BULK_ROW DE_ROW ENTRY BASE
  uv run python experiments/637-final-phase-factorial/readdef_factorial637.py --validate FACT X_ALONE Y_ALONE BASE
  (either form takes an optional `--leak-band v1|v2`; default v2)

A thin layer over 531's readdef.py, which is imported read-only and never modified, so its 74 verdicts are untouched. It
reuses readdef's own row rules. An optional one-line dispatch into readdef.main is in `readdef637.diff` (UNAPPLIED).

Why it works (637/PROPOSAL.md §1, measured on 8 of our board pairs in data/board-entries/):
  * a DE switch never moves the bulk members: Δpds raw exactly 0, |Δmse raw| ≤ 2.4e-8;
  * a bulk switch leaks into the DE members by a small amount: v1 (8 pairs) −0.00020 … −0.00039 overall; v2 (9 bulk-only
    board pairs, 658) mean −0.00017, sd 0.00021, range −0.00041 … +0.00019. See LEAK_BANDS.
So the bulk switch is read on pds + mse (its overall := bulk part + the leak allowance, charged to the bulk switch), and the
DE switch on fid + reach + jac + nmae after removing the leak allowance. A verdict counts only if it is the SAME at every
leak in the band LEAK_BAND; otherwise the read is "inconclusive_leak_band" and nothing is adopted from it.

Allowed: bulk rows 2, 3, 4 × DE rows 6, 7 (prior-only DE rows).
Refused:
  * row 5 (moves both channels);
  * row 8-9 (a ballast-moving bulk step leaked +0.0014 raw nmae on the board, 548 / 433k→446b, which is comparable to
    8-9's −0.002 raw nmae bin);
  * row 10-12 (a combination is read on its own).
In factorial mode the DE rows' "pds/mse pinned" falsifier cannot be read (the bulk switch moves pds/mse by design). It is
replaced by the bulk-channel check above, which holds by construction for DE switches (0 / 2.4e-8 on the board).

--validate (the A/B/C rehearsal, preregistered in PROPOSAL.md §1): PASS iff
  |Δpds_raw(FACT) − Δpds_raw(X)| ≤ 1e-6 and |Δmse_raw(FACT) − Δmse_raw(X)| ≤ 1e-6, and
  |ΔDE(FACT) − ΔDE(Y) − LEAK_CENTRE| ≤ 0.0006 overall,
all against the same BASE. Exit 0 on PASS, 3 on FAIL.
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "531-cloud-final-phase-tooling"))
import readdef as R  # noqa: E402  (read-only)

# The bulk switch's DE leak (overall units), versioned so a read registered under an older band reproduces exactly.
#   v1 (637, 2026-10-02): 8 of our board pairs, a consistently negative −0.00020 … −0.00039 -> centre −0.0003 ± 0.0002.
#   v2 (2026-10-03, [7772dc]'s 9 bulk-only board pairs in 658 / cdc1445c): mean −0.000174, sd 0.000213, rms 0.000266,
#       range −0.00041 … +0.00019 (a re-roll around a small negative mean, not a fixed offset) -> centre −0.0002 ± 0.0004.
# Default v2. `--leak-band v1` reproduces any read registered under v1 (e.g. 658's ALLSW R1).
LEAK_BANDS = {"v1": (-0.0003, (-0.0001, -0.0003, -0.0005)),
              "v2": (-0.0002, (0.0002, -0.0002, -0.0006))}
LEAK_VERSION = "v2"
LEAK_CENTRE, LEAK_BAND = LEAK_BANDS[LEAK_VERSION]


def set_leak_band(version: str) -> None:
    """Select the leak band (v1 / v2) for this process; the version is recorded in every verdict."""
    global LEAK_VERSION, LEAK_CENTRE, LEAK_BAND
    if version not in LEAK_BANDS:
        raise SystemExit(f"unknown leak band {version!r}; one of {sorted(LEAK_BANDS)}")
    LEAK_VERSION = version
    LEAK_CENTRE, LEAK_BAND = LEAK_BANDS[version]
VALIDATE_PIN = 1e-6
VALIDATE_DE_TOL = 0.0006
BULK_ROWS = ("2", "3", "4")
DE_ROWS = ("6", "7")
REFUSED = {
    "5": "row 5 moves both channels: read it as a single-switch entry",
    "8-9": "row 8-9 is refused in factorial mode: a ballast-moving bulk step leaked +0.0014 raw nmae (548), comparable to its −0.002 bin",
    "10-12": "a combination is read on its own, not factorially",
    "1": "row 1 is R0, not an arm",
    "13-14": "rows 13-14 are the reserve, not an arm",
}


def _key(row: str) -> str:
    k = R.ALIASES.get(str(row).strip())
    if k is None:
        raise SystemExit(f"unknown row {row!r}")
    return k


def bulk_view(d: dict, leak: float) -> dict:
    """The bulk switch's delta: bulk members as read, DE members replaced by the leak it is charged with."""
    v = copy.deepcopy(d)
    v["overall"] = d["bulk"] + leak
    v["de"] = leak
    v["fnr"] = leak
    return v


def de_view(d: dict, leak: float) -> dict:
    """The DE switch's delta: DE members minus the bulk switch's leak; bulk members zeroed (they belong to the bulk switch)."""
    v = copy.deepcopy(d)
    v["raw"]["pds"] = 0.0
    v["raw"]["mse"] = 0.0
    v["scaled"]["pds"] = 0.0
    v["scaled"]["mse"] = 0.0
    v["bulk"] = 0.0
    v["de"] = d["de"] - leak
    v["fnr"] = d["fnr"] - leak
    v["overall"] = d["de"] - leak
    return v


def _bulk_verdict(key: str, v: dict) -> tuple[str, bool | None, str, str]:
    if key == "2":
        code, adopt, _which, bin_s, adoption = R._row2(v, None)
    elif key == "3":
        code, adopt, _which, bin_s, adoption = R._row3(v, None)
    else:
        ok, bin_s = R._row4_one(v)
        code, adopt = ("adopt_arm" if ok else "base_stands"), ok
        adoption = "ADOPT this arm." if ok else "do not adopt. The base stands."
    return code, adopt, bin_s, adoption


def _de_verdict(key: str, v: dict) -> tuple[str, bool | None, str, str]:
    if key == "6":
        return R._row6(v)
    return R._row7(v)


def factorial_verdict(bulk_row: str, de_row: str, arm: dict, base: dict) -> dict:
    bk, dk = _key(bulk_row), _key(de_row)
    for k in (bk, dk):
        if k in REFUSED:
            raise SystemExit(f"REFUSED: {REFUSED[k]}")
    if bk not in BULK_ROWS:
        raise SystemExit(f"bulk row must be one of {BULK_ROWS}, got {bk}")
    if dk not in DE_ROWS:
        raise SystemExit(f"DE row must be one of {DE_ROWS}, got {dk}")
    hidden = R.members_hidden(arm) + R.members_hidden(base)
    if hidden:
        return {"code": "hidden", "bulk": None, "de": None, "lines": [f"MEMBERS HIDDEN: {hidden}"]}
    d = R.delta(arm, base)
    lines = R._fmt_block("factorial entry", arm, base, d)
    out = {"code": "read", "bulk_row": bk, "de_row": dk, "lines": lines, "leak_band": LEAK_VERSION}
    for chan, viewf, verd, key in (("bulk", bulk_view, _bulk_verdict, bk), ("de", de_view, _de_verdict, dk)):
        reads = [(leak,) + tuple(verd(key, viewf(d, leak))) for leak in LEAK_BAND]
        codes = {r[1] for r in reads}
        centre = [r for r in reads if r[0] == LEAK_CENTRE][0]
        if len(codes) == 1:
            out[chan] = dict(code=centre[1], adopt=centre[2], bin=centre[3], adoption=centre[4],
                             leak_band_codes=[r[1] for r in reads])
        else:
            out[chan] = dict(code="inconclusive_leak_band", adopt=False,
                             bin=f"verdict changes inside the leak band {LEAK_BAND}: {[(r[0], r[1]) for r in reads]}",
                             adoption="do not adopt from this factorial read; re-read this switch as a single-switch entry.",
                             leak_band_codes=[r[1] for r in reads])
    return out


def validate(fact: dict, x_alone: dict, y_alone: dict, base: dict) -> dict:
    hidden = R.members_hidden(fact) + R.members_hidden(x_alone) + R.members_hidden(y_alone) + R.members_hidden(base)
    if hidden:
        return {"pass": None, "lines": [f"MEMBERS HIDDEN: {hidden}"]}
    df, dx, dy = R.delta(fact, base), R.delta(x_alone, base), R.delta(y_alone, base)
    e_pds = abs(df["raw"]["pds"] - dx["raw"]["pds"])
    e_mse = abs(df["raw"]["mse"] - dx["raw"]["mse"])
    e_de = abs(df["de"] - dy["de"] - LEAK_CENTRE)
    ok = (e_pds <= VALIDATE_PIN + R.PIN_SLACK) and (e_mse <= VALIDATE_PIN + R.PIN_SLACK) and (e_de <= VALIDATE_DE_TOL)
    lines = [f"bulk: |dpds(F) - dpds(X)| {e_pds:.3e} (<= 1e-6), |dmse(F) - dmse(X)| {e_mse:.3e} (<= 1e-6)",
             f"DE: |dDE(F) - dDE(Y) - ({LEAK_CENTRE})| {e_de:.6f} overall (<= {VALIDATE_DE_TOL})",
             f"measured leak dDE(F) - dDE(Y) = {df['de'] - dy['de']:+.6f} (expected {LEAK_CENTRE} within the band {LEAK_BAND}; "
             f"leak band {LEAK_VERSION})"]
    return {"pass": ok, "lines": lines, "leak_band": LEAK_VERSION}


def render(v: dict) -> str:
    out = list(v["lines"])
    if v.get("code") == "hidden":
        return "\n".join(out + ["no factorial read: members hidden. The dark plan stands."])
    for chan in ("bulk", "de"):
        c = v[chan]
        row = v["bulk_row"] if chan == "bulk" else v["de_row"]
        out.append(f"[{chan} switch, row {row}] code {c['code']} (codes across leak band {c['leak_band_codes']})")
        out.append(f"  BIN: {c['bin']}")
        out.append(f"  ADOPTION: {c['adoption']}")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    a = list(sys.argv[1:] if argv is None else argv)
    if "--leak-band" in a:
        i = a.index("--leak-band")
        set_leak_band(a[i + 1])
        del a[i:i + 2]
    print(f"[637] leak band {LEAK_VERSION}: centre {LEAK_CENTRE:+.4f}, band {LEAK_BAND}", file=sys.stderr)
    if len(a) == 5 and a[0] == "--factorial":
        v = factorial_verdict(a[1], a[2], R.load_entry(a[3]), R.load_entry(a[4]))
        print(render(v))
        return 2 if v.get("code") == "hidden" else 0
    if len(a) == 5 and a[0] == "--validate":
        r = validate(*(R.load_entry(x) for x in a[1:5]))
        print("\n".join(r["lines"]))
        if r["pass"] is None:
            print("VALIDATION NOT READ: members hidden")
            return 2
        print("VALIDATION PASS: the factorial design is additive on this pair" if r["pass"] else
              "VALIDATION FAIL: keep one switch per D/E/F entry")
        return 0 if r["pass"] else 3
    print(__doc__.strip(), file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
