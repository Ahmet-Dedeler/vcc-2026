#!/usr/bin/env python3
"""604: UNAPPLIED day-of patch for the identity-augmented CLEAN4 block prior (council 600; 639/622 science).

What it does to a tree (--root, default the repo):
  1. installs scripts/def_identity_prior.py (a copy of experiments/604-identity-prior/def_identity_prior.py with the
     pinned reference profiles' sha256 written into REFS_SHA256);
  2. edits scripts/def_sameline.py (POST-BUNDLE: 538 -> 554 -> 590, sha 22c6f0f0...) so that
     `--programme-source identity` (lambda 0.5), `identity:<lam>` or `identity@<lam>` gives EVERY context its identity-augmented CLEAN4 as the block-prior default
     (programme_<c> / pdown_<c>; unmatched contexts get w 0 and no rows; a matched context keeps its same-line rows).
     Every other --programme-source value, and no --programme-source at all, runs the existing code byte for byte.
def_day.sh needs NO change:
  - identity BLOCK prior: DEF_DAY_PROGRAMME_SOURCE=identity (or identity:0.5) reaches def_sameline.py; def_day.sh already
    refuses noprog flags with a programme source; the hook can then read the same artifact (DEF_DAY_MDE_SIGN=programme:0).
  - identity HOOK only (block prior stays CLEAN4): DEF_DAY_MDE_SIGN=hookbank:identity:0 -- def_day's existing hookbank step
    runs def_sameline.py --programme-source identity into $W/hookprog.npz and exports VCC_DE_MDE_SIGN=file:...:0 (tau 0).
Order on the day: AFTER the 10-05 bundle (538's applier, 554, 590). 598 / 612 / 630 do not touch def_sameline.py.
usage: uv run python experiments/604-identity-prior/patch604.py [--root DIR] [--check]
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = HERE / "def_identity_prior.py"
REFS = HERE / "work/flex_refs604.npz"
POST_BUNDLE_SAMELINE = "22c6f0f08c1f02883887a34690bc5642f172026433d90b19275a578f24d86259"

ANCHOR = '''    if a.programme_source:
        if a.programme_source not in BANK_DIRS:
'''
NEW = '''    if a.programme_source == "identity" or a.programme_source.startswith(("identity:", "identity@")):
        # 604 (council 600; science 639 d4cd5f + 622 7772dc): every context's block-prior default is its OWN
        # identity-augmented CLEAN4 (scripts/def_identity_prior.py: Flex-only references = the bundle's other contexts +
        # the pinned A/B/C and H1 profiles, one per identity call, never the context's own call; score quantile-mapped
        # onto CLEAN4's P; lam 0 == CLEAN4 exactly). Same-line rows of matched contexts are kept; others get w 0.
        # source names: `identity` (= DIP.DEFAULT_LAM, the council's ship point 0.5), `identity:<lam>`, or the colon-free
        # `identity@<lam>` (def_day's hookbank form, DEF_DAY_MDE_SIGN=hookbank:identity:<tau>, splits on ':')
        import def_identity_prior as DIP  # noqa: PLC0415
        try:   # 607: an empty or non-numeric dose is a clean refusal, and a negative dose is refused
            ilam = DIP.DEFAULT_LAM if a.programme_source == "identity" else float(a.programme_source[len("identity:"):])
        except ValueError:
            raise SystemExit(f"--programme-source {a.programme_source!r}: the identity dose must be a number "
                             "(identity, identity:<lam> or identity@<lam>)") from None
        if not ilam >= 0.0:
            raise SystemExit(f"--programme-source {a.programme_source!r}: the identity dose must be >= 0")
        rgenes, rcpm, rcalls = DIP.load_refs()
        if rgenes != list(official):
            raise SystemExit("--programme-source identity: the pinned reference profiles' gene axis differs from the bundle's")
        icpm = {c: DIP.ctrl_cpm(bundle / f"context_{c}.h5ad", official) for c in sorted(call)}
        P_c4, ps_c4 = DIP.load_clean4(official)
        ires = DIP.identity_prior(official, icpm, {c: str(v) for c, v in call.items()}, rcpm, rcalls, ilam, P_c4, ps_c4)
        for c in sorted(call):
            out[f"programme_{c}"] = ires[c]["psign"]
            out[f"pdown_{c}"] = ires[c]["pdown"]
            out[f"programme_ok_{c}"] = np.asarray(True)
            if c not in matched:
                out[f"w_{c}"] = np.asarray(0.0)
                out[f"targets_{c}"] = np.zeros(0, dtype=np.int64)
                out[f"rows_{c}"] = np.zeros((0, len(official)), dtype=np.float32)
                out[f"signs_{c}"] = np.zeros((0, len(official)), dtype=np.int8)
                matched.append(c)
            prov["contexts"].setdefault(c, {"call": call[c]})["programme_source"] = {
                "identity_prior": ires[c]["prov"], "refs_sha256": DIP.file_sha256(DIP.REFS_PATH),
                "clean4_sha256": DIP.CLEAN4_SHA256}
        matched = sorted(matched)
        print(f"[sameline] --programme-source {a.programme_source}: every context's programme = its identity-augmented "
              f"CLEAN4 (lam {ilam:g}); contexts {matched} " + json.dumps({c: ires[c]["prov"] for c in ires}), flush=True)
    elif a.programme_source:
        if a.programme_source not in BANK_DIRS:
'''


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def installed_source() -> str:
    text = SRC.read_text()
    refs_sha = sha(REFS)
    old = "REFS_SHA256 = None   # pinned by patch604 when it installs this file (see experiments/604-identity-prior/RESULTS.md)"
    if text.count(old) != 1:
        raise SystemExit("def_identity_prior.py: the REFS_SHA256 placeholder is not there exactly once")
    return text.replace(old, f'REFS_SHA256 = "{refs_sha}"   # pinned by patch604 (flex_refs604.npz)')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    root = Path(a.root)
    same = root / "scripts/def_sameline.py"
    dst = root / "scripts/def_identity_prior.py"
    text = same.read_text()
    if "--programme-source identity:" in text or 'startswith(("identity:", "identity@"))' in text or 'startswith("identity:")' in text:
        raise SystemExit("def_sameline.py already carries 604's identity programme source (already applied)")
    if text.count(ANCHOR) != 1:
        raise SystemExit("def_sameline.py: the programme-source anchor is not there exactly once (apply the 10-05 bundle first?)")
    got = sha(same)
    if got != POST_BUNDLE_SAMELINE:
        print(f"[604] WARNING def_sameline.py is {got[:12]}, not the post-bundle {POST_BUNDLE_SAMELINE[:12]} "
              "(anchors still match; proceed only if that is intended)", flush=True)
    if dst.exists() and dst.read_text() != installed_source():
        raise SystemExit(f"{dst} exists and differs from 604's source")
    new_text = text.replace(ANCHOR, NEW)
    compile(new_text, "def_sameline.py (604)", "exec")
    compile(installed_source(), "def_identity_prior.py (604)", "exec")
    if a.check:
        print(f"[604] --check OK: def_sameline.py {got[:12]} -> {hashlib.sha256(new_text.encode()).hexdigest()[:12]}; "
              f"def_identity_prior.py would be installed ({hashlib.sha256(installed_source().encode()).hexdigest()[:12]})")
        return
    same.write_text(new_text)
    dst.write_text(installed_source())
    print(f"[604] applied: def_sameline.py {sha(same)}  def_identity_prior.py {sha(dst)}  (refs {sha(REFS)[:16]})")


if __name__ == "__main__":
    sys.exit(main())
