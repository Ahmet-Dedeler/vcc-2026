#!/usr/bin/env python3
"""Refuse to submit a package that has already been scored.

STANDING RULE (Ahmet, 2026-09-01; AGENTS.md section 6, docs/LOG.md):
    Never spend a submission slot on leaderboard rank, and never re-submit an already-scored
    package.  We get 2 slots/day and they are the only measurement of the real objective.
    Leaderboard rank is informational -- only the last entry after the Oct 22 test release is
    prize-eligible.  The goal is FIRST PLACE, not a better end-of-day row.

This exists because the rule was known and articulated in this repo hours before it was broken:
on 2026-08-31 a19 regressed, the displayed rank fell 12 -> 26, and the day's second slot went on
re-uploading a15 byte-identical -- buying rank and zero information.  Prose did not stop it, so
this is the version that can.

    # before `vcc submit`
    uv run python scripts/check_submission_novelty.py check submissions/a21-resample-x014.vcc

    # after a successful submit
    uv run python scripts/check_submission_novelty.py record submissions/a21-resample-x014.vcc \
        --entry-id QGpNVMwzzdSlhUoQtRdp --note "resampling emitter x=0.14"

`check` exits 1 if the package's sha256 is already in the ledger.  It is deliberately impossible
to pass with a flag -- if you genuinely need to re-submit identical bytes, delete the ledger
entry by hand and write down why in docs/LOG.md.

**The ledger is tracked by git and a missing ledger is a refusal, not a pass (242 #3; 247).**
Until 247 the ledger lived under the gitignored `submissions/*`, `load()` returned `{}` when it
was missing, and it had not been written to since a23 -- so on 2026-09-04 `check` cleared
`submissions/abc-knee.vcc`, scored 0.1740437 the day before, with the same `OK: novel package`
line a healthy run prints.  Now `.gitignore` carries `!submissions/SUBMITTED_LEDGER.json`, the
five scored packages since 2026-09-02 are backfilled, and a missing ledger exits 2 with a message.
`--ledger` exists so that refusal can be demonstrated without touching the real file.

**256 (2026-09-05): the ledger keys on the `.vcc` digest, and that digest is CLI-VERSION-DEPENDENT.**
`vcc` 0.2.0 (installed 22:35 local 09-04) writes `X/indptr` as int64 and adds a `meta.json` member, so a
byte-identical h5ad preps to a different `.vcc` sha than under the old CLI (`cd4svd5`: `6e8b41c6…` scored,
`7f8336fb…` re-prepped, every value inside equal). Keying on the upload alone would have printed
`OK: novel package` for a re-prep of a scored prediction. So `check` also resolves the PREDICTION's
identity -- `output.sha256` from the sibling `<stem>.h5ad.meta.json` (the builder's record), else the
sibling `<stem>.h5ad` hashed -- and refuses when it matches any record's `h5ad_sha256`; `record` stores
`h5ad_sha256` and `vcc_cli_version`.

**256b (2026-09-05 02:50Z): the prediction identity is REQUIRED, and bound to the container.** A package
whose prediction cannot be identified -- no sibling meta and no sibling h5ad, or a meta that is unreadable
or names no 64-hex `output.sha256` -- is REFUSED, not passed on the container digest: "cannot identify" is
not "novel". And when the `.vcc` carries a `meta.json` member (vcc >= 0.2.0 writes `n_obs`, `n_vars`,
`nnz`), those three must equal the sibling meta's `output` block where it records them, so a sibling meta
that belongs to a different prediction is caught rather than trusted by file name. Pass
`--allow-unidentified` ONLY for a package built outside the builder (a15-era) with a note in docs/LOG.md.

**256c (2026-09-05 13:00Z): CONTENT-BOUND provenance.** A sidecar can name the wrong prediction and still
agree in shape, so `check` now also computes the prediction's CONTENT digest from INSIDE the container --
`pred.h5ad.zst` streamed through `zstd -d` to a scratch file, then sha256 over the fields the scorer reads
(`X/data`, `X/indices`, `X/indptr` as int64, the obs index, `target_gene`, `context`, the var index), read
with h5py / anndata so the CLI's layout (int32 vs int64 indptr, chunking, extra members) does not matter --
and refuses when it equals any ledger record's `content_sha256`. `record` stores it. ~1-2 min and ~13 GB of
scratch per package (`work/novelty-scratch/`, cleaned). `--skip-content` exists for a dry test only and
prints that it was skipped; the arming scripts never pass it.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import sys
from pathlib import Path

LEDGER = Path(__file__).resolve().parent.parent / "submissions" / "SUBMITTED_LEDGER.json"


def sha256(path: Path, bufsize: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def load(ledger: Path = LEDGER, *, must_exist: bool = True) -> dict:
    """The ledger, or -- when it is missing and the caller is a CHECK -- a refusal.

    A guard with no state must refuse, not pass: the disarmed output used to be byte-identical
    to the healthy one (242 section 2.3).
    """
    if not ledger.exists():
        if must_exist:
            print("=" * 78, file=sys.stderr)
            print(f"REFUSED: the submission ledger is MISSING at {ledger}", file=sys.stderr)
            print("  This check cannot run without it, and 'cannot check' is not 'novel'.",
                  file=sys.stderr)
            print("  It is tracked by git (.gitignore: !submissions/SUBMITTED_LEDGER.json);",
                  file=sys.stderr)
            print("  restore it with `git checkout -- submissions/SUBMITTED_LEDGER.json`.",
                  file=sys.stderr)
            print("=" * 78, file=sys.stderr)
            raise SystemExit(2)
        return {}
    data = json.loads(ledger.read_text())
    if not isinstance(data, dict):
        raise SystemExit(f"REFUSED: {ledger} is not a JSON object")
    return data


def content_sha256(pkg: Path, scratch: Path | None = None) -> tuple[str | None, str]:
    """sha256 of the PREDICTION inside the container, layout-independent. (digest or None, diagnostic)."""
    import shutil
    import subprocess
    import tarfile
    import tempfile
    scratch = scratch or (Path(__file__).resolve().parent.parent / "work" / "novelty-scratch")
    scratch.mkdir(parents=True, exist_ok=True)
    if shutil.which("zstd") is None:
        return None, "the zstd command-line tool is not installed (brew install zstd)"
    try:
        with tarfile.open(pkg) as tar:
            try:
                m = tar.getmember("pred.h5ad.zst")
            except KeyError:
                return None, "container has no pred.h5ad.zst member"
            need = m.size * 6 + (2 << 30)
            free = shutil.disk_usage(scratch).free
            if free < need:
                return None, f"scratch has {free / 2**30:.1f} GiB free, the decompressed prediction needs ~{need / 2**30:.0f}"
            fh = tar.extractfile(m)
            fd, name = tempfile.mkstemp(dir=scratch, suffix=".h5ad")
            tmp = Path(name)
            try:
                os.close(fd)  # inside the try: a failure here must still unlink the scratch file; an open
                              # descriptor would keep an unlinked 13 GB file allocated until exit
                with open(tmp, "wb") as out:
                    proc = subprocess.Popen(["zstd", "-d", "-q", "-c"], stdin=subprocess.PIPE, stdout=out)
                    assert proc.stdin is not None
                    while chunk := fh.read(1 << 24):
                        proc.stdin.write(chunk)
                    proc.stdin.close()
                    if proc.wait() != 0:
                        return None, "zstd -d failed on the prediction member"
                return _hash_h5ad(tmp), "content of pred.h5ad.zst (X, obs index, target_gene, context, var index)"
            finally:
                tmp.unlink(missing_ok=True)
    except Exception as exc:
        return None, f"could not read the container ({exc})"


def _hash_h5ad(path: Path) -> str:
    import h5py
    import numpy as np
    from anndata.io import read_elem
    h = hashlib.sha256()
    with h5py.File(path, "r") as f:
        def ds(name: str, dtype=None) -> None:
            d = f[name]
            n = int(d.shape[0])
            step = max(1, (1 << 26) // max(1, d.dtype.itemsize))
            h.update(name.encode())
            for i in range(0, n, step):
                a = d[i:i + step]
                if dtype is not None:
                    a = a.astype(dtype)
                h.update(np.ascontiguousarray(a).tobytes())
        ds("X/data")
        ds("X/indices", np.int64)
        ds("X/indptr", np.int64)
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        h.update(b"obs_index"); h.update("\n".join(map(str, obs.index)).encode())
        for col in ("target_gene", "context"):
            h.update(col.encode()); h.update("\n".join(map(str, obs[col].astype(str))).encode() if col in obs else b"<absent>")
        h.update(b"var_index"); h.update("\n".join(map(str, var.index)).encode())
    return h.hexdigest()


def container_meta(pkg: Path) -> dict | None:
    """The `meta.json` member vcc >= 0.2.0 writes into the .vcc (n_obs, n_vars, nnz, cli_version); None when absent."""
    import tarfile
    try:
        with tarfile.open(pkg) as tar:
            try:
                m = tar.getmember("meta.json")
            except KeyError:
                return None
            fh = tar.extractfile(m)
            return json.loads(fh.read().decode("utf-8")) if fh else None
    except Exception as exc:
        print(f"note: could not read the container's meta.json ({exc})", file=sys.stderr)
        return None


def prediction_digest(pkg: Path) -> tuple[str | None, str]:
    """The h5ad the package was prepped from: (sha256 or None, how it was found / why not).

    The sibling meta is the builder's own record of the artifact it wrote; when the container also
    carries a meta.json (vcc >= 0.2.0) its shape must agree with the sibling's `output` block, so a
    sibling that belongs to another prediction is caught instead of trusted by file name."""
    stem = pkg.with_suffix("") if pkg.suffix == ".vcc" else pkg
    meta = Path(str(stem) + ".h5ad.meta.json")
    h5ad = Path(str(stem) + ".h5ad")
    if meta.exists():
        try:
            rec = json.loads(meta.read_text())
        except Exception as exc:  # a corrupt meta is a refusal, never a pass
            return None, f"sibling meta {meta.name} unreadable ({exc})"
        out = rec.get("output") or {}
        d = out.get("sha256")
        if not (isinstance(d, str) and len(d) == 64):
            return None, f"sibling meta {meta.name} names no 64-hex output.sha256"
        cm = container_meta(pkg)
        if cm:
            # the builder's meta writes `n_rows` and `nnz`; the container writes `n_obs`, `n_vars`, `nnz`
            for ckey, skeys in (("n_obs", ("n_obs", "n_rows")), ("n_vars", ("n_vars",)), ("nnz", ("nnz",))):
                got = cm.get(ckey)
                want = next((out[k] for k in skeys if out.get(k) is not None), None)
                if want is not None and got is not None and int(want) != int(got):
                    return None, (f"container meta.json {ckey}={got} disagrees with sibling meta output.{skeys[-1]}={want}: "
                                  f"{meta.name} does not describe this package")
        return d, f"output.sha256 of {meta.name}" + (" (container shape agrees)" if cm else "")
    if h5ad.exists():
        print(f"hashing {h5ad} ({h5ad.stat().st_size / 1e9:.2f} GB)...", file=sys.stderr)
        return sha256(h5ad), f"sha256 of {h5ad.name}"
    return None, "no sibling .h5ad.meta.json or .h5ad next to the package"


def cli_version() -> str:
    import subprocess
    try:
        return subprocess.run(["vcc", "--version"], capture_output=True, text=True, timeout=30).stdout.strip().splitlines()[0]
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ledger", type=Path, default=LEDGER,
                    help=f"the ledger file (default {LEDGER}); for demonstrations only")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="exit 1 if this package was already submitted")
    c.add_argument("package")
    c.add_argument("--skip-content", action="store_true",
                   help="skip the content digest (dry tests only; the arming scripts never pass it)")
    c.add_argument("--allow-unidentified", action="store_true",
                   help="check on the container digest alone when the prediction cannot be identified "
                        "(a package built outside the builder); requires a note in docs/LOG.md")
    r = sub.add_parser("record", help="record a package as submitted")
    r.add_argument("package")
    r.add_argument("--entry-id", required=True)
    r.add_argument("--note", default="")
    sub.add_parser("list", help="show the ledger")
    a = ap.parse_args()

    # `record` may create the ledger; `check` and `list` must find it.
    ledger = load(a.ledger, must_exist=(a.cmd != "record"))
    if a.cmd == "list":
        for h, rec in sorted(ledger.items(), key=lambda kv: kv[1].get("date", "")):
            print(f"{h[:16]}  {rec.get('date','?'):26} {rec.get('entry_id','?'):22} "
                  f"{rec.get('package','?')}  {rec.get('note','')}")
        return 0

    pkg = Path(a.package)
    if not pkg.exists():
        print(f"ERROR: no such package: {pkg}", file=sys.stderr)
        return 2
    print(f"hashing {pkg} ({pkg.stat().st_size / 1e9:.2f} GB)...", file=sys.stderr)
    digest = sha256(pkg)

    if a.cmd == "check":
        prev = ledger.get(digest)
        how = "the .vcc digest"
        if prev is None:
            pdig, src = prediction_digest(pkg)
            if pdig is None:
                if not a.allow_unidentified:
                    print("=" * 78, file=sys.stderr)
                    print("REFUSED: the package's PREDICTION cannot be identified, and 'cannot check' is not 'novel'.",
                          file=sys.stderr)
                    print(f"  reason   {src}", file=sys.stderr)
                    print("  The .vcc digest alone is CLI-version-dependent (256 section 4). Put the builder's", file=sys.stderr)
                    print("  <stem>.h5ad.meta.json (or the <stem>.h5ad) next to the package, or pass", file=sys.stderr)
                    print("  --allow-unidentified with a note in docs/LOG.md for a package built outside the builder.",
                          file=sys.stderr)
                    print("=" * 78, file=sys.stderr)
                    return 1
                print(f"note: {src}; --allow-unidentified: checked on the .vcc digest only", file=sys.stderr)
            else:
                print(f"prediction identity: {pdig[:16]}... ({src})", file=sys.stderr)
                hits = [(h, r) for h, r in ledger.items() if r.get("h5ad_sha256") == pdig]
                if hits:
                    prev = hits[0][1]
                    how = f"the PREDICTION digest {pdig[:16]}... ({src}) -- a re-prep of a scored h5ad"
        if prev is None and not a.skip_content:
            cdig, csrc = content_sha256(pkg)
            if cdig is None:
                print("=" * 78, file=sys.stderr)
                print("REFUSED: the package's CONTENT could not be digested, and 'cannot check' is not 'novel'.", file=sys.stderr)
                print(f"  reason   {csrc}", file=sys.stderr)
                print("=" * 78, file=sys.stderr)
                return 1
            print(f"content identity: {cdig[:16]}... ({csrc})", file=sys.stderr)
            hits = [(h, r) for h, r in ledger.items() if r.get("content_sha256") == cdig]
            if hits:
                prev = hits[0][1]
                how = f"the CONTENT digest {cdig[:16]}... -- the same prediction inside a different container"
        elif prev is None:
            print("note: --skip-content: the content digest was NOT computed", file=sys.stderr)
        if prev is None:
            print(f"OK: novel package, sha256 {digest[:16]}...")
            return 0
        print("=" * 78, file=sys.stderr)
        print("REFUSED: this package has ALREADY BEEN SCORED. Do not spend a slot on it.",
              file=sys.stderr)
        print(f"  matched  {how}", file=sys.stderr)
        print(f"  sha256   {digest}", file=sys.stderr)
        print(f"  entry    {prev.get('entry_id')}   submitted {prev.get('date')}", file=sys.stderr)
        print(f"  package  {prev.get('package')}   {prev.get('note','')}", file=sys.stderr)
        print("", file=sys.stderr)
        print("  A slot must buy information. An already-scored package provably cannot --",
              file=sys.stderr)
        print("  we know its result before we spend the slot. If a bad arm is showing on the",
              file=sys.stderr)
        print("  board, LEAVE IT: rank is informational, only the last entry after Oct 22 is",
              file=sys.stderr)
        print("  prize-eligible. See AGENTS.md section 6.", file=sys.stderr)
        print("=" * 78, file=sys.stderr)
        return 1

    pdig, src = prediction_digest(pkg)
    cdig, csrc = content_sha256(pkg)
    ledger[digest] = {
        "package": str(pkg),
        "entry_id": a.entry_id,
        "date": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "note": a.note,
        "bytes": pkg.stat().st_size,
        "h5ad_sha256": pdig,
        "h5ad_sha256_source": src,
        "content_sha256": cdig,
        "content_sha256_source": csrc,
        "vcc_cli_version": cli_version(),
    }
    a.ledger.parent.mkdir(parents=True, exist_ok=True)
    a.ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(f"recorded {digest[:16]}... -> {a.entry_id}  (now `git add {a.ledger}`; it is tracked)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
