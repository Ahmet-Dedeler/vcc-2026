"""303: one guarded R5 access; importing and checking never reads expression.

Coordinator contract (v1, deliberately incompatible with the old loose release)
--------------------------------------------------------------------------
Only the coordinator may create work/RELEASE.json. Exact JSON object:

    {
      "schema_version": 1, "experiment": "303", "claim": "C7", "reserve": "R5",
      "dataset": "Nakatake2020_hESC_TF_screen",
      "authorization": "ONE_SHOT_R5_EXPRESSION",
      "issued_by": "coordinator", "coordinator_id": "<nonblank identity>",
      "issued_at": "<YYYY-MM-DDTHH:MM:SSZ>",
      "manifest": {
        "files": {"<repository-relative path>": "<64 lowercase hex>", ...},
        "analysis_code": ["<additional final analysis .py path>", ...],
        "feature_files": ["<frozen feature artifact path>", ...],
        "r5": {
          "n_rows": <positive integer>,
          "row_index_sha256": "...", "symbols_sha256": "...",
          "obs_id_sha256": "...", "ensembl_sha256": "...",
          "conditions_sha256": "...", "sign_sha256": "...",
          "var_id_sha256": "...", "var_symbols_sha256": "..."
        }
      },
      "test_reports": {"<repository-relative receipt.json>": "<sha256>", ...}
    }

``required_files()`` enumerates mandatory file paths. Every analysis_code and
feature_files entry must ALSO occur in files; both lists must be nonempty.
Extra dependencies can be listed in files. No globbing or optional hashes.
The files include existing metadata / KA / synthetic outputs; they are frozen
inputs, not recursively interpreted manifests. work/test303.json is HISTORICAL:
its old gate tests expect the previous release schema/report["live"]. A recorded
pass there does not validate this version. Run the legacy pure synthetic portions
separately if desired; do not report the incompatible legacy gate as rerun/pass.
Fresh guard + analysis receipts below are mandatory regardless of that old report.
Old hashes inside meta303.json or INPUT-MANIFEST.json are NOT current code pins.

Freeze code + inputs + final preregistration first, then the manifest, then run
tests, then issue the release. Test receipts are OUTSIDE manifest.files:

    {"schema_version": 1, "manifest_sha256": "<manifest_sha256(manifest)>",
     "suite": "access_guard303", "passed": true, "tests_run": 1,
     "failures": 0, "errors": 0}

Require one receipt each for access_guard303 and frozen_analysis303. Counts are
strict integers, not booleans; passed is exactly true. The latter suite is the
coordinator's frozen feature/results-analysis validation, not a C7 outcome.
``test_access_guard303.py --manifest <core-manifest.json>`` prints a guard
receipt to stdout; the coordinator saves it. Never put release/receipts/access
outputs in the core manifest: manifest -> receipts -> release is acyclic.

All hashes are fresh, streamed through opened regular files with fstat/path
identity checks. The source stays open through HDF5 validation and extraction.
Live obs/var and cohort/mask identity are checked, never trusted from a sidecar.
run_access() is the ONLY expression-loading entrypoint, always validates again,
then O_EXCL-creates and fsyncs a marker AND its directory before any value slice.
Any marker, including empty/interrupted, permanently refuses another access.
Products also use O_EXCL and are never overwritten. There is no resume flag.

Limits: this enforces a cooperative local protocol, not OS isolation or a
cryptographic signature of the human coordinator. Someone able to modify this
code, forge receipts, monkeypatch Python, or delete markers can defeat it.
The coordinator must preserve the marker and audit the dependency list/receipts.
Hashing bytes is not expression interpretation; --dry-run without a release
does not even open the source. R5 remains locked until a separate release.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import resource
import stat
import sys
import time
from contextlib import ExitStack, contextmanager
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib303 as L  # noqa: E402

RELEASE = L.WORK / "RELEASE.json"
MARKER = L.WORK / "ACCESS-MARKER.json"
PRODUCT = L.WORK / "nak303.npz"
PREREG_FINAL = L.HERE / "PREREGISTRATION.md"
CODE_FILES = ("lib303.py", "stats303.py", "access303.py", "meta303.py",
              "ka303.py", "test303.py", "test_access_guard303.py")
_LOADED_HASHES = {p.name: L.sha256_file(p) for p in
                  (Path(__file__).resolve(), Path(L.__file__).resolve())}
_HEX = re.compile(r"[0-9a-f]{64}\Z")
_R5_HASHES = ("row_index_sha256", "symbols_sha256", "obs_id_sha256",
              "ensembl_sha256", "conditions_sha256", "sign_sha256",
              "var_id_sha256", "var_symbols_sha256")
REQUIRED_RELEASE_KEYS = (
    "schema_version", "experiment", "claim", "reserve", "dataset",
    "authorization", "issued_by", "coordinator_id", "issued_at", "manifest",
    "test_reports")


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise L.Refuse(message)


def _keys(obj, keys, label: str) -> None:
    _require(type(obj) is dict and set(obj) == set(keys),
             f"{label}: exact keys required: {sorted(keys)}")


def _hash(value, label: str) -> None:
    _require(type(value) is str and _HEX.fullmatch(value) is not None,
             f"{label}: expected 64 lowercase hexadecimal characters")


def _json(raw: bytes) -> dict:
    def pairs(items):
        result = {}
        for key, value in items:
            _require(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def invalid(value):
        raise L.Refuse(f"nonfinite JSON constant: {value}")

    try:
        result = json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid)
    except (ValueError, UnicodeError) as exc:
        raise L.Refuse(f"invalid JSON: {exc}") from exc
    _require(type(result) is dict, "JSON document must be an object")
    return result


def manifest_sha256(manifest: dict) -> str:
    """Stable core digest; receipts and release must never enter this object."""
    return hashlib.sha256(json.dumps(
        manifest, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("utf-8")).hexdigest()


def _relative(path: Path) -> str:
    return path.relative_to(L.ROOT).as_posix()


def required_files() -> tuple[str, ...]:
    """Explicit live dependencies; additions belong in the coordinator manifest."""
    paths = [
        L.DB, L.CANON, PREREG_FINAL, L.HELD279, L.PERSRC290, L.ZNF279,
        L.SYM302, L.SUMS302,
        *(L.HERE / name for name in CODE_FILES),
        *(L.H1BANK / name for name in
          ("delta.npy", "perts.txt", "genes.txt", "measured_mask.npy",
           "b_ctrl.npy", "meta.json", "rows.tsv")),
        *(L.WORK / name for name in
          ("COHORT.tsv", "AXIS.txt", "R5-ROWS.tsv", "meta303.json",
           "ka303.json", "test303.json")),
    ]
    return tuple(sorted(_relative(p) for p in paths))


def _path(name: str) -> Path:
    _require(type(name) is str and bool(name) and "\\" not in name
             and "\x00" not in name, "invalid manifest path")
    p = PurePosixPath(name)
    _require(not p.is_absolute() and p.as_posix() == name
             and not any(x in (".", "..") for x in p.parts),
             f"path must be canonical repository-relative: {name}")
    return L.ROOT / name


def _hash_map(obj, label: str) -> None:
    _require(type(obj) is dict and bool(obj), f"{label}: nonempty hash map required")
    for name, digest in obj.items():
        _path(name)
        _hash(digest, f"{label}[{name}]")


def _schema(rel: dict) -> None:
    _keys(rel, REQUIRED_RELEASE_KEYS, "release")
    _require(type(rel["schema_version"]) is int and rel["schema_version"] == 1,
             "release schema_version must be integer 1")
    for key, value in dict(experiment="303", claim="C7", reserve="R5",
                           dataset=L.RESERVE_DATASET, issued_by="coordinator",
                           authorization="ONE_SHOT_R5_EXPRESSION").items():
        _require(type(rel[key]) is str and rel[key] == value,
                 f"release {key} must be {value!r}")
    who = rel["coordinator_id"]
    _require(type(who) is str and bool(who.strip()) and who == who.strip()
             and all(ord(c) >= 32 for c in who), "invalid coordinator identity")
    issued = rel["issued_at"]
    _require(type(issued) is str and
             re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", issued) is not None,
             "issued_at must be a UTC timestamp YYYY-MM-DDTHH:MM:SSZ")
    try:
        dt = datetime.strptime(issued, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except ValueError as exc:
        raise L.Refuse("invalid issued_at calendar date") from exc
    _require(dt <= datetime.now(timezone.utc), "release issued_at is in the future")
    man = rel["manifest"]
    _keys(man, ("files", "analysis_code", "feature_files", "r5"), "manifest")
    _hash_map(man["files"], "manifest.files")
    _require(set(required_files()) <= set(man["files"]),
             f"manifest missing mandatory files: {sorted(set(required_files()) - set(man['files']))}")
    for key in ("analysis_code", "feature_files"):
        names = man[key]
        _require(type(names) is list and bool(names) and
                 all(type(n) is str for n in names), f"{key}: nonempty path list required")
        _require(names == sorted(set(names)), f"{key}: sorted unique paths required")
        _require(set(names) <= set(man["files"]), f"{key}: unbound file")
    _require(all(n.endswith(".py") for n in man["analysis_code"]),
             "analysis_code must enumerate Python source files")
    _require(not set(man["analysis_code"]) & set(man["feature_files"]),
             "code and feature artifacts must be separate")
    r5 = man["r5"]
    _keys(r5, ("n_rows", *_R5_HASHES), "manifest.r5")
    _require(type(r5["n_rows"]) is int and r5["n_rows"] > 0,
             "r5.n_rows must be a positive integer")
    for key in _R5_HASHES:
        _hash(r5[key], f"r5.{key}")
    _hash_map(rel["test_reports"], "test_reports")
    forbidden = {_relative(p) for p in
                 (RELEASE, MARKER, PRODUCT, L.WORK / "access303.json")}
    _require(not set(man["files"]) & (set(rel["test_reports"]) | forbidden),
             "manifest cycle: release, receipts and access outputs are not core inputs")
    _require(not set(rel["test_reports"]) & forbidden,
             "test receipts cannot be release/access outputs")


def _identity(st) -> tuple:
    return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns


class _PinnedFile:
    """Unbuffered opened regular file; never trusts a pathname after hashing."""

    def __init__(self, path: Path):
        self.path = path
        _require(path.resolve(strict=True) == path, f"symlink path refused: {path}")
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        self.fh = os.fdopen(fd, "rb", buffering=0)
        try:
            st = os.fstat(fd)
            _require(stat.S_ISREG(st.st_mode), f"not a regular file: {path}")
            self.identity = _identity(st)
            self.unchanged()
        except BaseException:
            self.fh.close()
            raise

    def unchanged(self) -> None:
        _require(self.path.resolve(strict=True) == self.path and
                 _identity(os.fstat(self.fh.fileno())) == self.identity and
                 _identity(self.path.stat(follow_symlinks=False)) == self.identity,
                 f"file changed since opening: {self.path}")

    def digest(self) -> str:
        self.unchanged()
        self.fh.seek(0)
        h = hashlib.sha256()
        while chunk := self.fh.read(1 << 22):
            h.update(chunk)
        self.unchanged()
        self.fh.seek(0)
        return h.hexdigest()

    def read(self) -> bytes:
        self.unchanged()
        self.fh.seek(0)
        result = self.fh.read()
        self.unchanged()
        return result


def _dataset(f: h5py.File, name: str) -> h5py.Dataset:
    # An external link, external storage, or VDS would escape the source hash.
    group = f
    parts = name.split("/")
    for part in parts:
        _require(isinstance(group.get(part, getlink=True), h5py.HardLink),
                 f"non-local or missing HDF5 object: {name}")
        group = group[part]
    _require(isinstance(group, h5py.Dataset), f"not an HDF5 dataset: {name}")
    _require(not group.is_virtual and not group.external,
             f"external/virtual storage refused: {name}")
    return group


def _strings(d: h5py.Dataset) -> np.ndarray:
    _require(d.ndim == 1 and d.dtype.kind in "OSU", f"invalid strings: {d.name}")
    raw = d[:]
    values = [v.decode("utf-8") if isinstance(v, bytes) else v for v in raw]
    _require(all(isinstance(v, str) and v and "\n" not in v and "\r" not in v
                 for v in values), f"invalid string metadata: {d.name}")
    return np.asarray(values, dtype=str)


def _metadata(f: h5py.File) -> L.DbMeta:
    """Read only obs/var; X/control are inspected for shape/storage, not values."""
    def cat(name):
        categories = _strings(_dataset(f, name + "/categories"))
        d = _dataset(f, name + "/codes")
        _require(d.ndim == 1 and d.dtype.kind in "iu", f"invalid categorical codes: {name}")
        codes = d[:]
        _require(np.all((codes >= 0) & (codes < len(categories))),
                 f"invalid categorical codes: {name}")
        return categories[codes]

    datasets = [_dataset(f, name) for name in
                ("X", "layers/control", "layers/measure_mask")]
    shape = datasets[0].shape
    _require(len(shape) == 2 and all(s > 0 for s in shape) and
             all(d.shape == shape and d.dtype.kind in "fiu" for d in datasets),
             "X/control/mask shape or dtype mismatch")
    _require(datasets[2].id != datasets[0].id and datasets[2].id != datasets[1].id,
             "measure_mask aliases a protected expression dataset")
    sign_d = _dataset(f, "obs/perturbation_sign")
    _require(sign_d.ndim == 1 and sign_d.dtype.kind in "iu",
             "perturbation_sign must have integer dtype")
    sign = sign_d[:]
    m = L.DbMeta(
        *shape, cat("obs/dataset"), _strings(_dataset(f, "obs/id")),
        cat("obs/perturbed_gene_name"), cat("obs/perturbed_gene"),
        cat("obs/condition"), sign.astype(np.int64),
        _strings(_dataset(f, "var/human_id")), cat("var/gene_name"))
    _require(all(a.shape == (shape[0],) for a in
                 (m.dataset, m.obs_id, m.pert_symbol, m.pert_gene, m.condition, m.sign))
             and all(a.shape == (shape[1],) for a in (m.var_id, m.var_symbol)),
             "obs/var lengths do not match source")
    return m


def _r5_identity(m: L.DbMeta) -> tuple[dict, np.ndarray]:
    rows = np.flatnonzero(m.dataset == L.RESERVE_DATASET).astype(np.int64)
    _require(rows.size > 0, "no live R5 rows")
    _require(len(set(m.pert_symbol[rows])) == rows.size and
             len(set(m.obs_id[rows])) == rows.size, "duplicate R5 row identities")
    _require(np.all(m.sign[rows] == 1), "R5 sign is not overexpression (+1)")
    out = dict(n_rows=int(rows.size), row_index_sha256=L.sha256_array(rows))
    for key, arr in (("symbols", m.pert_symbol[rows]), ("obs_id", m.obs_id[rows]),
                     ("ensembl", m.pert_gene[rows]), ("conditions", m.condition[rows]),
                     ("var_id", m.var_id), ("var_symbols", m.var_symbol)):
        out[key + "_sha256"] = L.sha256_strings(arr)
    out["sign_sha256"] = L.sha256_array(m.sign[rows])
    return out, rows


def _lines(raw: bytes, label: str) -> list[str]:
    lines = raw.decode("utf-8").splitlines()
    _require(bool(lines) and all(x and x == x.strip() for x in lines)
             and len(set(lines)) == len(lines), f"{label}: empty/duplicate/invalid lines")
    return lines


def _primary302(raw: bytes) -> set[str]:
    """Parse the actual saved 302 schema; no reserve or HDF5 I/O."""
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
    _require(reader.fieldnames is not None
             and len(set(reader.fieldnames)) == len(reader.fieldnames)
             and {"symbol", "cohort"} <= set(reader.fieldnames),
             "symmetry302.csv requires symbol/cohort columns (not tf)")
    rows = list(reader)
    _require(bool(rows) and all(type(r["symbol"]) is str and bool(r["symbol"])
                               and type(r["cohort"]) is str and bool(r["cohort"])
                               and None not in r for r in rows)
             and len({r["symbol"] for r in rows}) == len(rows),
             "symmetry302.csv has missing/duplicate row identities")
    return {r["symbol"] for r in rows if r["cohort"] == "primary"}


def _preflight(f, pins, m, live_r5, rows):
    """Validate frozen rows and measured axis with metadata and masks ONLY."""
    def raw(path):
        return pins[_relative(path)].read()

    canon = _lines(raw(L.CANON), "canonical genes")
    _require(canon.pop(0) == "gene_name", "canonical gene header must be gene_name")
    axis = _lines(raw(L.WORK / "AXIS.txt"), "axis")
    _require(set(axis) <= set(canon), "axis outside canonical genes")
    perts = _lines(raw(L.H1BANK / "perts.txt"), "KD perts")
    _require(_lines(raw(L.H1BANK / "genes.txt"), "KD genes") == canon,
             "KD gene order differs from canonical genes")
    row_of = {str(m.pert_symbol[i]): int(i) for i in rows}
    bank_of = {s: i for i, s in enumerate(perts)}
    cohort = sorted(set(row_of) & set(perts))
    _require(bool(cohort), "empty live OE/KD cohort")
    table = csv.DictReader(io.StringIO(raw(L.WORK / "COHORT.tsv").decode()), delimiter="\t")
    fields = ("order", "symbol", "r5_source_row", "bank_row", "own_column_on_axis",
              "in_302_cohort", "in_db_joung")
    _require(table.fieldnames == list(fields), "COHORT.tsv header mismatch")
    records = list(table)
    _require(len(records) == len(cohort), "COHORT.tsv differs from live cohort")
    primary302 = _primary302(raw(L.SYM302))
    joung = set(m.pert_symbol[m.dataset == L.JOUNG_DATASET])
    for k, (s, rec) in enumerate(zip(cohort, records)):
        expected = (str(k), s, str(row_of[s]), str(bank_of[s]), str(int(s in axis)),
                    str(int(s in primary302)), str(int(s in joung)))
        _require(tuple(rec.get(n) for n in fields) == expected and None not in rec,
                 f"COHORT.tsv live identity mismatch at order {k}")
    row_table = csv.DictReader(
        io.StringIO(raw(L.WORK / "R5-ROWS.tsv").decode()), delimiter="\t")
    row_records = list(row_table)
    _require(len(row_records) == len(rows), "R5-ROWS.tsv row count mismatch")
    for k, (i, rec) in enumerate(zip(rows, row_records)):
        expected = dict(order=str(k), source_row_index=str(i), obs_id=m.obs_id[i],
                        symbol=m.pert_symbol[i], ensembl=m.pert_gene[i], sign=str(m.sign[i]),
                        in_bank_cohort=str(int(m.pert_symbol[i] in cohort)))
        _require(all(rec.get(n) == v for n, v in expected.items()),
                 f"R5-ROWS.tsv live identity mismatch at order {k}")
    meta = _json(raw(L.WORK / "meta303.json"))
    for key in ("n_rows", "row_index_sha256", "symbols_sha256", "obs_id_sha256",
                "ensembl_sha256"):
        _require(type(meta["R5"][key]) is type(live_r5[key])
                 and meta["R5"][key] == live_r5[key], f"stale metadata R5.{key}")
    _require(meta["file"]["n_obs"] == m.n_obs and meta["file"]["n_var"] == m.n_var,
             "stale metadata dimensions")
    cohort_meta = meta["cohorts"]["nakatake_x_bank"]
    _require(cohort_meta["n"] == len(cohort)
             and cohort_meta["sha256"] == L.sha256_strings(cohort),
             "stale metadata cohort")
    _require(meta["axis"]["n"] == len(axis)
             and meta["axis"]["sha256"] == L.sha256_strings(axis), "stale metadata axis")
    proj = L.build_projection(m.var_symbol, canon)
    mm_pin = pins[_relative(L.H1BANK / "measured_mask.npy")]
    mm_pin.fh.seek(0)
    kd_mask = np.load(mm_pin.fh, allow_pickle=False)
    _require(kd_mask.dtype == np.bool_ and kd_mask.shape == (len(perts), len(canon)),
             "KD measured_mask shape/dtype mismatch")
    support = kd_mask[[bank_of[s] for s in cohort]].all(axis=0)
    del kd_mask
    selected = np.sort([row_of[s] for s in cohort])
    for lo in range(0, len(selected), L.BLOCK):
        mm = _dataset(f, "layers/measure_mask")[selected[lo:lo + L.BLOCK], :]
        _require(np.all((mm == 0) | (mm == 1)), "non-binary live R5 mask")
        support &= L.collapse_mask_only(mm.astype(bool), proj).all(axis=0)
    _require(axis == [g for g, yes in zip(canon, support) if yes],
             "AXIS.txt differs from live common measured axis")
    ka = _json(raw(L.WORK / "ka303.json"))
    _require(ka["summary"]["all_pass"] is True, "failed known-answer summary")
    for name in ("KA-C", "KA-L1", "KA-L2", "KA-L3", "KA-L4", "KA-L5", "KA-W"):
        _require(ka[name]["pass_"] is True, f"failed known answer {name}")
    tests = _json(raw(L.WORK / "test303.json"))
    _require(type(tests["n"]) is int and tests["n"] > 0
             and type(tests["n_fail"]) is int and tests["n_fail"] == 0
             and type(tests["results"]) is list
             and len(tests["results"]) == tests["n"]
             and all(r["pass_"] is True for r in tests["results"]),
             "failed/incomplete frozen synthetic tests")
    positions = {g: i for i, g in enumerate(canon)}
    return proj, axis, np.asarray([positions[g] for g in axis], np.int64)


def _unspent() -> None:
    # lexists detects even a dangling symlink or an interrupted empty marker.
    _require(not os.path.lexists(MARKER),
             f"ACCESS-MARKER present ({MARKER}): R5 SPENT; second access refused")
    _require(not os.path.lexists(PRODUCT) and
             not os.path.lexists(L.WORK / "access303.json"),
             "access output already exists without marker; refuse recovery/re-read")


@contextmanager
def _validated_source():
    """Keep all validated descriptors alive. This helper never reads expression."""
    try:
        _unspent()
        _require(PREREG_FINAL.is_file(), "final PREREGISTRATION.md does not exist")
        _require(RELEASE.is_file(), f"no coordinator release at {RELEASE}")
        with ExitStack() as stack:
            pins = {}

            def pin(path):
                p = _PinnedFile(path)
                stack.callback(p.fh.close)
                return p

            release_pin = pin(RELEASE)
            rel = _json(release_pin.read())
            _schema(rel)
            for name, expected in rel["manifest"]["files"].items():
                p = pin(_path(name))
                pins[name] = p
                _require(p.digest() == expected, f"fresh file digest mismatch: {name}")
            _require(bool(pins[_relative(PREREG_FINAL)].read().strip()),
                     "final PREREGISTRATION.md is empty")
            for name, loaded_hash in _LOADED_HASHES.items():
                _require(rel["manifest"]["files"][_relative(L.HERE / name)] == loaded_hash,
                         f"loaded {name} differs from frozen file; start a fresh process")
            suites = set()
            core_hash = manifest_sha256(rel["manifest"])
            for name, expected in rel["test_reports"].items():
                p = pin(_path(name))
                pins[name] = p
                _require(p.digest() == expected, f"test receipt hash mismatch: {name}")
                receipt = _json(p.read())
                _keys(receipt, ("schema_version", "manifest_sha256", "suite", "passed",
                                "tests_run", "failures", "errors"), "test receipt")
                _require(type(receipt["schema_version"]) is int
                         and receipt["schema_version"] == 1, "bad receipt version")
                _hash(receipt["manifest_sha256"], "receipt.manifest_sha256")
                _require(receipt["manifest_sha256"] == core_hash, "stale test receipt")
                _require(type(receipt["suite"]) is str
                         and receipt["suite"] in ("access_guard303", "frozen_analysis303")
                         and receipt["suite"] not in suites, "invalid/duplicate receipt suite")
                _require(receipt["passed"] is True
                         and type(receipt["tests_run"]) is int and receipt["tests_run"] > 0
                         and type(receipt["failures"]) is int and receipt["failures"] == 0
                         and type(receipt["errors"]) is int and receipt["errors"] == 0,
                         "failed/incomplete test receipt")
                suites.add(receipt["suite"])
            _require(suites == {"access_guard303", "frozen_analysis303"},
                     "guard and frozen-analysis test receipts both required")
            pins[_relative(RELEASE)] = release_pin
            source = pins[_relative(L.DB)]
            # Close HDF5 BEFORE its underlying file object (ExitStack is LIFO).
            f = stack.enter_context(h5py.File(source.fh, "r"))
            m = _metadata(f)
            r5, rows = _r5_identity(m)
            _require(r5 == rel["manifest"]["r5"], "live HDF5 R5 identity mismatch")
            proj, axis, axis_cols = _preflight(f, pins, m, r5, rows)
            for p in pins.values():
                p.unchanged()
            _unspent()
            report = dict(live=dict(source_sha256=rel["manifest"]["files"][_relative(L.DB)],
                                    r5=r5, manifest_sha256=core_hash),
                          release=rel, release_sha256=release_pin.digest(),
                          marker_absent=True)
            yield f, pins, m, rows, proj, axis, axis_cols, report
    except (OSError, ValueError, KeyError, TypeError, OverflowError, AttributeError, IndexError) as exc:
        raise L.Refuse(f"access validation/I/O failed: {exc}") from exc


def check_gate(strict_release: bool = True) -> dict:
    """Read-only validation; its return value is NOT an access capability."""
    _require(type(strict_release) is bool, "strict_release must be bool")
    _unspent()
    if not strict_release and not os.path.lexists(RELEASE):
        return dict(marker_absent=True, release=None, live={},
                    locked=True, reason="no coordinator release; no source opened")
    with _validated_source() as validated:
        return validated[-1]


def _sync_directory(directory: Path) -> None:
    fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _exclusive_json(path: Path, payload: dict) -> None:
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, sort_keys=True, indent=1, allow_nan=False)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    _sync_directory(path.parent)


def run_access() -> dict:
    """Validate, durably consume, then read exactly once. No ungated variant."""
    t0 = time.monotonic()
    with _validated_source() as (f, pins, m, rows, proj, axis, axis_cols, report):
        def unchanged():
            for p in pins.values():
                p.unchanged()

        unchanged()
        _exclusive_json(MARKER, dict(
            experiment="303", claim="C7", reserve="R5", dataset=L.RESERVE_DATASET,
            source=_relative(L.DB), started_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            pid=os.getpid(), **report,
            note="R5 SPENT, even if interrupted or no product. Never resume source reads."))
        # Marker fsync and directory fsync have both completed before this point.
        unchanged()
        delta_blocks, mask_blocks, csum, psum = [], [], [], []
        ctrl_first, ctrl_identical = None, True
        for lo in range(0, rows.size, L.BLOCK):
            sel = rows[lo:lo + L.BLOCK]
            unchanged()
            x = np.asarray(_dataset(f, "X")[sel, :], np.float64)
            unchanged()
            c = np.asarray(_dataset(f, "layers/control")[sel, :], np.float64)
            unchanged()
            mm = np.asarray(_dataset(f, "layers/measure_mask")[sel, :], np.float64)
            unchanged()
            if ctrl_first is None:
                ctrl_first = c[0].copy()
            ctrl_identical &= bool(np.array_equal(c, np.broadcast_to(ctrl_first, c.shape)))
            d, msk, _ = L.bank_transform(
                x, c, mm, proj, [m.obs_id[i] for i in sel], strict_sums=False)
            delta_blocks.append(d[:, axis_cols].astype(np.float32))
            mask_blocks.append(msk[:, axis_cols])
            csum.append(np.expm1(np.maximum(c, 0.0)).sum(axis=1))
            psum.append(np.expm1(np.maximum(c + x, 0.0)).sum(axis=1))
        unchanged()
        # A second full digest catches changes during the pass before publishing.
        _require(pins[_relative(L.DB)].digest() == report["live"]["source_sha256"],
                 "source changed during access")
        delta, mask = np.concatenate(delta_blocks), np.concatenate(mask_blocks)
        out = dict(rows=rows, symbols=m.pert_symbol[rows].astype("U"),
                   obs_id=m.obs_id[rows].astype("U"), axis=np.asarray(axis, dtype="U"),
                   delta=delta, mask=mask, b_ctrl=L.control_profile(ctrl_first, proj)[axis_cols],
                   control_cpm_sum=np.concatenate(csum), perturbed_cpm_sum=np.concatenate(psum),
                   control_rows_identical=np.asarray([ctrl_identical]))
        _require(all(np.isfinite(out[k]).all() for k in
                     ("delta", "b_ctrl", "control_cpm_sum", "perturbed_cpm_sum")),
                 "nonfinite extracted product")
        unchanged()
        fd = os.open(PRODUCT, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, "wb") as fh:
            np.savez_compressed(fh, **out)
            fh.flush()
            os.fsync(fh.fileno())
        _sync_directory(PRODUCT.parent)
        unchanged()
        rep = dict(
            n_rows=int(rows.size), n_axis=len(axis), control_rows_identical=ctrl_identical,
            control_cpm_sum=[float(out["control_cpm_sum"].min()),
                             float(out["control_cpm_sum"].max())],
            perturbed_cpm_sum=[float(out["perturbed_cpm_sum"].min()),
                               float(out["perturbed_cpm_sum"].max())],
            cpm_gate_pass=bool(
                np.all(np.abs(out["control_cpm_sum"] - L.SOURCE_TARGET_SUM) <= L.TARGET_SUM_ATOL)
                and np.all(np.abs(out["perturbed_cpm_sum"] - L.SOURCE_TARGET_SUM) <= L.TARGET_SUM_ATOL)),
            mask_row_constant=bool(np.all(mask == mask[0])), all_finite=True,
            product_sha256=L.sha256_file(PRODUCT), delta_sha256=L.sha256_array(delta),
            mask_sha256=L.sha256_array(mask), manifest_sha256=report["live"]["manifest_sha256"],
            release_sha256=report["release_sha256"], seconds=round(time.monotonic() - t0, 1),
            peak_rss_gib=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, 2))
        _exclusive_json(L.WORK / "access303.json", rep)
        return rep


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="validate only; never read reserve expression or create a marker")
    args = ap.parse_args()
    try:
        rep = check_gate(strict_release=False) if args.dry_run else run_access()
        print(json.dumps(rep, indent=1))
        return 0
    except L.Refuse as exc:
        print(f"REFUSED: {exc}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
