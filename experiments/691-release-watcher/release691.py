#!/usr/bin/env python3
"""691 — helpers for the Oct 22 release watcher (watch691.sh). Pure data handling: no build, no upload, no network.

  diff <baseline.json> <now.json>
      Prints one dataset id per line that is NEW or CHANGED (filename or description differs) AND available, compared
      with the baseline `vcc --json datasets list`. Exit 0 always; an empty print means nothing new.
  pick <now.json> <id> [<id> ...]
      Chooses the final bundle among the NEW ids (review 720 F5): one id -> it; several -> the single one whose id /
      filename / description mentions final, test or D/E/F; otherwise exit 2 (ambiguous: a person decides).
  unpack <zip> <dest_bundle_dir> <record.json>
      Treats the zip as untrusted data (AGENTS: downloaded archives are data). Refuses absolute paths, `..`, symlinks,
      anything but .h5ad / .csv / .json / .txt / .md, more than 64 members, or > 20 GB uncompressed. Extracts into a fresh
      sibling dir, finds the ONE directory holding manifest.json, and copies that bundle to <dest_bundle_dir> (which must
      not exist). Records every file's size and sha256.
  overrides <preflight.json>
      The DECISION-CARD §1 overrides that are pre-decided, for AUTO_R0 only, each behind its own arming flag (default off):
      AUTO_ALLOW_DEPTH691=1 -> DEF_DAY_ALLOW_DEPTH=1 iff the only failures are 641 depth rows and the three median ratios
      agree within 2 % and lie in [0.6, 1.8] (a uniform release depth; one context off alone stays a STOP);
      AUTO_ALLOW_WIDE691=1 -> DEF_DAY_ALLOW_WIDE_LIBS=1 iff the only other failures are 554 spread rows, every spread
      <= 3.5 (the card's own bound; 554's postbuild dose guard still checks the build). Prints the env assignments, exit 0;
      anything else prints why and exits 1.
  verify <bundle_dir> <record.json>
      Resume after a restart: exit 0 iff every file recorded by `unpack` is in <bundle_dir> with its size and sha256.
  summarize <bundle> <light_work_dir> <out.json>
      After the light stages: contexts, the panel against A/B/C's (SAME / NEW / OVERLAP n), the gene axis against A/B/C's
      (must be identical), per-context library median and q75/q25, identity calls, census gate. Prints the branch the day
      takes: SAME-PANEL (pre-place 453's operator 66610be7) or NEW-PANEL.

usage: uv run python -I experiments/691-release-watcher/release691.py <cmd> ...
"""
from __future__ import annotations

import hashlib
import json
import shutil
import stat
import sys
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
ABC = ROOT / "data/challenge"
ALLOWED_EXT = {".h5ad", ".csv", ".json", ".txt", ".md"}
MAX_MEMBERS, MAX_BYTES = 64, 20 * 2**30


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


def cmd_diff(base: Path, now: Path) -> int:
    b = {d["id"]: d for d in json.loads(base.read_text()).get("datasets", [])}
    n = json.loads(now.read_text()).get("datasets", [])
    for d in n:
        old = b.get(d["id"])
        changed = old is None or old.get("filename") != d.get("filename") or old.get("description") != d.get("description") \
            or (not old.get("available") and d.get("available"))
        if changed and d.get("available"):
            print(d["id"])
    return 0


FINAL_RE = __import__("re").compile(r"final|test|\bdef\b|d/e/f|\bd, ?e,? (and )?f\b", __import__("re").I)


def cmd_pick(now: Path, ids: list[str]) -> int:
    ds = {d["id"]: d for d in json.loads(now.read_text()).get("datasets", [])}
    ids = [i for i in ids if i in ds]
    if len(ids) == 1:
        d = ds[ids[0]]
        txt = f"{d['id']} {d.get('filename', '')} {d.get('description', '')}"
        if not FINAL_RE.search(txt):
            print(f"WARNING: the only new dataset does not say final/test/D-E-F: {txt[:200]!r}", file=sys.stderr)
        print(ids[0])
        return 0
    hits = [i for i in ids if FINAL_RE.search(f"{i} {ds[i].get('filename', '')} {ds[i].get('description', '')}")]
    if len(hits) == 1:
        print(hits[0])
        return 0
    print(f"AMBIGUOUS: {len(ids)} new datasets, {len(hits)} look final: "
          + "; ".join(f"{i} ({ds[i].get('filename')}: {str(ds[i].get('description'))[:80]})" for i in ids), file=sys.stderr)
    return 2


def cmd_unpack(zpath: Path, dest: Path, record: Path) -> int:
    if dest.exists():
        print(f"REFUSED: {dest} exists (a fresh bundle dir only)", file=sys.stderr)
        return 1
    with zipfile.ZipFile(zpath) as z:
        infos = z.infolist()
        if len(infos) > MAX_MEMBERS:
            print(f"REFUSED: {len(infos)} members > {MAX_MEMBERS}", file=sys.stderr)
            return 1
        total = 0
        for i in infos:
            pp = PurePosixPath(i.filename)
            mode = (i.external_attr >> 16) & 0o170000
            if pp.is_absolute() or ".." in pp.parts or i.filename.startswith(("/", "\\")) or ":" in pp.parts[0]:
                print(f"REFUSED: unsafe member path {i.filename!r}", file=sys.stderr)
                return 1
            if mode == stat.S_IFLNK:
                print(f"REFUSED: symlink member {i.filename!r}", file=sys.stderr)
                return 1
            if not i.is_dir() and pp.name and not pp.name.startswith(".") and pp.suffix.lower() not in ALLOWED_EXT \
                    and not i.filename.startswith("__MACOSX/"):
                print(f"REFUSED: unexpected member type {i.filename!r}", file=sys.stderr)
                return 1
            total += i.file_size
        if total > MAX_BYTES:
            print(f"REFUSED: {total} bytes uncompressed > {MAX_BYTES}", file=sys.stderr)
            return 1
        stage = dest.parent / (dest.name + ".unzipped")
        if stage.exists():
            print(f"REFUSED: {stage} exists", file=sys.stderr)
            return 1
        stage.mkdir(parents=True)
        z.extractall(stage)
    manifests = [m for m in stage.rglob("manifest.json") if "__MACOSX" not in m.parts]
    if len(manifests) != 1:
        print(f"REFUSED: {len(manifests)} manifest.json files in the archive (want 1): {manifests}", file=sys.stderr)
        return 1
    src = manifests[0].parent
    shutil.copytree(src, dest)
    files = {p.name: {"bytes": p.stat().st_size, "sha256": sha(p)} for p in sorted(dest.iterdir()) if p.is_file()}
    rec = {"zip": str(zpath), "zip_bytes": zpath.stat().st_size, "zip_sha256": sha(zpath), "bundle": str(dest),
           "staged_from": str(src), "files": files, "manifest": json.loads((dest / "manifest.json").read_text())}
    record.write_text(json.dumps(rec, indent=1) + "\n")
    print(f"unpacked {len(files)} files to {dest}; contexts {rec['manifest'].get('contexts')}")
    return 0


def cmd_overrides(pf: Path) -> int:
    import os
    d = json.loads(pf.read_text())
    failing = [c for c in d.get("checks", []) if not c.get("pass")]
    if not failing:
        print("NONE: preflight has no failing row", file=sys.stderr)
        return 1
    depth = [c for c in failing if "control library median within" in c.get("check", "")]
    wide = [c for c in failing if "control library q75/q25" in c.get("check", "")]
    other = [c["check"][:100] for c in failing if c not in depth and c not in wide]
    if other:
        print(f"NONE: failures no override covers: {other}", file=sys.stderr)
        return 1
    cf = d.get("control_files", {})
    env = []
    if depth:
        if os.environ.get("AUTO_ALLOW_DEPTH691") != "1":
            print("NONE: depth rows failed and AUTO_ALLOW_DEPTH691 is not armed", file=sys.stderr)
            return 1
        r = [float(v["library_median_ratio"]) for v in cf.values()]
        if not r or max(r) / min(r) > 1.02 or min(r) < 0.6 or max(r) > 1.8:
            print(f"NONE: depth not uniform or out of range: ratios {r}", file=sys.stderr)
            return 1
        env.append("DEF_DAY_ALLOW_DEPTH=1")
    if wide:
        if os.environ.get("AUTO_ALLOW_WIDE691") != "1":
            print("NONE: spread rows failed and AUTO_ALLOW_WIDE691 is not armed", file=sys.stderr)
            return 1
        sp = [float(v["library_spread_q75_q25"]) for v in cf.values()]
        if max(sp) > 3.5:
            print(f"NONE: a spread above 3.5 (no tested fallback, DECISION-CARD §1): {sp}", file=sys.stderr)
            return 1
        env.append("DEF_DAY_ALLOW_WIDE_LIBS=1")
    print(" ".join(env))
    return 0


def cmd_verify(bundle: Path, record: Path) -> int:
    try:
        rec = json.loads(record.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        print(f"no usable record {record}: {exc}", file=sys.stderr)
        return 1
    bad = []
    for name, f in rec.get("files", {}).items():
        p = bundle / name
        if not p.is_file() or p.stat().st_size != f["bytes"] or sha(p) != f["sha256"]:
            bad.append(name)
    if bad or not rec.get("files"):
        print(f"bundle {bundle} does not match {record}: {bad or 'no files recorded'}", file=sys.stderr)
        return 1
    print(f"bundle {bundle} matches its unpack record ({len(rec['files'])} files)")
    return 0


def read_csv_col(p: Path, col: str) -> list[str]:
    import pandas as pd  # noqa: PLC0415
    return pd.read_csv(p)[col].astype(str).tolist()


def cmd_summarize(bundle: Path, light: Path, out: Path) -> int:
    man = json.loads((bundle / "manifest.json").read_text())
    panel = set(read_csv_col(bundle / "pert_counts.csv", "target_gene"))
    abc_panel = set(read_csv_col(ABC / "pert_counts.csv", "target_gene"))
    genes = read_csv_col(bundle / "gene_names.csv", "gene_name")
    abc_genes = read_csv_col(ABC / "gene_names.csv", "gene_name")
    overlap = len(panel & abc_panel)
    branch = "SAME-PANEL" if panel == abc_panel else "NEW-PANEL"
    s = {"contexts": man.get("contexts"), "panel_id": man.get("panel_id"), "n_targets": len(panel),
         "panel_vs_abc": {"same": panel == abc_panel, "overlap": overlap},
         "gene_axis_identical_to_abc": genes == abc_genes, "n_genes": len(genes), "branch": branch}
    pf = light / "preflight.json"
    if pf.exists():
        s["preflight"] = json.loads(pf.read_text())
    idj = light / "identity.json"
    if idj.exists():
        d = json.loads(idj.read_text())
        s["identity_call"] = d.get("call")
        s["identity_top5"] = {c: [(r["sample"], round(r["score"], 2), r.get("lineage")) for r in v[:5]]
                              for c, v in (d.get("consensus") or {}).items()}
    cj = light / "census.json"
    if cj.exists():
        s["census"] = json.loads(cj.read_text())
    out.write_text(json.dumps(s, indent=1, default=str) + "\n")
    print(f"branch {branch}: contexts {s['contexts']}, targets {len(panel)}, overlap with A/B/C {overlap}, "
          f"gene axis identical to A/B/C: {s['gene_axis_identical_to_abc']}; identity {s.get('identity_call')}")
    if not s["gene_axis_identical_to_abc"]:
        print("STOP: the gene axis differs from A/B/C's (DEF-RUNBOOK 531 F6): do not ship, do not refit on the day")
        return 2
    if sorted(map(str, s["contexts"] or [])) == ["A", "B", "C"]:
        print("STOP: the bundle's contexts are A/B/C: this is not the final test set")
        return 3
    return 0


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    c, a = sys.argv[1], [Path(x) for x in sys.argv[2:]]
    if c == "diff" and len(a) == 2:
        return cmd_diff(*a)
    if c == "unpack" and len(a) == 3:
        return cmd_unpack(*a)
    if c == "pick" and len(a) >= 2:
        return cmd_pick(a[0], [str(x) for x in a[1:]])
    if c == "overrides" and len(a) == 1:
        return cmd_overrides(a[0])
    if c == "verify" and len(a) == 2:
        return cmd_verify(*a)
    if c == "summarize" and len(a) == 3:
        return cmd_summarize(*a)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
