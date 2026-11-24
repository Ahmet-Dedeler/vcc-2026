"""Rule 10: prove the A/B/C build is bit-identical under everything experiment 199 changed.

199 adds three NEW files under `scripts/` and touches no production module, so the claim is
not "the build still reproduces" but the stronger "the build's inputs are byte-unchanged and
nothing new is reachable from it". That is proven here four ways, all cheap enough to run
under the 10 GB RSS cap (a production build is 18.7–19.8 GiB, see 180 §11):

1. **The written meta's own hash block, recomputed from disk against the incumbent's.** The
   incumbent is `submissions/202-rc_callset.h5ad.meta.json`, the newest meta written by the
   builder at git HEAD (202 patched the builder after 191 and proved rc_hybrid byte-identical
   under the patch, so 191's meta records the pre-202 builder hash by design). Every one of the
   8 `hashes.implementation` and 42 `hashes.inputs` digests is recomputed from the working
   tree. This IS the meta diff: a build started now would write exactly these values.
2. **Working tree against git.** Each implementation file's blob is resolved out of the launch commit `ab8aa3e`
   and compared byte for byte, so an uncommitted edit cannot hide behind a matching meta.
3. **Reachability.** No production module imports, references, or shells out to any 199 file.
4. **The builder's own test suite**, which contains
   `test_default_emitter_reproduces_the_pre_change_builder_byte_for_byte` — it resolves the
   pre-`--emitter` builder revision out of git and compares the emitted `.h5ad` sha256
   through the REAL writer.

    uv run python experiments/199-def-runbook/rule10_199.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
INCUMBENT = ROOT / "submissions/202-rc_callset.h5ad.meta.json"   # newest meta written by the HEAD builder
NEW_199 = ["scripts/def_identity.py", "scripts/def_census.py", "scripts/def_preflight.py",
           "scripts/def_postbuild.py", "scripts/def_pin_panel_artifacts.py", "scripts/def_pipeline.py"]
FROZEN_TABLE_SHA = "a654d66bc6d093140b09cfe368eb6aa2c0ff591403876588e066f260f500b319"   # 202 meta frozen_bundle.sha256
BASE_REV = "ab8aa3e"   # the commit 199 was relaunched from; its builder is the incumbent's (6cc91501...)
PRODUCTION = ["scripts/make_submission_pm_cd4.py", "scripts/build_submission.py",
              "scripts/make_submission.py", "src/vcc2026/donor_transfer.py",
              "src/vcc2026/multi_donor.py", "src/vcc2026/whiten.py", "src/vcc2026/rc_emit.py",
              "experiments/061-perturbmatch-hct-count/build_candidate.py",
              "experiments/062-cd4-post-whiten/run.py"]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    meta = json.loads(INCUMBENT.read_text())
    result: dict = {"incumbent": str(INCUMBENT.relative_to(ROOT)),
                    "incumbent_output_sha256": meta["output"]["sha256"]}
    ok = True

    # --- 1. the meta's hash block, recomputed --------------------------------------------
    #
    # `hashes.inputs` keys are the builder's own logical names, not repo paths, and
    # `frozen_input_paths` maps them (bundle/ -> data/challenge, donors/ -> data/donors,
    # cd4/ -> the CD4 compact root, perturbmatch/ -> data/external/perturbmatch).
    #
    # THE WHITENER KEY IS A RULE-11 TRAP AND IT IS LIVE IN THE INCUMBENT'S META. Line 1403
    # of the builder hardcodes `paths["whitener/whitener_union.npz"] = operator_path`, so the
    # KEY says `whitener_union.npz` whatever operator was actually passed. The incumbent
    # shipped `whitener_union_pruned.npz`. The digest under that key is correct; the NAME is
    # not. Resolved here by searching disk for the recorded sha256 — never by the key.
    def resolve(section: str, key: str) -> Path | None:
        if section == "implementation":
            return ROOT / key
        head, _, rest = key.partition("/")
        base = {"bundle": ROOT / "data/challenge",
                "donors": ROOT / "data/donors",
                "perturbmatch": ROOT / "data/external/perturbmatch",
                "cd4": ROOT / "data/external/cd4t/cloud-compact/cd4_compact"}.get(head)
        if base is not None:
            return base / rest
        return None                                    # whitener: resolved by digest

    diffs, resolved = [], {}
    for section in ("implementation", "inputs"):
        for key, recorded in meta["hashes"][section].items():
            p = resolve(section, key)
            if p is None:                              # find the file whose digest matches
                cands = [q for q in sorted((ROOT / "data/donors").iterdir())
                         if q.is_file() and q.suffix == Path(key).suffix and sha256(q) == recorded]
                p = cands[0] if cands else None
            got = sha256(p) if p is not None and p.exists() else "<absent>"
            resolved[f"{section}/{key}"] = str(p.relative_to(ROOT)) if p else None
            if got != recorded:
                diffs.append({"section": section, "key": key,
                              "resolved_path": resolved[f"{section}/{key}"],
                              "recorded": recorded, "now": got})
    result["resolved_paths"] = resolved
    result["operator_key_is_misleading"] = {
        "meta_key": "whitener/whitener_union.npz",
        "actually": resolved.get("inputs/whitener/whitener_union.npz"),
        "sha256": meta["operator"]["sha256"],
        "note": "rule 11: identify the operator by sha256, never by the key or the filename",
    }
    n_fields = sum(len(meta["hashes"][s]) for s in ("implementation", "inputs"))
    result["meta_hash_fields_compared"] = n_fields
    result["meta_hash_differences"] = diffs
    ok &= not diffs
    print(f"[1] meta hashes recomputed: {n_fields} fields, {len(diffs)} differences")
    for d in diffs:
        print(f"    DIFF {d['section']}.{d['key']}  {d['recorded'][:16]}… -> {d['now'][:16]}…")

    # --- 2. working tree vs the launch commit --------------------------------------------
    # BASE_REV, not HEAD: during the rehearsal a sweep commit (04101c5) captured the pin
    # script's transient edit of the builder, so HEAD's builder is NOT the one the incumbent
    # was built with. ab8aa3e is the commit this experiment was launched from and its builder
    # is byte-identical to the one recorded in the incumbent's meta (section 1 proves that).
    tree = []
    for rel in PRODUCTION:
        p = ROOT / rel
        if not p.exists():
            continue
        blob = subprocess.run(["git", "show", f"{BASE_REV}:{rel}"], cwd=ROOT,
                              capture_output=True).stdout
        disk = p.read_bytes()
        same = blob == disk
        ok &= same
        tree.append({"path": rel, "matches_base_rev": same, "sha256": sha256_bytes(disk)})
        print(f"[2] {'==' if same else '!='} {BASE_REV}  {rel}  {sha256_bytes(disk)[:16]}…")
    result["base_rev"] = BASE_REV
    result["working_tree_vs_base_rev"] = tree

    # --- 3. reachability -----------------------------------------------------------------
    refs = []
    for new in NEW_199:
        stem = Path(new).stem
        hits = subprocess.run(
            ["git", "grep", "-n", "--", stem, "scripts/", "src/"], cwd=ROOT,
            capture_output=True, text=True).stdout.splitlines()
        # references FROM another 199 file (the driver calls the stages) are not production
        hits = [h for h in hits if not any(h.startswith(n + ":") for n in NEW_199)]
        refs.append({"file": new, "referenced_by_production": hits})
        ok &= not hits
        print(f"[3] {new}: {len(hits)} references from scripts/ or src/")
    result["reachability"] = refs

    # --- 4. the builder's own test suite -------------------------------------------------
    t1 = time.time()
    proc = subprocess.run(["uv", "run", "--with", "pytest", "python", "-m", "pytest", "-q",
                           "scripts/test_make_submission_pm_cd4.py"],
                          cwd=ROOT, capture_output=True, text=True)
    if "No module named pytest" in (proc.stdout + proc.stderr):
        proc = subprocess.run(["uv", "run", "--with", "pytest", "pytest", "-q",
                               "scripts/test_make_submission_pm_cd4.py"],
                              cwd=ROOT, capture_output=True, text=True)
    tail = (proc.stdout or proc.stderr).strip().splitlines()[-3:]
    result["pytest"] = {"returncode": proc.returncode, "tail": tail,
                        "seconds": round(time.time() - t1, 1)}
    ok &= proc.returncode == 0
    print(f"[4] python -m pytest scripts/test_make_submission_pm_cd4.py -> rc={proc.returncode} "
          f"({result['pytest']['seconds']}s)")
    for line in tail:
        print("    " + line)

    # --- 5. the production frozen table is byte-unchanged --------------------------------
    ft = sha256(ROOT / "scripts/frozen_inputs.json")
    result["frozen_inputs_json_sha256"] = ft
    result["frozen_inputs_json_unchanged"] = ft == FROZEN_TABLE_SHA
    ok &= ft == FROZEN_TABLE_SHA
    print(f"[5] scripts/frozen_inputs.json {ft[:16]}… {'==' if ft == FROZEN_TABLE_SHA else '!='} the incumbent's")

    result["pass"] = bool(ok)
    result["seconds"] = round(time.time() - t0, 1)
    (OUT / "rule10.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"\nRULE 10 {'PASS' if ok else 'FAIL'} — meta diff is "
          f"{len(diffs)} fields of {n_fields}   ({result['seconds']}s)")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
