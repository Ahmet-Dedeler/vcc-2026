#!/usr/bin/env python3
"""642 -- the HEK293 contingency, as an UNAPPLIED, default-OFF patch (apply ~10-06 or later, after 538 + 554 + 590 + 641).

Why: Arc's 2026-10-01 post lists "disease status" among the six lines' axes of diversity and shows immortalized lines as a
type (Fig. 3), so a D/E/F line may be absent from CCLE. Today such a line can never switch on same-line rows: CCLE has no
HEK293, RPE1, MCF10A, BJ or HaCaT, and 590's floor makes it abstain (real-Flex HEK293: SNU-398, 5/8, 0.647). HEK293(T) is
the one non-CCLE line with a genome-wide CRISPRi bank on disk (Orion HEK293T, 294/300 of the A/B/C panel).

What it changes (three files, every new behaviour off unless DEF_DAY_ALLOW_NONCCLE_REF names a reference):
1. `scripts/def_identity.py` (590-patched): a SECOND stage. With DEF_DAY_ALLOW_NONCCLE_REF=hek293 it re-scores each
   context against CCLE plus ONE pinned non-CCLE reference row (data/identity_refs/hek293_nonflex_ref642.npz, sha-pinned,
   built from NON-Flex HEK293 controls by experiments/642-hek293-contingency/build_ref642.py) with stage 1's method weights
   frozen. It writes `nonccle_ref` (per context: call, unanimity, min_top1_score, top 5; `call` = HEK293 only when the
   augmented consensus is the reference row, unanimous over the 8 methods and min_top1_score >= 0.75; stage 1 keeps
   590's 0.65). Every stage-1 key
   (call, unanimity, min_top1_score, methods, consensus) is computed exactly as before; with the knob unset the output is
   today's. A reference whose provenance lists a Flex source is refused (642: a Flex-built row pulls Satija's Flex HAP1 to
   a unanimous 8/8 call at 0.818).
2. `scripts/def_sameline.py` (530 + 590-patched): REGISTRY_NONCCLE {HEK293: orion_hek293t}. A context that the second
   stage called HEK293 never uses its CCLE call; it reads orion_hek293t ONLY with DEF_DAY_ALLOW_NONCCLE_REF=hek293 AND a
   non-empty DEF_DAY_NONCCLE_DECISION (the coordinator's decision, recorded in the provenance); otherwise it is unchanged.
   orion_hek293t is classed like orion_hct116 (PROGRAMME_OK, STRONG: own-target knockdown visible in 0.911 of rows, equal
   to orion_hct116's 0.911).
3. `scripts/def_day.sh` (538 + 554 + 641-patched): both knobs are logged at start and recorded in def_day.final.env.

The open question this does NOT answer: whether Orion HEK293T's rows transfer to an Arc-Flex HEK293 context (the 450
same-line mechanism, measured on H1 only). Switching on is the coordinator's explicit day-of decision.

    uv run python experiments/642-hek293-contingency/patch_hek293_642.py --root <repo> --check
    uv run python experiments/642-hek293-contingency/patch_hek293_642.py --root <repo>
"""
from __future__ import annotations

import argparse
import hashlib
import os
import tempfile
from pathlib import Path

MARK = "642: the HEK293 contingency"
FILES = (Path("scripts/def_identity.py"), Path("scripts/def_sameline.py"), Path("scripts/def_day.sh"))
NEED = {FILES[0]: "590: the lowest top-1 correlation", FILES[1]: "590: a registered, unanimous call",
        FILES[2]: "641: the depth row"}
REF_SHA256 = "1b829de027a3fa8ec0c7615b2388f75fef6dc9da46520d2c06bd272500d60704"

# ------------------------------------------------------------------------------------------------ def_identity.py
ID_IMPORT_A = "import json\nimport re\nimport time\n"
ID_IMPORT_B = "import json\nimport os  # 642: the HEK293 contingency (DEF_DAY_ALLOW_NONCCLE_REF)\nimport re\nimport time\n"
ID_FUNC_A = "def match(bundle: Path, contexts: list[str], genes_csv: Path | None = None) -> dict:\n"
ID_FUNC_B = '''# 642: the HEK293 contingency. A line ABSENT from CCLE (Arc's 2026-10-01 post: the six lines span "disease status"; CCLE
# has no HEK293, RPE1, MCF10A, BJ or HaCaT) can only be named by a reference row outside CCLE. The second stage below adds
# ONE pinned row and re-scores the contexts with stage 1's method weights frozen; it never changes a stage-1 key. Off
# unless DEF_DAY_ALLOW_NONCCLE_REF names the reference. Validated in experiments/642-hek293-contingency (real-Flex HEK293
# 8/8 at 0.820; 0 false calls on 423 labelled profiles).
NONCCLE_REFS = {
    "hek293": ("data/identity_refs/hek293_nonflex_ref642.npz", "''' + REF_SHA256 + '''"),
}
NONCCLE_FLOOR = 0.75   # the second stage's own floor (stage 1 keeps 590's 0.65): the true positive is 0.820, the
                       # nearest non-HEK293 profiles with the row top-1 in >= 5/8 methods sit at 0.46-0.69


def load_nonccle_ref(key: str):
    rel, want = NONCCLE_REFS[key]
    got = sha256(ROOT / rel)
    if got != want:
        raise SystemExit(f"642: {rel} sha256 {got[:16]}... is not the pinned {want[:16]}...: refusing an unpinned reference")
    z = np.load(ROOT / rel, allow_pickle=False)
    prov = json.loads(str(z["provenance"]))
    flex = [s.get("name") for s in prov.get("sources", []) if "flex" in str(s.get("chemistry", "")).lower()]
    if flex or not prov.get("sources"):
        # 642: a Flex-built HEK293 row pulled Satija's Flex HAP1 to a unanimous 8/8 call at 0.818 (same-platform references
        # attract same-platform queries; the contexts are Flex). Only non-Flex sources may build a reference.
        raise SystemExit(f"642: reference {key} lists Flex (or no) sources {flex}: refusing (never build one from Flex)")
    return z["genes"].astype(str), z["values"].astype(np.float64), prov


def nonccle_stage(key: str, X, genes, samples, feature_sets: dict, method_quality: dict, queries: dict, names) -> dict:
    """Score `names` (keys of `queries`) against CCLE + the pinned `key` row; method weights = stage 1's (frozen). A gene
    the reference lacks takes CCLE's per-gene median. `call[name]` is the line only when the augmented consensus is the
    row, unanimous over the methods and min_top1_score >= NONCCLE_FLOOR."""
    rgenes, rvals, prov = load_nonccle_ref(key)
    pos = {g: i for i, g in enumerate(rgenes)}
    med = np.median(X, axis=0)
    row = np.array([rvals[pos[g]] if g in pos and np.isfinite(rvals[pos[g]]) else med[j] for j, g in enumerate(genes)])
    rname = f"{key.upper()}_REF"
    Xa = np.vstack([X, row[None, :].astype(X.dtype)])
    sa = np.append(np.asarray(samples).astype(str), rname)
    top = {n: {} for n in names}
    for fname, fidx in feature_sets.items():
        for transform in ("pearson", "spearman"):
            k = f"{fname}:{transform}"
            Xm = Xa[:, fidx]
            if transform == "spearman":
                Xm = rankdata(Xm, axis=1).astype(np.float32)
            for n in names:
                qf = queries[n]
                valid = np.isfinite(qf[fidx])
                q = qf[fidx][valid]
                if transform == "spearman":
                    q = rankdata(q)
                sc = row_pearson(Xm[:, valid], q)
                order = np.argsort(-sc)[:25]
                top[n][k] = [(str(sa[i]), r + 1, float(sc[i])) for r, i in enumerate(order)]
            del Xm
    per = {}
    for n in names:
        acc, tw = {}, 0.0
        for k, rows in top[n].items():
            w = method_quality[k] ** 2
            tw += w
            for s, rank, _ in rows:
                acc[s] = acc.get(s, 0.0) + w * (26 - rank)
        ordered = sorted(acc, key=acc.get, reverse=True)[:5]
        top1 = {k: rows[0][0] for k, rows in top[n].items()}
        s1 = float(min(rows[0][2] for rows in top[n].values()))
        unanimous = len(set(top1.values())) == 1
        per[n] = {"call": ordered[0], "unanimous": unanimous, "n_methods_agree": sum(v == ordered[0] for v in top1.values()),
                  "min_top1_score": s1, "top5": [(s, acc[s] / max(tw, 1e-12)) for s in ordered],
                  "row_top1_methods": sum(v == rname for v in top1.values()),
                  "accepted": bool(ordered[0] == rname and unanimous and s1 >= NONCCLE_FLOOR)}
    return {"reference": key, "row_name": rname, "file": NONCCLE_REFS[key][0], "sha256": NONCCLE_REFS[key][1],
            "provenance": prov, "floor": NONCCLE_FLOOR, "weights": "stage 1 method_quality, frozen", "per_context": per,
            "call": {n: (key.upper() if per[n]["accepted"] else None) for n in names}}


def match(bundle: Path, contexts: list[str], genes_csv: Path | None = None) -> dict:
'''
ID_RET_A = (
    "        print(f\"  lowest top-1 correlation over the 8 methods: {out['min_top1_score'][ctx]:.3f} \"\n"
    "              f\"(590: same-line rows need >= 0.65; A/B/C 0.719-0.772)\")\n"
    "    return out\n"
)
ID_RET_B = (
    "        print(f\"  lowest top-1 correlation over the 8 methods: {out['min_top1_score'][ctx]:.3f} \"\n"
    "              f\"(590: same-line rows need >= 0.65; A/B/C 0.719-0.772)\")\n"
    "    # 642: the HEK293 contingency's second stage (off unless DEF_DAY_ALLOW_NONCCLE_REF names a reference). Every key\n"
    "    # above is final; this only adds `nonccle_ref`.\n"
    "    _nref = os.environ.get(\"DEF_DAY_ALLOW_NONCCLE_REF\", \"\")\n"
    "    if _nref:\n"
    "        if _nref not in NONCCLE_REFS:\n"
    "            raise SystemExit(f\"642: DEF_DAY_ALLOW_NONCCLE_REF={_nref!r} is not one of {sorted(NONCCLE_REFS)}\")\n"
    "        out[\"nonccle_ref\"] = nonccle_stage(_nref, X, genes, samples, feature_sets, method_quality, queries, list(contexts))\n"
    "        for ctx in contexts:\n"
    "            r = out[\"nonccle_ref\"][\"per_context\"][ctx]\n"
    "            print(f\"[642] {ctx}: with the {_nref} reference row: {r['call']} {r['n_methods_agree']}/8, min top-1 \"\n"
    "                  f\"{r['min_top1_score']:.3f}, row top-1 in {r['row_top1_methods']}/8 methods -> \"\n"
    "                  f\"{'CALLED ' + _nref.upper() if r['accepted'] else 'not called'}\", flush=True)\n"
    "    return out\n"
)

# ------------------------------------------------------------------------------------------------ def_sameline.py
SL_IMPORT_A = "import json\nimport time\nfrom pathlib import Path\n"
SL_IMPORT_B = "import json\nimport os  # 642: the HEK293 contingency (DEF_DAY_ALLOW_NONCCLE_REF, DEF_DAY_NONCCLE_DECISION)\nimport time\nfrom pathlib import Path\n"
SL_REG_A = '    "THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": ("wessels_thp1",),\n}\n'
SL_REG_B = SL_REG_A + (
    "# 642: the HEK293 contingency. Lines ABSENT from CCLE, named only by def_identity's second stage against a pinned\n"
    "# non-CCLE reference row (DEF_DAY_ALLOW_NONCCLE_REF=<key>, default off). Switch-on also needs DEF_DAY_NONCCLE_DECISION,\n"
    "# the coordinator's recorded decision. HEK293T (Orion) is a HEK293 derivative.\n"
    "REGISTRY_NONCCLE: dict[str, tuple[str, ...]] = {\n"
    "    \"HEK293\": (\"orion_hek293t\",),\n"
    "}\n"
)
SL_CLASS_A = 'STRONG_BANKS = {"orion_hct116", "k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2", "sunshine_calu3"}\n'
SL_CLASS_B = SL_CLASS_A + (
    "# 642: orion_hek293t is classed like orion_hct116 (same lab and platform, unstimulated genome-wide CRISPRi; own-target\n"
    "# knockdown visible in 0.911 of rows, equal to orion_hct116's 0.911). Reachable only through REGISTRY_NONCCLE.\n"
    "PROGRAMME_OK = PROGRAMME_OK | {\"orion_hek293t\"}\n"
    "STRONG_BANKS = STRONG_BANKS | {\"orion_hek293t\"}\n"
)
SL_DIRS_A = '    "h1_2025_full": ROOT / "data/donors/h1_2025_full",\n}\n'
SL_DIRS_B = ('    "h1_2025_full": ROOT / "data/donors/h1_2025_full",\n'
             '    "orion_hek293t": ROOT / "data/donors/orion_hek293t",   # 642: REGISTRY_NONCCLE only\n}\n')
SL_LOOP_A = (
    "    for c, sample in sorted(call.items()):\n"
    "        banks = REGISTRY.get(str(sample))\n"
    "        rec: dict[str, object] = {\"call\": sample, \"unanimous\": ident.get(\"call_unanimous_across_methods\", {}).get(c)}\n"
)
SL_LOOP_B = SL_LOOP_A + (
    "        # 642: def_identity's second stage named this context a line ABSENT from CCLE. Its CCLE call is then the nearest\n"
    "        # wrong line and never switches rows on. The line's own bank switches on only with DEF_DAY_ALLOW_NONCCLE_REF naming\n"
    "        # the reference AND a non-empty DEF_DAY_NONCCLE_DECISION (the coordinator's decision, recorded here).\n"
    "        _nonccle = None\n"
    "        _st2: dict = {}\n"
    "        _ncall = None if a.force else ((ident.get(\"nonccle_ref\") or {}).get(\"call\") or {}).get(c)\n"
    "        if _ncall:\n"
    "            _nb = REGISTRY_NONCCLE.get(str(_ncall))\n"
    "            _knob = os.environ.get(\"DEF_DAY_ALLOW_NONCCLE_REF\", \"\")\n"
    "            _decision = os.environ.get(\"DEF_DAY_NONCCLE_DECISION\", \"\").strip()\n"
    "            _st2 = ident[\"nonccle_ref\"][\"per_context\"][c]\n"
    "            rec.update(nonccle_ref_call=_ncall, nonccle_ref_sha256=ident[\"nonccle_ref\"].get(\"sha256\"),\n"
    "                       nonccle_unanimous=_st2.get(\"unanimous\"), nonccle_min_top1=_st2.get(\"min_top1_score\"))\n"
    "            if not _nb:\n"
    "                _why = \"no REGISTRY_NONCCLE bank\"\n"
    "            elif _knob != ident[\"nonccle_ref\"].get(\"reference\"):\n"
    "                _why = f\"DEF_DAY_ALLOW_NONCCLE_REF={_knob!r} does not name the reference\"\n"
    "            elif not _decision:\n"
    "                _why = \"no DEF_DAY_NONCCLE_DECISION recorded\"\n"
    "            else:\n"
    "                _why = None\n"
    "            if _why:\n"
    "                rec[\"mode\"] = (f\"none (642: the non-CCLE reference called {_ncall}, but {_why}; the CCLE call \"\n"
    "                               f\"{sample} is not used)\")\n"
    "                prov[\"contexts\"][c] = rec\n"
    "                print(f\"[sameline] {c}: the non-CCLE reference called {_ncall}; {_why}; context unchanged (642)\", flush=True)\n"
    "                continue\n"
    "            _nonccle, sample, banks = _ncall, _ncall, _nb\n"
    "            rec.update(unanimous=True, nonccle_decision=_decision)\n"
    "            print(f\"[sameline] {c}: the non-CCLE reference called {_ncall} (8/8, min top-1 {_st2['min_top1_score']:.3f}); \"\n"
    "                  f\"decision recorded: {_decision!r} (642)\", flush=True)\n"
)
SL_GATE_A = '            s1 = ident.get("min_top1_score", {}).get(c)\n'
SL_GATE_B = ('            # 642: a non-CCLE call is gated by its own (second-stage) score\n'
             '            s1 = _st2.get("min_top1_score") if _nonccle else ident.get("min_top1_score", {}).get(c)\n')

# ------------------------------------------------------------------------------------------------ def_day.sh
DAY_LOG_A = ("[ \"${DEF_DAY_ALLOW_DEPTH:-}\" = 1 ] && log \"641: the depth row: DEF_DAY_ALLOW_DEPTH=1: preflight's control-library "
             "median check is overridden\"\n")
DAY_LOG_B = DAY_LOG_A + (
    "[ -n \"${DEF_DAY_ALLOW_NONCCLE_REF:-}\" ] && log \"642: the HEK293 contingency: DEF_DAY_ALLOW_NONCCLE_REF=$DEF_DAY_ALLOW_NONCCLE_REF "
    "(def_identity's non-CCLE second stage); decision: ${DEF_DAY_NONCCLE_DECISION:-<none recorded>}\"\n"
)
DAY_ENV_A = "[ \"${DEF_DAY_ALLOW_DEPTH:-}\" = 1 ] && printf 'DEF_DAY_ALLOW_DEPTH=1\\n' >> \"$W/def_day.final.env\"\n"
DAY_ENV_B = DAY_ENV_A + (
    "# 642: the HEK293 contingency's knob and decision, recorded\n"
    "[ -n \"${DEF_DAY_ALLOW_NONCCLE_REF:-}\" ] && printf 'DEF_DAY_ALLOW_NONCCLE_REF=%q\\n' \"$DEF_DAY_ALLOW_NONCCLE_REF\" >> \"$W/def_day.final.env\"\n"
    "[ -n \"${DEF_DAY_NONCCLE_DECISION:-}\" ] && printf 'DEF_DAY_NONCCLE_DECISION=%q\\n' \"$DEF_DAY_NONCCLE_DECISION\" >> \"$W/def_day.final.env\"\n"
)


def edits() -> dict[Path, list[tuple[str, str]]]:
    return {FILES[0]: [(ID_IMPORT_A, ID_IMPORT_B), (ID_FUNC_A, ID_FUNC_B), (ID_RET_A, ID_RET_B)],
            FILES[1]: [(SL_IMPORT_A, SL_IMPORT_B), (SL_REG_A, SL_REG_B), (SL_CLASS_A, SL_CLASS_B), (SL_DIRS_A, SL_DIRS_B),
                       (SL_LOOP_A, SL_LOOP_B), (SL_GATE_A, SL_GATE_B)],
            FILES[2]: [(DAY_LOG_A, DAY_LOG_B), (DAY_ENV_A, DAY_ENV_B)]}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def plan(root: Path) -> dict[Path, str]:
    texts = {f: (root / f).read_text() for f in FILES}
    missing = [f"{f} (needs {NEED[f][:3]})" for f in FILES if NEED[f] not in texts[f]]
    if missing:
        raise SystemExit("642 applies after 538 + 554 + 590 + 641; missing: " + ", ".join(missing))
    marked = [str(f) for f in FILES if MARK in texts[f]]
    if marked:
        raise SystemExit(f"already patched ({', '.join(marked)}): refusing a second apply" if len(marked) == len(FILES)
                         else f"partially patched ({', '.join(marked)}): restore the three files and re-apply")
    ref = root / "data/identity_refs/hek293_nonflex_ref642.npz"
    if not ref.exists() or sha(ref) != REF_SHA256:
        raise SystemExit(f"{ref}: missing or not the pinned reference (sha256 {REF_SHA256[:16]}...)")
    new = {}
    for f, pairs in edits().items():
        t = texts[f]
        for old, rep in pairs:
            n = t.count(old)
            if n != 1:
                raise SystemExit(f"{f}: anchor found {n} times (need 1): {old.splitlines()[0][:90]!r}")
            t = t.replace(old, rep)
        if f.suffix == ".py":
            compile(t, str(f), "exec")
        new[f] = t
    return new


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    root = a.root.resolve()
    new = plan(root)
    for f in FILES:
        print(f"{f}\n  before {sha(root / f)}")
    if a.check:
        print("patches apply cleanly and compile (not written)")
        return
    for f in FILES:
        dst = root / f
        fd, tmp = tempfile.mkstemp(dir=dst.parent, prefix=f".{dst.name}.642-")
        with os.fdopen(fd, "w") as fh:
            fh.write(new[f])
        os.chmod(tmp, os.stat(dst).st_mode)
        os.replace(tmp, dst)
        print(f"{f}\n  after  {sha(dst)}")
    print("patched " + ", ".join(str(f) for f in FILES))


if __name__ == "__main__":
    main()
