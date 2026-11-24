#!/usr/bin/env python3
"""811: c0_l25's many-line consensus blend applied to a production stack dump (the C1 board-read hook's second step).

dumphook811.sh runs 707's swap first (swap707.py, unchanged) and then this, on the SWAPPED dump. This file changes only the
non-panel part of the bulk rows; carriers, perts, genes and env are copied unchanged, in swap707's key order and dtypes.

THE CONSTRUCTION (= composite811's C1 step on engine folds, the same functions):
  U(t) = base811.consensus over m811's DAY bank list (GROUPS minus KOLF2), equal weight per line, unit rows, non-panel genes;
  b'_t = |b_t| normalise((1 - lam) b_t/|b_t| + lam U_t/|U_t|) on the non-panel genes (base811.blend), lam 0.25, centring off,
  for every target carrying U; other rows and the panel columns are kept byte for byte.
  base811's `bank_rows`, `unit`, `consensus`, `blend` and its GROUPS / COMPOSITION literals are taken from base811.py by ast
  and exec'd unchanged (its sha256 is pinned), as 803's c0811h.py does: importing base811 would import engine801, whose write
  guard refuses writes outside 801 (the hook writes into def_day's work dir) and which is heavy.
  m811's day rule is restated in one line (DAY = GROUPS minus KOLF2) and asserted against m811.py's source text.
MAPPING engine fold <-> dump:
  - both live on the official 18,533-gene axis in the same effect space (the stack's composition rows); the dump's genes are
    asserted == data/challenge/gene_names.csv's order;
  - the engine blends the UNKEYED composition rows (fin, then key_scale, sharpening and render); the dump's bulk rows are the
    composition after the key (sharpening is a render knob applied to the loaded dump). The blend is norm-preserving and
    depends only on each row's direction, so a positive per-row scale commutes with it: blend(s*b) = s*blend(b) up to floating
    rounding, and the key's Frobenius constant is unchanged because no row norm changes. KA-scale measures this (811 RESULTS).
  - held group: on engine folds the held line's group is out of U; on the BOARD there is no held line (hold = none). The D/E/F
    identity rule (m811.identity_hold) would hold a matched line out; on A/B/C it is a no-op for the board read because the one
    matching bank line (Jurkat, context A) carries 0 of A/B/C's targets; a single dump cannot be per-context anyway.
BANKS: 811's compact banks (work/banks; manifest sha256 pinned below, every bank file's sha256 checked against it); the
composition donors from data/donors. The hook REFUSES a panel whose targets were not all in the banks' extraction panels.

usage:
  uv run python experiments/811-lolo-model/hook811.py --dump-in SWAPPED.npz --dump-out OUT.npz --bundle B [--lam 0.25]
env HOOK811_LAMBDA overrides --lam (dumphook811.sh passes it through; 0 = the KA-off path).
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE811 = HERE / "base811.py"
BASE811_SHA = "c29f9f976d39f76ab9acc2637b9cb9a66a8db79a786c37c4159d57a2b7a32652"
M811 = HERE / "m811.py"
M811_SHA = "e00192e644828a0d00b18ca5f927f93893ef7012661aac84a4736ca592a98394"
BANKS = HERE / "work" / "banks"
BANKS_MANIFEST_SHA = "4cc6d82a9e09f5cf6a0edad28a3f704e33db3af3e4ea06dca0ecf8ed12e7bba2"
LAM = 0.25


def sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def load811():
    """base811's functions and literals, exec'd from its pinned source; m811's day rule checked against m811's text."""
    assert sha256(BASE811) == BASE811_SHA, "base811.py changed"
    assert sha256(M811) == M811_SHA, "m811.py changed"
    assert 'DAY_GROUPS = {g: b for g, b in B.GROUPS.items() if g != "KOLF2"}' in M811.read_text(), "m811's day rule changed"
    man = json.loads((BANKS / "manifest.json").read_text())
    assert sha256(BANKS / "manifest.json") == BANKS_MANIFEST_SHA, "banks manifest changed"
    for b, rec in man["banks"].items():
        assert sha256(ROOT / rec["out"]) == rec["out_sha256"], f"bank {b} differs from its manifest"
    text = BASE811.read_text()
    genes = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    ns = {"np": np, "ROOT": ROOT, "GENES": genes, "BANKS": BANKS, "Path": Path}
    want, got = {"bank_rows", "unit", "consensus", "blend"}, set()
    for node in ast.parse(text).body:
        if isinstance(node, ast.FunctionDef) and node.name in want:
            exec(compile(ast.get_source_segment(text, node), "<base811>", "exec"), ns)
            got.add(node.name)
        elif isinstance(node, ast.Assign) and any(getattr(t, "id", None) in ("GROUPS", "COMPOSITION") for t in node.targets):
            ns[node.targets[0].id] = ast.literal_eval(node.value)
    assert got == want and "GROUPS" in ns and "COMPOSITION" in ns, got
    extraction = set()
    for p in man["panels"].values():
        with open(ROOT / p) as f:
            extraction |= {r["target_gene"] for r in csv.DictReader(f)}
    return ns, genes, extraction


def transform(bulk, perts, genes, panel, lam=LAM, hold=frozenset()):
    """bulk (n x G) rows for `perts` on `genes`; `panel` = the panel whose genes are excluded (nonpanel = not in panel)."""
    ns, off, extraction = load811()
    perts, genes, panel = [str(t) for t in perts], [str(g) for g in genes], [str(p) for p in panel]
    assert genes == off, "the rows' gene axis is not the official 18,533-gene order"
    missing = sorted(set(perts) - extraction)
    assert not missing, f"{len(missing)} targets were never extracted into 811's banks (e.g. {missing[:5]}): refuse"
    day = {g: b for g, b in ns["GROUPS"].items() if g != "KOLF2"}                  # m811.DAY_GROUPS
    groups = {g: b for g, b in day.items() if g not in set(hold)}                 # m811.day_consensus(hold=...)
    pset = set(panel)
    npn = np.array([g not in pset for g in genes])
    b = np.asarray(bulk, np.float64)
    if lam == 0:
        return b, dict(lam=0.0, note="lambda 0: rows returned unchanged")
    U, info = ns["consensus"](SimpleNamespace(n=len(perts), targets=perts, nonpanel=npn), groups)
    out = ns["blend"](b, U, float(lam), npn)
    carry = np.linalg.norm(U, axis=1) > 0
    lines = np.zeros(len(perts), int)
    for line, banks in groups.items():
        S = np.zeros((len(perts), len(genes)))
        for bk in banks:
            S += ns["unit"](ns["bank_rows"](bk, perts), npn)
        lines += np.linalg.norm(S, axis=1) > 0
    nb, no = np.linalg.norm(b[:, npn], axis=1), np.linalg.norm(out[:, npn], axis=1)
    moved = carry & (nb > 0)
    cos = [float(out[i, npn] @ b[i, npn] / (nb[i] * no[i])) for i in np.flatnonzero(moved)]
    ch = out - b
    return out, dict(lam=float(lam), groups=sorted(groups), hold=sorted(hold), n_carry=int(carry.sum()), n_moved=int(moved.sum()),
                     coverage=info["coverage"], lines_per_target_q=np.quantile(lines, [0, .25, .5, .75, 1]).tolist(),
                     lines_per_target_hist={int(k): int(v) for k, v in zip(*np.unique(lines, return_counts=True))},
                     h1_group_targets=int(info["coverage"].get("H1", 0)),
                     row_norm_max_rel_change=float(np.max(np.abs(no - nb) / np.maximum(nb, 1e-300))),
                     panel_cols_unchanged=bool(np.array_equal(out[:, ~npn], b[:, ~npn])),
                     unmoved_rows_unchanged=bool(np.array_equal(out[~moved], b[~moved])),
                     moved_row_cos_to_input_q=np.round(np.quantile(cos, [0, .25, .5, .75, 1]), 4).tolist() if cos else None,
                     change_energy_share_of_bulk=float((ch * ch).sum() / (b * b).sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dump-in", type=Path, required=True)
    ap.add_argument("--dump-out", type=Path, required=True)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--lam", type=float, default=LAM)
    a = ap.parse_args()
    import os
    lam = float(os.environ.get("HOOK811_LAMBDA", a.lam))
    t0 = time.time()
    D = np.load(a.dump_in, allow_pickle=False)
    perts = [str(x) for x in D["perts"]]
    lines = [ln for ln in (a.bundle / "pert_counts.csv").read_text().splitlines() if ln.strip()]
    j = [c.strip() for c in lines[0].split(",")].index("target_gene")
    panel = [ln.split(",")[j].strip() for ln in lines[1:]]
    assert sorted(perts) == sorted(panel) and len(perts) == 300, "dump perts != the bundle's panel"
    out, info = transform(D["bulk"], perts, D["genes"], panel, lam=lam, hold=frozenset())
    a.dump_out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.dump_out, bulk=out, carriers=np.asarray(D["carriers"]), perts=np.asarray(D["perts"]),
                        genes=np.asarray(D["genes"]), env=np.asarray(D["env"]))      # swap707's key order and dtypes
    side = dict(out=str(a.dump_out), sha256=sha256(a.dump_out), dump_in=str(a.dump_in), dump_in_sha256=sha256(a.dump_in),
                bundle=str(a.bundle), script_sha256=sha256(Path(__file__)), base811_sha256=BASE811_SHA, m811_sha256=M811_SHA,
                banks_manifest_sha256=BANKS_MANIFEST_SHA, seconds=round(time.time() - t0, 1), **info)
    Path(str(a.dump_out) + ".hook811.json").write_text(json.dumps(side, indent=1))
    print(json.dumps({k: v for k, v in side.items() if k not in ("lines_per_target_hist",)}, indent=1))


if __name__ == "__main__":
    main()
