"""803 §8 bulk_transform: 811's M811 arm c0_l<lam> (the many-line consensus U, no centring, blended norm-preservingly into each
target's composition row) on an H1 real-emitter surface.

811's code is reused VERBATIM, not imported: base811.py imports engine801 -> guard801, whose write guard would refuse the
emitter's own outputs. The function sources `bank_rows`, `unit`, `consensus`, `blend` and the `GROUPS` literal are taken from
base811.py by ast and exec'd unchanged (sha256 of the source file asserted); m811's day rule is restated in two lines:
DAY_GROUPS = GROUPS minus KOLF2 (m811.DAY_GROUPS), groups = DAY_GROUPS minus the held group (m811.day_consensus with
centre=False, which calls base811.consensus on SimpleNamespace(n, targets, nonpanel)). Held group: H1 (the truth line; VCC25
H1 and Pan H1 never enter). Banks: 803's H1-target extraction (bt/banks_h1, ext811.read_bank, KA'd against 811's banks);
the composition donors from data/donors (base811.bank_rows' own branch).
Axis: U lives on the official 18,533-gene axis (811). The H1 rows live on the surface's cap15 axis; U is mapped by gene name,
and the blend runs over npn = genes off the surface's panel AND on the official axis (genes outside the official axis keep
their values). blend is base811.blend (707's formula), applied to every target that carries U.
"""
from __future__ import annotations

import ast
import csv
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE811 = ROOT / "experiments/811-lolo-model/base811.py"
BASE811_SHA = "c29f9f976d39f76ab9acc2637b9cb9a66a8db79a786c37c4159d57a2b7a32652"
M811 = ROOT / "experiments/811-lolo-model/m811.py"
BANKS_H1 = ROOT / "experiments/803-de-engine/bt/banks_h1"


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _git_blob(p):
    b = Path(p).read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(b) + b).hexdigest()


def _load811():
    assert BASE811_SHA == "PINNED_AT_STAGING" or _sha(BASE811) == BASE811_SHA, "base811.py moved"
    text = BASE811.read_text()
    genes = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    ns = {"np": np, "ROOT": ROOT, "GENES": genes, "COMPOSITION": ("k562_gwps", "orion_hct116", "orion_hek293t"),
          "BANKS": BANKS_H1, "Path": Path}
    want = {"bank_rows", "unit", "consensus", "blend"}
    got = set()
    for node in ast.parse(text).body:
        if isinstance(node, ast.FunctionDef) and node.name in want:
            exec(compile(ast.get_source_segment(text, node), "<base811>", "exec"), ns)
            got.add(node.name)
        elif isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "GROUPS" for t in node.targets):
            ns["GROUPS"] = ast.literal_eval(node.value)
    assert got == want and "GROUPS" in ns, got
    assert ns["COMPOSITION"] == ast.literal_eval([n.value for n in ast.parse(text).body if isinstance(n, ast.Assign)
                                                  and any(getattr(t, "id", None) == "COMPOSITION" for t in n.targets)][0])
    return ns, genes


def transform(eff, genes, targets, surface, lam, hold=("H1",)):
    ns, off = _load811()
    tg = [str(t) for t in targets]
    gl = [str(g) for g in genes]
    day = {g: b for g, b in ns["GROUPS"].items() if g != "KOLF2"}                  # m811.DAY_GROUPS
    groups = {g: b for g, b in day.items() if g not in set(hold)}                 # m811.day_consensus(hold=...)
    pset = set(tg)
    nonpanel_off = np.array([g not in pset for g in off])
    U_off, info = ns["consensus"](SimpleNamespace(n=len(tg), targets=tg, nonpanel=nonpanel_off), groups)
    opos = {g: i for i, g in enumerate(off)}
    cols = np.array([opos.get(g, -1) for g in gl])
    U = np.zeros((len(tg), len(gl)))
    U[:, cols >= 0] = U_off[:, cols[cols >= 0]]
    npn = np.array([(g not in pset) for g in gl]) & (cols >= 0)
    fin = np.asarray(eff, np.float64)
    out = ns["blend"](fin, U, float(lam), npn)
    carry = np.linalg.norm(U[:, npn], axis=1) > 0
    ch = out - fin
    lines_per_target = np.zeros(len(tg), int)
    for line, banks in groups.items():
        S = np.zeros((len(tg), len(off)))
        for b in banks:
            S += ns["unit"](ns["bank_rows"](b, tg), nonpanel_off)
        lines_per_target += np.linalg.norm(S, axis=1) > 0
    info.update(lam=float(lam), hold=list(hold), groups={k: list(v) for k, v in groups.items()},
                n_carry=int(carry.sum()), n_targets=len(tg), lines_per_target=lines_per_target.tolist(),
                change_energy=float((ch * ch).sum() / (fin * fin).sum()),
                row_norm_npn_max_rel_change=float(np.max(np.abs(np.linalg.norm(out[:, npn], axis=1)
                                                                / np.maximum(np.linalg.norm(fin[:, npn], axis=1), 1e-300) - 1))),
                n_genes_off_official=int((cols < 0).sum()), base811_sha256=_sha(BASE811), base811_git_blob=_git_blob(BASE811),
                m811_sha256=_sha(M811), m811_git_blob=_git_blob(M811), banks_h1_manifest_sha256=_sha(BANKS_H1 / "manifest.json"))
    return out, info
