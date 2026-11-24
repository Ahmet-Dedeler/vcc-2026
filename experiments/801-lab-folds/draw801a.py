"""801 follow-up 5 -- the corrected default draw rule "A-only" (experiments/800-checker/DESIGN.md section 3.1 amendment,
2026-10-16 11:48Z): draw801.py's arc800 procedure (its Rule class, imported, unchanged) plus ONE change of rule, the pool also
removes the RPE1 and Jurkat essential screens (data/donors/rpe1/perts.txt U data/donors/nadig_jurkat/perts.txt). A HepG2
list: none on disk under data/donors and 814 names none (LOG l.200-212 lists Nadig HepG2 2,394 = the same DepMap set as
Jurkat's 2,394); recorded as absent.

Six INDEPENDENT draws A1..A6 (seeds 801301..801306) from the same full pool, as the lead asked: they may overlap each other and
L1-L6 (recorded), so
  * the exclusion is every pert_counts.csv on disk EXCEPT 801's own L1-L6 bundles (the explicit list
    work/manifest/panel_files_801a.txt, made on the Mac with `find experiments data -name pert_counts.csv`, sha256 checked
    here), which keeps L7/L8 (sealed), LABC (= A/B/C), A/B/C, R6, R7a and every other panel out;
  * capacity per cell is the whole pool cell (one draw's need), since the draws need not be disjoint (draw801's floor(pool/8)
    existed only for eight disjoint draws). Cells short even so move by draw801's frozen rule.
ASSERTED at draw time: 300 unique, 0 K562-essential, 0 essential-screen (RPE1 U Jurkat), 0 overlap with the exclusion
(incl. L7/L8), all on the axis, strata counts == quotas. Each draw is written to work/draws801a.json as it is drawn, with
its sha256, overlaps (vs L1-L6, A/B/C, each other) and strength-bin shares; the bundles (pert_counts.csv, gene_names.csv,
manifest.json, controls = copies of 199's relabelled files, sha256 checked) follow.

    uv run python experiments/801-lab-folds/draw801a.py          (EC2; light: < 1 GB, seconds)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import collections  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import draw801 as D  # noqa: E402
import diag801  # noqa: E402

WORK = HERE / "work"
OUT = WORK / "draws801a.json"
SEEDS = tuple((f"A{i}", 801300 + i) for i in range(1, 7))
SCREENS = (D.ROOT / "data/donors/rpe1/perts.txt", D.ROOT / "data/donors/nadig_jurkat/perts.txt")
LIST = WORK / "manifest/panel_files_801a.sha256"


class RuleA(D.Rule):
    def __init__(self):
        super().__init__()
        # the explicit exclusion list (not this machine's find): every panel on disk but 801's L1-L6
        self.files, self.excluded, self.panels = [], set(), {}
        for ln in LIST.read_text().splitlines():
            sha, rel = ln.split(None, 1)
            p = D.ROOT / rel.strip()
            assert D.sha256(p) == sha, f"{rel}: sha256 differs from the Mac's list"
            t = D.read_targets(p)
            self.files.append(p)
            self.excluded |= set(t)
            self.panels.setdefault(D.sha_list(set(t)), []).append(rel.strip())
        self.screen = set().union(*(D.lines(p) for p in SCREENS))

    def pool(self, earlier: set[str]) -> list[str]:
        return [t for t in super().pool(earlier) if t not in self.screen]


def main():
    t0 = time.time()
    assert not OUT.exists(), f"{OUT} exists: drawn once"
    D.N_DRAWS = 1                                            # independent draws: one draw's need per cell
    R = RuleA()
    pool = R.pool(set())
    qrec, q = R.quotas(pool)
    L = {k: set(json.loads((WORK / "draws.json").read_text())["draws"][k]["targets"]) for k in ("L1", "L2", "L3", "L4", "L5", "L6")}
    sealed = {k: set(D.read_targets(WORK / "sealed" / k / "pert_counts.csv")) for k in ("L7", "L8")}
    abc = set(R.abc)
    rec = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), rule="A-only (draw801a.py docstring; DESIGN 3.1 amendment)",
               script_sha256=D.sha256(Path(__file__)), draw801_sha256=D.sha256(HERE / "draw801.py"),
               screens={str(p.relative_to(D.ROOT)): D.sha256(p) for p in SCREENS}, screen_size=len(R.screen),
               hepg2_list="absent (none under data/donors; 814 names none)",
               n_panel_files=len(R.files), n_excluded_targets=len(R.excluded), panel_files=R.panels,
               essential_list_size=len(R.ess), basal_threshold_cpm=R.basal_thr, strength_edges=R.edges.tolist(),
               pool_size=len(pool), pool_without_screens_removed=len(D.Rule.pool(R, set())), strata=qrec, draws={})
    print(f"pool {len(pool)} (arc800 pool here {rec['pool_without_screens_removed']}); moves {qrec['moves']}", flush=True)
    byc = collections.defaultdict(list)
    for t in pool:
        byc[R.cell(t)].append(t)
    panels = {}
    for key, seed in SEEDS:
        rng = np.random.default_rng(seed)
        panel = []
        for k in sorted(q):
            if q[k]:
                panel += rng.choice(np.asarray(sorted(byc[k])), size=q[k], replace=False).tolist()
        panel = sorted(panel)
        s = set(panel)
        assert len(panel) == 300 and len(s) == 300
        assert not (s & R.ess), f"{key}: K562-essential"
        assert not (s & R.screen), f"{key}: essential-screen target"
        assert not (s & R.excluded), f"{key}: overlaps an excluded panel"
        assert not (s & sealed["L7"]) and not (s & sealed["L8"]), f"{key}: overlaps the sealed reserve"
        assert s <= R.gene_set
        got = collections.Counter(R.cell(t) for t in panel)
        assert all(got.get(k, 0) == q[k] for k in q), f"{key}: strata"
        d = dict(seed=seed, targets=panel, sha256=D.sha_list(panel), n=300, overlap_k562_essential=0, overlap_screens=0,
                 overlap_excluded=0, overlap_L={k: len(s & v) for k, v in L.items()}, overlap_abc=len(s & abc),
                 overlap_earlier_A={k: len(s & v) for k, v in panels.items()}, describe=diag801.describe(R, panel))
        rec["draws"][key] = d
        panels[key] = s
        OUT.write_text(json.dumps(rec, indent=1))
        print(f"[{key}] seed {seed} sha {d['sha256'][:16]}; overlap L {d['overlap_L']}; A/B/C {d['overlap_abc']}; "
              f"earlier A {d['overlap_earlier_A']}; bins {d['describe']['strength_bins']}", flush=True)
    for key, _ in SEEDS:
        out = WORK / key / "bundle"
        out.mkdir(parents=True, exist_ok=False)
        pd.DataFrame({"target_gene": rec["draws"][key]["targets"]}).to_csv(out / "pert_counts.csv", index=False)
        shutil.copyfile(D.SRC / "gene_names.csv", out / "gene_names.csv")
        man = json.loads((D.SRC / "manifest.json").read_text())
        man.update(partition="test", panel_id=f"vcc2026-synthetic-lab801a-{key}", contexts=[D.RENAME[c] for c in man["contexts"]],
                   per_context={D.RENAME[k]: v for k, v in man["per_context"].items()})
        (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
        for new in D.RENAME.values():
            ref = D.DEFB / f"context_{new}.h5ad"
            shutil.copyfile(ref, out / ref.name)
            assert D.sha256(out / ref.name) == D.sha256(ref)
        rec["draws"][key]["bundle"] = dict(path=str(out.relative_to(D.ROOT)),
                                           files={f: D.sha256(out / f) for f in ("pert_counts.csv", "gene_names.csv", "manifest.json")})
    rec["seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(rec, indent=1))
    print(f"-> {OUT.relative_to(D.ROOT)} sha256 {D.sha256(OUT)} ({rec['seconds']} s)")


if __name__ == "__main__":
    main()
