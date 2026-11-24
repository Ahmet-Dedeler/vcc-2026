"""552 -- draw reserve R7a (RESERVE.md section 2): draw534.py verbatim with R7's seeds and names. The exclusion scan
(`find experiments data -name pert_counts.csv`) now also finds R6's four bundles (534/work/partB/S1..S4), so R7a is
disjoint from R6 as well as from every tuned-on panel. Outputs under experiments/552-r7-reserve/work/. draw534's
docstring follows.

534 Part B -- draw reserve R6 (PREREGISTRATION-B section 2) and build its four synthetic bundles the way bundle238.py does.

bundle238's rules, verbatim, with the preregistered change: the eligible pool also excludes every target of every
`pert_counts.csv` under experiments/ and data/ at draw time (`find experiments data -name pert_counts.csv`, re-run here)
and each earlier R6 draw, so the four panels are disjoint from every panel on disk and from each other.

  S1  efr rule  K562 n HCT116 n HEK293T n CD4 atlas   seed 534101
  S2  efp rule  K562 u HCT116 u HEK293T               seed 534201
  S3  efr rule  (as S1)                               seed 534102
  S4  efp rule  (as S2)                               seed 534202
Drawn in this order. eligible = sorted(genes & cover - excluded - essential); panel = sorted(rng.choice(eligible, 300,
replace=False)) with rng = np.random.default_rng(seed) (bundle238's two lines).

ASSERTED AT DRAW TIME: 300 unique targets; 0 essential; 0 overlap with any excluded panel or earlier R6 draw; all on the
official gene axis. Each draw's list, pool size and sha256 go to work/partB/draws.json as it is drawn.

BUNDLES (bundle238's recipe): pert_counts.csv (target_gene), gene_names.csv (A/B/C's, copied), manifest.json (partition
test, a panel_id, contexts A/B/C -> D/E/F), and the three control files through 199's own `relabel`, asserted
byte-identical to 199's relabelled files (def-bundle for the efp rule, def-bundle-r for the efr rule, as bundle238 did).
Storage only: once a relabelled file is proven identical, it is replaced by an APFS clone (cp -c) of 199's file, so the
four bundles cost ~0 disk (the bytes are the same; sha256 re-checked after the clone).

    uv run python experiments/552-r7-reserve/draw552.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard552  # noqa: E402,F401
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent / "work"
ROOT = HERE.parents[2]
SRC = ROOT / "data/challenge"
E199 = ROOT / "experiments/199-def-runbook"
CD4_SOURCE = E199 / "work/GWCD4i.DE_stats.h5ad"
ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
DRAWS = HERE / "draws.json"
SPECS = (  # key, rule, pool, seed (RESERVE.md's table, in draw order)
    ("T1", "efr", "intersect", 552101),
    ("T2", "efp", "union", 552201),
    ("T3", "efr", "intersect", 552102),
    ("T4", "efp", "union", 552202),
)
spec = importlib.util.spec_from_file_location("defbundle199", E199 / "defbundle199.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def sha_list(xs) -> str:
    return hashlib.sha256("\n".join(sorted(xs)).encode()).hexdigest()


def essential() -> set[str]:
    return {ln.strip() for ln in ESSENTIAL_LIST.read_text().splitlines() if ln.strip()}


def cd4_targets() -> set[str]:
    """bundle238.cd4_targets (199r's reader), verbatim."""
    with h5py.File(CD4_SOURCE, "r") as h:
        node = h["obs"]["target_contrast_gene_name"]
        if isinstance(node, h5py.Group):
            return {c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]}
        return {v.decode() if isinstance(v, bytes) else str(v) for v in node[:]}


def donor_set(name: str) -> set[str]:
    return {ln.strip() for ln in (ROOT / f"data/donors/{name}/perts.txt").read_text().splitlines() if ln.strip()}


def panel_files() -> list[Path]:
    """`find experiments data -name pert_counts.csv`, re-run now (the preregistered exclusion)."""
    r = subprocess.run(["find", "experiments", "data", "-name", "pert_counts.csv"], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return sorted(ROOT / ln.strip() for ln in r.stdout.splitlines() if ln.strip())


def read_targets(p: Path) -> list[str]:
    df = pd.read_csv(p)
    col = "target_gene" if "target_gene" in df.columns else df.columns[0]      # conservative: exclude, never skip
    return [str(t) for t in df[col].dropna()]


def main() -> None:
    t0 = time.time()
    assert not DRAWS.exists(), f"{DRAWS} exists: R7a is drawn once"
    r6 = ROOT / "experiments/534-generalization-audit/work/partB/draws.json"
    r6_targets = set().union(*(set(d["targets"]) for d in json.loads(r6.read_text())["draws"].values()))
    genes = [str(g) for g in pd.read_csv(SRC / "gene_names.csv")["gene_name"]]
    gene_set = set(genes)
    ess = essential()
    files = panel_files()
    excluded, panels = set(), {}
    for p in files:
        t = read_targets(p)
        excluded |= set(t)
        panels.setdefault(sha_list(set(t)), []).append(str(p.relative_to(ROOT)))
    union_cover = set().union(*(donor_set(d) for d in DONORS))
    inter_cover = set.intersection(*(donor_set(d) for d in DONORS)) & cd4_targets()
    assert r6_targets <= excluded, "R6's bundles were not found by the exclusion scan"
    rec = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), rule="bundle238 + exclusion of every panel file on disk "
               "(R6's four included) and each earlier R7a draw (552 RESERVE.md section 2)", r6_excluded=len(r6_targets), n_panel_files=len(files), n_distinct_panels=len(panels),
               panel_files={k[:16]: v for k, v in panels.items()}, n_excluded_targets=len(excluded),
               essential_list_sha256=hashlib.sha256(ESSENTIAL_LIST.read_bytes()).hexdigest(), essential_list_size=len(ess),
               pools=dict(union=len(union_cover), intersect=len(inter_cover)), draws={})
    print(f"{len(files)} pert_counts.csv files, {len(panels)} distinct panels, {len(excluded)} excluded targets", flush=True)
    for k, v in panels.items():
        print(f"   {k[:16]}  {len(v)} file(s): {', '.join(v)}", flush=True)
    earlier: set[str] = set()
    for key, rule, pool, seed in SPECS:
        cover = union_cover if pool == "union" else inter_cover
        eligible_with = sorted(gene_set & cover - excluded - earlier)
        eligible = sorted(set(eligible_with) - ess)
        rng = np.random.default_rng(seed)
        panel = sorted(rng.choice(np.asarray(eligible), size=300, replace=False).tolist())
        assert len(panel) == 300 and len(set(panel)) == 300
        assert not (set(panel) & ess), f"{key}: essential targets"
        assert not (set(panel) & excluded), f"{key}: overlaps an excluded panel"
        assert not (set(panel) & earlier), f"{key}: overlaps an earlier R7a draw"
        assert not (set(panel) & r6_targets), f"{key}: overlaps R6"
        assert set(panel) <= gene_set, f"{key}: a target is off the gene axis"
        d = dict(rule=rule, pool=pool, seed=seed, pool_size_before_essential=len(eligible_with), pool_size=len(eligible),
                 essential_removed=len(eligible_with) - len(eligible), targets=panel, sha256=sha_list(panel),
                 overlap_excluded=0, overlap_essential=0, overlap_earlier_r7a=0, overlap_r6=0,
                 coverage={dn: len(set(panel) & donor_set(dn)) for dn in DONORS})
        rec["draws"][key] = d
        earlier |= set(panel)
        DRAWS.write_text(json.dumps(rec, indent=1))                # at draw time, one draw at a time
        print(f"[{key}] {rule} pool {len(eligible_with)} -> {len(eligible)} essential-free; seed {seed}; sha {d['sha256'][:16]}; "
              f"coverage {d['coverage']}; first {panel[:4]}", flush=True)

    # ------------------------------------------------------------------ bundles (bundle238's recipe)
    for key, rule, pool, seed in SPECS:
        out = HERE / key / "bundle"
        assert not out.exists(), out
        out.mkdir(parents=True)
        panel = rec["draws"][key]["targets"]
        pd.DataFrame({"target_gene": panel}).to_csv(out / "pert_counts.csv", index=False)
        shutil.copyfile(SRC / "gene_names.csv", out / "gene_names.csv")
        man = json.loads((SRC / "manifest.json").read_text())
        man["partition"] = "test"
        man["panel_id"] = f"vcc2026-synthetic-r7a-552-{key}"
        man["contexts"] = [base.RENAME[c] for c in man["contexts"]]
        man["per_context"] = {base.RENAME[k]: v for k, v in man["per_context"].items()}
        (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
        ref = E199 / "work" / ("def-bundle" if pool == "union" else "def-bundle-r")
        files_sha = {}
        for old, new in base.RENAME.items():
            f = out / f"context_{new}.h5ad"
            r = base.relabel(SRC / f"context_{old}.h5ad", f, new)
            want = base.sha256(ref / f.name)
            assert r["sha256"] == want, f"{key}: {f.name} differs from 199's relabelled file"
            os.remove(f)                                           # storage: the same bytes as an APFS clone
            subprocess.run(["cp", "-c", str(ref / f.name), str(f)], check=True)
            assert base.sha256(f) == want, f"{key}: clone of {f.name} differs"
            files_sha[f.name] = want
        assert base.sha256(out / "gene_names.csv") == base.sha256(ref / "gene_names.csv")
        assert [str(t) for t in pd.read_csv(out / "pert_counts.csv")["target_gene"]] == panel
        rec["draws"][key]["bundle"] = dict(path=str(out.relative_to(ROOT)), control_files_identical_to_199=files_sha,
                                           reference=str(ref.relative_to(ROOT)),
                                           files={f: base.sha256(out / f) for f in ("pert_counts.csv", "gene_names.csv", "manifest.json")})
        print(f"[{key}] bundle {out.relative_to(ROOT)}: controls byte-identical to {ref.name} (3/3), stored as clones", flush=True)
    after = panel_files()
    new = sorted(str(p.relative_to(ROOT)) for p in after if p not in files)
    rec["new_panel_files_after_draw"] = new
    if not all("552-r7-reserve" in p for p in new):                    # recorded, not fatal: the draws are already sealed
        print(f"NOTE: panel files other than R6's appeared during the draw: {new}", flush=True)
    rec["seconds"] = round(time.time() - t0, 1)
    DRAWS.write_text(json.dumps(rec, indent=1))
    print(f"-> {DRAWS.relative_to(ROOT)} ({rec['seconds']} s)")


if __name__ == "__main__":
    main()
