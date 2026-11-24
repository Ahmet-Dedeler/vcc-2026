"""801 -- draw the lab's eight Arc-like panels by the rule "arc800" (experiments/800-checker/DESIGN.md section 3.1), in order,
seeds 800101 ... 800108, and build bundles for the six WORKING panels (L1-L6). L7-L8 are the SEALED reserve: only their
target lists, seeds and sha256 are written (no bundle, no CD4 bank, no operator, no dump, no statistic).

THE RULE (frozen here, before the first draw; the dry run that fixed the fallback is recorded in RESULTS.md section 1):
  pool      official axis  n  (K562-GWPS u Orion HCT116 u Orion HEK293T perts)   [efp's union pool]
            - K562 essential list (data/donors/k562_essential/perts.txt)
            - every target of every `pert_counts.csv` under experiments/ and data/ at draw time (draw534's exclusion)
            - every earlier lab draw (L1 .. L(k-1); this script's own earlier draws)
            - basal expression below A/B/C's panel 5th percentile, where a target's basal expression is its pseudobulk
              control CPM (1e6 x gene total / all counts) in each of A/B/C's context files, averaged over A, B, C
  strata    cell = (coverage class) x (strength bin)
            coverage class: C3 all three cancer donors carry it; K2 K562 + at most one Orion; O2 both Orions, no K562;
                            O1 one Orion, no K562   (a target in none of the three is not in the pool)
            strength bin:   GWCD4i Rest `n_downstream` (data/external/cd4t/DE_stats.suppl_table.csv), in bins cut at the
                            A/B/C panel's quintiles of that value (np.quantile at .2/.4/.6/.8 over the A/B/C targets the
                            atlas carries; bin = np.digitize(x, edges) + 1 -> Q1..Q5); NA = the atlas lacks the target.
                            CD4 is the "pilot line" (Arc used an H1 pilot); it is never a fold.
  quota     each cell gets the A/B/C panel's own count in that cell (the panel has 300 targets, so these ARE A/B/C's
            shares). Capacity per draw = floor(pool cell / 8), so all eight disjoint draws fit. A cell whose A/B/C count
            exceeds its capacity is cut to capacity and the deficit moves, one target at a time, to (1) the same strength
            bin in the sibling class (C3<->K2: both K562-carried; O2<->O1: both K562-lacking), if it has spare capacity,
            else (2) the same class's next strength bin with spare capacity, stronger first (Q4: Q5, Q3, Q2, Q1; never into
            or out of NA). K562 coverage and CD4 coverage are preserved exactly by this order whenever (1) or (2) applies.
  draw      for draw k (seed s_k): rng = np.random.default_rng(s_k); for each cell in sorted order, take `quota` targets
            with rng.choice(sorted eligible cell members, quota, replace=False); panel = sorted(union).
ASSERTED AT DRAW TIME: 300 unique targets; 0 essential; 0 overlap with any excluded panel or earlier lab draw; all on the
official gene axis; every cell's count == its quota. Each draw is written to work/draws.json as it is drawn.

BUNDLES (L1-L6 only; draw534's recipe): pert_counts.csv (target_gene), gene_names.csv (A/B/C's, copied), manifest.json
(partition test, panel_id, contexts D/E/F), and the three relabelled A/B/C control files as APFS clones (cp -c) of 199's
relabelled files (experiments/199-def-runbook/work/def-bundle; 534 and 552 proved them byte-identical to relabel(A/B/C)),
sha256 re-checked after the clone.
SEALED (L7-L8): work/sealed/<L>/pert_counts.csv holds the target list only (one column, no bundle), so any later draw that
excludes "every pert_counts.csv on disk" excludes the sealed reserve too.

    uv run python experiments/801-lab-folds/draw801.py --dry-run    # pool cells, quotas; draws nothing, writes nothing
    uv run python experiments/801-lab-folds/draw801.py              # the draw (once: refuses if work/draws.json exists)
Light: < 1 GB, seconds (the context files are streamed in 8M-entry chunks).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import collections  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
SRC = ROOT / "data/challenge"
DEFB = ROOT / "experiments/199-def-runbook/work/def-bundle"
CD4_TABLE = ROOT / "data/external/cd4t/DE_stats.suppl_table.csv"
CD4_H5AD = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
RENAME = {"A": "D", "B": "E", "C": "F"}
DRAWS = WORK / "draws.json"
SPECS = tuple((f"L{i}", 800100 + i, "working" if i <= 6 else "sealed") for i in range(1, 9))
N_DRAWS = len(SPECS)
CLASSES = ("C3", "K2", "O2", "O1")
SIBLING = {"C3": "K2", "K2": "C3", "O2": "O1", "O1": "O2"}
QBINS = ("Q1", "Q2", "Q3", "Q4", "Q5")
BASAL_PCT = 5.0
CHUNK = 8 << 20


def sha_list(xs) -> str:
    return hashlib.sha256("\n".join(sorted(xs)).encode()).hexdigest()


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def lines(p: Path) -> set[str]:
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def panel_files() -> list[Path]:
    """`find experiments data -name pert_counts.csv`, run at draw time (draw534's exclusion)."""
    r = subprocess.run(["find", "experiments", "data", "-name", "pert_counts.csv"], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return sorted(ROOT / ln.strip() for ln in r.stdout.splitlines() if ln.strip())


def read_targets(p: Path) -> list[str]:
    df = pd.read_csv(p)
    col = "target_gene" if "target_gene" in df.columns else df.columns[0]      # conservative: exclude, never skip
    return [str(t) for t in df[col].dropna()]


def context_cpm(path: Path, n_genes: int) -> np.ndarray:
    """Pseudobulk CPM of one context file: 1e6 x per-gene total / all counts, streamed (never loads X)."""
    tot = np.zeros(n_genes)
    with h5py.File(path, "r") as h:
        idx, dat = h["X/indices"], h["X/data"]
        n = idx.shape[0]
        for a in range(0, n, CHUNK):
            tot += np.bincount(idx[a:a + CHUNK], weights=dat[a:a + CHUNK].astype(np.float64), minlength=n_genes)
    return 1e6 * tot / tot.sum()


def cd4_strength() -> dict[str, float]:
    df = pd.read_csv(CD4_TABLE)
    r = df[df["culture_condition"] == "Rest"]
    assert r["target_contrast_gene_name"].is_unique, "one Rest row per target"
    return {str(t): float(v) for t, v in zip(r["target_contrast_gene_name"], r["n_downstream"]) if np.isfinite(v)}


def cd4_h5ad_targets() -> set[str]:
    with h5py.File(CD4_H5AD, "r") as h:
        node = h["obs"]["target_contrast_gene_name"]
        if isinstance(node, h5py.Group):
            return {c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]}
        return {v.decode() if isinstance(v, bytes) else str(v) for v in node[:]}


class Rule:
    """Everything the rule needs, computed from the inputs; no panel draw involved."""

    def __init__(self):
        self.genes = [str(g) for g in pd.read_csv(SRC / "gene_names.csv")["gene_name"]]
        self.gene_set = set(self.genes)
        self.pos = {g: i for i, g in enumerate(self.genes)}
        self.donor = {d: lines(ROOT / f"data/donors/{d}/perts.txt") for d in DONORS}
        self.ess = lines(ESSENTIAL_LIST)
        self.nd = cd4_strength()
        self.abc = [str(t) for t in pd.read_csv(SRC / "pert_counts.csv")["target_gene"]]
        assert len(set(self.abc)) == 300
        self.cpm_ctx = {c: context_cpm(SRC / f"context_{c}.h5ad", len(self.genes)) for c in "ABC"}
        self.cpm = np.mean([self.cpm_ctx[c] for c in "ABC"], axis=0)
        abc_cpm = np.asarray([self.cpm[self.pos[t]] for t in self.abc])
        self.basal_thr = float(np.percentile(abc_cpm, BASAL_PCT))
        v = np.asarray([self.nd.get(t, np.nan) for t in self.abc])
        self.edges = np.quantile(v[np.isfinite(v)], [0.2, 0.4, 0.6, 0.8])
        self.files = panel_files()
        self.excluded, self.panels = set(), {}
        for p in self.files:
            t = read_targets(p)
            self.excluded |= set(t)
            self.panels.setdefault(sha_list(set(t)), []).append(str(p.relative_to(ROOT)))
        self.union = set().union(*self.donor.values())

    def cov(self, t: str) -> str | None:
        k, h, e = (t in self.donor[d] for d in DONORS)
        if k:
            return "C3" if (h and e) else "K2"
        if h and e:
            return "O2"
        return "O1" if (h or e) else None

    def sbin(self, t: str) -> str:
        x = self.nd.get(t)
        return "NA" if x is None else f"Q{int(np.digitize(x, self.edges)) + 1}"

    def cell(self, t: str):
        return (self.cov(t), self.sbin(t))

    def pool(self, earlier: set[str]) -> list[str]:
        return sorted(t for t in (self.gene_set & self.union) - self.ess - self.excluded - earlier
                      if self.cpm[self.pos[t]] >= self.basal_thr)

    def quotas(self, pool: list[str]):
        abc = collections.Counter(self.cell(t) for t in self.abc)
        assert None not in {c for c, _ in abc}, "an A/B/C target is carried by none of the three cancer donors"
        cap_n = collections.Counter(self.cell(t) for t in pool)
        cells = [(c, b) for c in CLASSES for b in QBINS + ("NA",)]
        cap = {k: cap_n.get(k, 0) // N_DRAWS for k in cells}
        q = {k: min(abc.get(k, 0), cap[k]) for k in cells}
        moves = []
        for k in cells:
            deficit = abc.get(k, 0) - q[k]
            c, b = k
            for _ in range(deficit):
                order = [(SIBLING[c], b)]
                if b != "NA":
                    i = QBINS.index(b)
                    order += [(c, QBINS[j]) for j in list(range(i + 1, 5)) + list(range(i - 1, -1, -1))]
                dest = next((d for d in order if q[d] < cap[d]), None)
                assert dest is not None, f"no spare capacity for a deficit in {k}"
                q[dest] += 1
                moves.append(dict(frm=list(k), to=list(dest)))
        assert sum(q.values()) == 300
        return dict(abc={f"{c}|{b}": abc.get((c, b), 0) for c, b in cells},
                    pool={f"{c}|{b}": cap_n.get((c, b), 0) for c, b in cells},
                    capacity={f"{c}|{b}": cap[(c, b)] for c, b in cells},
                    quota={f"{c}|{b}": q[(c, b)] for c, b in cells}, moves=moves), q


def bundle(key: str, panel: list[str], rec: dict) -> dict:
    out = WORK / key / "bundle"
    assert not out.exists(), out
    out.mkdir(parents=True)
    pd.DataFrame({"target_gene": panel}).to_csv(out / "pert_counts.csv", index=False)
    shutil.copyfile(SRC / "gene_names.csv", out / "gene_names.csv")
    man = json.loads((SRC / "manifest.json").read_text())
    man["partition"] = "test"
    man["panel_id"] = f"vcc2026-synthetic-lab801-{key}"
    man["contexts"] = [RENAME[c] for c in man["contexts"]]
    man["per_context"] = {RENAME[k]: v for k, v in man["per_context"].items()}
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    files = {}
    for new in RENAME.values():
        ref = DEFB / f"context_{new}.h5ad"
        want = sha256(ref)
        subprocess.run(["cp", "-c", str(ref), str(out / ref.name)], check=True)
        assert sha256(out / ref.name) == want, f"{key}: clone of {ref.name} differs"
        files[ref.name] = want
    assert sha256(out / "gene_names.csv") == sha256(DEFB / "gene_names.csv")
    assert [str(t) for t in pd.read_csv(out / "pert_counts.csv")["target_gene"]] == panel
    return dict(path=str(out.relative_to(ROOT)), control_files_clone_of=str(DEFB.relative_to(ROOT)), control_sha256=files,
                files={f: sha256(out / f) for f in ("pert_counts.csv", "gene_names.csv", "manifest.json")})


def main(dry: bool) -> None:
    t0 = time.time()
    R = Rule()
    h5t = cd4_h5ad_targets()
    base = dict(
        rule="arc800 (draw801.py docstring; experiments/800-checker/DESIGN.md 3.1)", script_sha256=sha256(Path(__file__)),
        n_panel_files=len(R.files), n_distinct_panels=len(R.panels), panel_files={k[:16]: v for k, v in R.panels.items()},
        n_excluded_targets=len(R.excluded), essential_list_sha256=sha256(ESSENTIAL_LIST), essential_list_size=len(R.ess),
        donor_perts_sha256={d: sha256(ROOT / f"data/donors/{d}/perts.txt") for d in DONORS},
        cd4_table_sha256=sha256(CD4_TABLE), cd4_rest_targets=len(R.nd), cd4_h5ad_targets=len(h5t),
        cd4_table_minus_h5ad=len(set(R.nd) - h5t), cd4_h5ad_minus_table=len(h5t - set(R.nd)),
        strength_edges=R.edges.tolist(), basal_threshold_cpm=R.basal_thr, basal_pct=BASAL_PCT,
        context_sha256={c: sha256(SRC / f"context_{c}.h5ad") for c in "ABC"},
        abc_panel_sha256=sha_list(R.abc), union_pool_before_exclusion=len(R.gene_set & R.union))
    pool0 = R.pool(set())
    qrec, q = R.quotas(pool0)
    base.update(pool_size_first_draw=len(pool0), strata=qrec)
    print(f"{len(R.files)} pert_counts.csv files, {len(R.panels)} distinct panels, {len(R.excluded)} excluded targets; "
          f"basal thr {R.basal_thr:.4f} CPM; strength edges {R.edges.tolist()}; pool {len(pool0)}", flush=True)
    for k in qrec["quota"]:
        if qrec["abc"][k] or qrec["quota"][k]:
            print(f"   {k:6s} abc {qrec['abc'][k]:3d}  pool {qrec['pool'][k]:4d}  cap {qrec['capacity'][k]:3d}  quota {qrec['quota'][k]:3d}")
    print(f"   moves: {qrec['moves']}", flush=True)
    if dry:
        print("dry run: nothing drawn, nothing written")
        return
    assert not DRAWS.exists(), f"{DRAWS} exists: the lab panels are drawn once"
    WORK.mkdir(parents=True, exist_ok=True)
    rec = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **base, draws={})
    earlier: set[str] = set()
    cells_sorted = sorted(q)
    for key, seed, status in SPECS:
        pool = R.pool(earlier)
        byc = collections.defaultdict(list)
        for t in pool:
            byc[R.cell(t)].append(t)
        rng = np.random.default_rng(seed)
        panel = []
        for k in cells_sorted:
            if q[k]:
                panel += rng.choice(np.asarray(sorted(byc[k])), size=q[k], replace=False).tolist()
        panel = sorted(panel)
        assert len(panel) == 300 and len(set(panel)) == 300
        assert not (set(panel) & R.ess), f"{key}: essential targets"
        assert not (set(panel) & R.excluded), f"{key}: overlaps an excluded panel"
        assert not (set(panel) & earlier), f"{key}: overlaps an earlier lab draw"
        assert set(panel) <= R.gene_set, f"{key}: a target is off the gene axis"
        got = collections.Counter(R.cell(t) for t in panel)
        assert all(got.get(k, 0) == q[k] for k in q), f"{key}: strata counts differ from the quotas"
        d = dict(seed=seed, status=status, targets=panel, sha256=sha_list(panel), n=300, overlap_excluded=0,
                 overlap_essential=0, overlap_earlier_lab=0)
        if status == "working":
            d["pool_size"] = len(pool)
        rec["draws"][key] = d
        earlier |= set(panel)
        DRAWS.write_text(json.dumps(rec, indent=1))                # at draw time, one draw at a time
        print(f"[{key}] {status} seed {seed}; sha {d['sha256'][:16]}", flush=True)
    for key, seed, status in SPECS:
        panel = rec["draws"][key]["targets"]
        if status == "working":
            rec["draws"][key]["bundle"] = bundle(key, panel, rec)
            print(f"[{key}] bundle {rec['draws'][key]['bundle']['path']}: controls cloned from 199's def-bundle (3/3)", flush=True)
        else:
            sd = WORK / "sealed" / key
            sd.mkdir(parents=True)
            pd.DataFrame({"target_gene": panel}).to_csv(sd / "pert_counts.csv", index=False)
            rec["draws"][key]["target_list"] = dict(path=str((sd / "pert_counts.csv").relative_to(ROOT)),
                                                    sha256=sha256(sd / "pert_counts.csv"))
    after = panel_files()
    rec["new_panel_files_after_draw"] = sorted(str(p.relative_to(ROOT)) for p in after if p not in R.files)
    rec["seconds"] = round(time.time() - t0, 1)
    DRAWS.write_text(json.dumps(rec, indent=1))
    print(f"-> {DRAWS.relative_to(ROOT)} sha256 {sha256(DRAWS)} ({rec['seconds']} s)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    main(ap.parse_args().dry_run)
