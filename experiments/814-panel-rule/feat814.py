"""814 -- target-level features for A/B/C's panel, the arc800 pool and L1-L6, and A/B/C's percentile among 1,000
redraws of a draw rule (arc800 = experiments/801-lab-folds/draw801.py, and the v2 variants named in RULES below).

EC2 only (reads the Orion banks, 2.7 GB each, and streams the A/B/C context files). Run from the repo root:
    uv run python experiments/814-panel-rule/feat814.py --rules R0,R1 --n 1000 --out experiments/814-panel-rule/work/run1

POOL (exactly draw801's rule, reconstructed and checked, never re-derived from a fresh `find`):
  excluded targets = the pert_counts.csv files draw801 recorded at its draw time (draws.json `panel_files`), read only to
  exclude; pool0 = (official axis n donor union) - K562 essential - excluded - below the basal floor.
  KNOWN ANSWERS (abort on any failure): n_excluded == draws.json's 4373; basal floor == 13.412203349132618 (bit-equal);
  strength edges, pool size 5125 and every stratum's quota == draws.json; the sequential redraw of L1..L6 with seeds
  800101..800106 reproduces each recorded sha256 (draw801's exact procedure). L7/L8 are never redrawn.
NULL POOL = pool0 minus the sealed L7/L8 target lists (draws.json, used as an exclusion set only; never written, printed or
  featurised). Redraws use seeds 814000 + i (never 800107/800108), the recorded quotas, draw801's per-cell procedure.

RULES (each a pool filter on the null pool; quotas as arc800; a cell short of its quota borrows by draw801's move order):
  R0  arc800 as built
  R1  R0 minus every target of the RPE1 / Jurkat essential-gene screens (data/donors/{rpe1,nadig_jurkat}/perts.txt)
  (later rules are appended here with their date before they are run; every rule run is a try, counted in RESULTS.md)

FEATURES (one row per gene of null pool u A/B/C u L1..L6; sources hashed in features814.json):
  banks k562 / hct / hek (production donors), h1f (VCC-2025 H1, all cells), pan (Pan 2025 H1 hESC, GSE295214 compact
  donor, 115): cov, basal CPM of the target in the bank's control, own-gene log effect (knockdown), off-target effect norm
  ||delta|| on genes > 10 CPM, robust DE count (518's statistic with a BANK-WIDE gene scale from a fixed 2,000-row sample,
  so it is comparable across panels);
  CD4 atlas (GWCD4i): n_downstream / n_total_de per condition, on-target effect, significance, baseMean, guide-level KD;
  DepMap Chronos (360's 17-line subset): mean, min, share < -0.5; A/B/C control CPM (mean, min over contexts);
  membership: RPE1, Jurkat (essential screens), VIPerturb K562, PerturbMatch, TeloHAEC, Pan H1, VCC-2025 H1 (full),
  Feng iPSC genome-wide fitness / non-fitness and targeted screens; n_public = number of those public screens.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import collections  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/801-lab-folds"))
import draw801 as D  # noqa: E402  (module-level code is definitions only)

DRAWS = ROOT / "experiments/801-lab-folds/work/draws.json"
SRC = ROOT / "data/challenge"
BANKS = {"k562": ROOT / "data/donors/k562_gwps", "hct": ROOT / "data/donors/orion_hct116",
         "hek": ROOT / "data/donors/orion_hek293t", "h1f": ROOT / "data/donors/h1_2025_full",
         "pan": ROOT / "experiments/115-gse295214-h1/work/compact-donor"}
LISTS = {"rpe1": ROOT / "data/donors/rpe1/perts.txt", "jurkat": ROOT / "data/donors/nadig_jurkat/perts.txt",
         "vip": ROOT / "data/donors/viperturb_k562/perts.txt", "pmatch": ROOT / "data/external/perturbmatch/perts.txt",
         "telo": ROOT / "data/external/igvf-telohaec/perts.txt"}
FENG = {"feng_fit": ROOT / "data/external/feng-ipsc/GenomeWideScreen_FitnessGenes_Cell-Metadata.tsv.gz",
        "feng_nonfit": ROOT / "data/external/feng-ipsc/GenomeWideScreen_NonFitnessGenes_Cell-Metadata.tsv.gz",
        "feng_tgt": ROOT / "data/external/feng-ipsc/TargetedScreen_Cell-Metadata.tsv.gz"}
CD4_TABLE = ROOT / "data/external/cd4t/DE_stats.suppl_table.csv"
CD4_KD = ROOT / "data/external/cd4t/guide_kd_efficiency.suppl_table.csv"
DEP = ROOT / "experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv"
SEED0 = 814000
FORBIDDEN_SEEDS = {800107, 800108}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def lines(p: Path) -> list[str]:
    return [ln.strip() for ln in Path(p).read_text().splitlines() if ln.strip()]


def say(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------------------------------------------ pool + KA
class Pool:
    def __init__(self):
        rec = json.loads(DRAWS.read_text())
        self.rec = rec
        self.genes = [str(g) for g in pd.read_csv(SRC / "gene_names.csv")["gene_name"]]
        self.gene_set, self.pos = set(self.genes), {g: i for i, g in enumerate(self.genes)}
        self.donor = {d: set(lines(ROOT / f"data/donors/{d}/perts.txt")) for d in D.DONORS}
        self.ess = set(lines(D.ESSENTIAL_LIST))
        self.nd = D.cd4_strength()
        self.abc = [str(t) for t in pd.read_csv(SRC / "pert_counts.csv")["target_gene"]]
        self.cpm_ctx = {c: D.context_cpm(SRC / f"context_{c}.h5ad", len(self.genes)) for c in "ABC"}
        self.cpm = np.mean([self.cpm_ctx[c] for c in "ABC"], axis=0)
        self.basal_thr = float(np.percentile(np.asarray([self.cpm[self.pos[t]] for t in self.abc]), D.BASAL_PCT))
        v = np.asarray([self.nd.get(t, np.nan) for t in self.abc])
        self.edges = np.quantile(v[np.isfinite(v)], [0.2, 0.4, 0.6, 0.8])
        excluded = set()
        for paths in rec["panel_files"].values():
            for p in paths:
                excluded |= set(D.read_targets(ROOT / p))
        self.excluded = excluded
        self.union = set().union(*self.donor.values())
        # KA rows on the rule's inputs
        ka = dict(n_excluded=(len(excluded), rec["n_excluded_targets"]),
                  basal_thr=(self.basal_thr, rec["basal_threshold_cpm"]),
                  edges=(self.edges.tolist(), rec["strength_edges"]),
                  context_sha=({c: sha256(SRC / f"context_{c}.h5ad") for c in "ABC"}, rec["context_sha256"]))
        for k, (a, b) in ka.items():
            assert a == b, f"KA FAIL {k}: {a} != {b}"
        self.pool0 = sorted(t for t in (self.gene_set & self.union) - self.ess - excluded
                            if self.cpm[self.pos[t]] >= self.basal_thr)
        assert len(self.pool0) == rec["pool_size_first_draw"], f"KA FAIL pool {len(self.pool0)}"
        q = rec["strata"]["quota"]
        self.quota = {tuple(k.split("|")): int(n) for k, n in q.items() if int(n)}
        # recompute the quotas from the rule and compare (draw801.Rule.quotas, without its find)
        cap_n = collections.Counter(self.cell(t) for t in self.pool0)
        assert {f"{c}|{b}": cap_n.get((c, b), 0) for c in D.CLASSES for b in D.QBINS + ("NA",)} == rec["strata"]["pool"], "KA FAIL pool cells"
        say(f"KA inputs PASS: excluded {len(excluded)}, floor {self.basal_thr!r}, pool0 {len(self.pool0)}, quotas {sum(self.quota.values())}")
        # KA: redraw L1..L6 sequentially (draw801's procedure) and match the recorded sha256
        earlier: set[str] = set()
        self.L = {}
        for key, seed, status in D.SPECS:
            if status != "working":
                break
            panel = self.draw([t for t in self.pool0 if t not in earlier], seed)
            assert D.sha_list(panel) == rec["draws"][key]["sha256"] == D.sha_list(rec["draws"][key]["targets"]), f"KA FAIL {key}"
            self.L[key] = panel
            earlier |= set(panel)
        say("KA redraw PASS: L1..L6 reproduce their recorded sha256")
        sealed = set()
        for key, seed, status in D.SPECS:
            if status == "sealed":
                sealed |= set(rec["draws"][key]["targets"])          # exclusion only
        self.null_pool = [t for t in self.pool0 if t not in sealed]
        self.n_sealed_removed = len(self.pool0) - len(self.null_pool)
        del sealed

    def cov(self, t):
        k, h, e = (t in self.donor[d] for d in D.DONORS)
        if k:
            return "C3" if (h and e) else "K2"
        if h and e:
            return "O2"
        return "O1" if (h or e) else None

    def sbin(self, t):
        x = self.nd.get(t)
        return "NA" if x is None else f"Q{int(np.digitize(x, self.edges)) + 1}"

    def cell(self, t):
        return (self.cov(t), self.sbin(t))

    def draw(self, pool, seed, quota=None):
        """draw801's per-cell draw: rng(seed); cells in sorted order; rng.choice(sorted members, quota, replace=False)."""
        assert seed not in FORBIDDEN_SEEDS
        quota = quota or self.quota
        byc = collections.defaultdict(list)
        for t in pool:
            byc[self.cell(t)].append(t)
        rng = np.random.default_rng(seed)
        panel = []
        for k in sorted(quota):
            if quota[k]:
                panel += rng.choice(np.asarray(sorted(byc[k])), size=quota[k], replace=False).tolist()
        assert len(set(panel)) == 300
        return sorted(panel)

    def quota_for(self, pool):
        """the recorded quotas; a cell short of members moves its deficit by draw801's order (capacity = cell size)."""
        n = collections.Counter(self.cell(t) for t in pool)
        q = dict(self.quota)
        cells = [(c, b) for c in D.CLASSES for b in D.QBINS + ("NA",)]
        moves = []
        for k in cells:
            have = q.get(k, 0)
            if have <= n.get(k, 0):
                continue
            q[k] = n.get(k, 0)
            c, b = k
            for _ in range(have - q[k]):
                order = [(D.SIBLING[c], b)]
                if b != "NA":
                    i = D.QBINS.index(b)
                    order += [(c, D.QBINS[j]) for j in list(range(i + 1, 5)) + list(range(i - 1, -1, -1))]
                dest = next(d for d in order if q.get(d, 0) < n.get(d, 0))
                q[dest] = q.get(dest, 0) + 1
                moves.append((k, dest))
        assert sum(q.values()) == 300
        return {k: v for k, v in q.items() if v}, moves


# ------------------------------------------------------------------------------------------------------ features
def bank_features(name: str, path: Path, genes_out: list[str]) -> pd.DataFrame:
    perts = lines(path / "perts.txt")
    bgenes = lines(path / "genes.txt")
    b = np.load(path / "b_ctrl.npy").astype(np.float64)
    ex = np.expm1(b)
    cpm = 1e6 * ex / ex.sum()
    gpos = {g: i for i, g in enumerate(bgenes)}
    keep = np.where(cpm > 10.0)[0]
    delta = np.load(path / "delta.npy", mmap_mode="r")
    pi = {p: i for i, p in enumerate(perts)}
    want = [t for t in genes_out if t in pi]
    rows = np.asarray(sorted(pi[t] for t in want))
    # bank-wide gene scale from a fixed sample of 2,000 rows (or all), independent of any panel
    rng = np.random.default_rng(814)
    samp = np.sort(rng.choice(len(perts), size=min(2000, len(perts)), replace=False))
    Rs = np.asarray(delta[samp][:, keep], np.float64)
    med = np.median(Rs, 0)
    scale = 1.4826 * np.median(np.abs(Rs - med), 0)
    ok = scale > 0
    del Rs
    out = {}
    for a in range(0, rows.size, 1000):
        blk = rows[a:a + 1000]
        R = np.asarray(delta[blk], np.float64)
        for j, r in enumerate(blk):
            t = perts[r]
            own = gpos.get(t)
            v = R[j, keep].copy()
            if own is not None:
                v[keep == own] = 0.0
            z = v[ok] / scale[ok]
            noise = np.median(np.abs(z)) / 0.6745
            out[t] = dict(basal_cpm=float(cpm[own]) if own is not None else np.nan,
                          own_lfc=float(R[j, own]) if own is not None else np.nan,
                          norm=float(np.sqrt(np.sum(v * v))),
                          cnt=float(np.sum(np.abs(z) > 4.0 * noise)) if noise > 0 else np.nan)
    df = pd.DataFrame.from_dict(out, orient="index")
    df.columns = [f"{name}_{c}" for c in df.columns]
    df[f"{name}_cov"] = 1.0
    say(f"bank {name}: {len(want)} of {len(genes_out)} rows, {keep.size} genes > 10 CPM, scale ok {int(ok.sum())}")
    return df


def feng_genes(p: Path) -> collections.Counter:
    c = collections.Counter()
    with gzip.open(p, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        j = hdr.index("Guide_Call")
        for ln in fh:
            g = ln.rstrip("\n").split("\t")[j]
            if "_" in g:
                s = g.rsplit("_", 1)[0]
                if s and s != "NonTarget" and not s.startswith("NonTarget"):
                    c[s] += 1
    return c


def build_features(P: Pool, genes_out: list[str]) -> tuple[pd.DataFrame, dict]:
    src = {}
    F = pd.DataFrame(index=pd.Index(genes_out, name="gene"))
    for c in "ABC":
        F[f"abc_cpm_{c}"] = [P.cpm_ctx[c][P.pos[t]] for t in genes_out]
    F["abc_cpm"] = F[[f"abc_cpm_{c}" for c in "ABC"]].mean(1)
    F["abc_cpm_min"] = F[[f"abc_cpm_{c}" for c in "ABC"]].min(1)
    F["abc_logcpm"] = np.log10(F["abc_cpm"] + 1)
    for name, path in BANKS.items():
        df = bank_features(name, path, genes_out)
        F = F.join(df, how="left")
        F[f"{name}_cov"] = F[f"{name}_cov"].fillna(0.0)
        src[name] = {f: sha256(path / f) for f in ("perts.txt", "genes.txt", "b_ctrl.npy", "delta.npy")}
    for name, p in LISTS.items():
        s = set(lines(p))
        F[f"in_{name}"] = [float(t in s) for t in genes_out]
        src[name] = sha256(p)
    for name, p in FENG.items():
        c = feng_genes(p)
        F[f"in_{name}"] = [float(c.get(t, 0) > 0) for t in genes_out]
        F[f"{name}_cells"] = [float(c.get(t, np.nan)) for t in genes_out]
        src[name] = sha256(p)
        say(f"{name}: {len(c)} genes")
    F["in_essscreen"] = ((F["in_rpe1"] + F["in_jurkat"]) > 0).astype(float)
    # CD4 atlas
    cd = pd.read_csv(CD4_TABLE)
    for cond in ("Rest", "Stim8hr", "Stim48hr"):
        r = cd[cd["culture_condition"] == cond].drop_duplicates("target_contrast_gene_name").set_index("target_contrast_gene_name")
        r = r.reindex(genes_out)
        F[f"cd4_{cond}_nds"] = r["n_downstream"].values
        F[f"cd4_{cond}_nde"] = r["n_total_de_genes"].values
        if cond == "Rest":
            F["cd4_kd_eff"] = r["ontarget_effect_size"].values
            F["cd4_kd_sig"] = r["ontarget_significant"].map({True: 1.0, False: 0.0}).values
            F["cd4_base"] = r["target_baseMean"].values
            F["cd4_cov"] = r["n_downstream"].notna().astype(float).values
    F["cd4_max_nds"] = F[["cd4_Rest_nds", "cd4_Stim8hr_nds", "cd4_Stim48hr_nds"]].max(1)
    F["cd4_log_nds"] = np.log1p(F["cd4_Rest_nds"])
    src["cd4_table"] = sha256(CD4_TABLE)
    kd = pd.read_csv(CD4_KD)
    kd = kd[kd["culture_condition"] == "Rest"].copy()
    kd["gene"] = kd.iloc[:, 0].astype(str).str.rsplit("-", n=1).str[0]
    g = kd.groupby("gene")
    F["cd4_guide_tmin"] = g["t_statistic"].min().reindex(genes_out).values
    F["cd4_guide_anysig"] = g["signif_knockdown"].any().astype(float).reindex(genes_out).values
    src["cd4_kd"] = sha256(CD4_KD)
    # DepMap
    dep = pd.read_csv(DEP, index_col=0)
    dep.columns = [re.sub(r"\s*\(\d+\)$", "", c) for c in dep.columns]
    dep = dep.loc[:, ~dep.columns.duplicated()]
    F["dep_mean"] = dep.mean().reindex(genes_out).values
    F["dep_min"] = dep.min().reindex(genes_out).values
    F["dep_frac05"] = (dep < -0.5).mean().reindex(genes_out).values
    F["dep_cov"] = dep.mean().reindex(genes_out).notna().astype(float).values
    src["depmap"] = sha256(DEP)
    for b in ("k562", "hct", "hek", "pan", "h1f"):
        F[f"{b}_lognorm"] = np.log(F[f"{b}_norm"])
        F[f"{b}_logcnt"] = np.log1p(F[f"{b}_cnt"])
    pub = ["k562_cov", "hct_cov", "hek_cov", "pan_cov", "h1f_cov", "cd4_cov", "in_vip", "in_pmatch", "in_telo",
           "in_essscreen", "in_feng_fit", "in_feng_nonfit", "in_feng_tgt"]
    F["n_public"] = F[pub].sum(1)
    F["prod_lognorm_mean"] = F[["k562_lognorm", "hct_lognorm", "hek_lognorm"]].mean(1)
    F["prod_logcnt_mean"] = F[["k562_logcnt", "hct_logcnt", "hek_logcnt"]].mean(1)
    return F, src


# ------------------------------------------------------------------------------------------------------ panel stats
def stat_defs(F: pd.DataFrame):
    flags = [c for c in F.columns if c.startswith("in_") or c.endswith("_cov") or c == "cd4_kd_sig" or c == "cd4_guide_anysig"]
    cont = [c for c in F.columns if c not in flags and not c.startswith("abc_cpm_")]
    return flags, cont


def panel_stats(M: np.ndarray, flags_idx, cont_idx, idx: np.ndarray) -> np.ndarray:
    """idx (k, 300) gene rows -> stats (k, n_stats): flag counts, then per continuous feature: median, q90, q10, mean."""
    V = M[idx]                                                    # (k, 300, n_feat)
    out = [np.nansum(V[:, :, flags_idx], 1)]
    C = V[:, :, cont_idx]
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out += [np.nanmedian(C, 1), np.nanquantile(C, 0.9, axis=1), np.nanquantile(C, 0.1, axis=1), np.nanmean(C, 1)]
    return np.concatenate(out, 1)


def main(rules: list[str], n: int, out: Path):
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    P = Pool()
    rule_pool = {"R0": list(P.null_pool)}
    ess_screen = set(lines(LISTS["rpe1"])) | set(lines(LISTS["jurkat"]))
    rule_pool["R1"] = [t for t in P.null_pool if t not in ess_screen]
    genes_out = sorted(set(P.null_pool) | set(P.abc) | set().union(*P.L.values()))
    F, src = build_features(P, genes_out)
    F["in_null_pool"] = [float(t in set(P.null_pool)) for t in genes_out]
    F["in_abc"] = [float(t in set(P.abc)) for t in genes_out]
    for k, v in P.L.items():
        F[f"in_{k}"] = [float(t in set(v)) for t in genes_out]
    F.to_csv(out / "features814.csv.gz")
    flags, cont = stat_defs(F.drop(columns=["in_null_pool", "in_abc"] + [f"in_{k}" for k in P.L]))
    names = [f"{c}:count" for c in flags] + [f"{c}:{s}" for s in ("med", "q90", "q10", "mean") for c in cont]
    cols = flags + cont
    M = F[cols].to_numpy(np.float64)
    row = {g: i for i, g in enumerate(genes_out)}
    fi, ci = np.arange(len(flags)), np.arange(len(flags), len(cols))
    ref = {"ABC": P.abc, **P.L}
    ref_stats = {k: panel_stats(M, fi, ci, np.asarray([[row[t] for t in v]]))[0] for k, v in ref.items()}
    res = dict(script_sha256=sha256(Path(__file__)), draws_json_sha256=sha256(DRAWS), sources=src, n_draws=n,
               seeds=[SEED0, SEED0 + n - 1], pool0=len(P.pool0), null_pool=len(P.null_pool),
               n_features=len(cols), stat_names=names, rules={})
    for R in rules:
        pool = rule_pool[R]
        quota, moves = P.quota_for(pool)
        idx = np.asarray([[row[t] for t in P.draw(pool, SEED0 + i, quota)] for i in range(n)])
        Snull = panel_stats(M, fi, ci, idx)
        pct = {}
        for k, v in ref_stats.items():
            lt = np.sum(Snull < v[None, :], 0)
            eq = np.sum(Snull == v[None, :], 0)
            ok = np.isfinite(v)
            p = np.where(ok, (lt + 0.5 * eq) / n, np.nan)
            pct[k] = p
        np.savez_compressed(out / f"null814-{R}.npz", stats=Snull, names=np.asarray(names), **{f"ref_{k}": v for k, v in ref_stats.items()})
        res["rules"][R] = dict(pool=len(pool), quota_moves=[[list(a), list(b)] for a, b in moves],
                               abc_pct=dict(zip(names, pct["ABC"].tolist())),
                               L_pct={k: dict(zip(names, pct[k].tolist())) for k in P.L},
                               null_mean=dict(zip(names, np.nanmean(Snull, 0).tolist())),
                               null_sd=dict(zip(names, np.nanstd(Snull, 0).tolist())),
                               abc_value=dict(zip(names, ref_stats["ABC"].tolist())))
        ext = sorted(((abs(p - 0.5), nm, p) for nm, p in zip(names, pct["ABC"]) if np.isfinite(p)), reverse=True)
        say(f"rule {R}: pool {len(pool)}, moves {len(moves)}; ABC most extreme stats:")
        for _, nm, p in ext[:40]:
            j = names.index(nm)
            print(f"   {nm:28s} ABC {ref_stats['ABC'][j]:10.4f}  null {np.nanmean(Snull[:, j]):10.4f} +- {np.nanstd(Snull[:, j]):8.4f}  "
                  f"pct {p:.3f}  | L: " + " ".join(f"{pct[k][j]:.2f}" for k in P.L))
    res["seconds"] = round(time.time() - t0, 1)
    (out / "feat814.json").write_text(json.dumps(res, indent=1, default=float))
    say(f"done {res['seconds']} s -> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", default="R0,R1")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--out", default=str(HERE / "work" / "run1"))
    a = ap.parse_args()
    main(a.rules.split(","), a.n, Path(a.out))
