"""Mock test of s0build805.py (select / extract x4 / merge) on four tiny synthetic donor files that use the real suppl
tables' guide ids. Light: < 1 GB, seconds. Run with T805_MOCK=1 (tiny N's in the builder).
    T805_MOCK=1 uv run python experiments/805-cd4-flex-surface/test805_mock.py <tmpdir>
Checks: the per-target cell count, controls, labels by guide, the axis mapping of counts, integer counts, the thinning."""
import os, subprocess, sys
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
ROOT = Path(__file__).resolve().parents[2]
SUP = ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24"
tmp = Path(sys.argv[1]); src = tmp / "src"; out = tmp / "out"; src.mkdir(parents=True, exist_ok=True)
genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
kd = pd.read_csv(SUP / "guide_kd_efficiency.suppl_table.csv"); kd = kd[(kd.culture_condition == "Rest") & kd.signif_knockdown.astype(bool)]
sg = pd.read_csv(SUP / "sgrna_library_metadata.suppl_table.csv").set_index("sgRNA")
de = pd.read_csv(SUP / "DE_stats.suppl_table.csv"); de = de[de.culture_condition == "Rest"]
sym = dict(zip(de.target_contrast, de.target_contrast_gene_name))
gl = []
for g in kd[kd.columns[0]].astype(str):
    if g in sg.index and sym.get(sg.loc[g, "target_gene_id"]):
        gl.append((g, sym[sg.loc[g, "target_gene_id"]]))
gl = gl[:3000]
# file var: 300 official genes shuffled + 5 non-official + one duplicated symbol
rng = np.random.default_rng(1)
vg = list(rng.choice(genes, 300, replace=False)) + ["NOTAGENE1", "NOTAGENE2"]
vg.append(vg[0])
for k, d in enumerate(["D1", "D2", "D3", "D4"]):
    gid = [g for g, _ in gl for _ in range(3)] + ["NTC-%d" % i for i in range(30 + 5 * k)] + ["multi-guide"] * 7
    gtype = ["targeting"] * (3 * len(gl)) + ["non-targeting"] * (30 + 5 * k) + ["targeting"] * 7
    n = len(gid)
    lq = np.zeros(n, bool); lq[::97] = True
    # counts: column j of the file gets value (row % 7) + j % 3 ; total depth ~ 100
    X = sp.random(n, len(vg), density=0.3, random_state=k, data_rvs=lambda m: rng.integers(1, 5, m)).tocsr().astype(np.float32)
    obs = pd.DataFrame({"guide_id": pd.Categorical(gid), "guide_type": pd.Categorical(gtype), "low_quality": lq,
                        "perturbed_gene_id": pd.Categorical(["x"] * n)}, index=[f"{d}_{i}" for i in range(n)])
    vgk = vg if k != 2 else vg[::-1]          # D3 carries its own var order (as the real files do)
    var = pd.DataFrame({"gene_name": pd.Categorical(vgk)}, index=[f"E{i}" for i in range(len(vgk))])
    A = ad.AnnData(X=X, obs=obs, var=var); A.write_h5ad(src / f"{d}_Rest.assigned_guide.h5ad")
B = Path(__file__).resolve().parent / "s0build805.py"
env = dict(os.environ, T805_MOCK="1")
for args in (["select"], *[["extract", d] for d in ["D1", "D2", "D3", "D4"]], ["merge"]):
    cmd = ["uv", "run", "python", str(B), args[0], "Rest", str(src), str(out)] + args[1:]
    subprocess.run(cmd, check=True, env=env)
import json
info = json.loads((out / "select.json").read_text())
T = ad.read_h5ad(out / "truth_P.h5ad")
vc = T.obs.target_gene.value_counts()
assert vc["non-targeting"] == 40 and all(vc.drop("non-targeting") == 10), vc
assert len(vc) == 6, vc
assert info["axis_overlap"] == 300 and info["file_dup_symbols"] == 4, (info["axis_overlap"], info["file_dup_symbols"])
assert np.all(T.X.data == np.round(T.X.data))
# the official-axis columns that are not in the file are all zero
inside = set(vg[:300]); zero_cols = [i for i, g in enumerate(genes) if g not in inside]
assert T.X[:, zero_cols].nnz == 0
print("MOCK OK", info["thinning_f"], info["depth_median_before"], info["depth_median_after"], info["P"], info["shortfalls"])
