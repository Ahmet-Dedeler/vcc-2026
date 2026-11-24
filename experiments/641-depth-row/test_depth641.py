#!/usr/bin/env python3
"""641 -- tests for the depth-row patch. Temp trees only: the live files are hashed at the start and the end and must not
change. Trees are 554's (`test_guard554.make_tree`: copies + 538's applier + 554's patch; data/ and submissions/ are
read-only symlinks), then 641 on top.

  A  patch mechanics: --check clean and writes nothing; apply; a second apply refused; a tree WITHOUT 554 refused (641
     anchors on 554's lines); a partial tree refused; an anchor miss writes nothing; `bash -n` / `py_compile` after.
  B  the row through the real CLI on small synthetic bundles (data/challenge's gene axis and panel, 240 control cells per
     context, controlled library medians): 20.0k / 17.5k / 23.0k -> no row (text identical to the 554 tree's); 29.5k
     (an un-downsampled release) and 16.0k -> one FAIL row per context; DEF_DAY_ALLOW_DEPTH=1 -> passing ALLOWED rows;
     the JSON always carries library_median_ratio.
  C  the two def_day.sh lines, executed in bash with a stub `log`: the override is logged and recorded only when set.
  D  real data (slow, ~2-4 min, <= 0.5 GB): def_preflight base (538+554) vs patched (538+554+641), knobs unset, on 554's six
     bundles: identical text and exit code; every context's ratio recorded. R7b's as-built Calu-3 controls (median 11,418)
     as every context: three FAIL depth rows (beside 554's spread rows), ALLOWED with DEF_DAY_ALLOW_DEPTH=1.

    uv run python experiments/641-depth-row/test_depth641.py -v        (D is skipped with TEST641_FAST=1)
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/554-dose-guard"))
sys.path.insert(0, str(HERE))
import test_guard554 as G  # noqa: E402
import patch_depth641 as D  # noqa: E402

PY = G.PY
PATCH = HERE / "patch_depth641.py"
LIVE = sorted({*G.LIVE, ROOT / "scripts/def_preflight.py", ROOT / "scripts/def_day.sh"})
BUNDLES = ["data/challenge", "experiments/180-def-readiness/work/def-bundle", "experiments/238-essential-free-panel/work/ef-bundle",
           "experiments/238-essential-free-panel/work/ef-bundle-r", "experiments/534-generalization-audit/work/partB/S1/bundle",
           "experiments/552-r7-reserve/work/T1/bundle"]
MASK = [re.compile(r"^\d+(\.\d+)?s$"), re.compile(r"^wrote ")]
RESULTS: dict = {}


def live_hashes() -> dict:
    return {str(p): G.sha(p) for p in LIVE}


def clean_env(**extra) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    env.update(extra)
    return env


def run_pre(tree: Path, bundle: Path, out: Path, env: dict) -> tuple[int, str, dict | None]:
    p = subprocess.run([str(PY), str(tree / "scripts/def_preflight.py"), "--bundle", str(bundle), "--out", str(out)],
                       cwd=tree, env=env, capture_output=True, text=True)
    text = "\n".join(ln for ln in p.stdout.splitlines() if not any(m.search(ln.strip()) for m in MASK))
    return p.returncode, text, (json.loads(out.read_text()) if out.exists() else None)


def apply641(tree: Path) -> subprocess.CompletedProcess:
    return subprocess.run([str(PY), str(PATCH), "--root", str(tree)], capture_output=True, text=True)


def depth_rows(js: dict) -> list[dict]:
    return [c for c in js["checks"] if "641" in c["check"]]


def synth_bundle(dest: Path, median: float, n_cells: int = 240, seed: int = 641) -> Path:
    """data/challenge's manifest / gene axis / panel, with three small synthetic control files whose per-cell library
    sizes have the given median (a log-normal spread like A/B/C's; integer counts on ~6k genes per cell)."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    dest.mkdir(parents=True)
    man = json.loads((ROOT / "data/challenge/manifest.json").read_text())
    genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    rng = np.random.default_rng(seed)
    for c in man["contexts"]:
        lib = np.exp(rng.normal(np.log(median), 0.45, n_cells))
        lib *= median / np.median(lib)                    # exact median before rounding
        rows, cols, vals = [], [], []
        p = rng.dirichlet(np.full(len(genes), 0.05))
        for i, L in enumerate(lib):
            x = rng.multinomial(int(round(L)), p)
            nz = np.flatnonzero(x)
            rows.append(np.full(nz.size, i)), cols.append(nz), vals.append(x[nz])
        X = sp.csr_matrix((np.concatenate(vals).astype(np.float32), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n_cells, len(genes)))
        X.sort_indices()
        n_ntc = man["per_context"][c]["n_ntc_ids"]
        obs = pd.DataFrame({"context": pd.Categorical([c] * n_cells), "target_gene": pd.Categorical(["non-targeting"] * n_cells),
                            "ntc_id": pd.Categorical([f"non-targeting-{i % n_ntc}" for i in range(n_cells)])},
                           index=[f"cell{i}" for i in range(n_cells)])
        ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes)).write_h5ad(dest / f"context_{c}.h5ad")
        man["per_context"][c]["control_cells"] = n_cells
    (dest / "manifest.json").write_text(json.dumps(man))
    for f in ("gene_names.csv", "pert_counts.csv"):
        shutil.copy(ROOT / "data/challenge" / f, dest / f)
    return dest


class Depth641(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if D.MARK in (ROOT / "scripts/def_preflight.py").read_text():
            raise unittest.SkipTest("641 is applied on HEAD (10-05 / 10-06 apply days): this suite validates the patch against a PRE-641 tree, which cannot be built from live files any more. Run it in a worktree at the pre-apply commit with the gitignored inputs linked (606 RESULTS, the P1/P2 worktree proof); post-apply regression = experiments/644-runbook-amendments/verify644.py")
        cls.live = live_hashes()
        cls.tmp = Path(tempfile.mkdtemp(prefix="test641-"))
        cls.base = G.make_tree(cls.tmp / "base", with554=True)          # 538 + 554
        cls.tree = G.make_tree(cls.tmp / "patched", with554=True)       # 538 + 554 + 641
        r = apply641(cls.tree)
        assert r.returncode == 0, r.stdout + r.stderr

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)
        assert live_hashes() == cls.live, "a live file changed"
        (HERE / "work").mkdir(exist_ok=True)
        (HERE / "work/test641.json").write_text(json.dumps(RESULTS, indent=1, default=str))

    # ---------------------------------------------------------------------------------------------- A: mechanics
    def test_A_mechanics(self):
        t = G.make_tree(self.tmp / "mech", with554=True)
        before = {f: G.sha(t / f) for f in D.FILES}
        r = subprocess.run([str(PY), str(PATCH), "--root", str(t), "--check"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("apply cleanly", r.stdout)
        self.assertEqual(before, {f: G.sha(t / f) for f in D.FILES})                     # --check writes nothing
        self.assertEqual(apply641(t).returncode, 0)
        after = {f: G.sha(t / f) for f in D.FILES}
        self.assertTrue(all(after[f] != before[f] for f in D.FILES))
        r2 = apply641(t)                                                                # second apply refused
        self.assertNotEqual(r2.returncode, 0)
        self.assertIn("refusing a second apply", r2.stderr + r2.stdout)
        self.assertEqual(after, {f: G.sha(t / f) for f in D.FILES})
        self.assertEqual(subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")]).returncode, 0)
        self.assertEqual(subprocess.run([str(PY), "-m", "py_compile", str(t / "scripts/def_preflight.py")]).returncode, 0)
        # without 554: refused, nothing written
        t0 = G.make_tree(self.tmp / "no554", with554=False)
        b0 = {f: G.sha(t0 / f) for f in D.FILES}
        r3 = apply641(t0)
        self.assertNotEqual(r3.returncode, 0)
        self.assertIn("554 is not applied", r3.stderr + r3.stdout)
        self.assertEqual(b0, {f: G.sha(t0 / f) for f in D.FILES})
        # partial: only def_preflight carries the mark
        tp = G.make_tree(self.tmp / "partial", with554=True)
        new = D.plan(tp)
        (tp / D.FILES[0]).write_text(new[D.FILES[0]])
        bp = {f: G.sha(tp / f) for f in D.FILES}
        r4 = apply641(tp)
        self.assertNotEqual(r4.returncode, 0)
        self.assertIn("partially patched", r4.stderr + r4.stdout)
        self.assertEqual(bp, {f: G.sha(tp / f) for f in D.FILES})
        # anchor miss (the 554 env line edited away): nothing written to either file
        tm = G.make_tree(self.tmp / "miss", with554=True)
        txt = (tm / D.FILES[1]).read_text().replace(D.DAY_ENV_A, D.DAY_ENV_A.replace("LOW_DOSE=1\\n", "LOW_DOSE=1 \\n"))
        (tm / D.FILES[1]).write_text(txt)
        bm = {f: G.sha(tm / f) for f in D.FILES}
        r5 = apply641(tm)
        self.assertNotEqual(r5.returncode, 0)
        self.assertIn("anchor found 0 times", r5.stderr + r5.stdout)
        self.assertEqual(bm, {f: G.sha(tm / f) for f in D.FILES})
        RESULTS["A"] = "PASS"

    # ---------------------------------------------------------------------------------------------- B: synthetic bundles
    def test_B_synthetic(self):
        rec = {}
        for tag, med, expect_rows in (("20.0k", 20030.0, 0), ("17.5k", 17500.0, 0), ("23.0k", 23000.0, 0),
                                      ("29.5k", 29482.0, 3), ("16.0k", 16000.0, 3)):
            b = synth_bundle(self.tmp / f"syn_{tag}", med)
            rb, tb, jb = run_pre(self.base, b, self.tmp / f"syn_{tag}_base.json", clean_env())
            rp, tp, jp = run_pre(self.tree, b, self.tmp / f"syn_{tag}_pat.json", clean_env())
            rows = depth_rows(jp)
            ratios = {c: r["library_median_ratio"] for c, r in jp["control_files"].items()}
            self.assertEqual(len(rows), expect_rows, (tag, rows))
            self.assertTrue(all(r["pass"] is False for r in rows))
            for c, r in jp["control_files"].items():
                self.assertAlmostEqual(r["library_median_ratio"], round(r["median_library"] / D_REF, 4), places=4)
            if expect_rows == 0:
                self.assertEqual((rb, tb), (rp, tp), tag)                                 # silent on PASS
            else:
                self.assertNotEqual(rp, 0)
                ra, ta, ja = run_pre(self.tree, b, self.tmp / f"syn_{tag}_allow.json", clean_env(DEF_DAY_ALLOW_DEPTH="1"))
                arows = depth_rows(ja)
                self.assertEqual(len(arows), 3)
                self.assertTrue(all(r["pass"] and "ALLOWED by DEF_DAY_ALLOW_DEPTH=1" in r["value"] for r in arows))
                self.assertTrue(all(r.get("library_median_allowed") is True for r in ja["control_files"].values()))
                self.assertEqual(ra, rb, tag)                                             # allowed == the 554 tree's rc
                rec[tag + "_allowed_rc"] = ra
            rec[tag] = dict(rc_base=rb, rc_patched=rp, ratios=ratios, rows=[r["value"] for r in rows])
        RESULTS["B"] = rec

    # ---------------------------------------------------------------------------------------------- C: def_day.sh lines
    def test_C_day_lines(self):
        text = (self.tree / "scripts/def_day.sh").read_text()
        add_log = D.DAY_LOG_B[len(D.DAY_LOG_A):]
        add_env = D.DAY_ENV_B[len(D.DAY_ENV_A):]
        self.assertIn(D.DAY_LOG_B, text)
        self.assertIn(D.DAY_ENV_B, text)
        # the 554 start-log block sits right after the first `log "start` line; ours follows 554's last line
        i_start = text.index('log "start')
        self.assertLess(i_start, text.index(add_log))
        rec = {}
        for allow in ("", "1"):
            w = self.tmp / f"dayW{allow or 0}"
            w.mkdir()
            script = f'set -euo pipefail\nW="{w}"\nlog() {{ printf "LOG %s\\n" "$*"; }}\n{add_log}{add_env}echo END\n'
            env = clean_env(**({"DEF_DAY_ALLOW_DEPTH": "1"} if allow else {}))
            p = subprocess.run(["/bin/bash", "-c", script], env=env, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(p.stdout.rstrip().endswith("END"))
            envf = w / "def_day.final.env"
            logged = "641: the depth row" in p.stdout
            recorded = envf.exists() and "DEF_DAY_ALLOW_DEPTH=1" in envf.read_text()
            self.assertEqual(logged, bool(allow))
            self.assertEqual(recorded, bool(allow))
            rec[f"allow{allow or 0}"] = dict(logged=logged, recorded=recorded)
        RESULTS["C"] = rec

    # ---------------------------------------------------------------------------------------------- D: real data
    @unittest.skipIf(os.environ.get("TEST641_FAST") == "1", "TEST641_FAST=1")
    def test_D_real(self):
        rec = {}
        for b in BUNDLES:
            rb, tb, _ = run_pre(self.base, ROOT / b, self.tmp / "rb.json", clean_env())
            rp, tp, jp = run_pre(self.tree, ROOT / b, self.tmp / "rp.json", clean_env())
            ratios = {c: r["library_median_ratio"] for c, r in jp["control_files"].items()}
            self.assertEqual((rb, tb), (rp, tp), b)
            self.assertEqual(depth_rows(jp), [])
            rec[b] = dict(rc=rp, identical=True, ratios=ratios)
        # R7b's as-built controls as every context (554's construction)
        import anndata as ad
        man = json.loads((ROOT / "data/challenge/manifest.json").read_text())
        wide = self.tmp / "r7b"
        wide.mkdir()
        a = ad.read_h5ad(ROOT / "experiments/552-r7-reserve/work/r7b/bundle/context_R.h5ad")
        for c in man["contexts"]:
            x = a.copy()
            x.obs["context"] = c
            x.obs["target_gene"] = "non-targeting"
            x.write_h5ad(wide / f"context_{c}.h5ad")
            man["per_context"][c]["control_cells"] = int(a.n_obs)
        (wide / "manifest.json").write_text(json.dumps(man))
        for f in ("gene_names.csv", "pert_counts.csv"):
            shutil.copy(ROOT / "data/challenge" / f, wide / f)
        r0, _, j0 = run_pre(self.tree, wide, self.tmp / "w0.json", clean_env())
        r1, _, j1 = run_pre(self.tree, wide, self.tmp / "w1.json", clean_env(DEF_DAY_ALLOW_DEPTH="1"))
        r2, _, j2 = run_pre(self.tree, wide, self.tmp / "w2.json", clean_env(DEF_DAY_ALLOW_DEPTH="1", DEF_DAY_ALLOW_WIDE_LIBS="1"))
        rb, _, jb = run_pre(self.base, wide, self.tmp / "wb.json", clean_env(DEF_DAY_ALLOW_WIDE_LIBS="1"))
        self.assertEqual(len(depth_rows(j0)), 3)
        self.assertTrue(all(not r["pass"] for r in depth_rows(j0)))
        self.assertTrue(all(r["pass"] for r in depth_rows(j1)))
        self.assertNotEqual(r1, 0)                       # 554's spread rows still fail without their own override
        self.assertEqual(r2, rb)                         # both overrides: the 554 tree's result with its override
        rec["r7b_as_every_context"] = dict(rc=r0, rows=[r["value"] for r in depth_rows(j0)],
                                           allowed_rows=[r["value"] for r in depth_rows(j1)], rc_allow_depth=r1,
                                           rc_allow_both=r2, rc_base_allow_wide=rb)
        RESULTS["D"] = rec


D_REF = D.PRE_CONST_B and 20030.0

if __name__ == "__main__":
    unittest.main()
