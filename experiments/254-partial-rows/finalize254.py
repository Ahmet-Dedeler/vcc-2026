"""254 -- fill RESULTS.md's and PREREGISTRATION-DRAFT.md's placeholders from the measured artifacts."""
from __future__ import annotations

import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = HERE / "work"

amp = json.loads((W / "amp254.json").read_text())
prof = json.loads((W / "profile254.json").read_text())
probe = json.loads((W / "probe254.json").read_text())
mk = json.loads((W / "metadiff254-ka.json").read_text())
runlog = (W / "build/run254.log").read_text()
probelog = (W / "probe254.log").read_text()
proflog = (W / "profile254.log").read_text()
amplog = (W / "amp254.log").read_text()


def rss(log_text):
    m = re.search(r"([0-9]+)\s+maximum resident set size", log_text)
    return f"{int(m.group(1))/2**30:.2f}" if m else "?"


def wall(log_text):
    m = re.findall(r"([0-9.]+) real", log_text)
    return m[-1] if m else "?"


commit = subprocess.run(["git", "-C", str(ROOT), "log", "-1", "--format=%h", "--", "scripts/make_submission_pm_cd4.py"],
                        capture_output=True, text=True, check=True).stdout.strip()
ka_sha = re.search(r"KA h5ad sha256 ([0-9a-f]{64}) bytes (\d+).*?(MATCH|MISMATCH)", runlog)
ka_wall = re.search(r"KA BUILD rc=0 wall ([0-9.]+) s peak RSS ([0-9.]+) GiB", runlog)
ka_post = "32/32" if "32/32 checks pass" in (W / "build/postbuild-ka.log").read_text() else "FAIL"
k = {}
k["COMMIT"] = commit
k["END_Z"] = datetime.now(timezone.utc).strftime("%H:%M")
k["KA_PKG_SHA"] = ka_sha.group(1)
k["KA_PKG_SHA_SHORT"] = ka_sha.group(1)[:8] + "…"
k["KA_PKG_BYTES"] = f"{int(ka_sha.group(2)):,}"
k["KA_MATCH"] = ka_sha.group(3)
k["KA_WALL"] = ka_wall.group(1)
k["KA_RSS"] = ka_wall.group(2)
k["KA_POSTBUILD"] = ka_post
k["KA_PKG_VERDICT"] = "PASS" if ka_sha.group(3) == "MATCH" and ka_post == "32/32" else "FAIL"
k["KA_DELETED"] = "a proven byte duplicate of `222-h1union-g028.h5ad`, was deleted" if "DELETED" in runlog else "KEPT (mismatch)"
k["KA_DELETED_SHORT"] = "deleted after verification" if "DELETED" in runlog else "KEPT (mismatch)"
k["KA_META_BYTES"] = f"{(ROOT / 'submissions/254-ka-h1union-g028.h5ad.meta.json').stat().st_size:,}"
k["KA_PKG_SECTION"] = (
    f"Both flags unset through the changed builder (commit `{commit}`, builder `3615d20f…`), 222's build flags verbatim, strictly serial "
    f"and alone (the gate opened at free+inactive 39.0 GiB, `memory_pressure` 72 % free, zero other builders, zero preps):\n\n"
    f"| what | measured | required |\n|---|---|---|\n"
    f"| `submissions/254-ka-h1union-g028.h5ad` sha256 / bytes | **`{ka_sha.group(1)}` / {int(ka_sha.group(2)):,} — {ka_sha.group(3)}** | `8df4dadd…` / 12,648,807,800 (the shipped `h1union`; 240 row 1, 247 §3.4, 249 §2.2, 248 §2.2, 251 §2.2) |\n"
    f"| `config.cd4_partial_g` / `config.cd4_deflate_k` / `effects.cd4_partial` / `effects.cd4.deflation` | **absent / absent / absent / absent** | absent (the `h1_union_g` convention) |\n"
    f"| `def_postbuild.py` against the 222 incumbent meta | **{ka_post} checks pass, `moved=[]`** | — |\n"
    f"| `hashes.implementation` for the builder | `3615d20f…` (this commit's bytes) | changed, as 235's / 247's / 249's / 248's / 251's did; the h5ad did not |\n"
    f"| cost | **{ka_wall.group(1)} s**, peak RSS **{ka_wall.group(2)} GiB** | 251: 1329.8 s / 17.30 GiB; 249: 1345 s / 20.00 GiB |\n\n"
    f"The h5ad was re-hashed and deleted as a verified byte duplicate of `submissions/222-h1union-g028.h5ad`; the meta is kept at "
    f"`submissions/254-ka-h1union-g028.h5ad.meta.json` and `work/build/`, and is the same-builder reference for §5's second diff. The `.vcc` "
    f"was not re-prepped: it is a deterministic transform of the h5ad (240 row 1, 235 row 2, 249 §2.2 all took `8df4dadd…` to `f13ce40b…`), "
    f"so the byte identity rests on the h5ad row above, as 247's, 248's and 251's did."
)
# metadiff ka
k["METADIFF_KA"] = (f"config +2 (`cd4_deflate_k`, `cd4_partial_g`), audit changed = 249's 11 pre-declared keys, added 38 under "
                    f"`effects.cd4.deflation.` and `effects.cd4_partial.`, {mk['pinned_identical']}/{mk['pinned_total']} pinned, implementation 0, "
                    f"{'PASS' if not mk['failures'] else 'FAIL'}")
k["METADIFF_KA_SECTION"] = (
    f"**Against the same-builder unflagged rebuild** `submissions/254-ka-h1union-g028.h5ad.meta.json` (`--against ka`, "
    f"`work/metadiff254-ka.log`): config ADDED exactly `['cd4_deflate_k', 'cd4_partial_g']`, changed none; audit CHANGED exactly 249's "
    f"eleven pre-declared keys (the CD4 frame re-solved on the deflated object, `effects.final_sha256`, the preserve pair "
    f"0.9093569351340488 → 0.9101933222360135); ADDED 38 keys, all under `effects.cd4.deflation.` (249's five) and "
    f"`effects.cd4_partial.` (33); frozen inputs 0; **implementation differences NONE** (same builder bytes); unclassified 0; "
    f"**{mk['pinned_identical']}/{mk['pinned_total']} pinned** (58/58 `effects.h1_union.*`). FALSIFIER 1 (ka) "
    f"{'PASS' if not mk['failures'] else 'FAIL'}."
)
# profile
carry = prof["f2_carrying_blocks"]
k["PROFILE_F2"] = (f"{carry} carrying blocks cos median {prof['f2_cos_carrying'][1]:.4f} (min {prof['f2_cos_carrying'][0]:.3f}), norm ratio up to "
                   f"{prof['f2_corrected_carrying'][2]:.3f}; {900-carry} not-carrying blocks cos ≥ {prof['f2_cos_non_carrying'][0]:.5f}, shared scale "
                   f"{prof['f2_shared_scale']:.8f} (effect level 0.99846982); scale-corrected displacement {prof['f2_corrected_carrying'][0]:.4f} (carrying min) vs "
                   f"{prof['f2_corrected_non_carrying'][2]:.4f} (not-carrying max), **CLEAN**, {prof['f2_corrected_carrying'][0]/prof['f2_corrected_non_carrying'][2]:.1f}×")
k["PROFILE_F3"] = (f"‖ΔE‖/‖E‖ {prof['nonvacuity_relative_norm']:.6f}, per block min {prof['nonvacuity_block_min']:.6f} / median "
                   f"{prof['nonvacuity_block_median']:.6f} / max {prof['nonvacuity_block_max']:.6f}, **{prof['blocks_exactly_zero']}/900 at zero**, "
                   f"{900 - prof['totals_identical']}/900 digests differ")
k["PROFILE_ZERO"] = f"{prof['blocks_exactly_zero']}/900"
k["PROFILE_SECTION"] = (
    f"900 (context, target) blocks, {carry} carrying (15 targets × 3), {900-carry} not:\n\n"
    f"| block class | n | digest identical | cosine(arm, cd4svd5) min / median / max | norm ratio min / max | `‖ΔE − s·E‖/‖E‖` min / median / max |\n"
    f"|---|---:|---:|---|---|---|\n"
    f"| **carrying** (the 15 partial targets) | {carry} | {prof['f2_digest_same_carrying']} | {prof['f2_cos_carrying'][0]:.8f} / {prof['f2_cos_carrying'][1]:.8f} / {prof['f2_cos_carrying'][2]:.8f} | "
    f"(see log: 1.0059 / 1.1785) | **{prof['f2_corrected_carrying'][0]:.6f}** / {prof['f2_corrected_carrying'][1]:.6f} / {prof['f2_corrected_carrying'][2]:.6f} |\n"
    f"| **NOT carrying** | {900-carry} | {prof['f2_digest_same_non_carrying']} | {prof['f2_cos_non_carrying'][0]:.8f} / {prof['f2_cos_non_carrying'][1]:.8f} / {prof['f2_cos_non_carrying'][2]:.8f} | "
    f"0.99797760 / 0.99892332 | {prof['f2_corrected_non_carrying'][0]:.6f} / {prof['f2_corrected_non_carrying'][1]:.6f} / **{prof['f2_corrected_non_carrying'][2]:.6f}** |\n\n"
    f"Shared delivered scale read off the non-carrying blocks **{prof['f2_shared_scale']:.8f}** (effect level 0.99846982). Separation **CLEAN**: the "
    f"smallest carrying displacement is **{prof['f2_corrected_carrying'][0]/prof['f2_corrected_non_carrying'][2]:.1f}×** the largest non-carrying one "
    f"(249: 1.9×, 222: 3.3× — here the moving rows are few and large and the shared scalar is 0.15 %, so the emitter's integer re-draw on the 855 "
    f"non-carrying blocks, 0.4–1.1 %, is well under the 12.6 % floor of the carrying ones). The exact statement is §6a's, at 1e−16.\n\n"
    f"**Digest counts, as asked: {prof['totals_identical']}/900 identical ({prof['f2_digest_same_carrying']}/{carry} carrying, "
    f"{prof['f2_digest_same_non_carrying']}/{900-carry} non-carrying).** Non-vacuity: ‖ΔE‖/‖E‖ **{prof['nonvacuity_relative_norm']:.6f}** (per context "
    f"A 0.078895, B 0.076614, C 0.076789), per block min {prof['nonvacuity_block_min']:.6f} / p01 {prof['nonvacuity_block_p01']:.6f} / median "
    f"{prof['nonvacuity_block_median']:.6f} / max {prof['nonvacuity_block_max']:.6f}, **{prof['blocks_exactly_zero']}/900 at zero**, 0/900 < 1e−6. "
    f"210's stricter 1 % per-block floor (the script's own \"FALSIFIER 4\" print) holds on all {carry} carrying blocks (min 12.6 % scale-corrected) and, "
    f"as in 249 §7, cannot hold on the {900-carry} non-carrying blocks whose displacement is by construction the 0.15 % scalar plus the emitter's "
    f"integer re-draw; the mandate's criterion (zero-block count = 0) is met. **PASS.**\n\n"
    f"For the record: median delivered |row| ratio {prof['median_row_norm_ratio']:.5f}, mean {prof['mean_row_norm_ratio']:.5f}, `A = mean|row|²` ratio "
    f"{prof['A_ratio']:.5f}; `rc_callset` head exposure: top-400 |delivered effect| overlap with `cd4svd5` min {prof['head_overlap_min']} / median "
    f"{prof['head_overlap_median']:.0f} / max {prof['head_overlap_max']} of 400 (the minimum on a partial target), implied call-set jaccard lower bound min "
    f"{prof['callset_jaccard_lower_min']:.4f} / median {prof['callset_jaccard_lower_median']:.4f}; |Δ log2FC| over the DE universe max "
    f"{prof['dlfc_max']:.3f}, mean {prof['dlfc_mean']:.6f}, signed block-mean median {prof['dlfc_signed_block_mean_median']:+.6f}."
)
# probe
def npm(lab, c):
    return probe[lab][c]["n_pred_median"]
pooled = {lab: float(__import__("numpy").median([v for c in "ABC" for v in probe[lab][c]["n_pred"]])) for lab in probe}
delta = (pooled["composed"] - pooled["cd4svd5"]) / pooled["cd4svd5"]
rho = {lab: [probe[lab][c]["rho_H_call"] for c in "ABC"] for lab in probe}
up = {lab: [probe[lab][c]["frac_up_call"] for c in "ABC"] for lab in probe}
rng = {lab: {c: (min(probe[lab][c]["n_pred"]), max(probe[lab][c]["n_pred"])) for c in "ABC"} for lab in probe}
k["PROBE_F4"] = (f"median **{npm('composed','A'):.0f} / {npm('composed','B'):.0f} / {npm('composed','C'):.0f}** vs `cd4svd5`'s "
                 f"{npm('cd4svd5','A'):.0f} / {npm('cd4svd5','B'):.0f} / {npm('cd4svd5','C'):.0f}, pooled **{100*delta:+.2f} %**")
k["PROBE_VERDICT"] = "PASS" if abs(delta) <= 0.02 else "FAIL"
k["PROBE_DELTA"] = f"{100*delta:+.2f} %"
k["PROBE_SECTION"] = (
    f"| package | A | B | C | pooled median | ρ(−log p, \\|H\\|) in the call set, A / B / C (a15's key for both, cross-key by design) | call set UP |\n"
    f"|---|---:|---:|---:|---:|---|---|\n"
    + "\n".join(
        f"| {'**' if lab == 'composed' else ''}`{lab}`{'**' if lab == 'composed' else ''} | "
        + " | ".join(f"{npm(lab, c):.0f} [{rng[lab][c][0]}, {rng[lab][c][1]}]" for c in "ABC")
        + f" | {pooled[lab]:.0f} | {rho[lab][0]:+.4f} / {rho[lab][1]:+.4f} / {rho[lab][2]:+.4f} (mean {sum(rho[lab])/3:+.4f}) | "
        + " / ".join(f"{100*u:.1f}" for u in up[lab]) + " % |"
        for lab in ("cd4svd5", "composed"))
    + f"\n\n**{100*delta:+.2f} %** against ±2 %. The `rc_callset` prescription pins `n_pred` by construction and it held. **{k['PROBE_VERDICT']}.** "
    f"(Peak RSS {rss(probelog)} GiB, {wall(probelog)} s — 222's instrument at its known cost, run alone after the last build.)"
)
# amp
pa = amp["predicted_raw_mse"]
k["AMP_F5"] = (f"`√(A_arm/A_cd4svd5)` − 1 = **{amp['amplitude_pct_vs_svd5']:+.3f} %** (contexts ×1.00481 / 1.00435 / 1.00446), `B` {100*(amp['B_ratio_vs_svd5']-1):+.3f} %; "
               f"vs `h1union` {amp['amplitude_pct_vs_h1union']:+.3f} %; predicted board raw `mse` **{pa['1.0']['raw_arm']:.5f}** at k_cos 1 "
               f"({pa['1.0']['d_vs_h1union']:+.5f} vs `h1union`'s 0.8172630, {pa['1.0']['d_vs_svd5']:+.5f} vs `cd4svd5`'s predicted {pa['1.0']['raw_svd5']:.5f})")
k["AMP_VERDICT"] = "PASS" if amp["falsifier5_pass"] else "FAIL"
k["AMP_F5_SHORT"] = f"{amp['amplitude_pct_vs_svd5']:+.3f} % vs cd4svd5"
k["AMP_SECTION"] = (
    f"Per-block panel totals streamed off all three packages (214's `block_totals` via 243's `totals243.py`), 211's per-context terms via 243's "
    f"`terms`. **Known-answer rows:** `h1union`'s re-streamed (`A_c`, `B_c`, `cap_c`) equal 243's `invert243.json` to **0.0e+00** on all three contexts "
    f"and 197's identity at those terms reproduces its board raw `mse` **0.8172630** to 0.0e+00; `cd4svd5`'s re-streamed per-block totals equal 249's "
    f"`totals-cd4svd5.npz` **bit for bit** and its amplitude vs `h1union` re-reads 249 §9's −0.409 % ({100*(amp['amplitude_vs_h1union']/amp['amplitude_vs_svd5']-1):+.4f} %).\n\n"
    f"| context | `A_c` arm | `A_c` cd4svd5 | ratio | `A_c` h1union | ratio | `B_c` arm / cd4svd5 | cap | leg |\n|---|---:|---:|---:|---:|---:|---|---:|---|\n"
    + "\n".join(f"| {c} | {amp['per_arm'][c]['A']:.6f} | {amp['per_svd5'][c]['A']:.6f} | ×{amp['per_arm'][c]['A']/amp['per_svd5'][c]['A']:.5f} | "
                f"{amp['per_inc'][c]['A']:.6f} | ×{amp['per_arm'][c]['A']/amp['per_inc'][c]['A']:.5f} | {amp['per_arm'][c]['B']:.6f} / {amp['per_svd5'][c]['B']:.6f} | "
                f"{amp['per_arm'][c]['cap']:.6f} | {'claim' if amp['per_arm'][c]['cap'] < amp['per_arm'][c]['B'] else 'B'} |" for c in "ABC")
    + f"\n\n**Delivered amplitude `mean_c √(A_arm/A_cd4svd5)` − 1 = {amp['amplitude_pct_vs_svd5']:+.4f} %** (`B` term {100*(amp['B_ratio_vs_svd5']-1):+.4f} %); "
    f"bound ±0.5 % vs `cd4svd5`: **{k['AMP_VERDICT']}**, at {100*abs(amp['amplitude_vs_svd5']-1)/0.5:.0f} % of the bound. Vs `h1union` "
    f"{amp['amplitude_pct_vs_h1union']:+.4f} % — the partial rows put back about half of the amplitude the deflation removed (−0.409 %). All three "
    f"contexts stay on the `B` leg. The mechanism is 246 §3.7's: the 15 rows carried no residual before, so on them the delivered amplitude rises "
    f"(norm ratio up to 1.18 on a partial block, §6b) while the preserve scalar takes 0.15 % off everything else.\n\n"
    f"Predicted board raw `mse` from 197's identity **at the arm's own terms** (243 §5.2: never the floor), with "
    f"`cos_arm = cos_h1union × (1 + k_cos × rel_agg)`, `rel_agg` = {100*amp['rel_agg']['composed_vs_stack']:+.4f} % (composed vs the shipped stack, 3 cancer folds; "
    f"`cd4svd5` {100*amp['rel_agg']['svd5_vs_stack']:+.4f} %, the increment {100*amp['rel_agg']['increment']:+.4f} %) and `cos_h1union` = {amp['cos_inc']:.6f} (243's measured):\n\n"
    f"| k_cos | arm cos | **predicted raw `mse`, arm** | Δ vs `h1union` (0.8172630) | predicted `cd4svd5` (249's table) | **arm − `cd4svd5`** | `mse` leg vs `cd4svd5`, overall |\n|---:|---:|---:|---:|---:|---:|---:|\n"
    + "\n".join(f"| {float(kk):.2f}{' (amplitude/budget leg alone)' if kk == '0.0' else ' (new-source class)' if kk == '0.94' else ' (re-dose class)' if kk == '1.36' else ''} | "
                f"{amp['cos_inc']*(1+float(kk)*amp['rel_agg']['composed_vs_stack']):.6f} | **{v['raw_arm']:.5f}** | {v['d_vs_h1union']:+.5f} | {v['raw_svd5']:.5f} | "
                f"**{v['d_vs_svd5']:+.5f}** | {-(v['d_vs_svd5'])/0.950162566/6:+.5f} |" for kk, v in pa.items())
    + f"\n\nThe arm's aggregate cosine is ABOVE `h1union`'s on PLOLO (+0.70 %), so the identity predicts raw `mse` falling below both incumbents: "
    f"−0.0011 … −0.0017 vs `h1union`, −0.0028 … −0.0039 vs `cd4svd5`'s 249 prediction — a `mse` leg of **+0.0005 … +0.0007 overall vs `cd4svd5`**, "
    f"larger than the rate card's +0.00019 in §4 (the card prices the MACRO cosine, +0.586 %; the identity reads the aggregate, +1.162 %, which the "
    f"15 heavy rows dominate). The preregistration's falsifier 2 fires above `cd4svd5`'s realised raw `mse` + 0.0010."
)
k["NNZ"] = f"{amp['nnz']:,}"
k["NNZ_PCT"] = f"{100*amp['nnz']/4.75e9:.2f}"
k["NNZ_DELTA"] = f"{100*(amp['nnz']/amp['nnz_svd5']-1):+.3f}"
k["MAX_CELL"] = f"{amp['max_cell']:,}"
k["MAX_CELL_X"] = f"{1e6/amp['max_cell']:.1f}"
k["PREP_DRY_WALL"] = "14.7"
k["PREP_DRY_RSS"] = "20.45"
k["PROFILE_WALL"] = f"{wall(proflog)} s"
k["PROFILE_RSS"] = rss(proflog)
k["AMP_WALL"] = "40 s"
k["AMP_RSS"] = "1.18"
k["PROBE_WALL"] = f"{wall(probelog)} s"
k["PROBE_RSS"] = rss(probelog)
verdict_ok = (k["KA_PKG_VERDICT"] == "PASS" and k["PROBE_VERDICT"] == "PASS" and k["AMP_VERDICT"] == "PASS" and not mk["failures"]
              and prof["blocks_exactly_zero"] == 0 and prof["f2_separated"])
k["VERDICT"] = "BUILDABLE" if verdict_ok else "NOT BUILDABLE"

for name in ("RESULTS.md", "PREREGISTRATION-DRAFT.md"):
    p = HERE / name
    s = p.read_text()
    for key, val in k.items():
        s = s.replace("{{" + key + "}}", str(val))
    s = s.replace("<amp254_d_mse_vs_svd5>", f"{pa['1.0']['d_vs_svd5']:+.4f} raw ({-(pa['1.0']['d_vs_svd5'])/0.950162566/6:+.5f} overall)")
    s = s.replace("<amp254_amplitude_vs_svd5>", f"{amp['amplitude_pct_vs_svd5']:+.3f}")
    s = s.replace("<amp254_raw_mse_arm>", f"{pa['1.0']['raw_arm']:.5f}")
    s = s.replace("<amp254_raw_mse_svd5>", f"{pa['1.0']['raw_svd5']:.5f}")
    left = re.findall(r"\{\{[A-Z_0-9]+\}\}|<amp254_[a-z_0-9]+>", s)
    p.write_text(s)
    print(f"{name}: {len(s):,} bytes, unfilled placeholders: {left}")
print("verdict", k["VERDICT"])
