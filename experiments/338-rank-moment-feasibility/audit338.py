"""Independent saved-histogram and primal/dual certificate audit.

Does not call the construction's Phi, bounds, or LP solver. This verifies
the recorded relaxation, not joint integer emission or response accuracy.
"""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mid_control_mass(control, values):
    """Group the combined support and accumulate CONTROL multiplicities."""
    _, inverse = np.unique(np.concatenate((control, values)), return_inverse=True)
    sizes = np.bincount(inverse[:len(control)], minlength=int(inverse.max())+1)
    lower = np.cumsum(sizes) - sizes
    return (2*lower + sizes)[inverse[len(control):]] / (2*len(control))


def baseline_rank(control, values):
    n, m = len(values), len(control)
    combined = np.concatenate((control, values))
    ranks = rankdata(combined, method="average")
    U = float(ranks[m:].sum() - n*(n+1)/2)
    _, counts = np.unique(combined, return_counts=True)
    tie = sum(int(t)**3-int(t) for t in counts)
    variance = n*m/12*((n+m+1) - tie/((n+m)*(n+m-1)))
    z = (U-n*m/2)/np.sqrt(variance) if variance else 0.
    return U, float(z)


def audit_bound(record, control, original_counts, original_depths, constrained):
    n = len(original_counts)
    total = int(original_counts.sum())
    zeros = int(np.count_nonzero(original_counts == 0))
    original_linear = original_counts.astype(float)/(original_depths.astype(float)/1e6)
    target_phi = float(mid_control_mass(control, np.log1p(original_linear)).sum())
    original_mean = float(original_linear.mean())
    np.testing.assert_allclose(record["baseline_phi_sum"], target_phi, rtol=0, atol=1e-10)
    np.testing.assert_allclose(record["baseline_mean_cpm"], original_mean, rtol=0, atol=1e-10)
    assert record["baseline_gene_total"] == total
    assert record["baseline_zero_cells"] == zeros
    groups = record["depth_groups"]
    K = min(total, min(g["depth"] for g in groups),
            max(64, 4*int(original_counts.max(initial=0))))
    assert record["count_support_max"] == K
    depth = np.repeat([g["depth"] for g in groups], K+1)
    level = np.tile(np.arange(K+1), len(groups))
    group_index = np.repeat(np.arange(len(groups)), K+1)
    coefficients = mid_control_mass(control, np.log1p(level/(depth/1e6)))
    A = [(group_index == j).astype(float) for j in range(len(groups))]
    b = [float(g["cells"]) for g in groups]
    A += [level/max(K, 1), (level == 0).astype(float)]
    b += [total/max(K, 1), float(zeros)]
    if constrained:
        A.append(coefficients)
        b.append(target_phi)
    A, b = np.array(A), np.array(b)
    scale = max(original_mean, 1.)
    objective = level/(depth/1e6)/n
    lookup = {(int(d), int(k)): i for i, (d, k) in enumerate(zip(depth, level))}
    checks = []
    for name, sign in (("minimum", 1), ("maximum", -1)):
        endpoint = record[name]
        if endpoint["status"] != 0:
            checks.append({"endpoint": name, "solver_status": endpoint["status"],
                           "certificate_available": False})
            continue
        x = np.zeros(len(level))
        seen = set()
        for item in endpoint["histogram"]:
            key = (item["depth"], item["count"])
            assert key not in seen and key in lookup
            seen.add(key)
            assert item["multiplicity"] >= 0.
            x[lookup[key]] = item["multiplicity"]
        count_error = abs(float(level.astype(np.longdouble)@x)-total)
        equality_error = float(np.max(np.abs(A@x-b)))
        phi_error = abs(float(coefficients@x)-target_phi)
        measured = float(objective@x)
        dual = np.asarray(endpoint["dual_equality_marginals"])
        reduced = sign*objective/scale-A.T@dual
        gap = float(sign*measured/scale-b@dual)
        assert count_error <= 1e-4, count_error
        assert equality_error <= 1e-6, equality_error
        assert not constrained or phi_error <= 1e-6, phi_error
        assert float(reduced.min()) >= -1e-7, reduced.min()
        assert abs(gap) <= 1e-6, gap
        np.testing.assert_allclose(measured, endpoint["mean_cpm"], rtol=0, atol=1e-7)
        np.testing.assert_allclose(b@dual, endpoint["dual_objective_scaled"],
                                   rtol=0, atol=1e-9)
        checks.append({
            "endpoint": name, "solver_status": 0, "certificate_available": True,
            "count_total_absolute_error": count_error,
            "equality_absolute_error": equality_error,
            "phi_sum_absolute_error": phi_error,
            "minimum_reduced_cost": float(reduced.min()),
            "scaled_primal_dual_gap": gap,
            "active_histogram_entries": len(seen),
        })
    return checks


def main():
    manifest = json.loads((HERE/"inputs.json").read_text())
    for name, expected in manifest["sha256"].items():
        assert sha(ROOT/name) == expected, name
    result = json.loads((WORK/"results.json").read_text())
    selection = json.loads((HERE/"selection.json").read_text())
    baseline = json.loads((WORK/"baseline.json").read_text())
    assert result["inputs_sha256"] == sha(HERE/"inputs.json")
    assert result["baseline_sha256"] == sha(WORK/"baseline.json")
    assert baseline["selected_blocks_sha256"] == sha(WORK/"selected-blocks.npz")
    expected_digest = json.loads((BASE/"work/metrics-cap15.json").read_text())[
        "emitted_count_digest_sha256"]
    assert baseline["emitted_count_digest_sha256"] == expected_digest
    with np.load(BASE/"work/reference.npz") as archive:
        # Deliberately load no real perturbed response members.
        control = archive["cs"]
        genes = archive["genes"].astype(str)
        keep = archive["keep"]
    with np.load(WORK/"selected-blocks.npz") as archive:
        blocks = archive["blocks"]
        targets = archive["targets"].astype(str)
        np.testing.assert_array_equal(archive["genes"].astype(str), genes)
    assert list(targets) == selection["targets"]
    with np.load(BASE/"work/stats-cap15.npz") as archive:
        expected_sums = archive["sums"]
    for block, index in zip(blocks, selection["target_indices"]):
        np.testing.assert_array_equal(block.sum(0), expected_sums[index])
    mapping = {str(target): i for i, target in enumerate(targets)}
    selected_genes = dict(zip(selection["genes"], selection["full_gene_indices"]))
    cases = result["cases"]
    assert len(cases) == 128
    assert {(c["target"], c["gene"]) for c in cases} == {
        (t, g) for t in targets for g in selection["genes"]}
    checks = []
    solved = 0
    rooms = {str(t): 0 for t in targets}
    max_U_error = max_z_error = 0.
    declared_checks = {(c["target"], c["gene"]): c for c in result["checks"]}
    for case in cases:
        target, gene = case["target"], case["gene"]
        block = blocks[mapping[target]]
        full = selected_genes[gene]
        j = case["control_gene_index"]
        assert genes[np.flatnonzero(keep)[j]] == gene
        counts = block[:, full]
        depths = block.sum(1)
        values = np.log1p(counts.astype(float)/(depths.astype(float)/1e6))
        U, z = baseline_rank(control[j], values)
        U_lookup = float(mid_control_mass(control[j], values).sum()*len(control[j]))
        assert abs(U-U_lookup) < 1e-7
        max_U_error = max(max_U_error, abs(U-U_lookup))
        # The recorded run compares its independent z to production. The
        # audit additionally reconstructs z from Python-integer tie sums.
        from vcc2026.rc_emit import wilcoxon_z, ctrl_tiesums
        cs = np.ascontiguousarray(control[j:j+1])
        native = float(wilcoxon_z(cs, ctrl_tiesums(cs),
                                 np.ascontiguousarray(values[None, :]))[0])
        max_z_error = max(max_z_error, abs(z-native))
        assert abs(z-native) < 1e-12
        assert declared_checks[(target, gene)]["baseline_lp_pass"]
        endpoints = {}
        for name in ("same_depth", "variable_depth", "without_rank"):
            endpoints[name] = audit_bound(case[name], control[j], counts, depths,
                                         name != "without_rank")
        same, changed, free = [case[k] for k in
                               ("same_depth", "variable_depth", "without_rank")]
        assert all(same[k]["status"] == 0 for k in ("minimum", "maximum"))
        mean = float(np.mean(counts/(depths/1e6)))
        tol = 1e-5*max(mean, 1)
        assert same["minimum"]["mean_cpm"] <= mean+tol
        assert same["maximum"]["mean_cpm"] >= mean-tol
        assert same["maximum"]["mean_cpm"]-same["minimum"]["mean_cpm"] <= (
            1e6*int(counts.sum())/400*(1/depths.min()-1/depths.max())+2*tol)
        feasible = all(changed[k]["status"] == 0 for k in ("minimum", "maximum"))
        assert case["constrained_extrema_solved"] == feasible
        solved += feasible
        room = 0.
        if feasible:
            room = min(mean-changed["minimum"]["mean_cpm"],
                       changed["maximum"]["mean_cpm"]-mean)/max(mean, 1.)
            if all(free[k]["status"] == 0 for k in ("minimum", "maximum")):
                assert free["minimum"]["mean_cpm"] <= changed["minimum"]["mean_cpm"]+1e-5
                assert free["maximum"]["mean_cpm"] >= changed["maximum"]["mean_cpm"]-1e-5
        np.testing.assert_allclose(room, case["symmetric_relative_room"], rtol=0, atol=1e-10)
        assert case["room_at_least_five_percent"] == bool(feasible and room >= .05)
        rooms[target] += bool(feasible and room >= .05)
        checks.append({"target": target, "gene": gene, "bounds": endpoints})
    fractions = {t: count/32 for t, count in rooms.items()}
    gates = {
        "baseline_and_rank_instruments": True,
        "at_least_95_percent_extrema_solved": solved >= .95*128,
        "each_target_at_least_quarter_with_five_percent_room":
            all(f >= .25 for f in fractions.values()),
    }
    assert gates == result["gates"]
    assert fractions == result["target_room_fractions"]
    status = ("EARNS_INTEGER_ROW_COUPLING_FEASIBILITY_ONLY" if all(gates.values()) else
              "DOES_NOT_EARN_INTEGER_CONSTRUCTION")
    assert status == result["status"]
    out = {
        "status": "PASS", "experiment_status": status, "gates_reproduced": gates,
        "cases": 128, "extrema_solved_pairs": solved, "checks": checks,
        "max_U_error": max_U_error, "max_z_error": max_z_error,
        "result_sha256": sha(WORK/"results.json"),
        "input_manifest_sha256": sha(HERE/"inputs.json"),
        "audit_code_sha256": sha(Path(__file__)),
        "reference_members_loaded": ["cs", "genes", "keep"],
        "reference_perturbation_values_loaded": False,
        "not_verified": ["independent whole-baseline re-emission",
                         "joint integer row construction", "biological accuracy"],
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    temporary = WORK/"independent-audit.json.tmp"
    temporary.write_text(json.dumps(out, indent=2, sort_keys=True, allow_nan=False)+"\n")
    temporary.replace(WORK/"independent-audit.json")
    print(json.dumps({k: v for k, v in out.items() if k != "checks"}, indent=2))


if __name__ == "__main__":
    main()
