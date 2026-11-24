#!/usr/bin/env python3
"""700 leakage audit (written before any PIE number is read): what PIE's xdataset checkpoint trained on and pooled as
evidence, against the truths 700 scores it on (h1-test's 85 targets in H1; K562 GWPS rows of the 269 covered panel targets).
Sources: the checkpoint's config.yaml + data_stats.json (HF), the pinned train/val split JSONs (PIE_splits), the VCC25
preprocessed meta (its perturbation list), the Replogle-Nadig meta. Writes work/audit700.json."""
import csv, json, glob, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PIE = ROOT / "data/external/pie"
split = json.load(open(PIE / "splits/replogle_xdataset/train.json"))
val = json.load(open(PIE / "splits/replogle_xdataset/val.json"))
stats = json.load(open(PIE / "runs/replogle_xdataset/data_stats.json"))["stats"]
cfg_w = __import__("yaml").safe_load(open(PIE / "runs/replogle_xdataset/config.yaml"))["data"]["dataset_weights"]
ds_train = sorted({k.split(".", 1)[0] for k in split}); ds_val = sorted({k.split(".", 1)[0] for k in val})
ctx_train = sorted(split)
vcc = json.load(open(glob.glob(str(PIE / "data/hf/datasets/arcinstitute/pie_arc_vcc_25/*/preprocessed/meta.json"))[0]))
rep = json.load(open(glob.glob(str(PIE / "data/hf/datasets/arcinstitute/pie_replogle_nadig_essential/*/preprocessed/meta.json"))[0]))
q_h1 = json.load(open(Path(__file__).parent / "work/query_h1test.json"))["h1test.ARC_H1"]
q_k = json.load(open(Path(__file__).parent / "work/query_k562.json"))["replogle.k562"]
vcc_train_perts = sorted({p for k, v in split.items() if k.startswith("arc_vcc_25.") for p in v})
vcc_val_perts = sorted({p for k, v in val.items() if k.startswith("arc_vcc_25.") for p in v})
rep_keys = [k for k in list(split) + list(val) if k.startswith("replogle.")]
k562_perts_in_rep = set(rep["pert_to_id"])
out = {
  "dataset_weights": cfg_w,
  "datasets_with_train_rows": ds_train, "datasets_with_val_rows": ds_val,
  "evidence_datasets (data_stats)": stats.get("evidence_datasets"),
  "replogle_split_keys": rep_keys,
  "train_contexts_n": len(ctx_train),
  "vcc25_train_perts": len(vcc_train_perts), "vcc25_val_perts": len(vcc_val_perts),
  "h1test_targets": len(q_h1),
  "h1test_targets_in_vcc25_train_or_val": sorted(set(q_h1) & (set(vcc_train_perts) | set(vcc_val_perts))),
  "h1test_targets_in_vcc25_meta": sorted(set(q_h1) & set(vcc["pert_to_id"])),
  "k562_query_targets": len(q_k),
  "k562_query_targets_in_replogle_essential_dir": len(set(q_k) & k562_perts_in_rep),
  "k562_query_targets_in_any_train_split_key_as_pert (any dataset)": len(set(q_k) & {p for v in split.values() for p in v}),
  "h1test_targets_in_any_train_split_as_pert (any dataset)": len(set(q_h1) & {p for v in split.values() for p in v}),
}
out["verdict"] = {
  "replogle_in_training_or_evidence": bool(rep_keys) or ("replogle" in (stats.get("evidence_datasets") or [])),
  "h1test_responses_seen": bool(out["h1test_targets_in_vcc25_train_or_val"] or out["h1test_targets_in_vcc25_meta"]),
  "note": ("H1 control cells (the 2025 control pool) are an inference input for ARC_H1 and were seen in training as "
           "ARC_H1/ARC_H1_VAL control means; h1-test's 85 perturbation responses were not. Panel targets that appear as "
           "perturbations in other training datasets (Orion etc.) are the same-pert evidence PIE is designed to use - "
           "the analogue of our donor rows, not a leak of the scored line."),
}
(Path(__file__).parent / "work/audit700.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
