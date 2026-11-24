#!/usr/bin/env python3
"""700 — PIE inference for new contexts with the pretrained xdataset checkpoint, on this Mac (MPS or CPU).

Wraps pie.predict.predict_loaded with three changes to the checkpoint's saved config, none of which touches the weights:
  1. data.source_dirs["context_text"] -> a local copy extended with the new contexts' Cellosaurus texts (700's
     data/external/pie/local_sources/context_text: + CVCL_0030 HeLa, CVCL_1108 CAL-33, embedded with
     text-embedding-3-large, as PIE's own builder does);
  2. trainer.precision -> "32-true" unless --bf16 (MPS autocast support for bf16 is not assumed);
  3. data.num_workers -> 0 (macOS spawn; the loader is cheap next to the forward pass).
Evidence comes from the checkpoint's own training config (evidence_source = the saved data config), so the first run
builds or reuses PIE_CACHE_DIR/evidence/<key> from the training datasets, as `pie infer` would.

usage: PIE_DATA_ROOT=... PIE_CACHE_DIR=... .venv/bin/python infer700.py --ckpt <best_auprc.ckpt> --prep <dir> [--prep <dir>]
         --query <query.json> --out <predictions.parquet> [--device mps|cpu] [--batch 8] [--context-text <dir>]
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import torch

from pie.predict import RowSource, load_checkpoint, predict_loaded, write_predictions_parquet

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
log = logging.getLogger("infer700")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--prep", action="append", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--device", default="mps" if torch.backends.mps.is_available() else "cpu")
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--context-text", default=None)
    ap.add_argument("--bf16", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    loaded = load_checkpoint(Path(a.ckpt))
    cfg = loaded.config
    data_upd = {"num_workers": 0}
    if a.context_text:
        sd = dict(cfg.data.source_dirs)
        sd["context_text"] = str(Path(a.context_text).resolve())
        data_upd["source_dirs"] = sd
    new_data = cfg.data.model_copy(update=data_upd)
    upd = {"data": new_data}
    if not a.bf16:
        upd["trainer"] = cfg.trainer.model_copy(update={"precision": "32-true"})
    loaded.config = cfg.model_copy(update=upd)
    log.info("checkpoint loaded in %.0fs; device %s; precision %s; context_text %s", time.time() - t0, a.device,
             loaded.config.trainer.precision, loaded.config.data.source_dirs.get("context_text"))
    preds = predict_loaded(loaded, RowSource("query", Path(a.query)), [Path(p) for p in a.prep], a.device, a.batch)
    out = write_predictions_parquet(preds, Path(a.out))
    n = sum(len(b.contexts) for b in preds.blocks)
    log.info("wrote %d rows to %s in %.0fs total", n, out, time.time() - t0)
    print(json.dumps({"rows": n, "out": str(out), "seconds": round(time.time() - t0, 1)}))


if __name__ == "__main__":
    main()
