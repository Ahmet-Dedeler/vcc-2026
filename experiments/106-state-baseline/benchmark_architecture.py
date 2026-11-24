#!/usr/bin/env python3
"""Run one synthetic STATE optimization step at the VCC feature dimensions."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import torch
import yaml

from state.tx.utils import get_lightning_module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("state", "state_sm"), default="state_sm")
    parser.add_argument("--device", choices=("cpu", "mps"), default="mps")
    parser.add_argument("--input-dim", type=int, default=18_080)
    parser.add_argument("--pert-dim", type=int, default=5_120)
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2] / "data/external/state-vcc/state"
    profile_path = repo / f"src/state/configs/model/{args.profile}.yaml"
    training_path = repo / "src/state/configs/training/default.yaml"

    profile = yaml.safe_load(profile_path.read_text())
    model_config = profile["kwargs"]
    backbone = model_config["transformer_backbone_kwargs"]
    backbone["max_position_embeddings"] = model_config["cell_set_len"]
    backbone["hidden_size"] = model_config["hidden_dim"]
    training_config = yaml.safe_load(training_path.read_text())
    training_config["batch_size"] = 1

    data_config = {
        "embed_key": None,
        "output_space": "all",
        "basal_mapping_strategy": "random",
        "control_pert": "non-targeting",
    }
    var_dims = {
        "input_dim": args.input_dim,
        "gene_dim": args.input_dim,
        "hvg_dim": args.input_dim,
        "output_dim": args.input_dim,
        "pert_dim": args.pert_dim,
        "batch_dim": 0,
        "gene_names": [f"G{i}" for i in range(args.input_dim)],
    }

    torch.manual_seed(args.seed)
    device = torch.device(args.device)
    if device.type == "mps" and not torch.backends.mps.is_available():
        raise SystemExit("MPS requested but unavailable")

    started = time.perf_counter()
    model = get_lightning_module(
        profile["name"],
        data_config,
        model_config,
        training_config,
        var_dims,
    ).to(device)
    build_seconds = time.perf_counter() - started

    sequence_length = int(model_config["cell_set_len"])
    ctrl = torch.rand(sequence_length, args.input_dim, device=device)
    pert_features = torch.zeros(sequence_length, args.pert_dim, device=device)
    pert_features[:, 0] = 1.0
    target = torch.rand(sequence_length, args.input_dim, device=device)

    optimizer = torch.optim.Adam(model.parameters(), lr=float(training_config["lr"]))
    optimizer.zero_grad(set_to_none=True)
    step_started = time.perf_counter()
    prediction = model(
        {"ctrl_cell_emb": ctrl, "pert_emb": pert_features},
        padded=True,
    ).reshape(1, sequence_length, args.input_dim)
    loss = model._compute_distribution_loss(
        prediction,
        target.reshape(1, sequence_length, args.input_dim),
    ).mean()
    loss.backward()
    optimizer.step()
    if device.type == "mps":
        torch.mps.synchronize()
    step_seconds = time.perf_counter() - step_started

    result = {
        "profile": args.profile,
        "device": str(device),
        "input_dim": args.input_dim,
        "pert_dim": args.pert_dim,
        "cell_set_len": sequence_length,
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
        "parameter_bytes": sum(
            parameter.numel() * parameter.element_size() for parameter in model.parameters()
        ),
        "build_seconds": build_seconds,
        "optimization_step_seconds": step_seconds,
        "loss": float(loss.detach().cpu()),
    }
    if device.type == "mps":
        result["mps_current_allocated_bytes"] = torch.mps.current_allocated_memory()
        result["mps_driver_allocated_bytes"] = torch.mps.driver_allocated_memory()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
