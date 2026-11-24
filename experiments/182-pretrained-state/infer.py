"""182 -- run the official `state tx infer` path on one prepared input.

Two compatibility shims, neither of which touches model maths:

1. The installed `transformers` validates `hidden_size % num_attention_heads == 0`.
   This checkpoint is hidden_size=328 / num_attention_heads=12 with an explicit
   head_dim=64, which is legal for the Llama attention implementation (q/k/v project
   328 -> 12*64) but trips the newer strict-dataclass validator. We disable that one
   validator. Lightning's `load_from_checkpoint` still loads the state dict strictly,
   so any real shape mismatch would still fail loudly.

2. Nothing else. Device stays CPU (no CUDA on this machine), which keeps the run
   deterministic.

Usage: infer.py <input.h5ad> <output.h5ad>
"""

import argparse
import sys

from transformers.models.llama.configuration_llama import LlamaConfig

# shim 1 -- see module docstring
LlamaConfig.validate = lambda self: None  # type: ignore[method-assign]

from state._cli._tx._infer import run_tx_infer  # noqa: E402

MODEL_DIR = "data/external/state-pretrained/ST-HVG-Replogle/zeroshot/hepg2"


def main():
    inp, out = sys.argv[1], sys.argv[2]
    args = argparse.Namespace(
        model_dir=MODEL_DIR,
        checkpoint=f"{MODEL_DIR}/checkpoints/best.ckpt",
        adata=inp,
        embed_key="X_hvg",
        pert_col="gene",
        output=out,
        celltype_col="cell_line",
        celltypes=None,
        batch_col="gem_group",
        control_pert=None,
        seed=20260903,
        max_set_len=None,
        quiet=False,
        tsv=None,
        all_perts=False,
        virtual_cells_per_pert=None,
        min_cells=None,
        max_cells=None,
    )
    run_tx_infer(args)


if __name__ == "__main__":
    main()
