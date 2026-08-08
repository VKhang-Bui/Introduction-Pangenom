#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$SCRIPT_DIR/pangenome_partitioning_pipeline.sh"
bash "$SCRIPT_DIR/run_pggb_params.sh"
bash "$SCRIPT_DIR/small_variants_eval.sh"
