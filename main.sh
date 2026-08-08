#!/usr/bin/env bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

PGGB_DIR="$ROOT_DIR/src/pggb"

bash "$PGGB_DIR/pangenome_partitioning_pipeline.sh"
bash "$PGGB_DIR/run_pggb_params.sh"
bash "$PGGB_DIR/small_variants_eval.sh"
