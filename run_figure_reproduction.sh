#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
python code/validate_figure_data.py
python code/regenerate_cluster_plot.py \
  --data data/AFPT_cluster_plot_data.csv \
  --output outputs
