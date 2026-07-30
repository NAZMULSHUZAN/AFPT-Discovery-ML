#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
python scripts/check_setup.py
python scripts/analyze_af2_results.py --config config.json
