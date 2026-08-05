#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config.json"

with CONFIG.open() as handle:
    cfg = json.load(handle)

required_packages = {
    "numpy": "numpy",
    "pandas": "pandas",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "Bio": "biopython",
}

print("=== Python packages ===")
for module, package in required_packages.items():
    status = "OK" if importlib.util.find_spec(module) else f"MISSING: pip install {package}"
    print(f"{package:12s} {status}")

print("\n=== DSSP ===")
dssp = cfg.get("dssp_executable", "mkdssp")
print(f"{dssp}: {'OK' if shutil.which(dssp) else 'NOT FOUND (optional)'}")

print("\n=== FASTA files ===")
fasta_dir = ROOT / cfg["fasta_dir"]
for cluster in cfg["clusters"]:
    path = fasta_dir / f"{cluster}_sequences.fasta"
    print(f"{path.relative_to(ROOT)}: {'OK' if path.exists() else 'MISSING'}")

print("\n=== Structure folders ===")
structures = ROOT / cfg["structures_dir"]
for cluster in cfg["clusters"]:
    path = structures / cluster
    print(f"{path.relative_to(ROOT)}: {'OK' if path.exists() else 'MISSING'}")
