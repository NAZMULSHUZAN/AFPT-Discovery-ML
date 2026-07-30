#!/usr/bin/env python3
"""Create cluster-level secondary-structure summary and figure from a peptide CSV."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(path: Path) -> dict:
    with path.open() as handle:
        config = json.load(handle)
    root = repo_root()
    for key in ("summary_csv", "secondary_structure_csv", "figures_dir"):
        config[key] = root / config[key]
    return config


def create_summary(dataframe: pd.DataFrame, cluster_order: list[str]) -> pd.DataFrame:
    required = {"cluster", "helix_frac", "sheet_frac", "turn_frac", "coil_frac"}
    missing = required - set(dataframe.columns)
    if missing:
        raise ValueError(f"Summary CSV is missing columns: {sorted(missing)}")

    valid = dataframe.dropna(
        subset=["helix_frac", "sheet_frac", "turn_frac", "coil_frac"]
    ).copy()
    if valid.empty:
        raise RuntimeError(
            "No DSSP-derived values are present. Install DSSP and rerun "
            "analyze_af2_results.py, or provide a completed summary CSV."
        )

    grouped = (
        valid.groupby("cluster", as_index=False)
        .agg(
            helix_pct=("helix_frac", lambda x: 100 * x.mean()),
            beta_pct=("sheet_frac", lambda x: 100 * x.mean()),
            turn_pct=("turn_frac", lambda x: 100 * x.mean()),
            coil_pct=("coil_frac", lambda x: 100 * x.mean()),
            n_peptides=("peptide_id", "count"),
        )
    )

    order_map = {name: index for index, name in enumerate(cluster_order)}
    grouped["_order"] = grouped["cluster"].map(order_map).fillna(len(order_map))
    return grouped.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def make_plot(summary: pd.DataFrame, output_path: Path) -> None:
    clusters = summary["cluster"].tolist()
    x = np.arange(len(clusters))
    width = 0.2

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x - 1.5 * width, summary["helix_pct"], width, label="Helix")
    ax.bar(x - 0.5 * width, summary["beta_pct"], width, label="Beta sheet")
    ax.bar(x + 0.5 * width, summary["turn_pct"], width, label="Turn")
    ax.bar(x + 1.5 * width, summary["coil_pct"], width, label="Coil")

    ax.set_xticks(x)
    ax.set_xticklabels(clusters, rotation=20, ha="right")
    ax.set_ylabel("Mean residue composition (%)")
    ax.set_title("Secondary Structure Composition Across AFPT Groups")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] {output_path}")


def main(config: dict) -> None:
    dataframe = pd.read_csv(config["summary_csv"])
    summary = create_summary(dataframe, config["clusters"])

    output_csv: Path = config["secondary_structure_csv"]
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_csv, index=False)
    print(f"[SAVED] {output_csv}")

    figure_path = config["figures_dir"] / "secondary_structure_clusters.png"
    make_plot(summary, figure_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=repo_root() / "config.json",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(load_config(args.config))
