from __future__ import annotations

import numpy as np
import pandas as pd


def align_ice_growth_rates(
    data: pd.DataFrame, experimental: pd.DataFrame
) -> pd.Series:
    """Align experimental ice-growth-rate values to the main feature table."""
    rate_map = dict(
        zip(
            experimental["Sequence_clean"],
            experimental["Ice_Growth_Rate"],
            strict=False,
        )
    )
    return data["Sequence_clean"].map(rate_map).rename("Ice_Growth_Rate")


def summarize_experimental_by_cluster(
    labels: np.ndarray,
    ice_growth_rate: pd.Series,
    favorable_threshold: float,
) -> pd.DataFrame:
    """Count matched and favorable experimental peptides in each cluster."""
    rates = ice_growth_rate.to_numpy(dtype=float)
    matched = ~np.isnan(rates)
    rows: list[dict[str, int]] = []

    for cluster_id in sorted(np.unique(labels)):
        cluster_mask = labels == cluster_id
        rows.append(
            {
                "raw_cluster": int(cluster_id),
                "size": int(cluster_mask.sum()),
                "n_experimental": int(np.sum(cluster_mask & matched)),
                "n_favorable": int(
                    np.sum(
                        cluster_mask
                        & matched
                        & (rates <= favorable_threshold)
                    )
                ),
            }
        )
    return pd.DataFrame(rows)


def select_experimentally_enriched_cluster(summary: pd.DataFrame) -> int:
    """Select the cluster with most favorable, then most matched, peptides."""
    if summary.empty:
        raise ValueError("Experimental cluster summary is empty.")
    ranked = summary.sort_values(
        ["n_favorable", "n_experimental", "size", "raw_cluster"],
        ascending=[False, False, False, True],
    )
    return int(ranked.iloc[0]["raw_cluster"])
