from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage


@dataclass(frozen=True)
class HierarchicalResult:
    target_indices: np.ndarray
    subcluster_labels: np.ndarray
    counts: dict[int, int]


def run_targeted_ward_clustering(
    matrix: np.ndarray,
    target_mask: np.ndarray,
    *,
    n_subclusters: int = 3,
    method: str = "ward",
    criterion: str = "maxclust",
) -> HierarchicalResult:
    """Apply hierarchical clustering only to the selected main cluster."""
    target_indices = np.flatnonzero(target_mask)
    if len(target_indices) < n_subclusters:
        raise ValueError("Target cluster is smaller than n_subclusters.")

    tree = linkage(matrix[target_indices], method=method)
    labels = fcluster(tree, t=n_subclusters, criterion=criterion)
    unique, counts = np.unique(labels, return_counts=True)
    return HierarchicalResult(
        target_indices=target_indices,
        subcluster_labels=labels,
        counts={int(key): int(value) for key, value in zip(unique, counts, strict=True)},
    )
