from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


@dataclass(frozen=True)
class MainClusteringResult:
    labels: np.ndarray
    k: int
    silhouette_cosine: float
    counts: dict[int, int]


def run_kmeans(
    matrix: np.ndarray,
    *,
    k: int,
    random_state: int = 42,
    n_init: int = 10,
) -> MainClusteringResult:
    """Run deterministic KMeans and calculate cosine silhouette."""
    labels = KMeans(
        n_clusters=k,
        n_init=n_init,
        random_state=random_state,
    ).fit_predict(matrix)

    if len(np.unique(labels)) < 2:
        raise RuntimeError("KMeans returned fewer than two distinct clusters.")

    counts = pd.Series(labels).value_counts().sort_index().to_dict()
    silhouette = silhouette_score(matrix, labels, metric="cosine")
    return MainClusteringResult(
        labels=labels,
        k=k,
        silhouette_cosine=float(silhouette),
        counts={int(key): int(value) for key, value in counts.items()},
    )


def select_best_k(
    matrix: np.ndarray,
    *,
    k_values: list[int],
    random_state: int = 42,
    n_init: int = 10,
) -> MainClusteringResult:
    """Select k with the largest cosine silhouette score."""
    valid = [k for k in k_values if 2 <= k < len(matrix)]
    if not valid:
        raise ValueError("No valid k values were provided.")

    results = [
        run_kmeans(
            matrix,
            k=k,
            random_state=random_state,
            n_init=n_init,
        )
        for k in valid
    ]
    return max(results, key=lambda result: (result.silhouette_cosine, -result.k))
