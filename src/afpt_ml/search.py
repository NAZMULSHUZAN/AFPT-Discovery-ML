from __future__ import annotations

import itertools
from collections.abc import Iterable

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize

from afpt_ml.preprocessing import prepare_search_frame


def _limited_combinations(
    candidates: list[str],
    subset_min: int,
    subset_max: int,
    cap: int,
) -> Iterable[tuple[str, ...]]:
    count = 0
    for size in range(subset_min, subset_max + 1):
        for combination in itertools.combinations(candidates, size):
            yield combination
            count += 1
            if count >= cap:
                return


def legacy_candidate_pool(
    scaled: pd.DataFrame, candidate_pool: int
) -> list[str]:
    """Reproduce the notebook's variance ranking with stable tie ordering.

    Because StandardScaler makes variances nearly equal, this ranking can be
    sensitive to floating-point and library-version differences. It is kept as
    an extension/legacy search and never overwrites the published feature set.
    """
    variances = scaled.var()
    ranked = variances.sort_values(ascending=False, kind="mergesort")
    return ranked.index[: min(candidate_pool, len(ranked))].tolist()


def _score_subset(
    scaled: pd.DataFrame,
    feature_tuple: tuple[str, ...],
    rates: np.ndarray,
    *,
    k_values: list[int],
    threshold: float,
    reward_min_count: int,
    silhouette_weight: float,
    experimental_weight: float,
    random_state: int,
    n_init: int,
) -> dict[str, object]:
    matrix = normalize(scaled.loc[:, list(feature_tuple)].to_numpy(), norm="l2")
    matched = ~np.isnan(rates)

    best: dict[str, object] = {
        "silhouette": -1.0,
        "jacs_score": 0.0,
        "combined": -1.0,
        "k": None,
        "best_cluster": None,
        "best_n_good": 0,
        "best_n_jacs": 0,
    }

    for k in k_values:
        labels = KMeans(
            n_clusters=k,
            n_init=n_init,
            random_state=random_state,
        ).fit_predict(matrix)
        if len(np.unique(labels)) < 2:
            continue

        silhouette = float(silhouette_score(matrix, labels, metric="cosine"))
        cluster_rows: list[tuple[int, int, int]] = []
        for cluster_id in sorted(np.unique(labels)):
            cluster_mask = labels == cluster_id
            n_jacs = int(np.sum(cluster_mask & matched))
            n_good = int(
                np.sum(cluster_mask & matched & (rates <= threshold))
            )
            cluster_rows.append((int(cluster_id), n_good, n_jacs))

        best_cluster, best_n_good, best_n_jacs = max(
            cluster_rows, key=lambda row: (row[1], row[2], -row[0])
        )
        if best_n_good >= reward_min_count:
            jacs_score = 1.0 + 0.01 * min(
                best_n_good - reward_min_count, 50
            )
        else:
            jacs_score = best_n_good / reward_min_count

        combined = (
            silhouette_weight * silhouette
            + experimental_weight * jacs_score
        )
        if combined > float(best["combined"]):
            best = {
                "silhouette": silhouette,
                "jacs_score": float(jacs_score),
                "combined": float(combined),
                "k": int(k),
                "best_cluster": best_cluster,
                "best_n_good": best_n_good,
                "best_n_jacs": best_n_jacs,
            }

    return {
        "n_features": len(feature_tuple),
        "features": feature_tuple,
        **best,
    }


def run_legacy_feature_search(
    data: pd.DataFrame,
    feature_columns: list[str],
    ice_growth_rate: pd.Series,
    config: dict[str, object],
    *,
    max_combinations: int | None = None,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """Re-run the legacy experimentally guided feature-subset search."""
    scaled = prepare_search_frame(data, feature_columns)
    candidate_pool = legacy_candidate_pool(
        scaled, int(config["candidate_pool"])
    )
    cap = int(max_combinations or config["max_combinations"])
    combinations = list(
        _limited_combinations(
            candidate_pool,
            int(config["subset_min"]),
            int(config["subset_max"]),
            cap,
        )
    )
    rates = ice_growth_rate.to_numpy(dtype=float)

    rows = Parallel(n_jobs=n_jobs)(
        delayed(_score_subset)(
            scaled,
            combination,
            rates,
            k_values=[int(value) for value in config["k_range"]],
            threshold=float(config["igr_threshold"]),
            reward_min_count=int(config["reward_min_count"]),
            silhouette_weight=float(config["silhouette_weight"]),
            experimental_weight=float(config["experimental_weight"]),
            random_state=int(config["random_state"]),
            n_init=int(config["kmeans_n_init"]),
        )
        for combination in combinations
    )
    return pd.DataFrame(rows).sort_values(
        ["combined", "silhouette"],
        ascending=[False, False],
    ).reset_index(drop=True)
