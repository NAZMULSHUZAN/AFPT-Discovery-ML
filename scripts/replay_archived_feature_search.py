"""Replay the archived 1,500 experimentally guided feature combinations."""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize

from afpt_ml.config import load_config, resolve_from_root
from afpt_ml.data import (
    load_experimental_data,
    load_feature_dataset,
    numeric_feature_columns,
)
from afpt_ml.experimental import align_ice_growth_rates
from afpt_ml.preprocessing import prepare_search_frame


def parse_features(value: object) -> tuple[str, ...]:
    """Parse the archived tuple stored in the CSV."""
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)

    parsed = ast.literal_eval(str(value))
    return tuple(str(item) for item in parsed)


def score_subset(
    scaled: pd.DataFrame,
    feature_tuple: tuple[str, ...],
    rates: np.ndarray,
    config: dict[str, object],
) -> dict[str, object]:
    """Reproduce the original notebook scoring for one feature subset."""
    matrix = normalize(
        scaled.loc[:, list(feature_tuple)].to_numpy(),
        norm="l2",
    )

    matched = ~np.isnan(rates)

    best = {
        "silhouette": -1.0,
        "jacs_score": 0.0,
        "combined": -1.0,
        "k": None,
        "best_cluster": None,
        "best_n_good": 0,
        "best_n_jacs": 0,
    }

    for k in config["k_range"]:
        labels = KMeans(
            n_clusters=int(k),
            n_init=int(config["kmeans_n_init"]),
            random_state=int(config["random_state"]),
        ).fit_predict(matrix)

        if len(np.unique(labels)) < 2:
            continue

        silhouette = float(
            silhouette_score(matrix, labels, metric="cosine")
        )

        # Exact legacy behavior: retain the first cluster having the
        # greatest number of favorable experimental peptides.
        best_cluster = None
        best_n_good = -1
        best_n_jacs = 0

        for cluster_id in np.unique(labels):
            cluster_mask = labels == cluster_id

            n_jacs = int(np.sum(cluster_mask & matched))
            n_good = int(
                np.sum(
                    cluster_mask
                    & matched
                    & (rates <= float(config["igr_threshold"]))
                )
            )

            if n_good > best_n_good:
                best_cluster = int(cluster_id)
                best_n_good = n_good
                best_n_jacs = n_jacs

        reward_minimum = int(config["reward_min_count"])

        if best_n_good >= reward_minimum:
            experimental_score = 1.0 + 0.01 * min(
                best_n_good - reward_minimum,
                50,
            )
        else:
            experimental_score = best_n_good / reward_minimum

        combined = (
            float(config["silhouette_weight"]) * silhouette
            + float(config["experimental_weight"]) * experimental_score
        )

        if combined > float(best["combined"]):
            best = {
                "silhouette": silhouette,
                "jacs_score": float(experimental_score),
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


def canonical_key(value: object) -> tuple[str, ...]:
    """Create an order-independent feature-set key."""
    return tuple(sorted(parse_features(value)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--n-jobs", type=int, default=8)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    config = load_config(root / "configs/published_search.json")
    search_config = config["feature_search"]

    data = load_feature_dataset(
        resolve_from_root(root, config["input_file"])
    )
    experimental = load_experimental_data(
        resolve_from_root(root, config["experimental_file"])
    )

    rates = align_ice_growth_rates(data, experimental)
    numeric_features = numeric_feature_columns(data)
    scaled = prepare_search_frame(data, numeric_features)

    archive_path = (
        root / "data/archive/feature_search_1500_archived.csv"
    )
    archived = pd.read_csv(archive_path)

    archived_combinations = [
        parse_features(value)
        for value in archived["features"]
    ]

    print(f"Archived combinations: {len(archived_combinations)}")
    print(f"Experimental matches: {rates.notna().sum()}")
    print("Recalculating all archived combinations...")

    recalculated_rows = Parallel(n_jobs=args.n_jobs)(
        delayed(score_subset)(
            scaled,
            feature_tuple,
            rates.to_numpy(dtype=float),
            search_config,
        )
        for feature_tuple in archived_combinations
    )

    recalculated = (
        pd.DataFrame(recalculated_rows)
        .sort_values(
            ["combined", "silhouette"],
            ascending=[False, False],
        )
        .reset_index(drop=True)
    )

    output_dir = root / "outputs/search_replay"
    output_dir.mkdir(parents=True, exist_ok=True)

    recalculated.to_csv(
        output_dir / "recalculated_feature_search_1500.csv",
        index=False,
    )

    archived = archived.copy()
    archived["feature_key"] = archived["features"].map(canonical_key)
    recalculated["feature_key"] = recalculated["features"].map(
        canonical_key
    )

    comparison = archived.merge(
        recalculated,
        on="feature_key",
        suffixes=("_archived", "_recalculated"),
        validate="one_to_one",
    )

    metric_columns = [
        "silhouette",
        "jacs_score",
        "combined",
    ]

    integer_columns = [
        "n_features",
        "k",
        "best_cluster",
        "best_n_good",
        "best_n_jacs",
    ]

    metric_differences = {
        column: float(
            np.max(
                np.abs(
                    comparison[f"{column}_archived"]
                    - comparison[f"{column}_recalculated"]
                )
            )
        )
        for column in metric_columns
    }

    exact_integer_match = {
        column: bool(
            (
                comparison[f"{column}_archived"]
                == comparison[f"{column}_recalculated"]
            ).all()
        )
        for column in integer_columns
    }

    archived_top = archived.sort_values(
        ["combined", "silhouette"],
        ascending=[False, False],
    ).iloc[0]

    recalculated_top = recalculated.iloc[0]

    same_best_features = (
        canonical_key(archived_top["features"])
        == canonical_key(recalculated_top["features"])
    )

    expected_features = tuple(
        sorted(
            line.strip()
            for line in (
                root / "configs/published_features.txt"
            ).read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    )

    reproduced_published_features = (
        canonical_key(recalculated_top["features"])
        == expected_features
    )

    top_metrics_match = all(
        np.isclose(
            float(archived_top[column]),
            float(recalculated_top[column]),
            atol=1e-10,
            rtol=0,
        )
        for column in metric_columns
    )

    top_counts_match = all(
        int(archived_top[column])
        == int(recalculated_top[column])
        for column in integer_columns
    )

    full_metric_match = all(
        difference <= 1e-10
        for difference in metric_differences.values()
    )

    full_integer_match = all(exact_integer_match.values())

    exact_replay = (
        len(comparison) == 1500
        and same_best_features
        and reproduced_published_features
        and top_metrics_match
        and top_counts_match
        and full_metric_match
        and full_integer_match
    )

    summary = {
        "archived_combinations": len(archived),
        "recalculated_combinations": len(recalculated),
        "matched_combinations": len(comparison),
        "same_best_features": same_best_features,
        "reproduced_published_features": reproduced_published_features,
        "top_metrics_match": top_metrics_match,
        "top_counts_match": top_counts_match,
        "maximum_metric_differences": metric_differences,
        "full_integer_match": full_integer_match,
        "exact_replay": exact_replay,
        "recalculated_top": {
            "features": list(parse_features(recalculated_top["features"])),
            "k": int(recalculated_top["k"]),
            "silhouette": float(recalculated_top["silhouette"]),
            "jacs_score": float(recalculated_top["jacs_score"]),
            "combined": float(recalculated_top["combined"]),
            "best_cluster": int(recalculated_top["best_cluster"]),
            "best_n_good": int(recalculated_top["best_n_good"]),
            "best_n_jacs": int(recalculated_top["best_n_jacs"]),
        },
    }

    with (
        output_dir / "archived_search_replay_summary.json"
    ).open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print("\n=== RECALCULATED TOP RESULT ===")
    print(json.dumps(summary["recalculated_top"], indent=2))

    print("\n=== REPLAY COMPARISON ===")
    print(json.dumps(
        {
            key: value
            for key, value in summary.items()
            if key != "recalculated_top"
        },
        indent=2,
    ))

    if exact_replay:
        print(
            "\nEXPERIMENTALLY GUIDED FEATURE SEARCH "
            "REPRODUCTION PASSED"
        )
    elif same_best_features and top_metrics_match and top_counts_match:
        print(
            "\nTOP PUBLISHED RESULT REPRODUCED, "
            "BUT SOME NON-TOP ARCHIVED ROWS DIFFER"
        )
    else:
        raise SystemExit(
            "\nEXPERIMENTALLY GUIDED FEATURE SEARCH "
            "REPRODUCTION FAILED"
        )


if __name__ == "__main__":
    main()
