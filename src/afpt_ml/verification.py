from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from afpt_ml.clustering import run_kmeans, select_best_k
from afpt_ml.config import load_config, resolve_from_root
from afpt_ml.data import (
    load_experimental_data,
    load_feature_dataset,
    load_published_features,
    numeric_feature_columns,
    validate_feature_columns,
)
from afpt_ml.experimental import (
    align_ice_growth_rates,
    select_experimentally_enriched_cluster,
    summarize_experimental_by_cluster,
)
from afpt_ml.hierarchy import run_targeted_ward_clustering
from afpt_ml.preprocessing import (
    prepare_baseline137_matrix,
    prepare_published_matrix,
)


@dataclass(frozen=True)
class PublishedVerificationResult:
    assignments: pd.DataFrame
    experimental_summary: pd.DataFrame
    summary: dict[str, object]


def _published_names(
    labels: np.ndarray, enriched_cluster: int
) -> tuple[np.ndarray, np.ndarray]:
    target_mask = labels == enriched_cluster
    main_names = np.where(target_mask, "Cluster-2", "Cluster-1")
    family_names = np.full(len(labels), "Cluster-1", dtype=object)
    return main_names, family_names


def run_published_pipeline(
    root: str | Path,
    *,
    config_path: str | Path = "configs/published_search.json",
    verify: bool = True,
) -> PublishedVerificationResult:
    """Run and optionally assert the locked published AFPT analysis."""
    root_path = Path(root)
    config = load_config(resolve_from_root(root_path, config_path))
    data = load_feature_dataset(resolve_from_root(root_path, config["input_file"]))
    features = load_published_features(
        resolve_from_root(root_path, config["published_features_file"])
    )
    validate_feature_columns(data, features)

    dataset_config = config["dataset"]
    numeric_columns = numeric_feature_columns(data)
    if len(data) != int(dataset_config["n_sequences"]):
        raise AssertionError(
            f"Expected {dataset_config['n_sequences']} rows; found {len(data)}."
        )
    if len(numeric_columns) != int(dataset_config["n_numeric_features"]):
        raise AssertionError(
            "Numeric feature count differs from the locked configuration: "
            f"{len(numeric_columns)}"
        )

    matrix = prepare_published_matrix(data, features)
    search_config = config["feature_search"]
    main_result = run_kmeans(
        matrix,
        k=int(config["published_result"]["best_k"]),
        random_state=int(search_config["random_state"]),
        n_init=int(search_config["kmeans_n_init"]),
    )

    experimental = load_experimental_data(
        resolve_from_root(root_path, config["experimental_file"])
    )
    ice_growth_rate = align_ice_growth_rates(data, experimental)
    experimental_summary = summarize_experimental_by_cluster(
        main_result.labels,
        ice_growth_rate,
        favorable_threshold=float(search_config["igr_threshold"]),
    )
    enriched_cluster = select_experimentally_enriched_cluster(
        experimental_summary
    )

    main_names, family_names = _published_names(
        main_result.labels, enriched_cluster
    )
    target_mask = main_result.labels == enriched_cluster
    hierarchy = run_targeted_ward_clustering(
        matrix,
        target_mask,
        n_subclusters=int(config["subclustering"]["n_subclusters"]),
        method=str(config["subclustering"]["method"]),
        criterion=str(config["subclustering"]["criterion"]),
    )
    for row_index, subcluster_id in zip(
        hierarchy.target_indices,
        hierarchy.subcluster_labels,
        strict=True,
    ):
        family_names[row_index] = f"C2-Sub-{int(subcluster_id)}"

    assignments = pd.DataFrame(
        {
            "Sequence": data["Sequence"],
            "Sequence_clean": data["Sequence_clean"],
            "Raw_Main_Label": main_result.labels,
            "Final_Main_Cluster": main_names,
            "Final_Family": family_names,
            "Ice_Growth_Rate": ice_growth_rate,
        }
    )
    main_counts = {
        str(key): int(value)
        for key, value in assignments["Final_Main_Cluster"]
        .value_counts()
        .sort_index()
        .items()
    }
    family_counts = {
        str(key): int(value)
        for key, value in assignments["Final_Family"]
        .value_counts()
        .sort_index()
        .items()
    }
    enriched_row = experimental_summary.loc[
        experimental_summary["raw_cluster"] == enriched_cluster
    ].iloc[0]

    summary: dict[str, object] = {
        "input_rows": len(data),
        "available_numeric_features": len(numeric_columns),
        "published_features": len(features),
        "main_silhouette_cosine": main_result.silhouette_cosine,
        "main_counts": main_counts,
        "family_counts": family_counts,
        "experimentally_enriched_raw_cluster": int(enriched_cluster),
        "favorable_experimental_peptides": int(enriched_row["n_favorable"]),
        "experimental_peptides_in_enriched_cluster": int(
            enriched_row["n_experimental"]
        ),
    }

    if verify:
        expected_main = {"Cluster-1": 407, "Cluster-2": 312}
        expected_families = {
            "C2-Sub-1": 134,
            "C2-Sub-2": 48,
            "C2-Sub-3": 130,
            "Cluster-1": 407,
        }
        published = config["published_result"]
        assert main_counts == expected_main, main_counts
        assert family_counts == expected_families, family_counts
        assert int(enriched_row["n_favorable"]) == int(
            published["favorable_experimental_peptides"]
        )
        assert int(enriched_row["n_experimental"]) == int(
            published["experimental_peptides_in_selected_cluster"]
        )
        assert np.isclose(
            main_result.silhouette_cosine,
            float(published["silhouette"]),
            atol=1e-10,
        )
        summary["verification"] = "PASSED"

    return PublishedVerificationResult(
        assignments=assignments,
        experimental_summary=experimental_summary,
        summary=summary,
    )


def save_published_result(
    result: PublishedVerificationResult, output_dir: str | Path
) -> None:
    """Save locked published assignments and summaries."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    result.assignments.to_csv(
        output_path / "published_cluster_assignments.csv", index=False
    )
    result.experimental_summary.to_csv(
        output_path / "experimental_cluster_summary.csv", index=False
    )
    (output_path / "published_verification_summary.json").write_text(
        json.dumps(result.summary, indent=2) + "\n",
        encoding="utf-8",
    )


def run_baseline137_pipeline(
    root: str | Path,
    *,
    config_path: str | Path = "configs/published_search.json",
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Run the original all-137-feature baseline clustering workflow."""
    root_path = Path(root)
    config = load_config(resolve_from_root(root_path, config_path))
    data = load_feature_dataset(resolve_from_root(root_path, config["input_file"]))
    features = numeric_feature_columns(data)
    matrix = prepare_baseline137_matrix(data, features)
    search_config = config["feature_search"]
    result = select_best_k(
        matrix,
        k_values=list(range(2, 7)),
        random_state=int(search_config["random_state"]),
        n_init=int(search_config["kmeans_n_init"]),
    )
    assignments = pd.DataFrame(
        {
            "Sequence": data["Sequence"],
            "Sequence_clean": data["Sequence_clean"],
            "Raw_Cluster": result.labels,
        }
    )
    summary: dict[str, object] = {
        "mode": "baseline137",
        "input_rows": len(data),
        "numeric_features": len(features),
        "selected_k": int(result.k),
        "silhouette_cosine": result.silhouette_cosine,
        "raw_counts": {str(key): value for key, value in result.counts.items()},
    }
    return assignments, summary
