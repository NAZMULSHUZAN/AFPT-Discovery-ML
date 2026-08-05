"""Verify the published AFPT clustering pipeline.

Workflow
--------
719 peptides × 137 available numeric features
→ archived published 10-feature subset
→ median imputation
→ StandardScaler
→ row-wise L2 normalization
→ KMeans, k=2
→ identify experimentally enriched cluster
→ Ward hierarchical clustering into 3 subclusters
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler, normalize


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "published_search.json"
OUTPUT_DIR = ROOT / "outputs" / "verification"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(value: object) -> str | None:
    """Return an uppercase sequence containing standard amino acids only."""
    if pd.isna(value):
        return None

    cleaned = "".join(
        character
        for character in str(value).strip().upper()
        if character in VALID_AA
    )
    return cleaned or None


def load_experimental_data(path: Path) -> pd.DataFrame:
    """Load the sheet containing Sequence and Ice_Growth_Rate."""
    workbook = pd.ExcelFile(path)

    for sheet in workbook.sheet_names:
        table = pd.read_excel(workbook, sheet_name=sheet)
        table.columns = [str(column).strip() for column in table.columns]

        if "IceGrowthRate" in table.columns:
            table = table.rename(
                columns={"IceGrowthRate": "Ice_Growth_Rate"}
            )

        required = {"Sequence", "Ice_Growth_Rate"}

        if required.issubset(table.columns):
            table["Sequence_clean"] = table["Sequence"].map(clean_sequence)
            table["Ice_Growth_Rate"] = pd.to_numeric(
                table["Ice_Growth_Rate"],
                errors="coerce",
            )

            return table.dropna(
                subset=["Sequence_clean", "Ice_Growth_Rate"]
            ).copy()

    raise ValueError(
        "No Excel sheet contained both Sequence and Ice_Growth_Rate."
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with CONFIG_PATH.open(encoding="utf-8") as handle:
        config = json.load(handle)

    data_path = ROOT / config["input_file"]
    experimental_path = ROOT / config["experimental_file"]
    feature_path = ROOT / config["published_features_file"]

    data = pd.read_csv(data_path)

    published_features = [
        line.strip()
        for line in feature_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    missing_features = [
        feature
        for feature in published_features
        if feature not in data.columns
    ]

    if missing_features:
        raise ValueError(
            f"Published features missing from input: {missing_features}"
        )

    if len(data) != config["dataset"]["n_sequences"]:
        raise ValueError(
            f"Expected {config['dataset']['n_sequences']} rows; "
            f"found {len(data)}."
        )

    data["Sequence_clean"] = data["Sequence"].map(clean_sequence)

    # Exact preprocessing used for the selected-feature clustering.
    feature_table = data[published_features].copy()

    for column in published_features:
        feature_table[column] = pd.to_numeric(
            feature_table[column],
            errors="coerce",
        )
        feature_table[column] = feature_table[column].fillna(
            feature_table[column].median()
        )

    scaled = StandardScaler().fit_transform(feature_table)
    normalized = normalize(scaled)

    search_config = config["feature_search"]

    kmeans = KMeans(
        n_clusters=config["published_result"]["best_k"],
        n_init=search_config["kmeans_n_init"],
        random_state=search_config["random_state"],
    )
    main_labels = kmeans.fit_predict(normalized)

    main_silhouette = silhouette_score(
        normalized,
        main_labels,
        metric=search_config["silhouette_metric"],
    )

    experimental = load_experimental_data(experimental_path)

    # This reproduces the notebook's sequence-to-IGR mapping behavior.
    igr_map = dict(
        zip(
            experimental["Sequence_clean"],
            experimental["Ice_Growth_Rate"],
        )
    )

    igr = data["Sequence_clean"].map(igr_map)
    experimental_mask = igr.notna().to_numpy()
    threshold = search_config["igr_threshold"]

    experimental_summary: list[dict[str, int]] = []

    for cluster_id in sorted(np.unique(main_labels)):
        cluster_mask = main_labels == cluster_id

        n_experimental = int(
            np.sum(cluster_mask & experimental_mask)
        )
        n_favorable = int(
            np.sum(
                cluster_mask
                & experimental_mask
                & (igr.to_numpy() <= threshold)
            )
        )

        experimental_summary.append(
            {
                "raw_cluster": int(cluster_id),
                "size": int(cluster_mask.sum()),
                "n_experimental": n_experimental,
                "n_favorable": n_favorable,
            }
        )

    experimental_summary_df = pd.DataFrame(experimental_summary)

    enriched_cluster = int(
        experimental_summary_df.sort_values(
            ["n_favorable", "n_experimental"],
            ascending=False,
        ).iloc[0]["raw_cluster"]
    )

    cluster_2_mask = main_labels == enriched_cluster
    cluster_1_mask = ~cluster_2_mask

    # Hierarchical clustering is applied only to the experimentally
    # enriched main cluster.
    hierarchy = linkage(
        normalized[cluster_2_mask],
        method=config["subclustering"]["method"],
    )

    subcluster_labels = fcluster(
        hierarchy,
        t=config["subclustering"]["n_subclusters"],
        criterion=config["subclustering"]["criterion"],
    )

    main_names = np.where(
        cluster_2_mask,
        "Cluster-2",
        "Cluster-1",
    )

    family_names = main_names.astype(object)

    cluster_2_indices = np.where(cluster_2_mask)[0]

    for row_index, subcluster_id in zip(
        cluster_2_indices,
        subcluster_labels,
    ):
        family_names[row_index] = f"C2-Sub-{int(subcluster_id)}"

    assignments = pd.DataFrame(
        {
            "Sequence": data["Sequence"],
            "Sequence_clean": data["Sequence_clean"],
            "Raw_Main_Label": main_labels,
            "Final_Main_Cluster": main_names,
            "Final_Family": family_names,
            "Ice_Growth_Rate": igr,
        }
    )

    main_counts = assignments[
        "Final_Main_Cluster"
    ].value_counts().sort_index()

    family_counts = assignments[
        "Final_Family"
    ].value_counts().sort_index()

    target_summary = experimental_summary_df[
        experimental_summary_df["raw_cluster"] == enriched_cluster
    ].iloc[0]

    print("\n=== PUBLISHED AFPT PIPELINE VERIFICATION ===")
    print(f"Input shape: {data.shape}")
    print(f"Published features: {len(published_features)}")
    print(f"Cosine silhouette: {main_silhouette:.6f}")

    print("\nMain clusters:")
    print(main_counts.to_string())

    print("\nExperimentally enriched cluster:")
    print(f"Raw label: {enriched_cluster}")
    print(
        "Favorable experimental peptides: "
        f"{int(target_summary['n_favorable'])}/"
        f"{int(target_summary['n_experimental'])}"
    )

    print("\nFinal families:")
    print(family_counts.to_string())

    # Membership counts expected from the archived published analysis.
    expected_main = {
        "Cluster-1": 407,
        "Cluster-2": 312,
    }

    expected_families = {
        "Cluster-1": 407,
        "C2-Sub-1": 134,
        "C2-Sub-2": 48,
        "C2-Sub-3": 130,
    }

    actual_main = main_counts.to_dict()
    actual_families = family_counts.to_dict()

    assert actual_main == expected_main, (
        f"Main cluster mismatch: {actual_main}"
    )

    assert actual_families == expected_families, (
        f"Family mismatch: {actual_families}"
    )

    assert int(target_summary["n_favorable"]) == 89
    assert int(target_summary["n_experimental"]) == 90

    assert np.isclose(
        main_silhouette,
        config["published_result"]["silhouette"],
        atol=1e-10,
    )

    assignments.to_csv(
        OUTPUT_DIR / "published_cluster_assignments.csv",
        index=False,
    )

    experimental_summary_df.to_csv(
        OUTPUT_DIR / "experimental_cluster_summary.csv",
        index=False,
    )

    summary = {
        "main_silhouette": float(main_silhouette),
        "main_counts": actual_main,
        "family_counts": actual_families,
        "experimentally_enriched_raw_cluster": enriched_cluster,
        "favorable_experimental_peptides": int(
            target_summary["n_favorable"]
        ),
        "experimental_peptides_in_enriched_cluster": int(
            target_summary["n_experimental"]
        ),
        "verification": "PASSED",
    }

    with (
        OUTPUT_DIR / "published_verification_summary.json"
    ).open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print("\nPUBLISHED PIPELINE VERIFICATION PASSED")
    print(f"Outputs: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
