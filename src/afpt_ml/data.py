from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

VALID_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
METADATA_CANDIDATES = {
    "Sequence",
    "Sequence_Clean",
    "Sequence_clean",
    "Group",
    "Accession",
    "ID",
    "Name",
    "Peptide_ID",
    "Published_Main_Cluster",
    "Published_Final_Family",
}


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


def load_feature_dataset(path: str | Path) -> pd.DataFrame:
    """Load the AFPT feature table and add a clean sequence column safely."""
    data = pd.read_csv(path)
    if "Sequence" not in data.columns:
        raise ValueError("The feature table must contain a 'Sequence' column.")

    sequence_clean = data["Sequence"].map(clean_sequence).rename("Sequence_clean")
    if "Sequence_clean" in data.columns:
        data = data.drop(columns=["Sequence_clean"])
    data = pd.concat([data.copy(), sequence_clean], axis=1)

    if data["Sequence_clean"].isna().any():
        raise ValueError("At least one sequence became empty after cleaning.")
    if data["Sequence_clean"].duplicated().any():
        duplicates = data.loc[
            data["Sequence_clean"].duplicated(keep=False), "Sequence_clean"
        ].tolist()
        raise ValueError(f"Cleaned sequences must be unique. Duplicates: {duplicates[:5]}")
    return data


def numeric_feature_columns(data: pd.DataFrame) -> list[str]:
    """Return numeric analysis columns while excluding known metadata columns."""
    return [
        column
        for column in data.select_dtypes(include=[np.number]).columns
        if column not in METADATA_CANDIDATES
    ]


def load_published_features(path: str | Path) -> list[str]:
    """Load one published feature name per line."""
    features = [
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not features:
        raise ValueError("The published feature list is empty.")
    if len(features) != len(set(features)):
        raise ValueError("The published feature list contains duplicates.")
    return features


def validate_feature_columns(data: pd.DataFrame, features: list[str]) -> None:
    """Confirm that requested features exist and are numeric-convertible."""
    missing = sorted(set(features) - set(data.columns))
    if missing:
        raise ValueError(f"Missing requested feature columns: {missing}")

    numeric = data[features].apply(pd.to_numeric, errors="coerce")
    all_missing = numeric.columns[numeric.isna().all()].tolist()
    if all_missing:
        raise ValueError(f"Features contain no numeric values: {all_missing}")


def load_experimental_data(path: str | Path) -> pd.DataFrame:
    """Load the workbook sheet containing sequence and ice-growth-rate data."""
    workbook = pd.ExcelFile(path)
    for sheet_name in workbook.sheet_names:
        table = pd.read_excel(workbook, sheet_name=sheet_name)
        table.columns = [str(column).strip() for column in table.columns]
        if "IceGrowthRate" in table.columns:
            table = table.rename(columns={"IceGrowthRate": "Ice_Growth_Rate"})

        required = {"Sequence", "Ice_Growth_Rate"}
        if not required.issubset(table.columns):
            continue

        table = table.copy()
        table["Sequence_clean"] = table["Sequence"].map(clean_sequence)
        table["Ice_Growth_Rate"] = pd.to_numeric(
            table["Ice_Growth_Rate"], errors="coerce"
        )
        return table.dropna(
            subset=["Sequence_clean", "Ice_Growth_Rate"]
        ).copy()

    raise ValueError(
        "No workbook sheet contained both Sequence and Ice_Growth_Rate."
    )
