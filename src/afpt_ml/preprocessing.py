from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, normalize


def median_imputed_frame(data: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    """Convert selected features to numeric and median-impute missing values."""
    frame = data[features].apply(pd.to_numeric, errors="coerce").copy()
    for column in features:
        median = frame[column].median(skipna=True)
        if pd.isna(median):
            raise ValueError(f"Feature '{column}' has no finite median.")
        frame[column] = frame[column].fillna(median)
    return frame


def prepare_published_matrix(
    data: pd.DataFrame, features: list[str]
) -> np.ndarray:
    """Published selected-feature preprocessing: scale, then row L2 normalize."""
    frame = median_imputed_frame(data, features)
    scaled = StandardScaler().fit_transform(frame)
    return normalize(scaled, norm="l2")


def prepare_baseline137_matrix(
    data: pd.DataFrame, features: list[str]
) -> np.ndarray:
    """Legacy all-feature baseline: median impute, then row L2 normalize."""
    frame = median_imputed_frame(data, features)
    return normalize(frame.to_numpy(dtype=float), norm="l2")


def prepare_search_frame(
    data: pd.DataFrame, features: list[str]
) -> pd.DataFrame:
    """Legacy feature-search preprocessing: median impute and StandardScaler."""
    frame = median_imputed_frame(data, features)
    scaled = StandardScaler().fit_transform(frame)
    return pd.DataFrame(scaled, columns=features, index=data.index)
