from pathlib import Path

from afpt_ml.config import load_config, resolve_from_root
from afpt_ml.data import (
    load_feature_dataset,
    load_published_features,
    numeric_feature_columns,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = load_config(ROOT / "configs" / "published_search.json")


def test_locked_input_integrity() -> None:
    data = load_feature_dataset(resolve_from_root(ROOT, CONFIG["input_file"]))
    features = load_published_features(
        resolve_from_root(ROOT, CONFIG["published_features_file"])
    )
    assert len(data) == int(CONFIG["dataset"]["n_sequences"])
    assert data["Sequence_clean"].nunique() == int(CONFIG["dataset"]["n_sequences"])
    assert len(numeric_feature_columns(data)) == int(CONFIG["dataset"]["n_numeric_features"])
    assert len(features) == int(CONFIG["published_result"]["n_selected_features"])
    assert set(features).issubset(data.columns)
