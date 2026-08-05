from pathlib import Path

from afpt_ml.config import load_config
from afpt_ml.verification import run_baseline137_pipeline

ROOT = Path(__file__).resolve().parents[1]
CONFIG = load_config(ROOT / "configs" / "published_search.json")


def test_baseline137_is_complete() -> None:
    assignments, summary = run_baseline137_pipeline(ROOT)
    assert len(assignments) == int(CONFIG["dataset"]["n_sequences"])
    assert summary["numeric_features"] == int(CONFIG["dataset"]["n_numeric_features"])
    assert 2 <= summary["selected_k"] <= 6
    assert sum(summary["raw_counts"].values()) == int(CONFIG["dataset"]["n_sequences"])
