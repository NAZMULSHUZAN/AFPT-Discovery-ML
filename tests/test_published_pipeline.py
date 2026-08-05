from pathlib import Path

import pytest

from afpt_ml.verification import run_published_pipeline

ROOT = Path(__file__).resolve().parents[1]


def test_exact_published_pipeline() -> None:
    result = run_published_pipeline(ROOT, verify=True)
    assert result.summary["verification"] == "PASSED"
    assert result.summary["main_counts"] == {
        "Cluster-1": 407,
        "Cluster-2": 312,
    }
    assert result.summary["family_counts"] == {
        "C2-Sub-1": 134,
        "C2-Sub-2": 48,
        "C2-Sub-3": 130,
        "Cluster-1": 407,
    }
    assert result.summary["main_silhouette_cosine"] == pytest.approx(
        0.3160492893162923, abs=1e-10
    )
    assert result.summary["favorable_experimental_peptides"] == 89
    assert result.summary["experimental_peptides_in_enriched_cluster"] == 90
