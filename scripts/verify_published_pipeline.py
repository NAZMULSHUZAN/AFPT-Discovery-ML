"""Backward-compatible wrapper for the modular AFPT verification package."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from afpt_ml.verification import run_published_pipeline, save_published_result


def main() -> None:
    result = run_published_pipeline(ROOT, verify=True)
    output = ROOT / "outputs" / "verification"
    save_published_result(result, output)

    print("\n=== PUBLISHED AFPT PIPELINE VERIFICATION ===")
    print(f"Input rows: {result.summary['input_rows']}")
    print(
        "Available numeric features: "
        f"{result.summary['available_numeric_features']}"
    )
    print(f"Published features: {result.summary['published_features']}")
    print(
        "Cosine silhouette: "
        f"{result.summary['main_silhouette_cosine']:.6f}"
    )
    print("\nMain clusters:")
    for name, count in result.summary["main_counts"].items():
        print(f"{name}: {count}")
    print("\nFinal families:")
    for name, count in result.summary["family_counts"].items():
        print(f"{name}: {count}")
    print("\nPUBLISHED PIPELINE VERIFICATION PASSED")
    print(f"Outputs: {output}")


if __name__ == "__main__":
    main()
