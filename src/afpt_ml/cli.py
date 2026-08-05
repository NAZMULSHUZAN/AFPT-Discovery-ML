from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from afpt_ml.config import load_config, resolve_from_root
from afpt_ml.data import (
    load_experimental_data,
    load_feature_dataset,
    numeric_feature_columns,
)
from afpt_ml.experimental import align_ice_growth_rates
from afpt_ml.search import run_legacy_feature_search
from afpt_ml.verification import (
    run_baseline137_pipeline,
    run_published_pipeline,
    save_published_result,
)


def _json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, AttributeError):
            pass
    return value


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(
        json.dumps(_json_safe(payload), indent=2) + "\n", encoding="utf-8"
    )


def _published(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    result = run_published_pipeline(root, config_path=args.config, verify=True)
    output = resolve_from_root(root, args.output)
    save_published_result(result, output)
    print(json.dumps(result.summary, indent=2))
    print("\nPUBLISHED PIPELINE VERIFICATION PASSED")


def _baseline137(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    assignments, summary = run_baseline137_pipeline(
        root, config_path=args.config
    )
    output = resolve_from_root(root, args.output)
    output.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(output / "baseline137_assignments.csv", index=False)
    _write_json(output / "baseline137_summary.json", summary)
    print(json.dumps(summary, indent=2))


def _search(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    config = load_config(resolve_from_root(root, args.config))
    data = load_feature_dataset(resolve_from_root(root, config["input_file"]))
    experimental = load_experimental_data(
        resolve_from_root(root, config["experimental_file"])
    )
    rates = align_ice_growth_rates(data, experimental)
    features = numeric_feature_columns(data)
    search = run_legacy_feature_search(
        data,
        features,
        rates,
        config["feature_search"],
        max_combinations=args.max_combinations,
        n_jobs=args.n_jobs,
    )

    run_name = args.run_name or datetime.now(timezone.utc).strftime(
        "%Y%m%dT%H%M%SZ"
    )
    output = resolve_from_root(root, args.output) / run_name
    output.mkdir(parents=True, exist_ok=False)
    search.to_csv(output / "feature_search_results.csv", index=False)
    best_features = list(search.iloc[0]["features"])
    (output / "best_feature_subset.txt").write_text(
        "\n".join(str(value) for value in best_features) + "\n",
        encoding="utf-8",
    )
    summary = {
        "mode": "legacy_extension_search",
        "tested_combinations": len(search),
        "best_result": search.iloc[0].to_dict(),
        "warning": (
            "This extension search does not overwrite the locked published "
            "feature set. The legacy candidate ranking can vary slightly "
            "across software environments."
        ),
    }
    _write_json(output / "search_summary.json", summary)
    print(search.head(10).to_string(index=False))
    print(f"\nSaved extension search to: {output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="afpt-cluster",
        description="AFPT analytical reproduction and extension CLI.",
    )
    parser.add_argument("--root", default=".", help="Repository root directory.")
    parser.add_argument(
        "--config", default="configs/published_search.json"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    published = subparsers.add_parser(
        "published", help="Run the locked exact published pipeline."
    )
    published.add_argument(
        "--output", default="outputs/published_modular"
    )
    published.set_defaults(function=_published)

    baseline = subparsers.add_parser(
        "baseline137", help="Run the original all-137-feature baseline."
    )
    baseline.add_argument("--output", default="outputs/baseline137")
    baseline.set_defaults(function=_baseline137)

    search = subparsers.add_parser(
        "search", help="Run a separate legacy extension feature search."
    )
    search.add_argument("--output", default="outputs/search_runs")
    search.add_argument("--run-name")
    search.add_argument("--max-combinations", type=int)
    search.add_argument("--n-jobs", type=int, default=1)
    search.set_defaults(function=_search)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
