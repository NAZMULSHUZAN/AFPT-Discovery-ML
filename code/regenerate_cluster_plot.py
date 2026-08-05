#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Regenerate the archived AFPT UMAP cluster figure."
    )
    parser.add_argument("--data", default="data/AFPT_cluster_plot_data.csv")
    parser.add_argument("--output", default="outputs")
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    required = [
        "UMAP_1", "UMAP_2", "Final_Family",
        "Experimental_IGR_available"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    family_colors = {
        "Cluster-1": "#E84A5F",
        "C2-Sub-1": "#17BECF",
        "C2-Sub-2": "#D84FD2",
        "C2-Sub-3": "#F4A261",
    }
    family_order = ["Cluster-1", "C2-Sub-1", "C2-Sub-2", "C2-Sub-3"]

    fig, ax = plt.subplots(figsize=(9.6, 5.5))
    for family in family_order:
        sub = df[df["Final_Family"] == family]
        ax.scatter(
            sub["UMAP_1"], sub["UMAP_2"],
            s=22, alpha=0.85,
            color=family_colors.get(family, "gray"),
            edgecolors="none",
            label=f"{family} (n={len(sub)})",
        )

    exp_mask = (
        df["Experimental_IGR_available"]
        .astype(str).str.lower()
        .map({"true": True, "false": False})
        .fillna(False)
    )
    exp = df[exp_mask]
    ax.scatter(
        exp["UMAP_1"], exp["UMAP_2"],
        s=24, alpha=0.95, color="#1F4E99",
        edgecolors="none",
        label=f"Experimental Ice Growth Rate (n={len(exp)})",
    )

    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()

    fig.savefig(out / "AFPT_cluster_plot.png", dpi=600, bbox_inches="tight")
    fig.savefig(out / "AFPT_cluster_plot.pdf", dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(df["Final_Family"].value_counts())
    print(f"Experimental peptides: {len(exp)}")

if __name__ == "__main__":
    main()
