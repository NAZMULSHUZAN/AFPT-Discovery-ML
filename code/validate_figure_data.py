#!/usr/bin/env python3
import pandas as pd

PATH = "data/AFPT_cluster_plot_data.csv"
EXPECTED = {
    "Cluster-1": 407,
    "C2-Sub-1": 134,
    "C2-Sub-2": 48,
    "C2-Sub-3": 130,
}

df = pd.read_csv(PATH)
required = {
    "Peptide_ID", "Sequence", "UMAP_1", "UMAP_2",
    "Final_Main_Cluster", "Final_Family",
    "Experimental_IGR_available", "Ice_Growth_Rate", "Low_IGR_hit",
}
missing = sorted(required - set(df.columns))
if missing:
    raise ValueError(f"Missing columns: {missing}")
if len(df) != 719:
    raise ValueError(f"Expected 719 rows, found {len(df)}")
if df["Peptide_ID"].nunique() != 719:
    raise ValueError("Peptide_ID values are not unique.")
if df["Sequence"].nunique() != 719:
    raise ValueError("Sequence values are not unique.")
counts = df["Final_Family"].value_counts().to_dict()
if counts != EXPECTED:
    raise ValueError(f"Family counts differ: {counts}")
exp = (
    df["Experimental_IGR_available"]
    .astype(str).str.lower()
    .map({"true": True, "false": False})
    .fillna(False).sum()
)
if int(exp) != 106:
    raise ValueError(f"Expected 106 experimental peptides, found {exp}")
print("Validation passed.")
print("Rows:", len(df))
print("Family counts:", counts)
print("Experimental peptides:", int(exp))
