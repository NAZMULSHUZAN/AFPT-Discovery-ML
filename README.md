This repository contains the reproducible computational workflow associated with the study **“Machine-Learning-Guided Design of Antifreezing Peptides (AFPTs).”**

The repository integrates experimentally guided feature selection, unsupervised clustering, hierarchical family analysis, reproducibility verification, and post-prediction structural analysis of antifreezing peptides.

## Main AFPT groups

The final analysis resolves four AFPT groups:

- Cluster-1
- C2-Sub-1
- C2-Sub-2
- C2-Sub-3

Cluster-2 was further analyzed by Ward hierarchical clustering because experimentally favorable low ice-growth-rate observations were predominantly concentrated in this main cluster.

## Reproducible analytical workflows

### 1. Published clustering reproduction

The locked published workflow uses:

- 719 AFPT sequences
- 137 available numeric features
- a selected 10-feature representation
- StandardScaler normalization
- row-wise L2 normalization
- K-means clustering with `k = 2`
- Ward hierarchical clustering within the experimentally enriched Cluster-2

The reproduced results are:

- Cluster-1: 407 peptides
- Cluster-2: 312 peptides
- C2-Sub-1: 134 peptides
- C2-Sub-2: 48 peptides
- C2-Sub-3: 130 peptides
- cosine silhouette score: 0.316049
- favorable experimental peptides in the enriched cluster: 89 of 90

Run:

```bash
afpt-cluster published
