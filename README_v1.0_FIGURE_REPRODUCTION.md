# AFPT Figure Reproduction — v1.0

This release reproduces the archived UMAP cluster visualization for the
machine-learning-guided antifreeze peptide study.

## Scope

Version 1.0 provides **figure reproduction only**. It redraws the figure from
archived UMAP coordinates and final cluster assignments. It does not refit UMAP
and does not rerun KMeans or Ward hierarchical clustering.

## Verified composition

| Group | Peptides |
|---|---:|
| Cluster-1 | 407 |
| C2-Sub-1 | 134 |
| C2-Sub-2 | 48 |
| C2-Sub-3 | 130 |
| Total | 719 |

The experimental overlay contains 106 unique peptides.

## Install

```bash
python -m pip install -r requirements-figure.txt
```

## Validate and reproduce

```bash
bash run_figure_reproduction.sh
```

Outputs:

```text
outputs/AFPT_cluster_plot.png
outputs/AFPT_cluster_plot.pdf
```

## Exact reproducibility statement

> The archived UMAP coordinates and final cluster assignments are provided to
> reproduce the reported cluster visualization. Version 1.0 does not refit UMAP
> or rerun the clustering algorithms.

A later release will add analytical clustering reproduction after the original
hierarchical-clustering workflow is fully verified.
