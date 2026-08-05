# AFPT Modular Analytical Pipeline

## Supported modes

### `published`

Locked analytical reproduction using the archived 10-feature set selected from
137 available numeric features. This mode must reproduce:

- cosine silhouette: `0.3160492893162923`
- Cluster-1: `407`
- Cluster-2: `312`
- C2-Sub-1: `134`
- C2-Sub-2: `48`
- C2-Sub-3: `130`
- favorable experimental peptides in the enriched cluster: `89/90`

The experimentally enriched main cluster is named Cluster-2. Ward hierarchical
clustering is applied only to this group to resolve three finer families.

### `baseline137`

Re-runs the original all-137-feature baseline:

1. numeric feature identification
2. median imputation
3. row-wise L2 normalization of the unscaled feature matrix
4. KMeans evaluation for k = 2..6
5. selection by cosine silhouette

This is a baseline analysis and is not expected to reproduce the locked final
10-feature families.

### `search`

Re-runs the legacy experimentally guided feature-subset search as an extension.
Results are written to a timestamped directory and never overwrite the locked
published configuration.

The legacy notebook ranked candidate features by variance after StandardScaler.
Because scaled variances are nearly equal, candidate order can be sensitive to
floating-point and software-version differences. Therefore, a new search may
identify a different near-optimal feature set. This does not change the
published reproduction.

## Commands

```bash
python -m pip install -e ".[dev]"
afpt-cluster published
afpt-cluster baseline137
afpt-cluster search --max-combinations 1500 --n-jobs 1
pytest
```

For a quick test of the extension search:

```bash
afpt-cluster search --max-combinations 10 --run-name smoke-test
```
