# Machine-Learning-Guided Design of Antifreezing Peptides

This repository contains the reproducible computational workflow associated
with the study:

**“Machine-Learning-Guided Design of Antifreezing Peptides.”**

The repository includes:

- the 719-sequence AFPT dataset with 137 calculated features;
- experimentally guided feature-subset selection;
- exact reproduction of the published 10-feature clustering solution;
- Cluster-2 hierarchical subclustering;
- replay of all 1,500 archived feature combinations;
- automated verification tests;
- paper-figure reproduction scripts;
- post-prediction AlphaFold/ColabFold and DSSP analysis tools.

---

## Final AFPT groups

The published analysis resolved four AFPT groups:

- **Cluster-1**
- **C2-Sub-1**
- **C2-Sub-2**
- **C2-Sub-3**

Ward hierarchical clustering was applied only within Cluster-2 because
experimentally favorable low ice-growth-rate observations were predominantly
concentrated in this main cluster.

---

# Quick reproduction guide

The following steps reproduce the main computational results without requiring
AlphaFold, ColabFold, or DSSP.

## 1. Download the repository

```bash
git clone https://github.com/NAZMULSHUZAN/AFPT-Discovery-ML.git
cd AFPT-Discovery-ML
```

## 2. Create a Python environment

Python 3.10 or newer is required.

### Option A: Using `venv`

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

### Option B: Using Conda

```bash
conda create -n afpt-ml python=3.11 -y
conda activate afpt-ml
```

## 3. Install the package

```bash
python -m pip install -e ".[dev]"
```

Confirm that the command-line tool is available:

```bash
afpt-cluster --help
```

---

# Reproduce the published clustering result

Run:

```bash
afpt-cluster published
```

This workflow performs:

```text
719 AFPT sequences
        ↓
137 available numeric features
        ↓
Published 10-feature subset
        ↓
Median imputation
        ↓
StandardScaler
        ↓
Row-wise L2 normalization
        ↓
K-means clustering, k = 2
        ↓
Experimentally enriched Cluster-2
        ↓
Ward hierarchical clustering
        ↓
Three Cluster-2 subfamilies
```

The expected output is:

```text
Input sequences: 719
Available numeric features: 137
Published selected features: 10
Cosine silhouette score: 0.3160492893

Cluster-1: 407
Cluster-2: 312

C2-Sub-1: 134
C2-Sub-2: 48
C2-Sub-3: 130

Favorable experimental peptides in Cluster-2: 89/90

PUBLISHED PIPELINE VERIFICATION PASSED
```

The generated tables are saved in:

```text
outputs/published_modular/
```

---

# Reproduce the experimentally guided feature selection

The published feature subset was selected using a combined objective:

```text
Combined score
= 70% cosine silhouette score
+ 30% experimental ice-growth-rate enrichment score
```

An experimentally characterized peptide was considered favorable for the
feature-selection objective when:

```text
Ice growth rate ≤ 16.520821
```

To replay the exact 1,500 archived candidate feature combinations, run:

```bash
python scripts/replay_archived_feature_search.py --n-jobs 8
```

Use a smaller number of jobs on computers with fewer CPU cores:

```bash
python scripts/replay_archived_feature_search.py --n-jobs 2
```

The expected best result is:

```text
Selected features: 10
Best k: 2
Cosine silhouette: 0.3160492893
Experimental enrichment score: 1.50
Combined score: 0.6712345025
Favorable experimental peptides: 89
Experimental peptides in the enriched cluster: 90
```

A successful run ends with:

```text
EXPERIMENTALLY GUIDED FEATURE SEARCH REPRODUCTION PASSED
```

The replay verifies that recalculation of all 1,500 archived combinations
produces:

- the same 10 selected features;
- the same `k = 2` solution;
- the same silhouette score;
- the same experimental enrichment;
- the same combined score;
- the same final cluster result.

Replay outputs are written to:

```text
outputs/search_replay/
```

## Important distinction

For exact reproduction of the published feature-selection result, use:

```bash
python scripts/replay_archived_feature_search.py --n-jobs 8
```

The general extension command:

```bash
afpt-cluster search
```

is intended for new exploratory searches. Because the legacy candidate ranking
used variance after standardization, a new search may produce a different
near-optimal feature ordering across software environments. It does not modify
or invalidate the locked published result.

---

# Reproduce the 137-feature baseline

Run:

```bash
afpt-cluster baseline137
```

This analysis uses all 137 numeric features and performs:

1. numeric-feature identification;
2. median imputation;
3. row-wise L2 normalization of the unscaled feature matrix;
4. K-means evaluation for `k = 2–6`;
5. selection using cosine silhouette score.

This baseline is included for comparison. It is not expected to generate the
same families as the selected 10-feature published solution.

Outputs are saved in:

```text
outputs/baseline137/
```

---

# Run a new extension search

A new exploratory search can be run without changing the locked published
configuration.

For a quick smoke test:

```bash
afpt-cluster search \
  --max-combinations 10 \
  --run-name smoke-test \
  --n-jobs 1
```

For a larger extension search:

```bash
afpt-cluster search \
  --max-combinations 1500 \
  --run-name extension-search-1500 \
  --n-jobs 8
```

Results are written to:

```text
outputs/search_runs/<run-name>/
```

Extension searches never overwrite the published configuration or archived
published results.

---

# Reproduce the paper figures

## Main clustering figures

Run:

```bash
bash run_figure_reproduction.sh
```

This script uses the prepared cluster-assignment and plotting tables to
regenerate the main clustering visualization.

Generated figures are written to:

```text
results/figures/
```

The numerical cluster assignments used for verification are also available in:

```text
outputs/verification/published_cluster_assignments.csv
```

The exact published cluster counts are:

```text
Cluster-1: 407
C2-Sub-1: 134
C2-Sub-2: 48
C2-Sub-3: 130
```

## Direct Python verification

The published result can also be reproduced without using the command-line
entry point:

```bash
python scripts/verify_published_pipeline.py
```

Expected final message:

```text
PUBLISHED PIPELINE VERIFICATION PASSED
```

---

# Automated tests

Run:

```bash
pytest
```

The tests verify:

- 719 input sequences;
- 719 unique sequences;
- 137 numeric features;
- absence of missing required data;
- presence of all 10 published features;
- exact Cluster-1 and Cluster-2 sizes;
- exact Cluster-2 subfamily sizes;
- the published cosine silhouette score;
- the experimental enrichment result of 89/90.

Expected result:

```text
3 passed
```

Code-style checking can be performed with:

```bash
ruff check src tests scripts
```

---

# Structural-analysis workflow

AlphaFold or ColabFold prediction is performed externally. The repository does
not run AlphaFold automatically.

After predicted PDB structures are available, this repository can:

- validate FASTA and structure-folder organization;
- detect representative PDB files;
- extract AlphaFold pLDDT values from PDB B-factors;
- calculate cluster-level pLDDT summaries;
- calculate DSSP secondary-structure composition;
- regenerate cluster-level structural figures.

## Required FASTA files

Place the following files in:

```text
data/fasta/
```

Required filenames:

```text
Cluster-1_sequences.fasta
C2-Sub-1_sequences.fasta
C2-Sub-2_sequences.fasta
C2-Sub-3_sequences.fasta
```

Example FASTA record:

```fasta
>peptide_001
ACDEFGHIKLMNP
```

The FASTA record ID should match the corresponding structure-directory name.

## Required structure layout

Place predicted structures in:

```text
structures/representative_models/
```

Example:

```text
structures/representative_models/
├── Cluster-1/
│   └── peptide_001/
│       └── peptide_001_rank_001.pdb
├── C2-Sub-1/
├── C2-Sub-2/
└── C2-Sub-3/
```

The analysis first searches for:

```text
*rank_001*.pdb
```

If that file is unavailable, it uses another `.pdb` file in the peptide
directory as a fallback.

## Check the structural-analysis setup

Run:

```bash
python scripts/check_setup.py
```

This reports:

- missing FASTA files;
- missing structure directories;
- missing Python dependencies;
- DSSP availability;
- configuration problems.

## Run AlphaFold-output analysis

```bash
bash run_analysis.sh
```

Equivalent Python command:

```bash
python scripts/analyze_af2_results.py --config config.json
```

The main output is:

```text
results/tables/summary_all_peptides.csv
```

This table can contain:

- peptide identifier;
- AFPT cluster or family;
- representative PDB path;
- mean pLDDT;
- minimum pLDDT;
- maximum pLDDT;
- helix fraction;
- beta-sheet fraction;
- turn fraction;
- coil fraction.

When DSSP is unavailable, pLDDT extraction can still run, but
secondary-structure columns are recorded as missing values.

---

# Install DSSP for secondary-structure analysis

DSSP is optional and is required only for secondary-structure assignment.

On Ubuntu or Debian:

```bash
sudo apt update
sudo apt install dssp
```

Check the installation:

```bash
mkdssp --version
```

Some systems provide the command as:

```bash
dssp
```

The executable can be changed in `config.json`:

```json
{
  "dssp_executable": "mkdssp"
}
```

DSSP assignments are summarized as:

- helix;
- beta sheet;
- turn;
- coil or random coil.

---

# Generate the secondary-structure figure

After the structure-analysis table has been generated, run:

```bash
bash run_secondary_structure_plot.sh
```

Equivalent Python command:

```bash
python scripts/secondary_structure_clusters.py --config config.json
```

Figures are written to:

```text
results/figures/
```

---

# Configuration

The structural-analysis paths are defined in `config.json`.

Example:

```json
{
  "fasta_dir": "data/fasta",
  "structures_dir": "structures/representative_models",
  "summary_csv": "results/tables/summary_all_peptides.csv",
  "secondary_structure_csv": "results/tables/secondary_structure_summary.csv",
  "figures_dir": "results/figures",
  "dssp_executable": "mkdssp"
}
```

All paths are interpreted relative to the repository root.

The locked clustering configuration is stored separately in:

```text
configs/published_search.json
```

The published 10-feature list is stored in:

```text
configs/published_features.txt
```

The archived 1,500-combination search is stored in:

```text
data/archive/feature_search_1500_archived.csv
```

---

# Important data files

```text
data/processed/AFPT_137_features.csv
```

Contains:

- 719 AFPT sequences;
- sequence metadata;
- 137 numeric features.

```text
data/raw/experimental/Reg_sequences.xlsx
```

Contains the experimental ice-growth-rate measurements used for experimental
guidance and enrichment calculations.

```text
outputs/verification/published_cluster_assignments.csv
```

Contains the reproduced final cluster and family assignments.

```text
docs/archived_search_replay_summary.json
```

Contains the compact verification summary for the exact 1,500-combination
experimental replay.

---

# Data and supplementary materials

The complete supplementary dataset and supporting materials are archived on
Zenodo:

**Zenodo DOI:** `https://doi.org/ZENODO_DOI_HERE`

The minimal files required for verification are included in this repository.
The Zenodo record serves as the archived supplementary-data record associated
with the study.

---

# Current reproducibility coverage

## Included in the current release

- published 10-feature clustering;
- 137-feature baseline clustering;
- exact experimental feature-search replay;
- main cluster and family assignments;
- figure-regeneration workflow;
- AlphaFold pLDDT post-processing;
- optional DSSP secondary-structure analysis;
- automated tests.

## Not yet included

The active-family classifier defined using:

```text
Positive:
C2-Sub-2 + C2-Sub-3

Negative:
Cluster-1 + C2-Sub-1
```

is not yet included in the current modular release. It should not be claimed as
reproducible from this repository until its original preprocessing,
cross-validation, model settings, and evaluation outputs have been separately
verified and added.

---

# Recommended reproduction order

For the fastest complete verification, run:

```bash
python -m pip install -e ".[dev]"

pytest

afpt-cluster published

python scripts/replay_archived_feature_search.py --n-jobs 8

bash run_figure_reproduction.sh
```

For structure-based figures, additionally run:

```bash
python scripts/check_setup.py

bash run_analysis.sh

bash run_secondary_structure_plot.sh
```

---

# Expected verification messages

A complete successful analytical reproduction should include:

```text
3 passed

PUBLISHED PIPELINE VERIFICATION PASSED

EXPERIMENTALLY GUIDED FEATURE SEARCH REPRODUCTION PASSED
```

---

# Citation

When using this repository, please cite the associated manuscript and the
Zenodo supplementary-data record.

The software version used for a study should also be identified by its GitHub
release or commit hash.
