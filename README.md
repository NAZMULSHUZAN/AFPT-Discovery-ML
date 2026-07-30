# Data-Driven Discovery of Antifreeze Peptides Using Machine Learning

This repository contains a reproducible structural-analysis workflow for comparing
machine-learning-defined antifreeze peptide groups:

- Cluster-1
- C2-Sub-1
- C2-Sub-2
- C2-Sub-3

The repository is designed so that all analysis after structure prediction can be
reproduced from local files. AlphaFold/ColabFold structure prediction and DSSP are
external software steps and are therefore treated as optional prerequisites.

## What is reproducible automatically?

Without AlphaFold or DSSP, the repository can:

- validate the folder and input-file structure;
- detect representative PDB files;
- extract AlphaFold pLDDT values from PDB B-factors;
- calculate sequence length;
- estimate simple N–O contact counts;
- calculate the exploratory amyloid-composition score used in the original script;
- create a combined peptide summary table;
- regenerate cluster-level figures from a previously prepared summary CSV.

With DSSP installed, the same analysis additionally calculates:

- helix fraction;
- beta-sheet fraction;
- turn fraction;
- coil fraction;
- cluster-level secondary-structure composition.

## Repository structure

```text
DataDriven_AFPT_Discovery_ML/
├── README.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── .gitignore
├── config.json
├── run_analysis.sh
├── run_secondary_structure_plot.sh
├── scripts/
│   ├── analyze_af2_results.py
│   ├── secondary_structure_clusters.py
│   └── check_setup.py
├── data/
│   ├── fasta/
│   └── example/
├── structures/
│   └── representative_models/
├── results/
│   ├── tables/
│   └── figures/
└── docs/
    ├── METHODS.md
    └── DATA_LAYOUT.md
```

## Installation

### Option 1: Conda

```bash
conda env create -f environment.yml
conda activate afpt-structure
```

### Option 2: pip

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Optional DSSP installation

DSSP is required only for secondary-structure assignment.

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install dssp
```

Check the executable:

```bash
mkdssp --version
```

Some systems expose it as `dssp` instead of `mkdssp`. You can define the executable
inside `config.json`.

## Required FASTA files

Place these files in `data/fasta/`:

```text
Cluster-1_sequences.fasta
C2-Sub-1_sequences.fasta
C2-Sub-2_sequences.fasta
C2-Sub-3_sequences.fasta
```

Each FASTA record ID should match the corresponding structure folder name.

Example:

```text
>peptide_001
ACDEFGHIKLMNP
```

## Required PDB folder layout

Place one structure directory per peptide:

```text
structures/representative_models/
├── Cluster-1/
│   └── peptide_001/
│       └── peptide_001_rank_001.pdb
├── C2-Sub-1/
├── C2-Sub-2/
└── C2-Sub-3/
```

The script first searches for `*rank_001*.pdb` and then uses any `.pdb` file as a
fallback.

## Configuration

Edit `config.json` only when necessary:

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

All paths are relative to the repository root.

## Check the setup

```bash
python scripts/check_setup.py
```

This reports missing FASTA files, structure directories, Python packages, and
whether DSSP is available.

## Run the full local analysis

```bash
bash run_analysis.sh
```

Equivalent Python command:

```bash
python scripts/analyze_af2_results.py --config config.json
```

If DSSP is unavailable, the script continues and writes `NaN` for the
secondary-structure columns.

Main output:

```text
results/tables/summary_all_peptides.csv
```

## Generate the cluster secondary-structure figure

After DSSP-derived values are available:

```bash
bash run_secondary_structure_plot.sh
```

Equivalent command:

```bash
python scripts/secondary_structure_clusters.py --config config.json
```

Outputs:

```text
results/tables/secondary_structure_summary.csv
results/figures/secondary_structure_clusters.png
```

## Reproducing the paper figure without rerunning DSSP

Commit the following derived table to the repository:

```text
results/tables/summary_all_peptides.csv
```

Then anyone can regenerate the cluster comparison figure without AlphaFold or
DSSP, provided the CSV already contains:

```text
cluster
helix_frac
sheet_frac
turn_frac
coil_frac
```

## Important methodological notes

1. AlphaFold/ColabFold is not installed automatically by this repository.
2. DSSP is not installed automatically because it is system software.
3. pLDDT values are read from the AlphaFold PDB B-factor field.
4. The N–O contact metric is a distance-based exploratory measure, not a
   geometry-validated hydrogen-bond assignment.
5. The included amyloid-composition score is an exploratory sequence score and
   should not be presented as a validated amyloid-prediction model.
6. Representative-structure selection criteria should be reported in the paper
   and README, for example rank-001 model, highest mean pLDDT, cluster medoid, or
   experimentally validated representative.

## Citation

Please cite the associated manuscript when it becomes available.
