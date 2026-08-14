Machine-Learning-Guided Design of Antifreezing Peptides

Selected for Journal Cover and 
published in the Journal of Chemical Information and Modeling (JCIM), 2026.
Link:(https://doi.org/10.1021/acs.jcim.6c01712)

* Machine Learning 
* Bioinformatics 
* Protein Design 
* Experimental Validation 
* Reproducible Research

My Contributions

* ML pipeline development
* Feature engineering/selection
* Clustering and validation
* Biological interpretation
* Reproducibility pipeline
* Peptide-design analysis

This repository contains the reproducible computational workflow associated
with the study **“Machine-Learning-Guided Design of Antifreezing Peptides.”**

The workflow integrates:

- experimentally guided feature selection;
- unsupervised clustering of 719 antifreezing peptides;
- hierarchical identification of AFPT families;
- exact reproduction of the published clustering results;
- MEME-based cluster-specific motif discovery;
- motif-informed de novo AFPT design;
- AlphaFold/ColabFold post-prediction analysis;
- optional DSSP secondary-structure analysis;
- automated verification tests.

## Final AFPT groups

The published analysis resolved four AFPT groups:

- Cluster-1
- C2-Sub-1
- C2-Sub-2
- C2-Sub-3

Cluster-2 was further divided using Ward hierarchical clustering because
experimentally favorable low ice-growth-rate observations were predominantly
concentrated in this main cluster.

---

## Quick installation

Clone the repository:

```bash
git clone https://github.com/NAZMULSHUZAN/AFPT-Discovery-ML.git
cd AFPT-Discovery-ML
```

Create and activate a Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

Install the package:

```bash
python -m pip install -e ".[dev]"
```

Confirm that the command-line tool is available:

```bash
afpt-cluster --help
```

---

## Reproduce the published clustering result

Run:

```bash
afpt-cluster published
```

The locked workflow performs:

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

Expected results:

```text
Cluster-1: 407
Cluster-2: 312

C2-Sub-1: 134
C2-Sub-2: 48
C2-Sub-3: 130

Cosine silhouette score: 0.316049
Favorable experimental peptides in the enriched cluster: 89/90
```

A successful run ends with:

```text
PUBLISHED PIPELINE VERIFICATION PASSED
```

Generated outputs are saved in:

```text
outputs/published_modular/
```

---

## Reproduce the experimentally guided feature selection

Candidate feature subsets were ranked using a combined objective:

```text
Combined score
= 70% cosine silhouette score
+ 30% experimental ice-growth-rate enrichment score
```

Peptides with an ice growth rate of:

```text
Ice growth rate ≤ 16.520821
```

were treated as favorable for the experimental-enrichment calculation.

To replay all 1,500 archived feature combinations, run:

```bash
python scripts/replay_archived_feature_search.py --n-jobs 8
```

Expected best result:

```text
Selected features: 10
Best k: 2
Cosine silhouette: 0.3160492893
Experimental enrichment score: 1.50
Combined score: 0.6712345025
Favorable experimental peptides: 89
Experimental peptides in the enriched cluster: 90
```

A successful replay ends with:

```text
EXPERIMENTALLY GUIDED FEATURE SEARCH REPRODUCTION PASSED
```

This confirms that recalculation of all 1,500 archived candidate combinations
reproduces the same selected 10 features, clustering solution, experimental
enrichment, and combined score.

---

## Reproduce the 137-feature baseline

Run:

```bash
afpt-cluster baseline137
```

This workflow uses all 137 numeric features and evaluates K-means solutions
using cosine silhouette score.

The 137-feature baseline is included for comparison and is not expected to
produce the same families as the selected 10-feature published solution.

Outputs are saved in:

```text
outputs/baseline137/
```

---

## Run a new extension search

A new exploratory feature search can be run without modifying the locked
published configuration.

Quick smoke test:

```bash
afpt-cluster search \
  --max-combinations 10 \
  --run-name smoke-test \
  --n-jobs 1
```

Larger extension search:

```bash
afpt-cluster search \
  --max-combinations 1500 \
  --run-name extension-search-1500 \
  --n-jobs 8
```

Extension-search outputs are saved separately in:

```text
outputs/search_runs/<run-name>/
```

For exact reproduction of the published feature-selection result, use the
archived-search replay script rather than the general extension-search command.

---

## MEME motif analysis and de novo AFPT design

Cluster-specific sequence motif discovery was performed using the MEME Suite.

Separate FASTA files were prepared for:

```text
Cluster-1_sequences.fasta
C2-Sub-1_sequences.fasta
C2-Sub-2_sequences.fasta
C2-Sub-3_sequences.fasta
```

These FASTA files are provided in:

```text
data/fasta/
```

Each cluster or subcluster FASTA file was analyzed separately using MEME to
identify conserved motifs and residue-position patterns.

The design workflow was:

```text
Final AFPT cluster assignments
        ↓
Cluster-specific FASTA files
        ↓
MEME motif discovery
        ↓
Identification of conserved sequence patterns
        ↓
Integration with experimental enrichment,
physicochemical properties, and structural interpretation
        ↓
De novo AFPT candidate design
```

MEME was used as a motif-discovery tool. The de novo peptides were not generated
automatically by MEME. Instead, conserved motifs and residue patterns identified
from favorable AFPT families were used to guide rational peptide design.

The cluster-specific FASTA files allow users to repeat the MEME analysis.
Exact motif reproduction additionally requires use of the same MEME version,
motif-number setting, motif-width range, background model, and sequence
distribution settings used in the original analysis.

When available, the original MEME settings and output files are archived in
the associated Zenodo supplementary-data record.

---

## Reproduce the paper figures

Run:

```bash
bash run_figure_reproduction.sh
```

Generated figures are written to:

```text
results/figures/
```

The reproduced numerical cluster assignments are available in:

```text
outputs/verification/published_cluster_assignments.csv
```

The published pipeline can also be verified directly using:

```bash
python scripts/verify_published_pipeline.py
```

Expected final message:

```text
PUBLISHED PIPELINE VERIFICATION PASSED
```

---

## Automated verification

Run:

```bash
pytest
```

The tests verify:

- 719 input and unique AFPT sequences;
- 137 available numeric features;
- presence of all 10 published features;
- exact Cluster-1 and Cluster-2 sizes;
- exact Cluster-2 subfamily sizes;
- published cosine silhouette score;
- experimental enrichment of 89 favorable peptides among 90 experimentally
  characterized peptides in the enriched cluster.

Expected result:

```text
3 passed
```

Code quality can be checked using:

```bash
ruff check src tests scripts
```

---

## AlphaFold/ColabFold structural analysis

AlphaFold or ColabFold structure prediction is performed externally and is not
run automatically by this repository.

After predicted PDB files are available, the repository can:

- validate the FASTA and structure-folder layout;
- detect representative PDB models;
- extract AlphaFold pLDDT values from PDB B-factors;
- calculate cluster-level confidence summaries;
- calculate DSSP secondary-structure composition;
- regenerate cluster-level structural figures.

Required FASTA files:

```text
data/fasta/
├── Cluster-1_sequences.fasta
├── C2-Sub-1_sequences.fasta
├── C2-Sub-2_sequences.fasta
└── C2-Sub-3_sequences.fasta
```

Example structure layout:

```text
structures/representative_models/
├── Cluster-1/
│   └── peptide_001/
│       └── peptide_001_rank_001.pdb
├── C2-Sub-1/
├── C2-Sub-2/
└── C2-Sub-3/
```

Check the local setup:

```bash
python scripts/check_setup.py
```

Run the structural analysis:

```bash
bash run_analysis.sh
```

Equivalent command:

```bash
python scripts/analyze_af2_results.py --config config.json
```

Main output:

```text
results/tables/summary_all_peptides.csv
```

---

## Optional DSSP analysis

DSSP is required only for secondary-structure assignment.

Ubuntu or Debian installation:

```bash
sudo apt update
sudo apt install dssp
```

Check the executable:

```bash
mkdssp --version
```

Generate the cluster-level secondary-structure figure:

```bash
bash run_secondary_structure_plot.sh
```

Equivalent command:

```bash
python scripts/secondary_structure_clusters.py --config config.json
```

The DSSP results summarize:

- helix fraction;
- beta-sheet fraction;
- turn fraction;
- coil or random-coil fraction.

---

## Recommended reproduction order

To reproduce the main paper analysis, run:

```bash
python -m pip install -e ".[dev]"

pytest

afpt-cluster published

python scripts/replay_archived_feature_search.py --n-jobs 8

bash run_figure_reproduction.sh
```

Expected successful messages:

```text
3 passed

PUBLISHED PIPELINE VERIFICATION PASSED

EXPERIMENTALLY GUIDED FEATURE SEARCH REPRODUCTION PASSED
```

For structural analysis, additionally run:

```bash
python scripts/check_setup.py

bash run_analysis.sh

bash run_secondary_structure_plot.sh
```

---

## Data and supplementary materials

The supplementary dataset, calculated features, experimental ice-growth-rate
data, final cluster and family labels, MEME input files, model settings, and
supporting materials are archived on Zenodo:

**Zenodo DOI:** `https://doi.org/10.5281/zenodo.21269941`

The reproducible analysis code is maintained in this GitHub repository:

**GitHub:**  
`https://github.com/NAZMULSHUZAN/AFPT-Discovery-ML`

---

## Current reproducibility coverage

Included in the current repository:

- exact published 10-feature clustering;
- all-137-feature baseline analysis;
- exact 1,500-combination experimental feature-search replay;
- final cluster and family assignments;
- paper-figure regeneration;
- cluster-specific FASTA files for MEME motif analysis;
- documentation of the MEME-informed de novo design workflow;
- AlphaFold pLDDT post-processing;
- optional DSSP secondary-structure analysis;
- automated verification tests.

The active-family classifier using:

```text
Positive: C2-Sub-2 + C2-Sub-3
Negative: Cluster-1 + C2-Sub-1
```

will be added after its original preprocessing, cross-validation settings,
model parameters, and evaluation outputs are independently verified.

---

## Citation

## Publication

**Machine-Learning-Guided Design of Antifreezing Peptides**  
Nazmul Shuzan, Jialun Wei, and Jie Zheng  
*Journal of Chemical Information and Modeling* (2026)

[Read the paper](https://doi.org/10.1021/acs.jcim.6c01712)
