# Data layout

## FASTA identifiers

The FASTA record identifier must exactly match the peptide structure-directory
name.

Example:

```text
data/fasta/C2-Sub-3_sequences.fasta
>pep_001
ACDEFGHIKLMNP
```

Corresponding structure:

```text
structures/representative_models/C2-Sub-3/pep_001/pep_001_rank_001.pdb
```

## Supported PDB naming

The analysis searches in this order:

1. `*rank_001*.pdb`
2. any `.pdb` file

## Derived files

`results/tables/summary_all_peptides.csv` is the main peptide-level output.

`results/tables/secondary_structure_summary.csv` is the cluster-level table used
for the paper-style secondary-structure plot.
