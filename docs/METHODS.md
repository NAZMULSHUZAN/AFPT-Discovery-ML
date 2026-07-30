# Methods: Structural Analysis

Representative AlphaFold/ColabFold models were organized according to the
machine-learning-defined peptide groups Cluster-1, C2-Sub-1, C2-Sub-2, and
C2-Sub-3. For each peptide, the rank-001 PDB model was selected when available.

AlphaFold confidence was summarized from the residue-level pLDDT values stored
in the PDB B-factor field. Secondary structure was assigned using DSSP and
classified as helix (H, G, or I), beta sheet (E or B), turn (T or S), or coil
(all remaining DSSP states). Cluster-level percentages were calculated as the
mean peptide-level fraction for each structural class.

An exploratory N–O proximity metric was calculated using a 3.5 Å distance
cutoff. Because donor–acceptor geometry and residue connectivity were not
evaluated, this measure should be described as an N–O contact metric rather
than a definitive hydrogen-bond count.

The structure-prediction and DSSP installation steps are external to this
repository. The supplied scripts reproduce all subsequent data processing,
summary generation, and plotting steps.
