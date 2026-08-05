#!/usr/bin/env python3
"""Analyze representative AlphaFold/ColabFold PDB files by AFPT cluster.

DSSP is optional. When unavailable, secondary-structure fields are written as NaN.
"""

from __future__ import annotations

import argparse
import json
import shutil
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser

try:
    from Bio.PDB.DSSP import DSSP
except ImportError:
    DSSP = None


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_config(path: Path) -> dict:
    with path.open() as handle:
        config = json.load(handle)
    root = repo_root()
    for key in ("fasta_dir", "structures_dir", "summary_csv",
                "secondary_structure_csv", "figures_dir"):
        config[key] = root / config[key]
    return config


def load_sequences(fasta_dir: Path, clusters: Iterable[str]) -> dict[tuple[str, str], str]:
    lookup: dict[tuple[str, str], str] = {}
    for cluster in clusters:
        fasta_path = fasta_dir / f"{cluster}_sequences.fasta"
        if not fasta_path.exists():
            print(f"[WARN] FASTA not found: {fasta_path}")
            continue
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            sequence = str(record.seq).replace("-", "").replace(" ", "").upper()
            lookup[(cluster, record.id)] = sequence
    print(f"[INFO] Loaded {len(lookup)} sequences.")
    return lookup


def find_best_pdb(job_dir: Path) -> Path | None:
    ranked = sorted(job_dir.glob("*rank_001*.pdb"))
    if ranked:
        return ranked[0]
    pdbs = sorted(job_dir.glob("*.pdb"))
    return pdbs[0] if pdbs else None


def plddt_stats(pdb_path: Path) -> tuple[float, float, float]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af2", str(pdb_path))
    values = [
        residue["CA"].get_bfactor()
        for model in structure
        for chain in model
        for residue in chain
        if "CA" in residue
    ]
    if not values:
        return np.nan, np.nan, np.nan
    array = np.asarray(values, dtype=float)
    return float(array.mean()), float(array.min()), float(array.max())


def secondary_structure(
    pdb_path: Path, dssp_executable: str
) -> tuple[float, float, float, float, int]:
    if DSSP is None or shutil.which(dssp_executable) is None:
        return np.nan, np.nan, np.nan, np.nan, 0

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af2", str(pdb_path))
    model = next(structure.get_models())

    try:
        dssp = DSSP(model, str(pdb_path), dssp=dssp_executable)
    except Exception as exc:  # noqa: BLE001
        print(f"[WARN] DSSP failed for {pdb_path}: {exc}")
        return np.nan, np.nan, np.nan, np.nan, 0

    codes = [dssp[key][2] for key in dssp]
    n = len(codes)
    if n == 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    helix = sum(code in {"H", "G", "I"} for code in codes) / n
    sheet = sum(code in {"E", "B"} for code in codes) / n
    turn = sum(code in {"T", "S"} for code in codes) / n
    coil = sum(code not in {"H", "G", "I", "E", "B", "T", "S"} for code in codes) / n
    return helix, sheet, turn, coil, n


def extract_atom_coordinates(pdb_path: Path, atom_name: str) -> list[np.ndarray]:
    coordinates = []
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != atom_name:
                continue
            try:
                coordinates.append(
                    np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        dtype=float,
                    )
                )
            except ValueError:
                continue
    return coordinates


def estimate_no_contacts(pdb_path: Path, cutoff: float = 3.5) -> int:
    nitrogens = extract_atom_coordinates(pdb_path, "N")
    oxygens = extract_atom_coordinates(pdb_path, "O")
    if not nitrogens or not oxygens:
        return 0
    oxygen_array = np.asarray(oxygens)
    return int(
        sum(np.sum(np.linalg.norm(oxygen_array - nitrogen, axis=1) < cutoff)
            for nitrogen in nitrogens)
    )


HYDROPHOBIC = set("IVLFWY")
BETA_FORMERS = set("VILFWY")


def exploratory_amyloid_score(sequence: str) -> float:
    if not sequence:
        return np.nan
    return float(
        sum(residue in HYDROPHOBIC for residue in sequence)
        + sum(residue in BETA_FORMERS for residue in sequence)
    )


def iter_jobs(structures_dir: Path, clusters: Iterable[str]):
    for cluster in clusters:
        cluster_dir = structures_dir / cluster
        if not cluster_dir.exists():
            print(f"[WARN] Structure folder not found: {cluster_dir}")
            continue
        for job_dir in sorted(path for path in cluster_dir.iterdir() if path.is_dir()):
            yield cluster, job_dir


def analyze(config: dict) -> pd.DataFrame:
    fasta_dir: Path = config["fasta_dir"]
    structures_dir: Path = config["structures_dir"]
    summary_csv: Path = config["summary_csv"]
    clusters = config["clusters"]
    dssp_executable = config.get("dssp_executable", "mkdssp")

    sequences = load_sequences(fasta_dir, clusters)
    dssp_available = DSSP is not None and shutil.which(dssp_executable) is not None
    print(f"[INFO] DSSP available: {dssp_available}")

    rows = []
    for cluster, job_dir in iter_jobs(structures_dir, clusters):
        peptide_id = job_dir.name
        pdb_path = find_best_pdb(job_dir)
        if pdb_path is None:
            print(f"[WARN] No PDB found in {job_dir}")
            continue

        sequence = sequences.get((cluster, peptide_id), "")
        mean_plddt, min_plddt, max_plddt = plddt_stats(pdb_path)
        helix, sheet, turn, coil, n_residues = secondary_structure(
            pdb_path, dssp_executable
        )
        no_contacts = estimate_no_contacts(pdb_path)
        denominator = n_residues or len(sequence)
        no_contacts_per_residue = (
            no_contacts / denominator if denominator else np.nan
        )

        rows.append(
            {
                "cluster": cluster,
                "peptide_id": peptide_id,
                "sequence": sequence,
                "length": len(sequence),
                "n_residues_dssp": n_residues,
                "mean_plddt": mean_plddt,
                "min_plddt": min_plddt,
                "max_plddt": max_plddt,
                "helix_frac": helix,
                "sheet_frac": sheet,
                "turn_frac": turn,
                "coil_frac": coil,
                "no_contacts_3p5A": no_contacts,
                "no_contacts_per_residue": no_contacts_per_residue,
                "exploratory_amyloid_score": exploratory_amyloid_score(sequence),
                "pdb_path": str(pdb_path.relative_to(repo_root())),
            }
        )

    if not rows:
        raise RuntimeError(
            "No structures were analyzed. Check the FASTA and structure folder layout."
        )

    dataframe = pd.DataFrame(rows)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(summary_csv, index=False)
    print(f"[SAVED] {summary_csv}")
    return dataframe


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=repo_root() / "config.json",
        help="Path to repository configuration JSON.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    analyze(load_config(arguments.config))
