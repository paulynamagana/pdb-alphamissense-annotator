"""
This module provides utility functions for handling PDB files, including loading structures,
extracting polypeptide sequences, modifying B-factors, and saving structures.

It uses the BioPython library for PDB parsing and manipulation.
"""

from Bio.PDB import PDBParser, PPBuilder, PDBIO, Structure
from Bio.Seq import Seq
import logging
from typing import Dict, Tuple
import numpy as np
import requests
from pathlib import Path


def load_structure(pdb_file_path: str) -> Structure:
    """
    Load a PDB structure from a file.

    Args:
        pdb_file_path: The absolute path to the PDB file.

    Returns:
        The loaded PDB structure as a Bio.PDB.Structure.Structure object.

    Raises:
        FileNotFoundError: If the PDB file cannot be found.
        Exception: For other errors during file parsing.
    """

    logging.info(f"Loading structure from {pdb_file_path}")
    parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure("protein", pdb_file_path)
        return structure
    except FileNotFoundError:
        logging.error(f"PDB file not found at: {pdb_file_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading PDB file: {e}")
        raise


def get_polypeptide_sequences(structure: Structure) -> Dict[str, Seq]:
    """
    Extracts polypeptide sequences from a PDB structure.

    This function iterates through the chains of a structure and uses the PPBuilder
    to identify and build polypeptide sequences.

    Args:
        structure: The Bio.PDB.Structure.Structure object.

    Returns:
        A dictionary mapping chain IDs to their corresponding polypeptide sequences
        as Bio.Seq.Seq objects.
    """

    logging.info("Extracting polypeptide sequences from structure")
    ppb = PPBuilder()
    sequences = {}

    for model in structure:
        for chain in model:
            # I need to create a list of all polypeptide fragments found in a chain, in the case there's a gap
            polypeptide_fragments = [
                pp.get_sequence() for pp in ppb.build_peptides(chain)
            ]

            # if there are fragments, then join them
            if polypeptide_fragments:
                full_sequence = sum(polypeptide_fragments, Seq(""))
                sequences[chain.id] = full_sequence
                logging.info(f"Found sequence for chain {chain.id}: {full_sequence}")

    return sequences


def clear_bfactors(structure: Structure):
    """
    Resets the B-factor of every atom in the structure to 0.0.

    Args:
        structure: The Bio.PDB.Structure.Structure object to modify.
    """
    # Iterate through every residue in the structure
    for residue in structure.get_residues():
        # Use get_unpacked_list() to get ALL atoms,to include all alternate locations
        for atom in residue.get_unpacked_list():
            atom.set_bfactor(-1.00)

    logging.info(f"Reset B-factors to -1.00")

    return structure


def map_and_update_bfactors(
    structure: Structure,
    chain_id: str,
    alignment_results: Dict[str, Tuple[object, str, str]],
    am_scores: Dict[int, float],
):
    """
    Combines mapping and B-factor updates into a single function.

    This function uses PPBuilder to ensure the list of residues matches the
    one used for alignment, resolving discrepancies caused by non-standard
    PDB records like ANISOU.
    """
    logging.info(f"Mapping and updating B-factors for chain {chain_id}...")

    # Get the definitive residue list using PPBuilder ---
    ppb = PPBuilder()
    chain_obj = structure[0][chain_id]
    pdb_residues = []
    for pp in ppb.build_peptides(chain_obj):
        pdb_residues.extend(list(pp))

    if not pdb_residues:
        logging.warning(
            f"PPBuilder found no residues for chain {chain_id}. Aborting update."
        )
        return

    # Build the residue map ---
    result_tuple = alignment_results.get(chain_id)
    if not result_tuple:
        logging.warning(f"No alignment result for chain {chain_id}. Aborting update.")
        return

    alignment, _, _ = result_tuple
    residue_map = {}
    afdb_coords = alignment.coordinates[0]
    pdb_coords = alignment.coordinates[1]

    for i in range(0, len(afdb_coords), 2):
        afdb_start_pos = afdb_coords[i]
        pdb_start_idx = pdb_coords[i]
        block_length = afdb_coords[i + 1] - afdb_start_pos

        for j in range(block_length):
            current_pdb_residue = pdb_residues[pdb_start_idx + j]
            pdb_residue_id = current_pdb_residue.get_id()
            afdb_position = (afdb_start_pos + j) + 1
            residue_map[pdb_residue_id] = afdb_position

    # Update B-factors using the new map
    residues_updated = 0
    # We iterate through the master list of all residues in the chain object
    for residue in chain_obj.get_residues():
        res_id = residue.get_id()
        score = None

        if res_id in residue_map:
            afdb_position = residue_map[res_id]
            score = am_scores.get(afdb_position)

        # Apply the score or the default value
        final_bfactor = score if score is not None else -1.00
        for atom in residue.get_atoms():
            atom.set_bfactor(final_bfactor)

        if score is not None:
            residues_updated += 1
            res_name = residue.get_resname()
            res_num = res_id[1]
            logging.debug(
                f"Chain {chain_id}: Updated PDB Res {res_name}{res_num} "
                f"-> AFDB Pos {afdb_position} with score {score:.3f}"
            )


def save_pdb(structure: Structure, original_file: str, output_path: str):
    """
    Saves a PDB structure to a file.

    Args:
        structure: The Bio.PDB.Structure.Structure object to save.
        output_path: The absolute path where the PDB file will be saved.
    """

    io = PDBIO()
    io.set_structure(structure)

    output = output_path + original_file + "_AM_modified.pdb"

    io.save(output)
    logging.info(f"Structure saved to {output}")


def modify_pdb_with_am_scores(
    pdb_url: str, average_scores_array: np.ndarray, output_path: str
):
    """Modifies a PDB file with AlphaMissense pathogenicity scores.

    Args:
        pdb_url (str): URL to the PDB file.
        average_scores_file (numpy.ndarray): Array of average pathogenicity scores.
        output_path: The file path to save the modified PDB.
    """
    if not pdb_url:
        logging.error("PDB URL is not provided.")
        return

    try:
        logging.info("Retrieving PDB file")
        response = requests.get(pdb_url, timeout=10)
        response.raise_for_status()
        pdb_content = response.text.splitlines()
    except requests.exceptions.RequestException as e:
        logging.error(f"Failed to retrieve the PDB file: {e}")
        return

    try:
        logging.info(f"Writing modified PDB data to: {output_path}")
        with open(output_path, "w", encoding="utf-8") as out_file:
            for line in pdb_content:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    residue_number = int(line[22:26].strip())

                    if residue_number < len(average_scores_array) and not np.isnan(
                        average_scores_array[residue_number]
                    ):
                        score = average_scores_array[residue_number]
                        score_str = f"{score:.2f}"
                        while len(score_str) < 6:
                            score_str = " " + score_str
                        modified_line = line[:60] + score_str + line[66:]
                        out_file.write(modified_line + "\n")
                    else:
                        out_file.write(line + "\n")
                else:
                    out_file.write(line + "\n")
    except IOError as e:
        logging.error(f"Error writing to file: {e}")
