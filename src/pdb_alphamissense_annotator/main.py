"""
Main workflow to annotate a PDB with AlphaMissense scores.

stitched everything together

Args:
    pdb_file (str): Path to input PDB file
    uniprot_id (str): UniProt accession
    chains (list of str): Chains to annotate
    output_dir(str): Annotated PDB output
"""

import logging
import sys
import os
from typing import Tuple, List, Dict, Any
from pathlib import Path
import numpy as np

# Add the root directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from src.pdb_alphamissense_annotator.am_utils import (
    get_alphamissense_info,
    process_alphamissense_data,
    calculate_average_pathogenicity,
    extract_plddt_scores_from_url,
    extract_am_data,
)
from src.pdb_alphamissense_annotator.pdb_utils import (
    load_structure,
    get_polypeptide_sequences,
    clear_bfactors,
    map_and_update_bfactors,
    save_pdb,
    modify_pdb_with_am_scores,
)
from src.pdb_alphamissense_annotator.sequence_utils import (
    align_sequences,
    save_alignment_report,
)

from src.pdb_alphamissense_annotator.plot_utils import plot_scores, plot_am_heatmap


# configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def annotate_pdb(
    pdb_file_path: str, uniprot_id: str, chains: List[str], output_path: str
):
    """
    Main workflow to annotate a PDB with AlphaMissense scores.

    Args:
        pdb_file (str): Path to input PDB file
        uniprot_id (str): UniProt accession
        chains (list of str): Chains to annotate
        output_path (str): Annotated PDB path
    """

    structure = load_structure(pdb_file_path)  # Loading PDB structure

    uniprot_description, afdb_sequence, am_url, pdbUrl = get_alphamissense_info(
        uniprot_id
    )  # Loading AlphaMissense scores for UniProt accession
    am_data = process_alphamissense_data(am_url)
    am_scores_dict = calculate_average_pathogenicity(am_data)
    logging.info(f"Obtained AlphaMissense scores {uniprot_id}: {uniprot_description}")

    logging.info(f"Annotating chains with AM scores...")
    pdb_sequences = get_polypeptide_sequences(structure)

    all_alignments = align_sequences(pdb_sequences, afdb_sequence, chains)
    save_alignment_report(
        all_alignments, f"{output_path}report_{Path(pdb_file_path).stem}.txt"
    )

    structure_clear = clear_bfactors(structure)

    for chain in chains:
        map_and_update_bfactors(structure_clear, chain, all_alignments, am_scores_dict)

    save_pdb(structure_clear, Path(pdb_file_path).stem, output_path)


def afdb_structure_am_coloured(uniprot_ids: List[str], output_path: str):
    """Modifies a PDB file with AlphaMissense pathogenicity scores.

    Args:
        pdb_url (str): URL to the PDB file.
        average_scores_file (numpy.ndarray): Array of average pathogenicity scores.
    """

    for uniprot_id in uniprot_ids:
        try:
            _, _, am_url, pdb_url = get_alphamissense_info(uniprot_id)
            if not am_url or not pdb_url:
                logging.error(
                    f"Could not retrieve valid URLs for {uniprot_id}. Exiting."
                )
                continue

        except Exception as e:
            logging.error(f"Failed to get info for {uniprot_id}: {e}. Skipping.")
            continue

        am_data = process_alphamissense_data(am_url)
        avg_scores_dict = calculate_average_pathogenicity(am_data)
        plddt_scores_dict = extract_plddt_scores_from_url(pdb_url)
        am_data_complete = extract_am_data(am_url)

        if not avg_scores_dict or not plddt_scores_dict:
            logging.error(f"Could not retrieve valid scores for {uniprot_id}. Exiting.")
            continue

        max_res = max(max(avg_scores_dict.keys()), max(plddt_scores_dict.keys()))
        residues = list(range(1, max_res + 1))

        am_scores_list = [avg_scores_dict.get(res) for res in residues]
        plddt_scores_list = [plddt_scores_dict.get(res) for res in residues]

        # Create a NumPy array for the modified PDB file, handling 1-based indexing
        am_scores_array = np.full(max_res + 1, np.nan)
        for res, score in avg_scores_dict.items():
            am_scores_array[res] = score

        # 4. Generate outputs
        logging.info("Generating line plot of AM vs. pLDDT scores...")
        plot_scores(am_scores_list, plddt_scores_list, uniprot_id, output_path)

        logging.info("Generating AlphaMissense heatmap...")
        plot_am_heatmap(am_data_complete, output_path, uniprot_id)

        logging.info("Creating modified PDB file with AM scores...")
        modified_pdb_path = f"{output_path}/{uniprot_id}_AM_scores.pdb"
        modify_pdb_with_am_scores(pdb_url, am_scores_array, modified_pdb_path)

        logging.info(
            f"Successfully completed report for {uniprot_id}. Outputs are in '{output_path}'."
        )
