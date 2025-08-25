"""
sequence alignment, reporting
"""

from Bio import Align
from Bio.Seq import Seq
import logging
from typing import Tuple, List, Dict, Any
import os


def align_sequences(
    pdb_sequences: Dict[str, Seq], afdb_sequence: str, chain_id_list: List[str] | str
) -> Dict[str, Tuple[object, str, str]]:
    """
    Aligns multiple PDB-derived sequences to a single reference sequence, one by one.

    This function takes a dictionary of polypeptide sequences (keyed by chain ID)
    and performs a local alignment for each one against a single canonical sequence.

    Args:
        pdb_sequences: A dictionary mapping chain IDs to their Bio.Seq.Seq objects.
        afdb_sequence: The canonical sequence (e.g., from AlphaFold DB) to align against.
        chain_id_list: A list of chain IDs to align.

    Returns:
        A dictionary where keys are the original chain IDs and values are tuples
        containing the best alignment object, the aligned PDB sequence (with gaps),
        and the aligned reference sequence (with gaps).
    """

    if not pdb_sequences or not afdb_sequence:
        logging.warning("Empty sequences provided for alignment")
        return {}

    aligner = Align.PairwiseAligner()  # initialise the PairwiseAligner

    aligner.match_score = 1.0
    aligner.mismatch_score = (
        -100.0
    )  # very harsh penalty, as residues have to be with 100% identity to transfer AM annotatrions
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -1.0
    aligner.mode = "local"

    all_alignments = {}
    logging.info(
        f"Aligning {len(pdb_sequences)} PDB sequence(s) against the reference sequence."
    )

    for chain_id, pdb_seq in pdb_sequences.items():
        if chain_id in chain_id_list:
            try:
                alignments = aligner.align(str(afdb_sequence), str(pdb_seq))

                if not alignments:
                    logging.warning(
                        f"No alignment could be generated for chain {chain_id}."
                    )
                    continue

                best_alignment = alignments[0]
                afdb_aligned = best_alignment.aligned[0]
                pdb_aligned = best_alignment.aligned[0]

                # Store the results for this chain
                all_alignments[chain_id] = (best_alignment, pdb_aligned, afdb_aligned)

            except Exception as e:
                logging.error(
                    f"An error occurred during alignment for chain {chain_id}: {e}"
                )

    return all_alignments


def save_alignment_report(
    alignment_results: Dict[str, Tuple[object, str, str]], output_path: str
):
    """
    Saves a consolidated alignment report for one or more chains to a single file.

    Args:
        alignment_results (List[Dict[str, Any]]): A list of dictionaries, where each
            dictionary contains the results for one chain. Expected keys:
            'chain_id', 'alignment'.
        output_path (str): The path to save the report file.
    """
    if not alignment_results:
        logging.warning("No alignment results provided to save.")
        return

    # Ensure the output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w") as out_file:
        # About section
        out_file.write("=== AlphaMissense Alignment Report ===\n")
        out_file.write(
            "This report describes the alignment between a PDB-derived protein chain\n"
        )
        out_file.write("sequence and its canonical AlphaFold Database sequence.\n\n")
        out_file.write(
            "The alignment ensures residue-level correspondence so that AlphaMissense\n"
        )
        out_file.write(
            "pathogenicity scores can be accurately mapped to the structure (B-factors).\n"
        )

        for chain, result_tuple in alignment_results.items():
            if not result_tuple:
                logging.warning(
                    f"Skipping report for chain {chain} due to empty result."
                )
                continue

            alignment, pdb_aligned, afdb_aligned = result_tuple

            out_file.write(f"\n=== Alignment Report for Chain {chain} ===\n\n")
            out_file.write("=== Alignment Snippet ===\n")
            # Use the format() method for a cleaner, standard alignment view
            out_file.write(alignment.format())
            out_file.write("\n=== PDB Aligned ===\n")
            out_file.write(f"{str(pdb_aligned)}\n\n")

            out_file.write("\n=== AFDB Aligned ===\n")
            out_file.write(f"{str(afdb_aligned)}\n\n")

            out_file.write("\n=== Alignment coordinates ===\n")
            out_file.write(f"{alignment.coordinates}\n\n")

        out_file.write("\n=== End of Report ===\n\n")

    logging.info(f"Consolidated alignment report saved to {output_path}")
