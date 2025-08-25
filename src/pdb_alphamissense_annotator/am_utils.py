"""
This module provides utility functions for fetching and processing AlphaMissense data.

It includes functions to retrieve data from the AlphaFold DB API, process the
AlphaMissense CSV file, and calculate per-residue average pathogenicity scores.
"""

import requests
import logging
import pandas as pd
import numpy as np
from typing import Tuple, Optional, Dict


def get_alphamissense_info(
    uniprot_id: str,
) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Retrieves the AlphaMissense data URL and protein description from the AlphaFold DB API.

    Args:
        uniprot_id: The UniProt accession ID.

    Returns:
        A tuple containing the UniProt description and the URL to the AlphaMissense CSV file.
        Returns (None, None) on failure.
    """
    afdb_api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"
    logging.info(f"Retrieving AlphaMissense data for {uniprot_id} from {afdb_api_url}")

    try:
        response = requests.get(afdb_api_url)
        response.raise_for_status()
        data = response.json()

        if not data or not isinstance(data, list):
            logging.warning(f"API returned empty or unexpected data for {uniprot_id}")
            return None, None, None, None

        uniprot_description = data[0].get("uniprotDescription")
        am_url = data[0].get("amAnnotationsUrl")
        pdbUrl = data[0].get("pdbUrl")
        afdb_sequence = data[0].get("sequence")

        if am_url is None:
            logging.warning(
                f"{uniprot_id} does not have AlphaMissense data in the AlphaFold database"
            )

        return uniprot_description, afdb_sequence, am_url, pdbUrl

    except requests.exceptions.RequestException as e:
        logging.error(f"Network or HTTP error for {uniprot_id}: {e}")
        return None, None, None, None
    except (ValueError, IndexError) as e:
        logging.error(f"Error parsing API response for {uniprot_id}: {e}")
        return None, None, None, None


def process_alphamissense_data(am_url: str) -> Optional[pd.DataFrame]:
    """
    Downloads and processes an AlphaMissense CSV file into a structured DataFrame.

    Args:
        am_url: The URL to the AlphaMissense CSV file.

    Returns:
        A pandas DataFrame with parsed and cleaned AlphaMissense data, or None on failure.
    """
    if not am_url:
        return None

    logging.info(f"Processing AlphaMissense data from {am_url}")

    try:
        am_df = pd.read_csv(am_url)
    except Exception as e:
        logging.error(
            f"Error reading or processing AlphaMissense CSV from {am_url}: {e}"
        )
        return None

    # Efficiently parse columns
    parsed_df = pd.DataFrame(
        {
            "residue_number": pd.to_numeric(
                am_df["protein_variant"].str.extract(r"(\d+)", expand=False),
                errors="coerce",
            ),
            "reference_aa": am_df["protein_variant"].str.extract(
                r"^([A-Z])", expand=False
            ),
            "alternative_aa": am_df["protein_variant"].str.extract(
                r"([A-Z])", expand=False
            ),
            "pathogenicity_score": pd.to_numeric(
                am_df["am_pathogenicity"], errors="coerce"
            ),
        }
    )

    parsed_df.dropna(inplace=True)
    parsed_df["residue_number"] = parsed_df["residue_number"].astype(int)

    return parsed_df


def calculate_average_pathogenicity(
    am_data: Optional[pd.DataFrame],
) -> Optional[Dict[int, float]]:
    """
    Calculates the average pathogenicity score for each residue.

    Args:
        am_data: A DataFrame containing AlphaMissense data with 'residue_number' and 'pathogenicity_score'.

    Returns:
        A dictionary mapping residue numbers to their average pathogenicity score.
        Returns None if the input data is invalid.
    """
    if am_data is None or am_data.empty:
        logging.warning(
            "Cannot calculate averages: No AlphaMissense data provided or data is empty."
        )
        return None

    # Group by residue number and calculate the mean pathogenicity score
    average_scores = (
        am_data.groupby("residue_number")["pathogenicity_score"]
        .mean()
        .round(3)
        .to_dict()
    )

    return average_scores


def extract_am_data(am_url):
    """Extracts AlphaMissense data from the URL and saves it as a CSV file.

    Args:
        am_url (str): URL to the AlphaMissense data CSV.

    Returns:
        pandas.DataFrame or None: The AlphaMissense data, or None on error.
    """
    try:
        am_file = pd.read_csv(am_url)
    except Exception as e:
        logging.error(f"Error reading AlphaMissense file: {e}")
        return None

    reference_aa = am_file["protein_variant"].str.extract(r"^([A-Z])")[0]
    alternative_aa = am_file["protein_variant"].str.extract("([A-Z])$")[0]
    residue_number = pd.to_numeric(
        am_file["protein_variant"].str.extract(r"([0-9]+)")[0]
    )
    pathogenicity_score = pd.to_numeric(am_file["am_pathogenicity"])

    am_data = pd.DataFrame(
        {
            "reference_aa": reference_aa,
            "residue_number": residue_number,
            "alternative_aa": alternative_aa,
            "pathogenicity_score": pathogenicity_score,
        }
    )

    return am_data


def extract_plddt_scores_from_url(pdb_file_url: str) -> Dict[int, float]:
    """
    Extracts pLDDT scores from an AlphaFold PDB file URL.
    The pLDDT score is stored in the B-factor column.

    Returns:
        A dictionary mapping residue numbers to their pLDDT score.
    """
    plddt_scores = {}
    if not pdb_file_url:
        return {}

    try:
        response = requests.get(pdb_file_url, timeout=10)
        response.raise_for_status()
        pdb_content = response.text.splitlines()

        for line in pdb_content:
            if line.startswith(("ATOM", "HETATM")):
                residue_number = int(line[22:26].strip())
                plddt_score = float(line[60:66].strip())
                plddt_scores[residue_number] = plddt_score

    except requests.exceptions.RequestException as e:
        logging.error(f"Error fetching PDB data: {e}")
    except ValueError:
        logging.error("Error parsing pLDDT scores in PDB file.")

    return plddt_scores
