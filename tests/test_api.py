"""
Tests for the AlphaFold DB API integration.
"""

import pytest
import sys
import os

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.pdb_alphamissense_annotator.am_utils import get_alphamissense_info


@pytest.mark.api
def test_get_alphamissense_info_api_call():
    """
    Tests the real API call to the AlphaFold DB to ensure the contract is met.
    This test requires an internet connection and may be slower.
    """
    # A well-known, stable UniProt ID (TP53)
    uniprot_id = "P04637"

    description, sequence, amurl, pdburl = get_alphamissense_info(uniprot_id)

    # Check that the API returns the expected data types
    assert isinstance(description, str)
    assert isinstance(sequence, str)
    assert isinstance(amurl, str)
    assert isinstance(pdburl, str)

    # Check that the returned values are not empty
    assert description
    assert sequence
    assert amurl
    assert pdburl

    # Check for expected content in the URL
    assert "AF-P04637-F1-aa-substitutions.csv" in amurl
