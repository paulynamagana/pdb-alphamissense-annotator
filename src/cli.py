"""
CLI to run the annotations workflow
python -m src.cli am2pdb --pdb_file_path --uniprot_id --chains, --output_path"""

import os
import sys
import logging
import typer
from typing import Tuple, List, Dict, Any
from typing import Annotated

# -- Path setup for local imports
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

# create typer application
app = typer.Typer(help="A cli for adding AlphaMissense data")


@app.command()
def am2pdb(
    pdb_file_path: Annotated[str, typer.Argument(help="PDB file to annotate")],
    uniprot_id: Annotated[
        str,
        typer.Argument(
            help="UniProt ID from which you want to extract the AlphaMissense data"
        ),
    ],
    chains: Annotated[str, typer.Argument(help="chains or chain to annotate")],
    output_path: Annotated[str, typer.Argument(help="Path to the output files")],
):
    from src.pdb_alphamissense_annotator.main import annotate_pdb

    # Split the comma-separated string into a clean list of chain IDs
    chains_list = [chain.strip() for chain in chains.split(",")]

    annotate_pdb(pdb_file_path, uniprot_id, chains_list, output_path)
    logging.info(f"Finalised transfering AlphaMissense annotations")


@app.command()
def afdb_pdb_plots(
    uniprot_ids: Annotated[
        str,
        typer.Argument(help="uniprot Ids to get the pdb with AlphaMissense and pliots"),
    ],
    output_path: Annotated[str, typer.Argument(help="Path to the output files")],
):

    from src.pdb_alphamissense_annotator.main import afdb_structure_am_coloured

    uniprot_ids_list = [id.strip() for id in uniprot_ids.split(",")]
    afdb_structure_am_coloured(uniprot_ids_list, output_path)


if __name__ == "__main__":
    app()
