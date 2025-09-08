# AlphaMissense PDB Annotator
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-311/)
[![CI](https://github.com/paulynamagana/pdb-alphamissense-annotator/actions/workflows/ci.yml/badge.svg)](https://github.com/paulynamagana/pdb-alphamissense-annotator/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



Annotate any PDB structure with AlphaMissense pathogenicity scores at the residue level.

* A "tool" that a friend of mine requested to add AlphaMissense average scores into the experimental PDB file obatained from X-ray crystallography. The tricky part is that the experimental structure was not complete and had missing regions.


## Features
- Input: Any PDB file (experimental or modeled).
- Input: UniProt accession (e.g. P04637 for TP53).
- Retrieves AlphaMissense scores via the AlphaFold Database (AFDB) API.
- Performs sequence alignment between the PDB chain(s) and the UniProt sequence (from AFDB).
- Assigns AlphaMissense scores to residues where alignment is valid.
- Adds scores in the B-factor 

Handles:
- Multiple chains in one PDB.
- Multiple copies of the same chain (you can annotate one UniProt against multiple chain IDs).
- Missing residues or mismatches (set to -1.00 by default).

Outputs:
- Annotated PDB file with AM scores in B-factors.
    If the input file contains multiple conformations, only the Bfactors in conformation A will be changes to AlphaMissense scores. 
- Alignment + coverage report.
- Modified PDB file from AFDB with AlphaMissense scores
- Plot of AM scores vd plddt
- AlphaMissense Heatmap from AFDB

## Installation
```bash
git clone https://github.com/paulynamagana/pdb-alphamissense-annotator.git
cd pdb-alphamissense-annotator
pip install -r requirements.txt
```

## Usage

### Annotate a PDB file with AlphaMissense scores

To run the annotation, use the `am2pdb` command. You need to provide the path to your PDB file, the UniProt ID for the protein, the chain(s) you want to annotate, and the output directory.

#### Command

```bash
python -m src.cli am2pdb <PDB_FILE_PATH> <UNIPROT_ID> <CHAINS> <OUTPUT_PATH>
```

#### Arguments

*   `<PDB_FILE_PATH>`: The path to the input PDB file you want to annotate.
*   `<UNIPROT_ID>`: The UniProt accession number for the protein (e.g., `Q92794`). The script will fetch AlphaMissense data for this ID.
*   `<CHAINS>`: The specific chain ID(s) from the PDB file to annotate. For multiple chains, provide them as a single comma-separated string (e.g., `"A,B"`).
*   `<OUTPUT_PATH>`: The directory where the annotated PDB file and the alignment report will be saved.

#### Example

```bash
python -m src.cli am2pdb examples/0411-01_DPFhx_2025_pdbredo-fixed_zn_charge-delete_oxt.pdb Q92794 "A" results/
```

This command will:
1.  Load the structure from `examples/example1.pdb`.
2.  Fetch AlphaMissense scores for UniProt ID `Q92794`.
3.  Align the UniProt sequence with chain `A` of the PDB file.
4.  Write a new PDB file with AlphaMissense scores in the B-factor column to the `results/` directory.
5.  Generate an alignment report in the `results/` directory.




### Color AFDB structure by AlphaMissense score

You can also generate a PDB file of an AlphaFold structure colored by its AlphaMissense score and create plots.

#### Command

```bash
python -m src.cli afdb-pdb-plots <UNIPROT_IDS> <OUTPUT_PATH>
```

#### Arguments

*   `<UNIPROT_IDS>`: A comma-separated list of UniProt IDs for which to generate colored PDBs and plots.
*   `<OUTPUT_PATH>`: The directory where the output files will be saved.

#### Example

```bash
python -m src.cli afdb-pdb-plots "P04637,Q92794" results/
```
This command will generate colored PDB files and plots for TP53 (`P04637`) and another protein (`Q92794`) and save them in the `results/` directory.

## Project structure
```
pdb-alphamissense-annotator/
│
├── src/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py                              # CLI entrypoint
│   ├── pdb_alphamissense_annotator/
│   │   ├── am_utils.py                     # AlphaMissense CSV parsing, per-residue aggregation
│   │   ├── main.py                         # orchestrates workflow
│   │   ├── pdb_utils.py                    # PDB parsing, B-factor resetting, writing
│   │   ├── plot_utils.py                   # Plotting utils
│   │   └── sequence_utils.py               # sequence extraction, alignment, chain reporting
│   └── pymol_visualisation/                # PyMOL script
│
├── tests/
│   ├── test_api.py
│
├── .github/workflows/
│   └── ci.yml                              # GitHub Actions: run pytest, linting
│
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt                        #  dependencies
```


## 📚 Citation
If you use AlphaMissense data, please cite:
Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science (2023).
