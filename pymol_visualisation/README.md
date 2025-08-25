# AlphaMissense colouring Script for PyMOL
This script visualises AlphaMissense scores in PyMOL by colouring structures using a continuous gradient:

- Blue for scores near 0 (benign)
- Grey around 0.5 (uncertain)
- Red for scores near 1 (pathogenic)

The scores must be pre-loaded into the B-factor column of your PDB or structure file.


## Step-by-Step Instructions

### 1. Load Your Structure into PyMOL
Launch PyMOL and load your structure


### 2. Load the colouring Script from GitHub (no download needed)
In PyMOL, run this command to load the script directly:

```
run https://raw.githubusercontent.com/paulynamagana/pdb-alphamissense-annotator/main/src/pymol_visualisation/coloram.py
```

### 4. Run the colouring Function
To apply the gradient colouring based on B-factors:

```
coloram
```

## That's it!
