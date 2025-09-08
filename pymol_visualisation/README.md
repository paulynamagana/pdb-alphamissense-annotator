# AlphaMissense Coloring Script for PyMOL
This script visualizes AlphaMissense scores in PyMOL by coloring structures using a continuous gradient:

- Blue for scores near 0 (benign)
- Grey around 0.5 (uncertain)
- Red for scores near 1 (pathogenic)

The scores must be pre-loaded into the B-factor column of your PDB or structure file.


## Step-by-Step Instructions

### 1. Load Your Structure into PyMOL
Launch PyMOL and load your structure


### 2. Load the Coloring Script from GitHub (no download needed)
In PyMOL, run this command to load the script directly:

```
run https://raw.githubusercontent.com/paulynamagana//pdb-alphamissense-annotator/visualisation/coloram.py
```

### 4. Run the Coloring Function
To apply the gradient coloring based on B-factors:

```
coloram
```

## That's it!