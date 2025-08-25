import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import pandas as pd
import matplotlib.ticker as ticker


def plot_scores(am_average_scores, plddt_scores, uniprot_id, output_path):
    """
    Plots pathogenicity and rescaled pLDDT scores against residue number using Seaborn.

    Args:
        pathogenicity_scores (list): List of pathogenicity scores.
        plddt_scores (list): List of pLDDT scores.
    """

    # Rescale pLDDT scores from 0-100 to 0-1
    plddt_rescaled = [score / 100 for score in plddt_scores]

    # Create a DataFrame for easy plotting with Seaborn
    data = {
        "Residue number": range(1, len(am_average_scores) + 1),
        "Average pathogenicity": am_average_scores,
        "pLDDT (rescaled)": plddt_rescaled,
    }
    df = pd.DataFrame(data)

    # Create line plots with Seaborn
    plt.figure(figsize=(12, 6))
    ax = plt.subplot(111)
    sns.lineplot(
        x="Residue number", y="Average pathogenicity", data=df, label="Average AM score"
    )
    sns.lineplot(
        x="Residue number", y="pLDDT (rescaled)", data=df, label="pLDDT (rescaled)"
    )

    # Add labels, title, and legend
    plt.xlabel("Residue number")
    plt.ylabel("Score")
    plt.title(f"Average AM and pLDDT scores per position ({uniprot_id})")
    plt.legend()
    plt.ylim(0, 1)
    plt.xlim(1, len(am_average_scores))

    if len(am_average_scores) > 200:
        ax.xaxis.set_major_locator(
            ticker.MultipleLocator(100)
        )  # Set marks every 100 residues
    else:
        ax.xaxis.set_major_locator(
            ticker.MultipleLocator(50)
        )  # Set marks every 100 residues

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        fancybox=True,
        shadow=True,
        ncol=5,
    )

    plt.grid(axis="y", linestyle="--")  # Optional grid

    # Save the plot to a file
    plt.tight_layout(rect=[0, 0.1, 1, 1])  # Leaves space below for the legend

    # Save the plot to a file
    output_file = os.path.join(output_path, f"graph_plDDT-AM-score_{uniprot_id}.png")
    plt.savefig(output_file)  # save file
    plt.close()


def plot_am_heatmap(am_data, output_path, uniprot_id):
    """
    Plots a heatmap for the AlphaMissense data with residue number on the primary x-axis
    and reference residues on the secondary x-axis, alternative_aa on y axis.
    Colour coded following the original colours from the AlphaFold Database, see: https://alphafold.ebi.ac.uk/entry/Q5VSL9
    """

    # custom palette
    def create_custom_colormap():
        cdict = {
            "red": [
                (0.0, 56 / 255, 56 / 255),
                (0.34, 204 / 255, 204 / 255),
                (0.464, 204 / 255, 204 / 255),
                (1.0, 165 / 255, 165 / 255),
            ],
            "green": [
                (0.0, 83 / 255, 83 / 255),
                (0.34, 204 / 255, 204 / 255),
                (0.464, 204 / 255, 204 / 255),
                (1.0, 13 / 255, 13 / 255),
            ],
            "blue": [
                (0.0, 163 / 255, 163 / 255),
                (0.34, 204 / 255, 204 / 255),
                (0.464, 204 / 255, 204 / 255),
                (1.0, 18 / 255, 18 / 255),
            ],
        }
        return LinearSegmentedColormap("CustomMap", segmentdata=cdict)

    # custom colormap
    custom_cmap = create_custom_colormap()

    # pivot table
    pivot_table = pd.pivot_table(
        am_data,
        values="pathogenicity_score",
        index="alternative_aa",
        columns="residue_number",
    )

    # Determine the length of the protein sequence
    sequence_length = pivot_table.shape[1]

    # Dynamically adjust the figure size.  You can adjust the scaling factors as needed.
    fig_width = max(10, sequence_length * 0.15)  # Minimum width of 10 inches
    fig_height = 6  # Keep height constant, or adjust as needed
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.heatmap(
        pivot_table,
        cmap=custom_cmap,
        vmin=0,
        vmax=1,
        cbar_kws={"label": "AlphaMissense score"},
    )  # Limits for the color scale

    ax.set_xlabel("Residue Number")
    ax.set_ylabel("Alternative Amino Acid")
    plt.title(f"AlphaMissense Pathogenicity Heatmap ({uniprot_id})")

    primary_xticks = range(0, pivot_table.shape[1], 20)  # Ticks for the primary x-axis
    ax.set_xticks(primary_xticks)
    ax.set_xticklabels(pivot_table.columns[primary_xticks])
    ax.set_facecolor("black")  # Set background black for matching AA
    plt.yticks(rotation=0)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks([i / 10.0 for i in range(11)])
    cbar.set_ticklabels([f"{i / 10.0:.1f}" for i in range(11)])

    # Create the second x-axis
    ax2 = ax.twiny()  # Create a twin Axes sharing the y-axis

    # Get unique residue numbers and corresponding reference amino acids, preserving order
    unique_residues = am_data[["residue_number", "reference_aa"]].drop_duplicates(
        subset=["residue_number"]
    )
    unique_residues = unique_residues.sort_values("residue_number")

    # Set the tick locations and labels for the top x-axis
    secondary_xticks = range(0, pivot_table.shape[1], 1)  # Show all ticks
    ax2.set_xticks(secondary_xticks)

    ref_aa_labels = [
        (
            unique_residues[unique_residues["residue_number"] == i][
                "reference_aa"
            ].values[0]
            if i in unique_residues["residue_number"].values
            else ""
        )
        for i in secondary_xticks
    ]

    ax2.set_xticklabels(ref_aa_labels)
    ax2.set_xlabel(" ")  # Label the top x-axis
    # Make the labels smaller
    ax2.tick_params(axis="x", labelsize="small")

    # Adjust layout to prevent labels from overlapping
    plt.tight_layout()

    output_file = os.path.join(output_path, f"AM_heatmap_{uniprot_id}.png")
    plt.savefig(output_file)
    plt.close()
