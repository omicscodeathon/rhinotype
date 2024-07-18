import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# read fasta function to read sequences from a FASTA file
def read_fasta(fasta_file):
    sequences = []
    headers = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        headers.append(record.id)
    return {"sequences": sequences, "headers": headers}

# Compare sequence function to return substitution type
def compare_sequences(seqA, seqB):
    differences = [i for i in range(len(seqA)) if seqA[i] != seqB[i]]
    subsType = [seqB[i] for i in differences]
    return pd.DataFrame({"position": differences, "subsType": subsType})

# SNPeek function to plot differences between sequences
def SNPeek(fastaData, showLegend=False):
    sequences = fastaData["sequences"]
    seqNames = fastaData["headers"]
    genomeLength = max([len(seq) for seq in sequences])

    # Define color map for nucleotides substitutions
    colorMap = {"A": "green", "T": "red", "C": "blue", "G": "yellow"}

    diffList = []
    # Map substitutions types to colors
    for i in range(1, len(sequences)):
        diff = compare_sequences(sequences[0], sequences[i])
        diff["color"] = diff["subsType"].map(colorMap).fillna("black")
        diffList.append(diff)

    # Get the maximum length of the sequences
    genomeLength = max(max([len(seq) for seq in sequences]), 2)

    # Plot the differences
    plt.figure(figsize=(15, 10))  # Increase figure size
    plt.xlim(1, genomeLength)
    plt.ylim(0.5, len(sequences))
    plt.xlabel(f"Genome Position of {seqNames[-1]}, acting as reference", fontsize=10)  # Adjust font size if necessary
    plt.yticks(ticks=np.arange(1, len(sequences) + 1), labels=seqNames, fontsize=10)  # Adjust tick font size
    plt.gca().yaxis.set_tick_params(labelsize=8)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    # Plot vertical lines for each difference
    for i, diff in enumerate(diffList, start=1):
        for _, row in diff.iterrows():
            plt.plot([row["position"], row["position"]], [i - 0.25, i + 0.25], color=row["color"], marker="|")

    # Add legend
    if showLegend:
        plt.legend(["A", "T", "C", "G", "Other"], 
                    handlelength=0.8, markerfirst=False, 
                    loc="upper left", bbox_to_anchor=(1, 1),
                    facecolor="white", framealpha=0.7)

    # Save the plot as a PNG file
    path = os.path.join(os.path.dirname(__file__), '../figures/SNPeek.png')
    plt.savefig(path)
    # Show the plot
    plt.show()
