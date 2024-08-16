import os
# Import scripts
from getprototypeseqs import getprototypeseqs, user_input
from readfasta import read_fasta
from SNPeek import SNPeek
from assign_types import assign_types
from pairwise_distances import pairwise_distances
from overall_mean_distance import overall_mean_distance
from count_SNPs import count_snp
from plot_frequency import plot_frequency
from plot_distances import plot_distances
from plot_tree import plot_tree
from plot_AA import plot_AA

# Function 1: getprototypeseqs
getprototypeseqs()

# Function 2: read_fasta

# Change mode either if user has VP1 sequence or VP4/2 sequence
def user_read_sequence():
    global fasta_data
    if user_input.capitalize() == "Vp1":
        fasta_file = os.path.join(os.path.dirname(__file__), "../data/vp1_test.fasta")
        fasta_data = read_fasta(fasta_file=fasta_file)
        print(fasta_data)
    elif user_input.capitalize() == "Vp4/2":
        fasta_file = os.path.join(os.path.dirname(__file__), "../data/test.fasta")
        fasta_data = read_fasta(fasta_file=fasta_file)
        print(fasta_data)
    else:
        print("Please provide either VP1 or VP4/2 sequence")


# Call user read sequence function
user_read_sequence()

# Function 3: SNPeek
SNPeek(fasta_data)

# Function 4: assign_types

# Change input_align to either vp1 or vp4/2 align sequences
def user_assign_types():
    global fasta_align
    if user_input.capitalize() == "Vp1":
        input_align = os.path.join(os.path.dirname(__file__), "../data/vp1_align.fasta")
        fasta_align = read_fasta(fasta_file=input_align)
    elif user_input.capitalize() == "Vp4/2":
        input_align = os.path.join(os.path.dirname(__file__), "../data/input_aln.fasta")
        fasta_align = read_fasta(fasta_file=input_align)
    else:
        print("Please provide either VP1 or VP4/2 sequence")


# Call user assign types function
user_assign_types()

# Assign genotypes according to the p-distance model and vp region provided
genotypes = assign_types(fasta_align, model="p-distance")

# Print first 5 genotypes
print(genotypes.head())

# Function 5: pairwise_distances
output = pairwise_distances(fasta_align, "p-distance", gap_deletion=True)

# Function 6: overall_mean_distance
overall_mean_distance(fasta_align, model="p-distance", gap_deletion=True)

# Function 7: count_snp
print(count_snp(fasta_data))

# Function 8: plot_frequency
genotypes = assign_types(fasta_align, model="p-distance", gap_deletion=True)
plot_frequency(genotypes)

# Function 9: plot_distances
distance_to_prototypes = pairwise_distances(fasta_align, "p-distance")
plot_distances(distance_to_prototypes)

# Function 10: plot_tree
plot_tree(distance_to_prototypes)

# Function 11: plot_AA

# User input for VP1 or VP4/2 sequence to pick translated sequence
def user_plot_AA():
    if user_input.capitalize() == "Vp1":
        test = os.path.join(os.path.dirname(__file__), "../data/vp1_test_translated.fasta")
        plot_AA(test)
    elif user_input.capitalize() == "Vp4/2":
        test = os.path.join(os.path.dirname(__file__), "../data/test.translated.fasta")
        plot_AA(test)
    else:
        print("Please provide either VP1 or VP4/2 translated sequence")

# Call user plot AA function
user_plot_AA()
