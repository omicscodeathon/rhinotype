import os

# Import scripts
from getprototypeseqs import getprototypeseqs
from readfasta import read_fasta
from SNPeek import SNPeek
# from assign_types import assign_types
from pairwise_distances import pairwise_distances
# from overall_mean_distance import overall_mean_distance
from count_SNPs import count_snp
# from plot_frequency import plot_frequency
from plot_distances import plot_distances
from plot_tree import plot_tree
from plot_AA import plot_AA

# Function 1: getprototypeseqs
getprototypeseqs()

# Function 2: read_fasta
test = read_fasta(fasta_file = "data/test.fasta")
print(test)

# Function 3: SNPeek
fasta_file = read_fasta(fasta_file="data/test.fasta")
SNPeek(fasta_file)

# Function 4: assign_types
# assign_types()

# Function 5: pairwise_distances
test = os.path.join("data", "input_aln.fasta")
fastaD = read_fasta(test)
output = pairwise_distances(fastaD, "p-distance", gap_deletion=True)

# Function 6: overall_mean_distance
# overall_mean_distance(output)

# Function 7: count_snp
test = os.path.join("data", "test.fasta")
fastaData = read_fasta(test)
count_snp(fastaData)

# Function 8: plot_frequency

# Function 9: plot_distances
test = os.path.join("data", "input_aln.fasta")
fastaData = read_fasta(test)
distance_to_prototypes = pairwise_distances(fastaData, "p-distance")
plot_distances(distance_to_prototypes)

# Function 10: plot_tree
plot_tree(distance_to_prototypes)

# Function 11: plot_AA
test = os.path.join("data", "test.translated.fasta")
plot_AA(test)
