# Import scripts
from getprototypeseqs import getprototypeseqs
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
test = os.path.join("data", "vp1_test.fasta")
fasta_data = read_fasta(fasta_file = test)

print(fasta_data)

# Function 3: SNPeek
SNPeek(fasta_data)

# Function 4: assign_types
input_align = os.path.join("data", "vp1_align.fasta")
fasta_align = read_fasta(fasta_file=input_align)

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
test = os.path.join("data", "vp1_test_translated.fasta")
plot_AA(test)
>>>>>>> 29a0b85e853ced9a9812ec16314495c209164cab
