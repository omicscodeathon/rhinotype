import os
from getprototypeseqs import getprototypeseqs
from readfasta import read_fasta
from SNPeek import SNPeek
from genetic_distances import calc_p_distance, calc_jukes_cantor_distance, calc_kimura_2p_distance, calc_tamura_3p_distance
from pairwise_distances import pairwise_distances
from assign_types import assign_types
from plot_distances import plot_distances
from plot_frequency import plot_frequency

#getprototypeseqs()

path = os.path.join(os.path.dirname(__file__), '../data/test.fasta')

#test = read_fasta(fasta_file = path)

#print(test)

#fasta_file = read_fasta(fasta_file = path)

#SNPeek(fasta_file)

fasta_data = read_fasta(path)

# Calculate different distances
#p_distances = calc_p_distance(fasta_data)
#jc_distances = calc_jukes_cantor_distance(fasta_data)
#k2p_distances = calc_kimura_2p_distance(fasta_data)
#tamura_distances = calc_tamura_3p_distance(fasta_data)
assigned_types_df = assign_types(fasta_data, model="JC")
#result = pairwise_distances(fasta_data, gap_deletion=True)

#print("p-Distances:\n", p_distances)
#print("Jukes-Cantor Distances:\n", jc_distances)
#print("Kimura 2-Parameter Distances:\n", k2p_distances)
#print("Tamura 3-Parameter Distances:\n", tamura_distances)
print("assigned_type_df:\n", assigned_types_df)
plot_frequency(assigned_types_df, show_legend=True)
#print(result)
#plot_distances(result)