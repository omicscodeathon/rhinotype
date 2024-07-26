import numpy as np
import math
from Bio import SeqIO
from genetic_distances import delete_missing_data_sites, count_snps_helper, calc_p_distance, calc_jukes_cantor_distance, calc_transition_transversions, calc_kimura_2p_distance, calc_tamura_3p_distance



# Calculate overall mean distance using p-distance
def overall_p_distance(fasta_data, gap_deletion=True):
    sequences = fasta_data['sequences']

    if gap_deletion:
        sequences = delete_missing_data_sites(sequences)

    num_sequences = len(sequences)
    total_distance = 0
    num_comparisons = 0

    for i in range(num_sequences - 1):
        for j in range(i + 1, num_sequences):
            seq_i_chars = np.array(list(sequences[i]))
            seq_j_chars = np.array(list(sequences[j]))
            distance = np.sum(seq_i_chars != seq_j_chars) / len(seq_i_chars)
            total_distance += distance
            num_comparisons += 1

    overall_p_distance = total_distance / num_comparisons
    return overall_p_distance

# Calculate overall mean distance using Jukes-Cantor model
def overall_jc_distance(fasta_data, gap_deletion=True):
    sequences = fasta_data['sequences']

    if gap_deletion:
        sequences = delete_missing_data_sites(sequences)

    num_sequences = len(sequences)
    total_jc_distance = 0
    num_comparisons = 0

    for i in range(num_sequences - 1):
        for j in range(i + 1, num_sequences):
            seq_i_chars = np.array(list(sequences[i]))
            seq_j_chars = np.array(list(sequences[j]))
            p_distance = np.sum(seq_i_chars != seq_j_chars) / len(seq_i_chars)

            if p_distance < 0.75:
                jc_distance = -3/4 * math.log(1 - 4/3 * p_distance)
                total_jc_distance += jc_distance
                num_comparisons += 1

    overall_mean_jc_distance = total_jc_distance / num_comparisons if num_comparisons > 0 else np.nan
    return overall_mean_jc_distance

# Calculate overall mean distance using Kimura 2 parameter model
def overall_k2p_distance(fasta_data, gap_deletion=True):
    sequences = fasta_data['sequences']

    if gap_deletion:
        sequences = delete_missing_data_sites(sequences)

    num_sequences = len(sequences)
    total_distance = 0
    num_comparisons = 0

    for i in range(num_sequences - 1):
        for j in range(i + 1, num_sequences):
            seq_i_chars = np.array(list(sequences[i]))
            seq_j_chars = np.array(list(sequences[j]))

            transitions = np.sum((seq_i_chars != seq_j_chars) & (
                ((seq_i_chars == 'A') & (seq_j_chars == 'G')) |
                ((seq_i_chars == 'G') & (seq_j_chars == 'A')) |
                ((seq_i_chars == 'C') & (seq_j_chars == 'T')) |
                ((seq_i_chars == 'T') & (seq_j_chars == 'C'))
            ))

            transversions = np.sum((seq_i_chars != seq_j_chars) & (
                ~((seq_i_chars == 'A') & (seq_j_chars == 'G')) &
                ~((seq_i_chars == 'G') & (seq_j_chars == 'A')) &
                ~((seq_i_chars == 'C') & (seq_j_chars == 'T')) &
                ~((seq_i_chars == 'T') & (seq_j_chars == 'C'))
            ))

            P = transitions / len(seq_i_chars)
            Q = transversions / len(seq_i_chars)

            if (1 - 2*P - Q) > 0 and (1 - 2*Q) > 0:
                k2p_distance = 0.5 * math.log(1 / (1 - 2*P - Q)) + 0.25 * math.log(1 / (1 - 2*Q))
                total_distance += k2p_distance
                num_comparisons += 1

    overall_mean_distance = total_distance / num_comparisons if num_comparisons > 0 else np.nan
    return overall_mean_distance

# Calculate overall mean distance using Tamura 3 parameter model
def overall_t3p_distance(fasta_data, gap_deletion=True):
    sequences = fasta_data['sequences']

    if gap_deletion:
        sequences = delete_missing_data_sites(sequences)

    num_sequences = len(sequences)
    total_distance = 0
    num_comparisons = 0

    for i in range(num_sequences - 1):
        for j in range(i + 1, num_sequences):
            seq_i_chars = np.array(list(sequences[i]))
            seq_j_chars = np.array(list(sequences[j]))

            transitions = np.sum((seq_i_chars != seq_j_chars) & (
                ((seq_i_chars == 'A') & (seq_j_chars == 'G')) |
                ((seq_i_chars == 'G') & (seq_j_chars == 'A')) |
                ((seq_i_chars == 'C') & (seq_j_chars == 'T')) |
                ((seq_i_chars == 'T') & (seq_j_chars == 'C'))
            ))

            transversions = np.sum((seq_i_chars != seq_j_chars) & (
                ~((seq_i_chars == 'A') & (seq_j_chars == 'G')) &
                ~((seq_i_chars == 'G') & (seq_j_chars == 'A')) &
                ~((seq_i_chars == 'C') & (seq_j_chars == 'T')) &
                ~((seq_i_chars == 'T') & (seq_j_chars == 'C'))
            ))

            P = transitions / len(seq_i_chars)
            Q = transversions / len(seq_i_chars)
            gc_content = np.sum(np.isin(seq_i_chars, ['G', 'C'])) / len(seq_i_chars)

            if (1 - 2*P - Q) > 0 and (1 - 2*Q) > 0:
                G = 1 / (1 + (1 - 2*Q) * (gc_content / (1 - gc_content)))
                t3p_distance = (0.5 * math.log(1 / (1 - 2*P - Q))) + (G * 0.5 * math.log(1 / (1 - 2*Q)))
                total_distance += t3p_distance
                num_comparisons += 1

    overall_mean_distance = total_distance / num_comparisons if num_comparisons > 0 else np.nan
    return overall_mean_distance

# Main function to calculate overall mean distance based on the chosen model
def overall_mean_distance(fasta_data, model='p-distance', gap_deletion=True):
    if model == "p-distance":
        result = overall_p_distance(fasta_data, gap_deletion)
    elif model == "JC":
        result = overall_jc_distance(fasta_data, gap_deletion)
    elif model == "Kimura2p":
        result = overall_k2p_distance(fasta_data, gap_deletion)
    elif model == "Tamura3p":
        result = overall_t3p_distance(fasta_data, gap_deletion)
    else:
        raise ValueError("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p'")

    return result