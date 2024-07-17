import numpy as np
import pandas as pd
from readfasta import read_fasta

def delete_missing_data_sites(seqs):
    # Convert list of sequences into a matrix where each sequence is a row
    seq_matrix = np.array([list(seq) for seq in seqs])
    
    # Find columns that do not contain gaps (assuming "-" as the gap symbol)
    valid_columns = np.all(seq_matrix != '-', axis=0)
    
    # Extract these columns to create a cleaned matrix
    cleaned_matrix = seq_matrix[:, valid_columns]
    
    # Convert cleaned matrix back to a list of sequences
    cleaned_seqs = [''.join(row) for row in cleaned_matrix]
    
    return cleaned_seqs

def count_snps_helper(fasta_data, gap_deletion=True):
    refs = fasta_data['sequences']
    ref_headers = fasta_data['headers']
    
    # Optionally remove sites with missing data
    if gap_deletion:
        refs = delete_missing_data_sites(refs)
    
    # Convert all cleaned sequences to a matrix of character vectors
    seq_matrix = np.array([list(seq) for seq in refs])
    
    # Initialize the SNP matrix
    n = len(seq_matrix)
    snp_matrix = np.full((n, n), np.nan)
    
    # Loop over all pairs of sequences using matrix indices
    for i in range(n):
        for j in range(i, n):
            seq1 = seq_matrix[i, :]
            seq2 = seq_matrix[j, :]
            
            # Calculate SNPs
            if len(seq1) == len(seq2) and len(seq1) > 0:
                snps = np.sum(seq1 != seq2)
                snp_matrix[i, j] = snp_matrix[j, i] = snps
            else:
                snp_matrix[i, j] = snp_matrix[j, i] = np.nan
    
    return snp_matrix

def calc_p_distance(fasta_data, gap_deletion=True):
    # Count the SNPs between each query and each reference sequence
    snp_counts = count_snps_helper(fasta_data, gap_deletion=gap_deletion)
    
    # Extract query headers
    query_headers = ref_headers = fasta_data['headers']
    
    # Directly calculate reference sequence lengths
    refs = fasta_data['sequences']
    
    # Optionally remove sites with missing data
    if gap_deletion:
        refs = delete_missing_data_sites(refs)
    
    ref_lengths = np.array([len(seq) for seq in refs])
    
    # Prepare a matrix for p-distances with appropriate dimensions and names
    p_distances_matrix = np.empty((len(query_headers), len(ref_lengths)))
    p_distances_matrix[:] = np.nan
    p_distances_matrix = pd.DataFrame(p_distances_matrix, index=query_headers, columns=ref_headers)
    
    # Calculate p-distance for each query-reference pair
    for q in range(len(snp_counts)):
        for i in range(len(ref_lengths)):
            if not np.isnan(snp_counts[q, i]):
                p_distances_matrix.iloc[q, i] = snp_counts[q, i] / ref_lengths[i]
    
    return p_distances_matrix

def calc_jukes_cantor_distance(fasta_data, gap_deletion=True):
    # Calculate p-distance for multiple queries
    p_dist = calc_p_distance(fasta_data, gap_deletion=gap_deletion)
    
    # Initialize a matrix to store Jukes-Cantor distances
    jc_dist = np.zeros_like(p_dist)
    
    # Apply the Jukes-Cantor formula to each element in the p-distance matrix
    jc_dist = -3/4 * np.log(1 - 4/3 * p_dist)
    
    # Handling cases where p_dist >= 0.75, setting JC distance to Inf
    # JC assumes that all nucleotide substitutions are equally probable and independent,
    # which might not always hold true in real data.
    # The handling of cases where p_dist >= 0.75 with Inf is to indicate that the Jukes-Cantor model
    # might not be valid for these high levels of divergence due to the assumption of the model being violated.
    jc_dist[p_dist >= 0.75] = np.inf
    
    # Return the Jukes-Cantor genetic distance matrix
    return jc_dist

# Function to help calculate proportions of transitions and transversions in Kimura2P and Tamura3p
def calc_transition_transversions(ref_chars, query_chars):
    transitions = sum((ref_chars == 'A') & (query_chars == 'G') |
                      (ref_chars == 'G') & (query_chars == 'A') |
                      (ref_chars == 'C') & (query_chars == 'T') |
                      (ref_chars == 'T') & (query_chars == 'C'))

    transversions = sum((ref_chars == 'A') & (query_chars == 'C') |
                        (ref_chars == 'A') & (query_chars == 'T') |
                        (ref_chars == 'G') & (query_chars == 'C') |
                        (ref_chars == 'G') & (query_chars == 'T') |
                        (ref_chars == 'C') & (query_chars == 'A') |
                        (ref_chars == 'C') & (query_chars == 'G') |
                        (ref_chars == 'T') & (query_chars == 'A') |
                        (ref_chars == 'T') & (query_chars == 'G'))

    total_sites = len(ref_chars)
    P = transitions / total_sites
    Q = transversions / total_sites

    return {'P': P, 'Q': Q}

def calc_kimura_2p_distance(fasta_data, gap_deletion=True):
    """Calculates Kimura 2-parameter distance between all pairs of sequences."""
    refs = queries = fasta_data['sequences']
    ref_headers = query_headers = fasta_data['headers']
    
    if gap_deletion:
        refs = queries = delete_missing_data_sites(refs)
        
    k2p_matrix = np.full((len(queries), len(refs)), np.nan)
    
    for q in range(len(queries)):
        query_chars = np.array(list(queries[q]))
        
        for r in range(len(refs)):
            ref_chars = np.array(list(refs[r]))
            
            if len(query_chars) == len(ref_chars):
                tt = calc_transition_transversions(ref_chars, query_chars)
                P = tt['P']
                Q = tt['Q']
                
                if (1 - 2 * P - Q) > 0 and (1 - 2 * Q) > 0:
                    K2P_distance = -0.5 * np.log((1 - 2 * P - Q) * np.sqrt(1 - 2 * Q))
                else:
                    K2P_distance = np.nan
                
                k2p_matrix[q, r] = K2P_distance
            else:
                k2p_matrix[q, r] = np.nan
    
    return k2p_matrix

def calc_tamura_3p_distance(fasta_data, gap_deletion=True):
    """Calculates Tamura 3-parameter distance between all pairs of sequences."""
    refs = queries = fasta_data['sequences']
    ref_headers = query_headers = fasta_data['headers']
    
    if gap_deletion:
        refs = queries = delete_missing_data_sites(refs)
        
    tamura_matrix = np.full((len(queries), len(refs)), np.nan)
    
    for q in range(len(queries)):
        query_chars = np.array(list(queries[q]))
        
        for r in range(len(refs)):
            ref_chars = np.array(list(refs[r]))
            
            if len(query_chars) == len(ref_chars):
                gc_content_ref = np.sum(np.isin(ref_chars, ['G', 'C'])) / len(ref_chars)
                gc_content_query = np.sum(np.isin(query_chars, ['G', 'C'])) / len(query_chars)
                tt_results = calc_transition_transversions(ref_chars, query_chars)
                P = tt_results['P']
                Q = tt_results['Q']
                theta1 = gc_content_ref
                theta2 = gc_content_query
                C = theta1 + theta2 - 2 * theta1 * theta2
                
                if (1 - P / C - Q) > 0 and (1 - 2 * Q) > 0:
                    distance = -C * np.log(1 - P / C - Q) - 0.5 * (1 - C) * np.log(1 - 2 * Q)
                else:
                    distance = np.nan
                
                tamura_matrix[q, r] = distance
            else:
                tamura_matrix[q, r] = np.nan
    
    return tamura_matrix
