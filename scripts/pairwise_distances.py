from genetic_distances import calc_p_distance, calc_jukes_cantor_distance, calc_kimura_2p_distance, calc_tamura_3p_distance

def pairwise_distances(fasta_data, model="p-distance", gap_deletion=True):
    """Calculates pairwise genetic distances using the specified model."""
    if model == "p-distance":
        result = calc_p_distance(fasta_data, gap_deletion=gap_deletion)
    elif model == "JC":
        result = calc_jukes_cantor_distance(fasta_data, gap_deletion=gap_deletion)
    elif model == "Kimura2p":
        result = calc_kimura_2p_distance(fasta_data, gap_deletion=gap_deletion)
    elif model == "Tamura3p":
        result = calc_tamura_3p_distance(fasta_data, gap_deletion=gap_deletion)
    else:
        raise ValueError("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p'")
    
    return result