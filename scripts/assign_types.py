import os
import pandas as pd
import numpy as np
from pairwise_distances import pairwise_distances
from getprototypeseqs import user_input


def assign_types(fasta_data, model='p-distance', gap_deletion=True, threshold=0.105):
    # Read prototype sequences
    try:
        # Read prototype sequences
        if user_input.capitalize() == 'Vp1':
            path = os.path.join(os.path.dirname(__file__), '../data/vp1_test.csv')
        elif user_input.capitalize() == 'Vp4/2':
            path = os.path.join(os.path.dirname(__file__), '../data/prototypes.csv')
        else:
            raise Exception("Please provide either VP1 or VP4/2 prototypes csv")
        prototypes_df = pd.read_csv(path)
    except FileNotFoundError:
        raise Exception("Prototypes file not found. Please check the file path.")
    
    names_to_keep = prototypes_df['Accession'].tolist()
    
    
    # Run pairwiseDistances to calculate distances
    distances = pairwise_distances(fasta_data, model=model, gap_deletion=gap_deletion)
    
    # Filter columns based on the prototypes
    distances = distances.loc[:, distances.columns.isin(names_to_keep)]
    
    # Filter rows to remove the same names as in the list
    distances = distances.loc[~distances.index.isin(names_to_keep), :]
    
    # Initialize lists to store output data
    query_vec = []
    assigned_type_vec = []
    distance_vec = []
    ref_seq_vec = []
    
    # Iterate over each row (query) in the distances DataFrame
    for i, row in distances.iterrows():
        query_header = i
        valid_cols = row[row < threshold].index
        
        if len(valid_cols) == 0:
            # If no valid columns found, mark as "unassigned"
            assigned_type_vec.append("unassigned")
            distance_vec.append(np.nan)
            # Find the column with the minimum distance
            min_dist_col = row.idxmin()
            ref_seq_vec.append(min_dist_col)
        else:
            # Choose the one with the minimum distance
            min_distance_col = row[valid_cols].idxmin()
            assigned_type = min_distance_col
            assigned_type_vec.append(assigned_type.replace("RV", "").split("_")[-1])
            distance_vec.append(row[min_distance_col])
            ref_seq_vec.append(assigned_type)
        
        query_vec.append(query_header)
    
    # Create a DataFrame from the results
    output_df = pd.DataFrame({
        'query': query_vec,
        'assignedType': assigned_type_vec,
        'distance': distance_vec,
        'reference': ref_seq_vec
    })
    
    return output_df