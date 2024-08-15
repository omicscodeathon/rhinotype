import os
import pandas as pd
from genetic_distances import calc_p_distance, calc_jukes_cantor_distance, calc_kimura_2p_distance, calc_tamura_3p_distance
from pairwise_distances import pairwise_distances


def assign_types(fasta_data, model='p-distance', gap_deletion=True, threshold=0.105):
    # Read prototype sequences
    try:
        # Read prototype sequences
        path = os.path.join(os.path.dirname(__file__), '../data/vp1_test.csv')
        prototypes_df = pd.read_csv(path)
    except FileNotFoundError:
        raise Exception("Prototypes file not found. Please check the file path.")
    
    names_to_keep = prototypes_df['Accession'].tolist()
    
    
    # Calculate pairwise distances
    distance_matrix = pairwise_distances(fasta_data, model=model, gap_deletion=gap_deletion)

    # Get accession names from fasta data headers
    headers = fasta_data["headers"]
    
    # Convert to pandas DataFrame for easier manipulation
    distance_df = pd.DataFrame(distance_matrix, index=headers, columns=headers)

    # Filter columns based on the prototypes
    columns = distance_df.columns.to_list()
    trans_columns = [column.replace('RVB', 'B') for column in columns]
    trans_col_to_drop = [column for column in trans_columns if column not in names_to_keep]
    col_to_drop = [column.replace('B', 'RVB') for column in trans_col_to_drop]
    distance_df = distance_df.drop(columns=col_to_drop)
    
    
    # Filter rows to remove the same names as in the list
    rows = distance_df.index.to_list()
    trans_rows = [row.replace('RVB', 'B') for row in rows]
    trans_row_to_drop = [row for row in trans_rows if row in names_to_keep]
    row_to_drop = [row.replace('B', 'RVB') for row in trans_row_to_drop]
    distance_df = distance_df.drop(row_to_drop, axis=0, errors='ignore')

    
    # Initialize lists to store output data
    query_list = []
    assigned_type_list = []
    distance_list = []
    ref_seq_list = []
    
    # Iterate over each row (query) in the distances DataFrame
    for query in distance_df.index:
        valid_cols = distance_df.loc[query] < threshold
        
        if valid_cols.sum() == 0:
            # If no valid columns found, mark as "unassigned"
            assigned_type_list.append('unassigned')
            distance_list.append(None)

            # Find the column with the minimum distance
            min_dist_col = distance_df.loc[query].idxmin()
            ref_seq_list.append(min_dist_col)
        else:
            # Choose the one with the minimum distance
            min_distance_col = distance_df.loc[query, valid_cols].idxmin()
            assigned_type = min_distance_col
            assigned_type_list.append(assigned_type.split('_')[-1].replace('RV', ''))
            distance_list.append(distance_df.loc[query, min_distance_col])
            ref_seq_list.append(assigned_type)
        
        query_list.append(query)
    
    # Create output DataFrame
    output_df = pd.DataFrame({
        'query': query_list,
        'assignedType': assigned_type_list,
        'distance': distance_list,
        'reference': ref_seq_list
    })
    
    return output_df