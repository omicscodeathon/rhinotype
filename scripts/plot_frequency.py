import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_frequency(assigned_types_df, show_legend=False):
    # Add 'species' column based on the first letter of 'assignedType'
    assigned_types_df['species'] = assigned_types_df['assignedType'].str[0]

    # Save to csv file assigned types
    path = os.path.join(os.path.dirname(__file__), '../data/vp1_assigned_types.csv')
    assigned_types_df.to_csv(path, index=False)
    
    # Aggregate counts by assignedType
    types_counts = assigned_types_df.groupby('assignedType').size().reset_index(name='query')
    types_counts['label'] = types_counts['assignedType'] + ', ' + types_counts['query'].astype(str)
    types_counts['species'] = types_counts['assignedType'].str[0]
    
    # Transform the data frame
    types_counts = types_counts.assign(
        end_y=types_counts['query'].cumsum(),
        start_y=types_counts['query'].cumsum().shift(fill_value=0)
    )
    
    # Replace species "u" with "Other"
    types_counts.loc[types_counts['species'] == 'u', 'species'] = 'Other'
    
    # Define colors for each species
    color_map = {"A": "blue", "B": "red", "C": "green", "Other": "grey"}
    colors = types_counts['species'].map(color_map)
    
    # Plot the bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(types_counts['assignedType'], types_counts['query'], color=colors)
    plt.title("Frequency of Types")
    plt.xlabel("RV Type")
    plt.ylabel("Count")
    
    if show_legend:
        # Add legend
        plt.legend(handles=[plt.Rectangle((0,0),1,1, color=color_map[key]) for key in color_map],
                   labels=[key for key in color_map], title="Species", loc='upper right')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    path = os.path.join(os.path.dirname(__file__), '../figures/frequency.png')
    plt.savefig(path)
    plt.show()
