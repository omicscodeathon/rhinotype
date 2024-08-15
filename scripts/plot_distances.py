import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pandas as pd

def plot_distances(distances_matrix):
    # Ensure the input is a DataFrame to leverage labels
    if not isinstance(distances_matrix, pd.DataFrame):
        distances_matrix = pd.DataFrame(distances_matrix)
    
    # Plotting the heatmap
    plt.figure(figsize=(18, 14))
    sns.heatmap(distances_matrix, cmap='YlOrRd', cbar=True, xticklabels=True, yticklabels=True)
    plt.title("Genetic distances between sequences")
    
    # Rotate the labels and adjust the font size for better readability
    plt.xticks(rotation=90, fontsize=4)
    plt.yticks(rotation=0, fontsize=4)
    
    path = os.path.join(os.path.dirname(__file__), '../figures/distances.png')
    plt.savefig(path)
    plt.show()
