import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import os
from getprototypeseqs import user_input

def plot_tree(distance_matrix):
    # Convert the data to a numpy array if it's not already
    distance_matrix = np.array(distance_matrix)
    
    # Convert the symmetric distance matrix to a condensed distance matrix
    condensed_distance_matrix = squareform(distance_matrix)
    
    # Perform hierarchical clustering using complete linkage
    hc = linkage(condensed_distance_matrix, method='complete')
    
    # Plot the dendrogram
    plt.figure(figsize=(15, 8))
    dendrogram(hc, leaf_rotation=90, leaf_font_size=8)
    plt.title(f"{user_input} simple tree")
    plt.xlabel("")
    plt.ylabel("Genetic distance")
    path = os.path.join(os.path.dirname(__file__), '../figures/tree.png')
    plt.savefig(path)
    plt.show()
