import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

def plot_tree(distance_matrix):
    # Convert the data to a numpy array if it's not already
    distance_matrix = np.array(distance_matrix)
    
    # Perform hierarchical clustering using complete linkage
    hc = linkage(distance_matrix, method='complete')
    
    # Plot the dendrogram
    plt.figure()
    dendrogram(hc, leaf_rotation=90, leaf_font_size=8)
    plt.title("A simple tree")
    plt.xlabel("")
    plt.ylabel("Genetic distance")
    plt.show()
