import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_distances(distances_matrix):
    # Convert the data to a matrix
    distances_matrix = np.array(distances_matrix)
    
    # Plotting the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(distances_matrix, cmap='YlOrRd', cbar=True)
    plt.show()
