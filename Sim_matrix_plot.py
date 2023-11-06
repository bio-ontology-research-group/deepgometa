import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

domain = 'mf'

##############

###Run this if Sim output doesn't have protein names
##Create empty protein name matrix
protein_file = domain + '_specific_combined.txt'
protein_names = [line.split()[0] for line in open(protein_file, 'r') if line.strip()]
protein_df = pd.DataFrame(index = protein_names, columns = protein_names)

##Read Sim.groovy results
sim_file = domain + '_specific_SimOut.txt'
sim_df = pd.read_csv(sim_file, sep = " ", header = None)

##Add protein names to Sim
sim_df.index = protein_names
sim_df.columns = protein_names
sim_df.to_csv(domain + '_specific_SimOut.txt', sep = '\t')

##############

###Run this if Sim output already has protein names
# sim_file = domain + '_specific_SimOut.txt'
# sim_df = pd.read_csv(sim_file, sep = "\t", index_col = 0)

##############

# Apply PCA
pca = PCA(n_components = 2)
PCs = pca.fit_transform(sim_df)
subPC = pd.DataFrame(data = PCs, columns = ['PC1', 'PC2'], index = sim_df.index)
explained_variance = pca.explained_variance_ratio_

# Create a color mapping for 'aquatic' and 'terrestrial' labels
with open('protein_names_combined.txt', 'r') as file:
    protein_names = [line.strip().split('\t') for line in file]
color_map = {'aquatic': 'blue', 'terrestrial': 'brown'}
protein_labels_dict = {name: label for name, label in protein_names}
colors = [color_map.get(protein_labels_dict.get(protein, 'unknown'), 'gray') for protein in subPC.index]

# Plot PCA with colored points
plt.scatter(subPC['PC1'], subPC['PC2'], c = colors)
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f}%)')
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f}%)')
plt.title('PCA Plot of Protein Sequence Similarity')

legend_labels = ['Aquatic', 'Terrestrial'] #Legend
legend_colors = ['blue', 'brown']
legend_handles = [plt.Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = color, markersize = 10, label = label) for color, label in zip(legend_colors, legend_labels)]
plt.legend(handles = legend_handles, title = 'Protein Types', loc = 'upper right')

plt.savefig(domain + '_Sim_PCA.png') #Save plot
plt.close()
