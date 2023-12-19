import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
import numpy as np

dataset = 'India'
pheno = 'diet' #Phenotype to be tested, make sure to add all other phenotypes in features variable
k = 2

# Load the function abundance data
function_abundance_path = 'binary_function_df.tsv'
function_abundance_data = pd.read_csv(function_abundance_path, sep='\t', index_col=0)

# Load the metadata
metadata_path = 'sample_mapping.txt'
metadata = pd.read_csv(metadata_path, sep='\t')

# Set the index of the function abundance data to match the 'WGS' column in the metadata
function_abundance_data.index.name = 'WGS'
merged_data = metadata.merge(function_abundance_data, on='WGS', how='left')

# Extract the 'Sample' names for labeling
sample_labels = merged_data.reset_index()['Sample']

# Extracting the features for PCA (excluding '16S', 'WGS', and 'Host' columns)
features = merged_data.drop(columns=['16S', 'WGS', pheno, 'Sample', 'age', 'location'])

# Apply PCA
pca = PCA(n_components=2)
pca_results = pca.fit_transform(features)

# Create a DataFrame for the PCA results
pca_df = pd.DataFrame(pca_results, columns=['PCA1', 'PCA2'])
pca_df[pheno] = merged_data[pheno]
pca_df['Sample'] = sample_labels

# K-means
def cluster_purity(labels, cluster_assignments):
    unique_clusters = np.unique(cluster_assignments)
    total_purity = 0
    for cluster in unique_clusters:
        cluster_labels = labels[cluster_assignments == cluster]
        most_common = np.bincount(cluster_labels).argmax()
        purity = np.sum(cluster_labels == most_common) / len(cluster_labels)
        total_purity += purity * len(cluster_labels)
    return total_purity / len(labels)

num_clusters = k  # Modify this as needed
kmeans = KMeans(n_clusters=num_clusters, random_state=0)
cluster_assignments = kmeans.fit_predict(pca_results)

pheno_labels = pd.factorize(merged_data[pheno])[0]
purity = cluster_purity(pheno_labels, cluster_assignments)
print("Overall Cluster Purity:", purity)

# Plotting
plt.figure(figsize=(12, 8))
scatter = sns.scatterplot(x='PCA1', y='PCA2', hue=pheno, data=pca_df, legend="full", palette="tab10")
plt.title('WGS PCA Based on Relative Function Abundance (' + dataset + ', ' + pheno + ')')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')

# Adding sample labels to the points
for line in range(0,pca_df.shape[0]):
     scatter.text(pca_df.PCA1[line]+0.02, pca_df.PCA2[line], pca_df.Sample[line], horizontalalignment='left', size='medium', color='black')

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Save the plot to a file
plt.savefig(pheno + '_WGS_plot.png', bbox_inches='tight')
