import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Sample data
data = {
    "ACGT-AAAG": [1.04, 1.16, 1.66, 1.02, 1.09, 2.52],
    "ACGT-ACGT": [1.25, 1.28, 2.37, 1.16, 1.3, 2.6],
    "AAAG-AAAG": [1.18, 1.22, 1.78, 1.02, 0.84, 2.55]
}
crops = ["Sorghum", "Fonio", "Chickpea", "Cassava", "Quinoa", "Pigeon Pea"]
df = pd.DataFrame(data, index=crops)

scaler = StandardScaler()
normalized_data = scaler.fit_transform(df)
df_normalized = pd.DataFrame(normalized_data, index=crops, columns=df.columns)

# Perform hierarchical clustering
linked = linkage(df_normalized, method='ward')

# Plot dendrogram
plt.figure(figsize=(10, 7))
dendrogram(linked, labels=df.index, leaf_rotation=45)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Crops")
plt.ylabel("Distance")
plt.show()

# Perform K-means clustering
kmeans = KMeans(n_clusters=3, n_init=10,random_state=0)
df['Cluster'] = kmeans.fit_predict(df_normalized)

# Print the clusters
print(df)

# Plot a heatmap with clustering
sns.clustermap(df_normalized, metric="euclidean", method="ward", cmap="viridis")
plt.title("Heatmap with Clustering")
plt.show()
